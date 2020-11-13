############################################################################
############################################################################
###########                                                      ###########
###########          Data preparation (GSE104812)                ###########
###########             Author: Benjamin Planterose              ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########             b.planterose@erasmusmc.nl                ###########
###########                                                      ###########
############################################################################
############################################################################

## Load libraries ##

library(data.table)
library(minfi)
library(Biobase)
library(GEOquery)
library(FlowSorted.Blood.450k)
library(genefilter)
library(gplots)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(genefilter)
library(scales)

## Load functions ##

# It will be used to perform discovery of whole-blood cell composition sensitive probes
# The original function is derived minfi:::pickCompProbes. It has been slightly modified 
# so that it can accept custom normalisations.
pickCompProbes_modified <- function (colData, beta_ref, cellTypes = NULL, numProbes = 50, compositeCellType = compositeCellType, 
                                     probeSelect = probeSelect) 
{
  splitit <- function(x) 
  {
    split(seq(along = x), x)
  }
  p <- beta_ref
  pd <- as.data.frame(colData)
  if (!is.null(cellTypes)) 
  {
    if (!all(cellTypes %in% pd$CellType)) stop("elements of argument 'cellTypes' is not part of 'mSet$CellType'")
    keep <- which(pd$CellType %in% cellTypes)
    pd <- pd[keep, ]
    p <- p[, keep]
  }
  pd$CellType <- factor(pd$CellType, levels = cellTypes)
  ffComp <- rowFtests(p, pd$CellType)
  prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[, i]))
  r <- matrixStats::rowRanges(p)
  compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
  names(compTable)[1] <- "Fstat"
  names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", "high", "range")
  tIndexes <- splitit(pd$CellType) # Indices for different cell types
  print(tIndexes)
  tstatList <- lapply(tIndexes, function(i) 
  {
    x <- rep(0, ncol(p))
    x[i] <- 1
    return(rowttests(p, factor(x)))
  })
  if (probeSelect == "any") 
  {
    probeList <- lapply(tstatList, function(x) 
    {
      y <- x[x[, "p.value"] < 1e-08, ]
      yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
      
      c(rownames(yAny)[1:(numProbes * 2)])
    })
  }
  
  else 
  {
    print('it when in')
    probeList <- lapply(tstatList, function(x) 
    {
      y <- x[x[, "p.value"] < 1e-08, ]
      yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
      yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
      
      c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
    })
  }
  trainingProbes <- unique(unlist(probeList))
  p <- p[trainingProbes, ]
  pMeans <- colMeans(p)
  names(pMeans) <- pd$CellType
  form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse = "+")))
  phenoDF <- as.data.frame(model.matrix( ~ pd$CellType - 1))
  colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
  if (ncol(phenoDF) == 2) 
  {
    X <- as.matrix(phenoDF)
    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
  }
  else 
  {
    tmp <- minfi:::validationCellType(Y = p, pheno = phenoDF, modelFix = form)
    coefEsts <- tmp$coefEsts
  }
  out <- list(coefEsts = coefEsts, compTable = compTable, sampleMeans = pMeans)
  return(out)
}

# It will perform cell composition correction of the beta value matrix
cell.comp.correction <- function(delta.beta, delta.cell.counts, sig.cpg, cell.comp)
{
  beta_comp <- matrix(rep(0, times = nrow(delta.beta)*nrow(delta.cell.counts)), nrow = nrow(delta.beta))
  rownames(beta_comp) <- rownames(delta.beta)
  colnames(beta_comp) <- rownames(delta.cell.counts)
  matching <- na.omit(match(sig.cpg, rownames(delta.beta)))
  beta_comp[matching,] <- cell.comp[matching,-1] # also erase p-value column
  beta.values.corrected <- delta.beta - beta_comp%*%delta.cell.counts
  
  return(beta.values.corrected)
}

#############################################   QC   #############################################

# Read pheno
setwd('/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/')
phenotype <- getGEO('GSE104812', destdir=".")
pheno <- phenotype[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)
pheno <- pheno[, c(33:34)]
colnames(pheno) = c("age", "sex")
pheno$age = as.numeric(pheno$age)
pheno$sex = as.factor(pheno$sex)

# Bad Samples/probes
setwd('/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/')
rgSet <- read.metharray.exp(getwd(), extended = T)
qc <- ENmix::QCinfo(rgSet)
qc$badsample # character(0)
bad.probes <- qc$badCpG
write.table(file = 'bad_probes.txt', x = bad.probes, quote = F, sep = '\t', col.names = F, row.names = F)

nbthre = 3
detPthre = 1e-06
qcmat <- qc$nbead < nbthre | qc$detP > detPthre
badValuePerSample <- apply(qcmat, 2, sum)/nrow(qcmat)
qc_df = data.frame(badFreq = badValuePerSample, bisul = qc$bisul)
thr = mean(qc_df$bisul) - 3*sd(qc_df$bisul)
dim(qc_df) # 48 2

setwd('/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/QC/')
tiff(filename = "children_aging.tiff", height = 427, width = 550, units = "px")
plot(qc_df$badFreq, qc_df$bisul, ylim = c(0, 35000), xlim = c(0, 0.055), pch = 19, col = alpha("black", 0.2), cex = 0.4,
     xlab = "Percent of low-quality data", ylab = "Average bisulfite conversion intensity", main = "Children population", las = 1, cex.axis = 0.8)
abline(v = 0.05, lty= 2, col = alpha("red3", 0.5))
abline(h = thr, lty = 2, col = alpha("red3", 0.5))
dev.off()


# Check Sex
sqn <- preprocessQuantile(rgSet)
Sex <- getSex(sqn)

setwd('/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/QC/')
tiff(filename = "sex1_children_aging.tiff", height = 762, width = 656, units = "px")
plot(x = Sex$xMed, y = Sex$yMed, type = "n", xlab = "X chr, median total intensity (log2)", 
     ylab = "Y chr, median total intensity (log2)")
text(x = Sex$xMed, y = Sex$yMed, labels = Sex$predictedSex, col = ifelse(Sex$predictedSex ==  "M", "deepskyblue", "deeppink3"))
legend("bottomleft", c("M", "F"), col = c("deepskyblue", "deeppink3"), pch = 16)
dev.off()


sex.mat <- table(Sex$predictedSex, pheno$sex) 
#     Female Male
# F     36    0
# M      0   21


rgSamples <- unlist(lapply(X = strsplit(colnames(sqn), split = '_'), function(X) { paste(X[[1]])}))
pheno <- pheno[rgSamples,]
pheno <- as.data.frame(pheno)
levels(pheno$sex) <- c('F', 'M')
sex.mat <- table(Sex$predictedSex, pheno$sex)

#   F   M
# F 19   0
# M   0 29

# Perform balloon plot
setwd('/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/QC/')
tiff(filename = "sex2_children_aging.tiff", height = 762, width = 656, units = "px")
balloonplot(as.table(t(sex.mat)), xlab = 'Predicted', ylab = 'Registered', main = '')
dev.off()



#############################################   Preparation   #############################################


# Read pheno
setwd('/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/')
phenotype <- getGEO('GSE104812', destdir=".")
pheno <- phenotype[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)
pheno <- pheno[, c(33:34)]
colnames(pheno) = c("age", "sex")
pheno$age = as.numeric(pheno$age)
pheno$sex = as.factor(pheno$sex)


# Probes to remove
data(SNPs.147CommonSingle)
f.SNP <- c(rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$Probe_maf >= 0.01],
           rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$CpG_maf > 0],
           rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$SBE_maf > 0])
SNP_probes <- na.omit(unique(f.SNP))
length(SNP_probes) # 99337
setwd("/media/ben/DATA/Ben/1_evCpGs/annotation/")
CR_1 <- as.vector(read.table('crossreactive_Chen.txt', header = F)$TargetID) # Chen YA et al. Epigenetics. 2013 Feb;8(2):203-9. doi: 10.4161/epi.23470. Epub 2013 Jan 11.
kobor <- fread('GPL16304-47833.txt') # Price ME et al. Epigenetics Chromatin. 2013 Mar 3;6(1):4. doi: 10.1186/1756-8935-6-4.
CR_2 <- unique(c(kobor$ID[kobor$Autosomal_Hits == 'A_YES'], kobor$ID[kobor$XY_Hits == 'XY_YES']))
CR_probes <- unique(c(CR_1, CR_2))
length(CR_probes) # 41937

setwd('/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/QC/')
bad_probes = as.vector(read.table("bad_probes.txt", header = F)$V1)
length(bad_probes) # 

# Read IDAT
setwd('/media/ben/DATA/Ben/1_evCpGs/data/aging_children/GSE104812_RAW/')
rgSet <- read.metharray.exp(getwd())

annotation = getAnnotation(rgSet)
Y_probes = rownames(annotation)[annotation$chr == "chrY"]
X_probes = rownames(annotation)[annotation$chr == "chrX"]
annotation = as.data.frame(annotation)
probes2remove <- unique(c(Y_probes, X_probes, SNP_probes, CR_probes, bad_probes))
length(probes2remove) # 138248


library(FlowSorted.Blood.450k)
# Read isolated cell type profiles
data(FlowSorted.Blood.450k)
pheno2 <- colData(FlowSorted.Blood.450k)
indices <- which(is.na(match(FlowSorted.Blood.450k$CellType, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", 'WBC'))))
sum(is.na(indices))
FlowSorted.Blood.450k <- FlowSorted.Blood.450k[,-indices]

sample1 <- sampleNames(rgSet)
length(sample1) # 48
sample2 <- sampleNames(FlowSorted.Blood.450k)
pheno2 <- pheno2[sample2, ]

RGSET <- combineArrays(rgSet, FlowSorted.Blood.450k)
coldata <- colData(RGSET)
rm(rgSet, FlowSorted.Blood.450k); gc()

######################### Preprocessing - SQN ##############################

# Perform normalisation
um.sqn <- preprocessQuantile(RGSET)

indices <- match(probes2remove, rownames(um.sqn))
sum(is.na(indices))
um.sqn <- um.sqn[-indices,]
dim(um.sqn)

nocombat.beta.sqn <- minfi::getBeta(um.sqn)
nocombat.CN.sqn <- 2^getCN(um.sqn)
nocombat.beta.sqn <- nocombat.beta.sqn*nocombat.CN.sqn/(nocombat.CN.sqn+100) # Apply offset of 100
rm(nocombat.CN.sqn, um.sqn) ; gc()


# Cell composition correction
compData <- pickCompProbes_modified(coldata, nocombat.beta.sqn, 
                                    cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"), compositeCellType = 'Blood', 
                                    probeSelect = 'both', numProbes = 50)
coefs <- compData$coefEsts
#rm(nocombat.beta.sqn_cell) ; gc()
cell.counts <- minfi:::projectCellType(nocombat.beta.sqn[rownames(coefs), sample1], coefs, lessThanOne = F, nonnegative = T)
cell.comp <- compData$compTable # Obtain isolated cell profile
cell.comp <- cell.comp[,c(-1, -9, -10, -11)]
cell.comp <- as.matrix(cell.comp)
sig.comp.cpg <- names(which(cell.comp[,'p.value'] < 1e-08)) # Cell composition significant CpGs
length(sig.comp.cpg) # 51490
std_comp <- colMeans(cell.counts)
std_mat <- matrix(rep(std_comp, times = nrow(cell.counts)), byrow = T, nrow = nrow(cell.counts))
std.delta.cell.counts <- t(cell.counts - std_mat) # create standard profile
beta.sqn <- cell.comp.correction(nocombat.beta.sqn[, sample1], std.delta.cell.counts, sig.comp.cpg, cell.comp)
dim(beta.sqn) # 347264     48

rm(compData, coefs, cell.counts, cell.comp, sig.comp.cpg, std_comp, std_mat, std.delta.cell.counts); gc()



########## Export ###############
setwd('/media/ben/DATA/Ben/1_evCpGs/data/aging_children/')
fwrite(data.table(beta.sqn, keep.rownames = T), paste(Sys.Date(), 'SQN_nocombat_cellcomp.txt', sep = '_'), quote = F, 
       row.names = T, col.names = T, sep = '\t', nThread = 4)
########## Export ###############

densityPlot(beta.sqn)
