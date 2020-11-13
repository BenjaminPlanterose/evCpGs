############################################################################
############################################################################
###########                                                      ###########
###########             Data preparation (GSE61496)              ###########
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

library(FlowSorted.Blood.450k)
library(minfi)
library(GEOquery)
library(data.table)
library(ENmix)
library(gplots)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(genefilter)

## Load functions ##

# It assists to prepare batch ids from the chip sentrix_id
create.batch <- function(batch.f)
{
  batch <- rep(1, times = batch.f[1])
  for(i in 2:length(batch.f))
  {
    batch <- c(batch, rep(i, times = batch.f[i]))
  }
  return(batch)
}

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

## Technical QC ##

# Read IDAT in extended format (included nBead, detP and sdGreen/Red)
setwd('/media/ultron/2tb_disk2/RAW_DATA/2018/twins_project/GSE61496_external_validation/')
rgSet <- read.metharray.exp(getwd(), extended = T)

# Perform QC
setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/data/qc/Danish/')
qc <- QCinfo(rgSet)
bad_cpgs <- qc$badCpG
write.table(bad_cpgs, file = 'bad_cpgs.txt', row.names = F, col.names = F, quote = F)


## Check that sex has been properly assign ##

# Preprocess with StrQN
sqn <- preprocessQuantile(rgSet)

# Predict sex and plot
Sex <- getSex(sqn); plotSex(Sex)

# Check for sex
rgSamples <- unlist(lapply(X = strsplit(colnames(sqn), split = '_'), function(X) { paste(X[[1]])}))
pheno <- pheno[rgSamples,]
pheno <- as.data.frame(pheno)
pheno$Sex <- as.factor(pheno$Sex)
levels(pheno$Sex)[1] <- c('M')
levels(pheno$Sex)[2] <- c('F')
sex.mat <- table(Sex$predictedSex, pheno$Sex)
sex.mat <- sex.mat[,c(2,1)]
#     F   M
# F 142   0
# M   0 160
balloonplot(as.table(t(sex.mat)), xlab = 'Predicted', ylab = 'Registered', main = '')


#############################################   Preparation   #############################################

# Read IDATs
setwd('/media/ultron/2tb_disk2/RAW_DATA/2018/twins_project/GSE61496_external_validation/')
rgSet <- read.metharray.exp(getwd())

# Extract GEO ID
samples <- colnames(rgSet)
rgSamples <- unlist(lapply(X = strsplit(samples, split = '_'), function(X) { paste(X[[1]])}))
chips <- unlist(lapply(X = strsplit(samples, split = '_'), function(X) { paste(X[[2]])}))

#Download GDS file, put it in the current directory, and load it:
setwd('/media/ultron/2tb_disk2/RAW_DATA/2018/twins_project/GSE61496_external_validation/')
phenotype <- getGEO('GSE61496', destdir=".")
pheno <- phenotype[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)
dim(pheno) # 312 40
pheno <- pheno[,c(36, 38, 39)]
colnames(pheno) <- c('Age', 'Pair_id', 'Sex')
head(pheno)
pheno <- data.matrix(as.data.frame(pheno))

# Get rid of low quality + NAs in phenotype
indices <- match(c("GSM1506315", "GSM1506549", 'GSM1506321', 'GSM1506543', 
                   'GSM1506317', 'GSM1506323', 'GSM1506424', 'GSM1506430',
                   'GSM1506529', 'GSM1506535'), rownames(pheno))
pheno <- pheno[-indices, ]

# Probes to remove
setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/data/')
# Read CpGs to remove
# XY-probes
Y_probes <- unique(as.vector(read.table(file = "y_chromosome_probes.txt", header = F)$V1))
length(Y_probes) # 416
X_probes <- unique(as.vector(read.table(file = "x_chromosome_probes.txt", header = F)$V1))
length(X_probes) # 11232

# Bad QC - Obtained via QCinfo function from ENmix r-package (detPthre=0.000001, nbthre=3, samplethre=0.05, CpGthre=0.05, bisulthre=NULL, outlier=TRUE)
setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/data/qc/Danish/')
bad_probes <- unique(as.vector(read.table(file = "bad_cpgs.txt", header = F)$V1))
length(bad_probes) # 4545

# SNP probes - dbSNP v.147
data(SNPs.147CommonSingle)
f.SNP <- c(rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$Probe_maf >= 0.01],
           rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$CpG_maf > 0],
           rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$SBE_maf > 0])
SNP_probes <- na.omit(unique(f.SNP))
length(SNP_probes) # 99337

# CR probes
setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/data/')
CR_1 <- as.vector(read.table('crossreactive_Chen.txt', header = T)$TargetID) # Chen YA et al. Epigenetics. 2013 Feb;8(2):203-9. doi: 10.4161/epi.23470. Epub 2013 Jan 11.
kobor <- fread('GPL16304-47833.txt') # Price ME et al. Epigenetics Chromatin. 2013 Mar 3;6(1):4. doi: 10.1186/1756-8935-6-4.
CR_2 <- unique(c(kobor$ID[kobor$Autosomal_Hits == 'A_YES'], kobor$ID[kobor$XY_Hits == 'XY_YES']))
CR_probes <- unique(c(CR_1, CR_2))
length(CR_probes) # 41993
probes2remove <- unique(c(Y_probes, X_probes, SNP_probes, CR_probes, bad_probes))
length(probes2remove) # 139755

# Read and process isolated cell type references
data(FlowSorted.Blood.450k)
pheno2 <- colData(FlowSorted.Blood.450k)
indices <- which(is.na(match(FlowSorted.Blood.450k$CellType, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", 'WBC'))))
sum(is.na(indices))
FlowSorted.Blood.450k <- FlowSorted.Blood.450k[,-indices]
sample1 <- sampleNames(rgSet)
length(sample1) # 306
sample2 <- sampleNames(FlowSorted.Blood.450k)
pheno2 <- pheno2[sample2, ]

# Combine IDATs for combined normalization
RGSET <- combineArrays(rgSet, FlowSorted.Blood.450k)

coldata <- colData(RGSET)
rm(rgSet, FlowSorted.Blood.450k); gc()

######################### Preprocessing - SQN ##############################

# Perform normalisation
um.sqn <- preprocessQuantile(RGSET)

# Get rid of probes2remove
indices <- match(probes2remove, rownames(um.sqn))
sum(is.na(indices))
um.sqn <- um.sqn[-indices,]

# Apply offset of 100
nocombat.beta.sqn <- minfi::getBeta(um.sqn)
nocombat.CN.sqn <- 2^getCN(um.sqn)
nocombat.beta.sqn <- nocombat.beta.sqn*nocombat.CN.sqn/(nocombat.CN.sqn+100) # Apply offset of 100
rm(nocombat.CN.sqn, um.sqn) ; gc()

# Separate in two matrices (ComBat cannot be applied to the cell references because the chips are confounded with cell types)
nocombat.beta.sqn_cell <- nocombat.beta.sqn[,sample2]
nocombat.beta.sqn <- nocombat.beta.sqn[, sample1]

# Construct model matrix
index <- match(rgSamples, rownames(pheno))
age1 <- pheno[index, 'Age']
gender1 <- pheno[index, 'Sex']
df <- as.data.frame(t(nocombat.beta.sqn))
df$age <- age1
df$gender <- as.factor(gender1)
modcombat <- model.matrix(~age + gender, data=df)
rm(df) ; gc()

# Batch
# Create batch vector
batch <- create.batch(table(chips)[match(names(table(chips)), unique(chips))])

# ComBat
combat.beta.sqn <- sva::ComBat(dat=na.omit(nocombat.beta.sqn), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=F)

rm(nocombat.beta.sqn) ; gc()

########## QC ###############
sum(combat.beta.sqn < 0) # 379
sum(combat.beta.sqn > 1) # 112
sum(is.na(combat.beta.sqn)) # 0

########## QC ###############

# Cell composition correction
compData <- pickCompProbes_modified(coldata, cbind(combat.beta.sqn, nocombat.beta.sqn_cell), 
                                    cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"), compositeCellType = 'Blood', 
                                    probeSelect = 'both', numProbes = 50)
coefs <- compData$coefEsts
rm(nocombat.beta.sqn_cell) ; gc()
cell.counts <- minfi:::projectCellType(combat.beta.sqn[rownames(coefs),], coefs, lessThanOne = F, nonnegative = T)
cell.comp <- compData$compTable # Obtain isolated cell profile
cell.comp <- cell.comp[,c(-1, -9, -10, -11)]
cell.comp <- as.matrix(cell.comp)
sig.comp.cpg <- names(which(cell.comp[,'p.value'] < 1e-08)) # Cell composition significant CpGs
length(sig.comp.cpg) # 
std_comp <- colMeans(cell.counts)
std_mat <- matrix(rep(std_comp, times = nrow(cell.counts)), byrow = T, nrow = nrow(cell.counts))
std.delta.cell.counts <- t(cell.counts - std_mat) # create standard profile
beta.sqn <- cell.comp.correction(combat.beta.sqn, std.delta.cell.counts, sig.comp.cpg, cell.comp)

########## QC ###############
densityPlot(cell.counts)
sum(beta.sqn < 0) # 386
sum(beta.sqn > 1) # 129
########## QC ###############


########## Export ###############
setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/data/beta/Danish_cohort/')
fwrite(data.table(beta.sqn, keep.rownames = T), paste(Sys.Date(), 'sqn_combat_cellcomp.txt', sep = '_'), quote = F, row.names = T, col.names = T, sep = '\t', nThread = 4)
gc()
########## Export ###############