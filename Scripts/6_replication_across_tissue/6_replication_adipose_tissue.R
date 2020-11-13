############################################################################
############################################################################
###########                                                      ###########
###########             Replication (adipose tissue)             ###########
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
library(genefilter)
library(data.table)
library(minfi)
library(parallel)
library(gplots)
library(minfi)
library(equivalence)
library(DescTools)
library(scales)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## Load functions ##

# Transforms data as imported from data.table into a data.frame
process.beta.fread <- function(beta)
{
  cpgs <- as.vector(beta$rn)
  beta <- beta[,-(1:2)]
  beta <- as.matrix(beta)
  rownames(beta) <- cpgs
  return(beta)
}

# It arranges the beta value matrix based on the twin design. The resulting matrix will have 
# columns such that:  {Twin 11, Twin 12, Twin 21, Twin 22, ...}
arrange.beta <- function(beta, matching)
{
  indices <- numeric(length = 2*nrow(matching))
  indices[seq(1, 2*nrow(matching), 2)] <- match(as.vector(matching$V1), colnames(beta))
  indices[seq(2, 2*nrow(matching), 2)] <- match(as.vector(matching$V2), colnames(beta))
  return(beta[,indices])
}

# Computes |deltaBeta| for all combinations of non-twin pairs
compute.delta.beta.notwin <- function(row)
{
  # Columns are assumed to be arranged by matching - t11 t12 t21 t22 t31 t32
  result <- numeric()
  l <- length(row)
  odd <- seq(1, l-2, by = 2)
  even <- seq(2, l-2, by = 2)
  
  result <- c(unlist(lapply(X = odd, FUN = function(X, row, l) {row[(X+2):l] - row[X]}, row, l)),
              unlist(lapply(X = even, FUN = function(X, row, l) {row[(X+1):l] - row[X]}, row, l)))
  return(result)
}

# Parallellized function that computes trimmed mean of |deltaBeta| for twins and non-twins
# It is employed to decide on the epsilon employed for equivalence testing
hvCpG.epsilon <- function(beta, matching)
{
  beta <- arrange.beta(beta, matching) # arrange so that twin L1 - Twin R1 - Twin L2 - Twin R2 - etc
  
  #np <- detectCores(logical = FALSE)
  np <- detectCores(logical = FALSE) - 8L
  cl <- makeCluster(np)
  clusterExport(cl, c('process.CpG2', 'compute.delta.beta.notwin', 'Trim'), envir=environment()) 
  
  r <- parLapply(cl = cl, X = 1:nrow(beta), fun = function(X) process.CpG2(beta[X,])) # process per row
  names(r) <- rownames(beta)
  
  stopCluster(cl)
  return(r)
}

# Subfunction employed by hvCpG.epsilon that actually computes trimmed mean of |deltaBeta| 
# for twins and non-twins
process.CpG2 <- function(row)
{
  twin <- row[seq(1, length(row), 2)] - row[seq(2, length(row), 2)]
  notwin <- compute.delta.beta.notwin(row)
  
  twin <- Trim(abs(twin), trim = 0.2)
  notwin <- Trim(abs(notwin), trim = 0.2)
  
  res <- c(mean(twin), mean(notwin))
  names(res) <- c('Twin', 'No twin')
  return(res)
}

# Subfunction employed by hvCpG.epsilon that actually performs evCpG discovery 
hvCpG.discovery <- function(beta, matching, epsilon)
{
  beta <- arrange.beta(beta, matching) # arrange so that twin L1 - Twin R1 - Twin L2 - Twin R2 - etc
  
  #np <- detectCores(logical = FALSE)
  np <- detectCores(logical = FALSE) - 8L
  cl <- makeCluster(np)
  clusterExport(cl, c('process.CpG', 'rtost', 'compute.delta.beta.notwin', 'epsilon'), envir=environment()) 
  
  r <- parLapply(cl = cl, X = 1:nrow(beta), fun = function(X) process.CpG(beta[X,])) # process per row
  names(r) <- rownames(beta)
  
  stopCluster(cl)
  return(r)
}

# Subfunction employed by hvCpG.discovery that actually performs equivalence testing based on
# a TOST that employs Yuen´s t-test
process.CpG <- function(row)
{
  twin <- row[seq(1, length(row), 2)] - row[seq(2, length(row), 2)]
  notwin <- compute.delta.beta.notwin(row)
  TOST <- rtost(x = abs(twin), y = abs(notwin), paired = F, epsilon = epsilon, 
                var.equal = FALSE, alpha = 0.05, tr = 0.2)
  res <- TOST$p.value
  
  return(res)
}


# 1) Prepare matching

setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
pheno = fread("E-MTAB-1866.sdrf.txt")
pheno = as.data.frame(pheno[,c(1,6,7,8,9,10, 11, 28)])
rownames(pheno) = pheno$`Factor Value[individual]`
dim(pheno) # 648  8
zig = pheno$`Characteristics[twin zygosity]`
names(zig) = pheno$`Factor Value[individual]`
MZ = names(which(zig == "monozygotic"))
length(MZ) # 240
pheno = pheno[(pheno$`Factor Value[individual]` %in% MZ) & (pheno$`Characteristics[co-twin]` != "not applicable"),]
rownames(pheno) = pheno$`Factor Value[individual]`
dim(pheno) # 196 8
colnames(matching) = c("V1", "V2")
sort(table(c(matching$V1, matching$V2)))
# TWPID17664 TWPID17665  TWPID4936  TWPID4937
#          1          1          1          1
# TWPID17665"  "TWPID4936" cotwins are not in pheno. Should have been labelled not applicable
pheno = pheno[!(rownames(pheno) %in% c("TWPID17665",  "TWPID4936")),]
matching = pheno[,4:5]
colnames(matching) = c("V1", "V2")
dim(matching) # 194   2
rownames(matching) = NULL
for(i in 1:nrow(matching))
{
  matching[i,] <- sort(matching[i,])
}
matching <- unique(matching)
dim(matching) # 97  2
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
write.table(matching, "matching.txt", row.names = F, col.names = F, quote = F)


# 2) Prepare dataset and compute IQRs

setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
data = fread("MuTHER_Fat_450K_norm_AE_030913.txt", header = T)
data = data[-1,]; gc()
CpGs = data$`Hybridization REF`
data = data[, -1]
data = data.matrix(data); gc()
rownames(data) = CpGs
dim(data) # 485577   648

iqrs_whole = rowIQRs(data, na.rm = T)
names(iqrs_whole) = CpGs
sum(iqrs_whole > 0.07)
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
write.table(iqrs_whole, "iqrs_whole.txt", col.names = F, row.names = T, quote = F)
data = data[, rownames(pheno)]; gc()
dim(data) # 485577    194

setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
fwrite(data.table(data, keep.rownames = T), paste(Sys.Date(), 'muther_fat_MZ.txt', sep = '_'), quote = F, 
       row.names = T, col.names = T, sep = '\t', nThread = 4)


# 3) Prepare probes2remove

setwd("/media/ben/DATA/Ben/1_evCpGs/data/test/")
rgSet = read.metharray.exp(getwd())
annotation = getAnnotation(rgSet)
Y_probes = rownames(annotation)[annotation$chr == "chrY"]
X_probes = rownames(annotation)[annotation$chr == "chrX"]
annotation = as.data.frame(annotation)
data(SNPs.147CommonSingle)
f.SNP <- c(rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$Probe_maf >= 0.01],
           rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$CpG_maf > 0],
           rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$SBE_maf > 0])
SNP_probes <- na.omit(unique(f.SNP))
length(SNP_probes) # 99337
setwd("/media/ben/DATA/Ben/1_evCpGs/annotation/")
CR_1 <- as.vector(read.table('crossreactive_Chen.txt', header = T)$TargetID) # Chen YA et al. Epigenetics. 2013 Feb;8(2):203-9. doi: 10.4161/epi.23470. Epub 2013 Jan 11.
kobor <- fread('GPL16304-47833.txt') # Price ME et al. Epigenetics Chromatin. 2013 Mar 3;6(1):4. doi: 10.1186/1756-8935-6-4.
CR_2 <- unique(c(kobor$ID[kobor$Autosomal_Hits == 'A_YES'], kobor$ID[kobor$XY_Hits == 'XY_YES']))
CR_probes <- unique(c(CR_1, CR_2))
length(CR_probes) # 41993
dataset <- fread('ICC_values.csv')
lowICC <- as.vector(dataset$`ilmnid(CpG site)`[dataset$`ICC value` < 0.37])
length(lowICC) # 270527
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
iqrs_whole <- read.table("iqrs_whole.txt")
name = iqrs_whole$V1
iqrs_whole = iqrs_whole$V2
names(iqrs_whole) = name
lowIQR = names(which(iqrs_whole < 0.07))
length(lowIQR) # 461859
probes2remove <- unique(c(Y_probes, X_probes, SNP_probes, CR_probes, lowICC, lowIQR))
length(probes2remove) # 480362

plot(density(iqrs_whole))
abline(v = 0.07)
sum(iqrs_whole > 0.07)

# 4) Read and prepare data
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
data = fread("2020-04-03_muther_fat_MZ.txt"); gc()
data = process.beta.fread(data); gc()
data = na.omit(data); gc()
SNP_beta = data[startsWith(rownames(data), "rs"),]
dim(SNP_beta) # 65 194
data = data[-which(startsWith(rownames(data), "rs")),]
dim(data) # 480320    194
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
matching = read.table("matching.txt", header = F)
dim(matching) # 97  2
probes2remove2 <- unique(c(Y_probes, X_probes, SNP_probes, CR_probes))
length(probes2remove2) # 137912
data = arrange.beta(data, matching)
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
pheno = fread("E-MTAB-1866.sdrf.txt")
pheno = as.data.frame(pheno[,c(1,6,7,8,9,10, 11, 28)])
rownames(pheno) = pheno$`Factor Value[individual]`
pheno = pheno[colnames(data),]

## QC
densityPlot(data)
densityPlot(delta)
mdsPlot(SNP_beta)
##

# Establishing an epsilon
X = data[!(rownames(data) %in% probes2remove),]; dim(X) # 8142   196
means1 <- hvCpG.epsilon(X, matching); gc()
twin <- unlist(lapply(means1, function(x) x[1]))
notwin <- unlist(lapply(means1, function(x) x[2]))
plot(density(twin), xlim = c(0,0.4), lwd = 2, col = 'chocolate', main = paste('SQN, mean(delta.beta.twin) =', round(mean(abs(twin)), 4), ', mean(delta.beta.unrelated) =', round(mean(abs(notwin)), 4)))
lines(density(notwin), lwd = 2, col = 'chocolate2')
abline(v = median(twin), lty = 2)
epsilon = median(twin); print(epsilon) # 0.04196169

# Replication
setwd("/media/ben/DATA/Ben/1_evCpGs/")
stochCpG <- as.vector(read.table(file = 'evCpGs.txt')$V1)
cross_stoch = stochCpG[stochCpG %in% rownames(data)]
length(cross_stoch) # 332
X = data[cross_stoch,]; dim(X) # 332   194
epsilon = 0.04196169
pvals <- hvCpG.discovery(X, matching, epsilon); gc()
pvals = unlist(pvals)

setwd("~/Desktop/")
x = pvals[stochCpG]
names(x) = stochCpG
write.table(x, "pvals_rep.txt", col.names = T, quote = F)

# Multiple testing correction
a.pvals = p.adjust(pvals, method = "bonferroni")
sum(a.pvals < 0.05) # 308
sig = names(which(a.pvals < 0.05))
sig2 = sig[!(sig %in% lowIQR)]; length(sig2) # 154
154/333 # 46.2%

# Visualize
setwd("~/Desktop/")
tiff(filename = "IQR_log10pval", res = 500, width = 2.7, height = 2.7, units = "in")
par(mar=c(2, 2, 1, 1), mgp=c(1, 0.3, 0), las=0)
plot(iqrs_whole[cross_stoch], -log10(pvals), col = alpha("black", 0.3), pch =19, 
     xlim = c(0, 0.155), cex = 0.5, axes = F, cex.lab = 0.5)
axis(1, seq(0, 0.15, 0.05), tck=-0.02, cex.axis = 0.5)
axis(2, seq(0, 50, 10), tck=-0.02, cex.axis = 0.5, las = 2)
abline(h = -log10(0.05/length(cross_stoch)), lty = 2, col = "red2")
abline(v = 0.07, lty = 2, col = "red2")
dev.off()

# Compute counts per sector
sum(iqrs_whole[cross_stoch] < 0.07) # 154
sum(iqrs_whole[cross_stoch] > 0.07) # 178
sum(iqrs_whole[cross_stoch] > 0.07 & pvals > 0.05/332) # 24
sum(iqrs_whole[cross_stoch] > 0.07 & pvals <= 0.05/332) # 154

setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
write.table(x = sig2, file = "replicated.txt",row.names = F, col.names = F, quote = F)

###################################  Compare with E-risk  ###################################

# Compare behaviour of E-risk markers in Adipose
setwd('/media/ben/DATA/Ben/1_evCpGs/discovery/')
pvals1 <- readRDS('y_rTOST_SQN_ComBat_cellcomp_pvalues.rds')
tested0 = names(pvals1)
tested = tested0[!(tested0 %in% stochCpG)]
length(tested0) # 4652
length(tested) # 4319
cross_tested = tested[tested %in% rownames(data)]
length(cross_tested) # 4302
setwd('/media/ben/DATA/Ben/1_evCpGs/data/mQTL/')
control <- as.vector(read.table(file = 'control.txt')$V1)
length(control) # 998
cross_control = control[control %in% rownames(data)][1:length(cross_stoch)]
length(cross_control) # 332
plot(ecdf(delta[cross_stoch,]))
lines(ecdf(delta[cross_control,]), col = "red")
lines(ecdf(delta[cross_tested,]), col = "blue")
#

###################################  Prepare .bed files  ###################################

# Prepare bed file for visualization
data = arrange.beta(data, matching)
delta = abs(data[,seq(1,ncol(data),2)] - data[,seq(2,ncol(data),2)])
delta = delta[!(rownames(delta) %in% probes2remove2),]
dim(delta) #  343844     97
means = rowMedians(delta)
names(means) = rownames(delta)
variable = names(which(means > 0.04))
length(variable) #  12227
head(annotation)
jkl = annotation[variable, ]
bed = data.frame(chr = jkl$chr, start = as.integer(jkl$pos), end = as.integer(jkl$pos)+1L, ID = rownames(jkl), value = 1000)
setwd("/media/ben/DATA/Ben/1_evCpGs/tracks/450K/")
write.table(x = "track name=pairedReads description=Clone Paired Reads useScore=1", file = "MuTHER_deltabeta.bed", quote = F, row.names = F, col.names = F)
fwrite(bed, "MuTHER_deltabeta.bed", quote = F, row.names = F, col.names = F, sep = "\t", append = T)

# Read data
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
sig2 <- as.vector(read.table(file = 'replicated.txt', header = F)$V1); length(sig2) # 154
jkl = annotation[sig2, ]
bed = data.frame(chr = jkl$chr, start = as.integer(jkl$pos), end = as.integer(jkl$pos)+1L, ID = rownames(jkl), value = 1000)
setwd("/media/ben/DATA/Ben/1_evCpGs/tracks/")
write.table(x = "track name=pairedReads description=Clone Paired Reads useScore=1", file = "replicated.bed", quote = F, row.names = F, col.names = F)
fwrite(bed, "replicated.bed", quote = F, row.names = F, col.names = F, sep = "\t", append = T)

###################################  Annotation  ###################################


# 
setwd("/media/ben/DATA/Ben/1_evCpGs/data/erisk/")
erisk_bad_qc = as.vector(read.table("bad_cpgs.txt", header = F)$V1)
bg0 = rownames(annotation)[!(rownames(annotation) %in% c(X_probes, Y_probes, CR_probes, SNP_probes, erisk_bad_qc))]
length(bg0) # 346555
bg02 = bg0[bg0 %in% rownames(data)]
length(bg02) # 342868


## Protocadherin enrichment
bg = bg02[!(bg0 %in% sig2)]
bg_genes <- annotation[bg,]$UCSC_RefGene_Name
length(bg_genes) # 343717
target_genes <- annotation[sig2,]$UCSC_RefGene_Name
length(target_genes) # 154
a <- sum(startsWith(target_genes, prefix = 'PCDHA')) + sum(startsWith(target_genes, prefix = 'PCDHB')) + sum(startsWith(target_genes, prefix = 'PCDHG'))
b <- sum(startsWith(bg_genes, prefix = 'PCDHA')) + sum(startsWith(bg_genes, prefix = 'PCDHB')) + sum(startsWith(bg_genes, prefix = 'PCDHG'))
c <- length(target_genes) - a
d <- length(bg_genes) - b
contrast <- matrix(c(a,c,b,d), nrow = 2)
colnames(contrast) <- c('target', 'bg')
rownames(contrast) <- c('PCDH_ABG', 'Not PCDH_ABG')
contrast
#              target     bg
# PCDH_ABG         12    526
# Not PCDH_ABG    142 343191

fisher.test(contrast)$p.value # 3.723307e-17
# data:  contrast
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   27.63680 99.85289
# sample estimates:
#   odds ratio 
# 55.14952 


##  GO term enrichment
GOterms <- gometh(sig.cpg = sig2, all.cpg = bg02, collection = 'GO', array.type = '450K', 
                  plot.bias = T, prior.prob = T)
sum(GOterms$FDR < 0.05)
GOterms <- GOterms[order(GOterms$FDR, decreasing = F),]
setwd("/media/ben/DATA/Ben/1_evCpGs/data/GO/")
ev = as.data.frame(fread("GO_ev.txt", header = T, sep = ",", fill = T)[1:19])
replic = fread("GO_replicated.txt", header = T, sep = ",", fill = T)
replic = as.data.frame(replic[replic$Term %in% ev$Term,])
GOmat = cbind(ev[, c(1,2,6,7)], replic[match(replic$Term, ev$Term), c(6,7)])
ev_fdr = -log10(GOmat[,4])
rep_fdr = -log10(GOmat[,6])
#
setwd("~/Desktop/")
tiff(filename = "GO.tiff", width = 4, height = 8, units = "in", res = 300)
par(mar=c(2, 8, 0, 0.5), mgp=c(3, 1, 0), las=0)
barplot(rbind(rep_fdr[length(rep_fdr):1], ev_fdr[length(ev_fdr):1]), horiz = T, names.arg = rev(GOmat$Ont), las = 1,
        col = c('darkgoldenrod2', "red1"), border = NA, beside = T, cex.names = 0.4, xaxt = "n")
axis(1, seq(0, 30, 10))
abline(v = -log10(0.05), lty = 2, col = "red")
dev.off()
#

# KEGG enrichment
KEGGterms <- gometh(sig.cpg = sig, all.cpg = bg0, collection = 'KEGG', array.type = '450K', plot.bias = T, prior.prob = T)
sum(KEGGterms$FDR < 0.1) # 0
head(KEGGterms[order(KEGGterms$P.DE, decreasing = F),], 15)

