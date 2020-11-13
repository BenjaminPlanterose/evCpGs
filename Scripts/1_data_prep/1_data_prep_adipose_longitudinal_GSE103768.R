############################################################################
############################################################################
###########                                                      ###########
###########          Data preparation (GSE103768)                ###########
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
library(minfi)
library(GEOquery)
library(ENmix)
library(data.table)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(gplots)
library(scales)

## Technical QC ##
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/GSE103768_RAW/")
rgSet <- read.metharray.exp(getwd(), extended = T); gc()
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/QC/")
qc <- QCinfo(rgSet)
bad_cpgs <- qc$badCpG
write.table(bad_cpgs, file = 'bad_cpgs.txt', row.names = F, col.names = F, quote = F)

nbthre = 3
detPthre = 1e-06
qcmat <- qc$nbead < nbthre | qc$detP > detPthre
badValuePerSample <- apply(qcmat, 2, sum)/nrow(qcmat)
thr = mean(qc_df$bisul) - 3*sd(qc_df$bisul)
qc_df = data.frame(badFreq = badValuePerSample, bisul = qc$bisul)
dim(qc_df) # 57 2

setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/QC/")
tiff(filename = "longitudinal_adipose.tiff", height = 427, width = 550, units = "px")
plot(qc_df$badFreq, qc_df$bisul, ylim = c(0, 35000), xlim = c(0, 0.055), pch = 19, col = alpha("black", 0.2), cex = 0.4,
     xlab = "Percent of low-quality data", ylab = "Average bisulfite conversion intensity", main = "Adipose Longitudinal", las = 1, cex.axis = 0.8)
abline(v = 0.05, lty= 2, col = alpha("red3", 0.5))
abline(h = thr, lty = 2, col = alpha("red3", 0.5))
dev.off()


## Sex QC ##
sqn <- preprocessQuantile(rgSet)
Sex <- minfi::getSex(sqn)

setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/QC/")
tiff(filename = "sex1_longitudinal_adipose.tiff", height = 762, width = 656, units = "px")
plot(x = Sex$xMed, y = Sex$yMed, type = "n", xlab = "X chr, median total intensity (log2)", 
     ylab = "Y chr, median total intensity (log2)")
text(x = Sex$xMed, y = Sex$yMed, labels = Sex$predictedSex, col = ifelse(Sex$predictedSex ==  "M", "deepskyblue", "deeppink3"))
legend("bottomleft", c("M", "F"), col = c("deepskyblue", "deeppink3"), pch = 16)
dev.off()


sex.mat <- table(Sex$predictedSex, pheno$`gender:ch1`) 
#     Female Male
# F     36    0
# M      0   21

# Perform balloon plot
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/QC/")
tiff(filename = "sex2_longitudinal_adipose.tiff", height = 762, width = 656, units = "px")
balloonplot(as.table(t(sex.mat)), xlab = 'Predicted', ylab = 'Registered', main = '')
dev.off()



## 450K 65 SNP probe QC ##
snps <- getSnpBeta(rgSet)
Names = colnames(snps)
colnames(snps) = unlist(lapply(strsplit(Names, split = "_"), function(x) x[2]))
heatmap.2(snps, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 0.5, scale = 'none',
          breaks=seq(0, 1, 0.05), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)



###############################################################################################


setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/GSE103768_RAW/")
rgSet = read.metharray.exp(getwd()); gc()
samples <- colnames(rgSet)
rgSamples <- unlist(lapply(X = strsplit(samples, split = '_'), function(X) { paste(X[[1]])}))

setwd('/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/')
phenotype <- getGEO('GSE103768', destdir=".")
pheno <- phenotype[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)
pheno = pheno[, c(1, 35:38)]
title = as.vector(pheno$title)
individuals = unlist(lapply(strsplit(title, split = "_"), function(x) x[1]))
table(individuals) # 3 of each


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

setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/QC/")
bad_probes = as.vector(read.table("bad_cpgs.txt", header = F)$V1)
length(bad_probes)

annotation = getAnnotation(rgSet)
Y_probes = rownames(annotation)[annotation$chr == "chrY"]
X_probes = rownames(annotation)[annotation$chr == "chrX"]
annotation = as.data.frame(annotation)
probes2remove <- unique(c(Y_probes, X_probes, SNP_probes, CR_probes, bad_probes))
length(probes2remove) # 138525

rm(SNPs.147CommonSingle, bad_probes, CR_1, CR_2, CR_probes, f.SNP, SNP_probes, X_probes, Y_probes, kobor); gc()


# Exclude probes2remove
um.sqn <- preprocessQuantile(rgSet)
um.sqn <- um.sqn[!(rownames(um.sqn) %in% probes2remove),]
dim(um.sqn) # 346987    894

# Apply offset of alpha = 100 (alpha argument of getBeta does not work)
nocombat.beta.sqn <- minfi::getBeta(um.sqn)
nocombat.CN.sqn <- 2^getCN(um.sqn)
nocombat.beta.sqn <- nocombat.beta.sqn*nocombat.CN.sqn/(nocombat.CN.sqn+100) 
rm(nocombat.CN.sqn, um.sqn) ; gc()


########## Export ###############
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/")
fwrite(data.table(nocombat.beta.sqn, keep.rownames = T), paste(Sys.Date(), 'SQN_nocombat_nocellcomp.txt', sep = '_'), quote = F, 
       row.names = T, col.names = T, sep = '\t', nThread = 4)
########## Export ###############
