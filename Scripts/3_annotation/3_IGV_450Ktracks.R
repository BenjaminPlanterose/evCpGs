############################################################################
############################################################################
###########                                                      ###########
###########              Prepare 450K IGV tracks                 ###########
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
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

## Load functions ##

# Annotate CpGs based on Illumina annotation (hg19)
functional.annotation.hg19 <- function(hvCpGs)
{
  old.dir <- getwd()
  # setwd("where")
  RGSET <- read.metharray.exp(getwd())
  
  annotation <- getAnnotation(RGSET)
  #annotation2 <- PAMES::illumina450k_hg38
  dim(annotation)
  #names <- annotation2$`Composite Element REF`
  names <- rownames(annotation)
  
  manifest <- getManifest(RGSET)
  
  
  # Chr
  chr <- annotation$chr[match(hvCpGs, names)]
  pos <- annotation$pos[match(hvCpGs, names)]
  
  # strand
  strand <- annotation$strand[match(hvCpGs, as.vector(annotation$Name))]
  
  # Type
  I <- getProbeInfo(manifest, type = c('I'))
  I_annot <- as.vector(I$nCpG)
  names(I_annot) <- I$Name
  
  II <- getProbeInfo(manifest, type = c('II'))
  II_annot <- as.vector(II$nCpG)
  names(II_annot) <- II$Name
  
  # Types
  types <- character(length = length(hvCpGs))
  types[na.omit(match(names(I_annot), hvCpGs ))] <- 'I'
  types[na.omit(match(names(II_annot), hvCpGs ))] <- 'II'
  
  # Number of CpGs assessed by probes
  n.CpGs <- numeric(length = length(hvCpGs))
  n.CpGs[na.omit(match(names(I_annot), hvCpGs ))] <- as.vector(na.omit(I_annot[hvCpGs]))
  n.CpGs[na.omit(match(names(II_annot), hvCpGs ))] <- as.vector(na.omit(II_annot[hvCpGs]))
  
  
  # Island status
  island <- annotation$Relation_to_Island[match(hvCpGs, as.vector(annotation$Name))]
  
  
  # Feature group
  feature <- annotation$Regulatory_Feature_Group[match(hvCpGs, as.vector(annotation$Name))]
  feature[which(feature == '')] <- NA
  
  # Seqs
  #seqs <- annotation$Forward_Sequence[match(hvCpGs, as.vector(annotation$Name))]
  
  
  # Enhancer
  enhancer <- annotation$Enhancer[match(hvCpGs, as.vector(annotation$Name))]
  enhancer[which(enhancer == '')] <- 'FALSE'
  enhancer <- as.logical(enhancer)
  
  # Annot
  UCSC_name <- annotation$UCSC_RefGene_Name[match(hvCpGs, as.vector(annotation$Name))]
  #UCSC_accession <- annotation$UCSC_RefGene_Accession[match(hvCpGs, as.vector(annotation$Name))]
  UCSC_group <- annotation$UCSC_RefGene_Group[match(hvCpGs, as.vector(annotation$Name))]
  
  # DMR
  dmr <- annotation$DMR[match(hvCpGs, as.vector(annotation$Name))]
  
  # Probe MAF
  MAF <- annotation$Probe_maf[match(hvCpGs, as.vector(annotation$Name))]
  
  
  # DNAse-I hypersensitivity
  DHS <- annotation$DHS[match(hvCpGs, as.vector(annotation$Name))]
  DHS <- as.logical(DHS)
  DHS[which(DHS == '')] <- NA
  
  # In 27 K
  #annotation$Methyl27_Loci[match(hvCpGs, as.vector(annotation$Name))]
  
  #Gene_type <- annotation2$Gene_Type[match(hvCpGs, names)]
  
  
  
  final.annot <- data.frame(chr = chr, pos= pos, strand = strand, type = types, nCpGs = n.CpGs,
                            island.status = island, feature.group = feature,
                            enhancer = enhancer, UCSC.name = UCSC_name, 
                            UCSC.group = UCSC_group, DMR = dmr, MAF = MAF,
                            DHS = DHS)
  #, Gene_type = Gene_type)
  rownames(final.annot) <- hvCpGs
  setwd(old.dir)
  return(final.annot)                    
  
  
}



#############################################   Read data   #############################################

# evCpGs
# setwd("where")
stochCpG <- as.vector(read.table(file = 'stochCpG.txt')$V1)
sig = stochCpG


# # E-risk
# setwd("where")
beta1 <- fread('2019-08-21_SQN_combat_cellcomp.txt', nThread = 4)
beta1 <- process.beta.fread(beta1); dim(beta1) # 346555    852
# setwd("where")
matching_Erisk <- read.table(file = 'matching_MZ.txt', header = T)
beta1 <- arrange.beta(beta1, matching_Erisk)

# TwinsUK
# setwd("where")
TwinsUK <- fread('2019-09-09_sqn_combat_cellcomp.txt', nThread = 4)
TwinsUK <- process.beta.fread(TwinsUK)
dim(TwinsUK) # 346705    656
matching_UK <- read.table('twin_matching.txt', header = F)
TwinsUK <- arrange.beta(TwinsUK, matching_UK)

# Danish Cohort
# setwd("where")
danish <- fread('2019-09-09_sqn_combat_cellcomp.txt')
danish <- process.beta.fread(danish)
replicates_excluded <- c('GSM1506278', 'GSM1506587',
                         'GSM1506580', 'GSM1506342',
                         'GSM1506333', 'GSM1506438', 'GSM1506500',
                         'GSM1506327', 'GSM1506494', 'GSM1506554')
rgSamples <- unlist(lapply(strsplit(colnames(danish), split = '_'), function(x) x[1]))
danish <- danish[,!(rgSamples %in% replicates_excluded)]
matching_danish <- read.table('matching_danish.txt', header = T)
danish <- arrange.beta(danish, matching_danish)


#############################################   evCpG track   #############################################

annotation <- functional.annotation.hg19(sig)
# StochCpGs
igv_stochCpG <- data.frame(seqname = annotation$chr,
                           start = annotation$pos, end = annotation$pos+1,
                           name = rownames(annotation),
                           score = 1000)
# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1', file = 'stochCpG_IGV.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = igv_stochCpG, file = 'stochCpG_IGV.bed', nThread = 4, sep = '\t', col.names = F, append = T)


#############################################   High IQR in E-risk track   #############################################

# setwd("where")
pvals1 <- readRDS('y_rTOST_SQN_ComBat_cellcomp_pvalues.rds')
high_IQRs <- names(pvals1); length(high_IQRs) # 4652
annotation <- functional.annotation.hg19(high_IQRs)
igv_high_IQRs <- data.frame(seqname = annotation$chr,
                            start = annotation$pos, end = annotation$pos+1,
                            name = rownames(annotation),
                            score = 1000)

# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1', file = 'highIQR_Erisk_IGV.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = igv_high_IQRs, file = 'highIQR_Erisk_IGV.bed', nThread = 4, sep = '\t', col.names = F, append = T)


#############################################   Whole 450K track   #############################################

# Whole 450K
old.dir <- getwd()
# setwd("where")
RGSET <- read.metharray.exp(getwd())
annotation <- getAnnotation(RGSET)
setwd(old.dir)
dim(annotation) # 485512    33
igv_450K <- data.frame(seqname = annotation$chr,
                       start = annotation$pos, end = annotation$pos+1,
                       name = rownames(annotation),
                       score = 1000)

# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1', file = '450K.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = igv_450K, file = '450K.bed', nThread = 4, sep = '\t', col.names = F, append = T)

#############################################   Whole 450K_binned coverage track (only chr5)   #############################################

# Whole 450K binned by coverage (chr5)
annotation <- annotation[annotation$chr == 'chr5',]
bins <- seq(from = 1, to = 180915260, by = 1000) # http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
counts <- binCounts(x = annotation$pos, bx = bins)
coverage_track <- data.frame(seqname = 'chr5',
                             start = as.integer(bins[-length(bins)]), end = as.integer(bins[2:length(bins)] - 1),
                             name = paste('chr5', paste(bins[-length(bins)], bins[2:length(bins)] - 1, sep = '-'), sep = ':'),
                             score = 1000*counts/max(counts))
dim(coverage_track)
# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1 graphType=line viewLimits=0:1000', file = 'coverage_bg_450K_chr5.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = coverage_track, file = 'coverage_bg_450K_chr5.bed', nThread = 4, sep = '\t', col.names = F, append = T)

#############################################   Whole 450K track untrustable probes   #############################################

# Whole 450K without cross-reactives
# Probes to remove
# setwd("where")
# Read CpGs to remove
# XY-probes
Y_probes <- unique(as.vector(read.table(file = "y_chromosome_probes.txt", header = F)$V1))
length(Y_probes) # 416
X_probes <- unique(as.vector(read.table(file = "x_chromosome_probes.txt", header = F)$V1))
length(X_probes) # 11232

# SNP probes - dbSNP v.147
data(SNPs.147CommonSingle)
f.SNP <- c(rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$Probe_maf >= 0.01],
           rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$CpG_maf > 0],
           rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$SBE_maf > 0])
SNP_probes <- na.omit(unique(f.SNP))
length(SNP_probes) # 99337

# CR probes
CR_1 <- as.vector(read.table('crossreactive_Chen.txt', header = T)$TargetID) # Chen YA et al. Epigenetics. 2013 Feb;8(2):203-9. doi: 10.4161/epi.23470. Epub 2013 Jan 11.
kobor <- fread('GPL16304-47833.txt') # Price ME et al. Epigenetics Chromatin. 2013 Mar 3;6(1):4. doi: 10.1186/1756-8935-6-4.
CR_2 <- unique(c(kobor$ID[kobor$Autosomal_Hits == 'A_YES'], kobor$ID[kobor$XY_Hits == 'XY_YES']))
CR_probes <- unique(c(CR_1, CR_2))
length(CR_probes) # 41993
probes2remove <- unique(c(Y_probes, X_probes, SNP_probes, CR_probes))
length(probes2remove) # 137912

#
annotation <- annotation[!(rownames(annotation) %in% probes2remove),]
dim(annotation) # 347600     33

igv_450K_filtered <- data.frame(seqname = annotation$chr,
                                start = annotation$pos, end = annotation$pos+1,
                                name = rownames(annotation),
                                score = 1000)
# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1', file = '450K_filtered.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = igv_450K_filtered, file = '450K_filtered.bed', nThread = 4, sep = '\t', col.names = F, append = T)



#############################################   Delta beta tracks   #############################################

# E-risk
delta.beta <- abs(beta1[, seq(1, ncol(beta1),2)] - beta1[, seq(2, ncol(beta1),2)])
dim(delta.beta) # 346555    426
medians <- rowMedians(delta.beta)
names(medians) <- rownames(delta.beta)
bigMeds <- names(which(medians > 0.04))
annotation <- functional.annotation.hg19(bigMeds)
delta.beta.igv <- data.frame(seqname = annotation$chr,
                             start = annotation$pos, end = annotation$pos+1,
                             name = rownames(annotation),
                             score = 1000)

# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1', file = 'E_risk_deltabeta.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = delta.beta.igv, file = 'E_risk_deltabeta.bed', nThread = 4, sep = '\t', col.names = F, append = T)

# TwinsUK
delta.beta <- abs(TwinsUK[, seq(1, ncol(TwinsUK),2)] - TwinsUK[, seq(2, ncol(TwinsUK),2)])
dim(delta.beta) # 346705    328
medians <- rowMedians(delta.beta)
names(medians) <- rownames(delta.beta)
bigMeds <- names(which(medians > 0.04))
annotation <- functional.annotation.hg19(bigMeds)
delta.beta.igv <- data.frame(seqname = annotation$chr,
                             start = annotation$pos, end = annotation$pos+1,
                             name = rownames(annotation),
                             score = 1000)

# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1', file = 'TwinsUK_deltabeta.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = delta.beta.igv, file = 'TwinsUK_deltabeta.bed', nThread = 4, sep = '\t', col.names = F, append = T)

# Danish
delta.beta <- abs(danish[, seq(1, ncol(danish),2)] - danish[, seq(2, ncol(danish),2)])
dim(delta.beta) # 345757    146
medians <- rowMedians(delta.beta)
names(medians) <- rownames(delta.beta)
bigMeds <- names(which(medians > 0.04))
annotation <- functional.annotation.hg19(bigMeds)
delta.beta.igv <- data.frame(seqname = annotation$chr,
                             start = annotation$pos, end = annotation$pos+1,
                             name = rownames(annotation),
                             score = 1000)

# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1', file = 'Danish_deltabeta.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = delta.beta.igv, file = 'Danish_deltabeta.bed', nThread = 4, sep = '\t', col.names = F, append = T)