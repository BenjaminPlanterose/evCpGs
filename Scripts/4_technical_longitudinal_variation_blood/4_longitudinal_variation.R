############################################################################
############################################################################
###########                                                      ###########
###########           Addressing longitudinal variation          ###########
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

library(GEOquery)
library(data.table)
library(minfi)
library(gplots)
library(ICC)
library(ggplot2)

# computes ICC
compute_ICC = function(i, beta, individual, time)
{
  
  df = data.frame(value = beta[i,], individual = individual, time = time)
  icc = ICCest(x = individual, y = value, data = df, CI.type = "Smith", alpha = 0.05)
  icc = c(icc[[1]], icc[[2]], icc[[3]])
  names(icc) = c("estimate", "low95", "high95")
  
  if(i / 100 == i %/% 100)
  {
    print(i)
  }
  
  
  return(icc)
}

# Read data
setwd("/media/ben/DATA/Ben/1_evCpGs/")
stochCpG <- as.vector(read.table(file = 'evCpGs.txt')$V1)
length(stochCpG) # 333

setwd('/media/ben/DATA/Ben/1_evCpGs/discovery/')
pvals1 <- readRDS('y_rTOST_SQN_ComBat_cellcomp_pvalues.rds')
tested0 = names(pvals1)
tested = tested0[!(tested0 %in% stochCpG)]
length(tested0) # 4652
length(tested) # 4319


setwd('/media/ben/DATA/Ben/1_evCpGs/data/mQTL/')
control <- as.vector(read.table(file = 'control.txt')$V1)
length(control) # 998


#####################################################################

# Retrieve data
setwd('/media/ben/DATA/Ben/1_evCpGs/data/long/')
phenotype <- getGEO('GSE51388', destdir=".")
pheno <- phenotype[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)
#pheno$data_processing[1]
pheno <- pheno[, c(1, 36:39)]

title = as.vector(pheno$title)
group <- unlist(lapply(strsplit(title, split = '_'), function(x) x[1]))
eset <- na.omit(as.matrix(phenotype$GSE51388_series_matrix.txt.gz))
dim(eset) # 362822     60
sum(startsWith(rownames(eset), "rs")) # 0

pheno_A = pheno[group == "Group A",]
pheno_B = pheno[group == "Group B",]
eset_A = eset[, rownames(pheno_A)]
eset_B = eset[, rownames(pheno_B)]
dim(eset_A) # 362822     24
dim(eset_B) # 362822     36

title = as.vector(pheno_A$title)
individuals <- unlist(lapply(strsplit(title, split = 'Group A_MZ pair_ '), function(x) x[2]))
colnames(eset_A) = individuals

title = as.vector(pheno_B$title)
individuals <- unlist(lapply(strsplit(title, split = 'Group B_'), function(x) x[2]))
individuals[1:8] <- unlist(lapply(strsplit(individuals[1:8], split = '#11_'), function(x) x[2]))
colnames(eset_B) = individuals


#
cross_stoch = stochCpG[which(stochCpG %in% rownames(eset))]
length(cross_stoch) # 296
cross_control = control[which(control %in% rownames(eset))]
length(cross_control) # 869
cross_tested = tested[tested %in% rownames(eset)]
length(cross_tested) # 3651

#
L_A_ev <- eset_A[cross_stoch,]
L_A_c <- eset_A[cross_control[1:length(cross_stoch)],]
L_B_ev <- eset_B[cross_stoch,]
L_B_c <- eset_B[cross_control[1:length(cross_stoch)],]
L_B_c2 <- eset_B[cross_tested,]




# Control
heatmap.2(L_A_c, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 1,
          breaks=seq(0, 1, 0.04), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)

heatmap.2(L_B_c, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 1,
          breaks=seq(0, 1, 0.04), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)

# 
heatmap.2(L_A_ev, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 1,
          breaks=seq(0, 1, 0.04), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)
heatmap.2(L_B_ev, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 1,
          breaks=seq(0, 1, 0.04), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)


##########################################################################

pair = sapply(strsplit(colnames(L_B_ev), split = "_"), function(x) x[1])
pair = as.factor(pair)
levels(pair) = brewer.pal(n = 8, name = "Dark2")
pair = as.character(pair)


heatmap.2(L_B_c, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 1,
          breaks=seq(0, 1, 0.04), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90,
          ColSideColors = pair)

heatmap.2(L_B_ev, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 1,
          breaks=seq(0, 1, 0.04), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90,
          ColSideColors = pair)



###########################################################################



individual = unlist(lapply(strsplit(individuals, split = "_"), function(x) x[1]))
time = unlist(lapply(strsplit(individuals, split = "_"), function(x) x[2]))
time = unlist(lapply(strsplit(time, split = " m"), function(x) x[1]))
time = as.numeric(time)

ICC_c = sapply(1:nrow(L_B_c), function(x) compute_ICC(x, L_B_c, individual, time)[1])
ICC_c2 = sapply(1:nrow(L_B_c2), function(x) compute_ICC(x, L_B_c2, individual, time)[1])
ICC_ev = sapply(1:nrow(L_B_ev), function(x) compute_ICC(x, L_B_ev, individual, time)[1])


data = data.frame(value = c(ICC_c, ICC_c2, ICC_ev), group = c(rep("c1", times = length(ICC_c)), 
                                                              rep("non-ev", times = length(ICC_c2)),
                                                              rep("ev", times = length(ICC_ev))))
ggplot(data, aes(x=group, y=value, fill=group)) + geom_violin()


setwd("/media/ben/DATA/Ben/1_evCpGs/data/long/")
ICC = fread("GSE61151_Summary_icc_M.txt")
dim(ICC) # 484949      5

head(ICC)
icc = ICC$ICC
names(icc) = ICC$ID_REF

cross_stoch = stochCpG[stochCpG %in% ICC$ID_REF]; length(cross_stoch) # 333
cross_tested = tested[tested %in% ICC$ID_REF]; length(cross_tested) # 4319
cross_control = (control[control %in% ICC$ID_REF])[1:length(cross_stoch)]; length(cross_control) # 333


data = data.frame(value = c(icc[cross_control], icc[cross_tested], icc[cross_stoch]), group = c(rep("c1", times = length(cross_control)),
                                                                                 rep("non-ev", times = length(cross_tested)),
                                                                                 rep("ev", times = length(cross_stoch))))
ggplot(data, aes(x=group, y=value, fill=group)) + geom_violin()





