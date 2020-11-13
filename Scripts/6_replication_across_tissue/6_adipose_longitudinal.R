############################################################################
############################################################################
###########                                                      ###########
###########                Adipose longitudinal                  ###########
###########             Author: Benjamin Planterose              ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########             b.planterose@erasmusmc.nl                ###########
###########                                                      ###########
############################################################################
############################################################################
# 

## Load libraries ##

library(minfi)
library(GEOquery)
library(ENmix)
library(data.table)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(gplots)
library(scales)
library(data.table)
library(ICC)
library(ggplot2)


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

# Computes ICC
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

# Read longitudinal adipose tissue data
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/")
nocombat.beta.sqn = fread("2020-04-06_SQN_nocombat_nocellcomp.txt")
nocombat.beta.sqn = process.beta.fread(nocombat.beta.sqn)
dim(nocombat.beta.sqn) # 346987     57

# Obtain phenotypes
setwd('/media/ben/DATA/Ben/1_evCpGs/data/adipose_long/')
phenotype <- getGEO('GSE103768', destdir=".")
pheno <- phenotype[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)
pheno = pheno[, c(1, 35:38)]
title = as.vector(pheno$title)
individuals = unlist(lapply(strsplit(title, split = "_"), function(x) x[1]))
time = unlist(lapply(strsplit(title, split = "_"), function(x) x[2]))
table(individuals) # 3 of each
table(time)

## Quick QC
densityPlot(nocombat.beta.sqn)
mdsPlot(nocombat.beta.sqn, sampGroups = individuals, pch = 19)
##

# Obtain replicated evCpGs in adipose tissue
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
stochCpG <- as.vector(read.table(file = 'replicated.txt', header = F)$V1); length(stochCpG) # 154
cross_stoch = stochCpG[stochCpG %in% rownames(nocombat.beta.sqn)]; length(cross_stoch) # 154

setwd('/media/ben/DATA/Ben/1_evCpGs/discovery/')
pvals1 <- readRDS('y_rTOST_SQN_ComBat_cellcomp_pvalues.rds')
tested0 = names(pvals1); length(tested0) # 4652
tested = tested0[!(tested0 %in% stochCpG)]; length(tested) # 4319
cross_tested = tested[tested %in% rownames(nocombat.beta.sqn)]; length(cross_tested) # 4306

setwd('/media/ben/DATA/Ben/1_evCpGs/data/mQTL/')
control <- as.vector(read.table(file = 'control.txt')$V1); length(control) # 998
cross_control = control[control %in% rownames(nocombat.beta.sqn)][1:length(cross_stoch)]; length(cross_control) # 154

X = nocombat.beta.sqn[cross_stoch,]
Y = nocombat.beta.sqn[cross_control,]
Z = nocombat.beta.sqn[cross_tested,]
ID = unlist(lapply(strsplit(title, split = "ers"), function(x) x[2]))
colnames(X) = individuals
colnames(Y) = individuals
colnames(Z) = individuals

heatmap.2(Y, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 0.5, scale = 'none',
          breaks=seq(0, 1, 0.05), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)
heatmap.2(X, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 0.5, scale = 'none',
          breaks=seq(0, 1, 0.05), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)

# Exclude individual 3
X2 = X[, !(individuals %in% "WeightRegainers3")]
Y2 = Y[, !(individuals %in% "WeightRegainers3")]
Z2 = Z[, !(individuals %in% "WeightRegainers3")]
individuals2 = individuals[!(individuals %in% "WeightRegainers3")]
time2 = time[!(individuals %in% "WeightRegainers3")]


heatmap.2(Y2, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 0.5, scale = 'none',
          breaks=seq(0, 1, 0.05), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)
heatmap.2(X2, tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 0.5, scale = 'none',
          breaks=seq(0, 1, 0.05), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)


ICC_c = sapply(1:nrow(Y2), function(x) compute_ICC(x, Y2, individuals2, time2)[1])
ICC_c2 = sapply(1:nrow(Z2), function(x) compute_ICC(x, Z2, individuals2, time2)[1])
ICC_ev = sapply(1:nrow(X2), function(x) compute_ICC(x, X2, individuals2, time2)[1])

data = data.frame(value = c(ICC_c, ICC_c2, ICC_ev), group = c(rep("c1", times = length(ICC_c)), 
                                                              rep("non-ev", times = length(ICC_c2)),
                                                              rep("ev", times = length(ICC_ev))))

ggplot(data, aes(x=group, y=value, fill=group)) + geom_violin()

