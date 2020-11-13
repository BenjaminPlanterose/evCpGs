############################################################################
############################################################################
###########                                                      ###########
###########          Visualizing evCpGs epigenome-wide           ###########
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
library(hexbin)
library(matrixStats)
library(lattice)

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

# Perform hexbin plot for concordance/range
concordance_range_plot <- function(hvCpG, beta, main, lab_range)
{
  MAE_max <- 1/6
  main <- paste(main, paste('(n = ', ncol(beta)/2, ' twin pairs)', sep =''))
  cross_hvCpGs <- hvCpG[hvCpG %in% rownames(beta)]; length(cross_hvCpGs) # 100
  MAE_vec_complete <- lapply(1:nrow(beta), function(x) 1/2*mean(abs(beta[x, seq(1, ncol(beta), 2)] - 
                                                                      beta[x, seq(2, ncol(beta), 2)])))
  MAE_vec_complete <- unlist(MAE_vec_complete)
  names(MAE_vec_complete) <- rownames(beta)
  MAE_vec <- MAE_vec_complete[cross_hvCpGs]
  range <- rowRanges((beta[, seq(1, ncol(beta), 2)] + beta[, seq(2, ncol(beta), 2)])/2)
  range <- range[,2] - range[,1]
  names(range) <- rownames(beta)
  
  cex = 1
  hexbinplot(1 - MAE_vec_complete/MAE_max ~ range, aspect = 1, xbin = 150, 
             cex.title = cex, # legend title size
             cex.labels = cex/2, # legend label size
             xlab = list(label = 'Twin-averaged CpG Methylation Range', cex = cex), 
             ylab = list(label = 'Concordance (1 - MAE/MAEmax)', cex = cex), 
             main = list(label = main, cex = cex),
             scales=list(cex=cex, tck=c(0.2,0), tick.number = 4),
             panel = function(x, y, ...){
               panel.hexbinplot(x,y,...)
               panel.text(lab_range[1], lab_range[2], 'stochCpGs', col ='dodgerblue3', cex = cex)
               panel.points(x = range[cross_hvCpGs], y = 1 - MAE_vec/MAE_max, pch=1, cex=cex/5)})
}


#############################################   Visualize stochCpGs epigenome wide   #############################################

# Read data
# setwd("where")
stochCpG <- as.vector(read.table(file = 'stochCpG.txt')$V1)
sig = stochCpG

# E-risk
# setwd("where")
beta1 <- fread('2019-08-21_SQN_combat_cellcomp.txt', nThread = 4)
beta1 <- process.beta.fread(beta1); dim(beta1) # 346555    852
# setwd("where")
matching_Erisk <- read.table(file = 'matching_MZ.txt', header = T)
beta1 <- arrange.beta(beta1, matching_Erisk)
# setwd("where")
tiff(filename = paste('hexbin_twinsuk', 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
concordance_range_plot(sig, beta1, 'E-risk', c(0.23,0.93))
dev.off()

# TwinsUK
# setwd("where")
TwinsUK <- fread('2019-09-09_sqn_combat_cellcomp.txt', nThread = 4)
TwinsUK <- process.beta.fread(TwinsUK)
dim(TwinsUK) # 346705    656
matching_UK <- read.table('twin_matching.txt', header = F)
TwinsUK <- arrange.beta(TwinsUK, matching_UK)
dim(TwinsUK) # 346705    656
# setwd("where")
tiff(filename = paste('hexbin_twinsuk', 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
concordance_range_plot(sig, TwinsUK, 'TwinsUK', c(0.4,0.95))
dev.off()

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
dim(danish) # 345757    292

# setwd("where")
tiff(filename = paste('hexbin_danish', 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
concordance_range_plot(sig, danish, 'Danish Twin Registry', c(0.3,0.95))
dev.off()
