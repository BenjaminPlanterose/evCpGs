############################################################################
############################################################################
###########                                                      ###########
###########              Addressing technical noise              ###########
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

library(scales)
library(data.table)
library(minfi)
library(equivalence)
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

# It arranges the beta value matrix based on the twin design. The resulting matrix will have 
# columns such that:  {Twin 11, Twin 12, Twin 21, Twin 22, ...}
arrange.beta <- function(beta, matching)
{
  indices <- numeric(length = 2*nrow(matching))
  indices[seq(1, 2*nrow(matching), 2)] <- match(as.vector(matching$V1), colnames(beta))
  indices[seq(2, 2*nrow(matching), 2)] <- match(as.vector(matching$V2), colnames(beta))
  return(beta[,indices])
}

# It will help extract concordance/CpG range in Twin/Replicate data
concordance_range <- function(hvCpG, beta)
{
  cross_hvCpGs <- hvCpG[hvCpG %in% rownames(beta)]; length(cross_hvCpGs) # 100
  beta <- beta[cross_hvCpGs,]
  MAE_vec_complete <- lapply(1:nrow(beta), function(x) 1/2*mean(abs(beta[x, seq(1, ncol(beta), 2)] - 
                                                                      beta[x, seq(2, ncol(beta), 2)])))
  MAE_vec_complete <- unlist(MAE_vec_complete)
  names(MAE_vec_complete) <- rownames(beta)
  range <- rowRanges((beta[, seq(1, ncol(beta), 2)] + beta[, seq(2, ncol(beta), 2)])/2)
  range <- range[,2] - range[,1]
  names(range) <- rownames(beta)
  return(list(MAE_vec_complete[cross_hvCpGs], range[cross_hvCpGs]))
}

# It performs a vector plot with starting points at replicate data towards twin data
# plus segregating the trends into 4 categories.
plot.arrow.field <- function(c_r_rep, c_r_twin, xlim, ylim)
{
  MAE_max <- 1/6
  x0 = c_r_rep[[2]]
  x1 = c_r_twin[[2]]
  y0 = 1 - c_r_rep[[1]]/MAE_max
  y1 = 1 - c_r_twin[[1]]/MAE_max
  deltax = x1 - x0
  deltay = y1 - y0
  trend1 <- which(deltax >= 0 & deltay > 0)
  trend2 <- which(deltax >= 0 & deltay <= 0)
  trend3 <- which(deltax < 0 & deltay > 0)
  trend4 <- which(deltax < 0 & deltay <= 0)
  trend.list <- list(trend1, trend2, trend3, trend4)
  mains <- c('DeltaX >= 0 & DeltaY > 0', 'DeltaX >= 0 & DeltaY <= 0', 
             'DeltaX < 0 & DeltaY > 0', 'DeltaX < 0 & DeltaY <= 0')
  par(mfrow = c(2,2))
  for(i in 1:length(trend.list))
  {
    plot(x = x0[trend.list[[i]]], y = y0[trend.list[[i]]], ylim = ylim, xlim = xlim,
         main = mains[i], type = 'n', xlab = 'Twin-averaged CpG Methylation Range', 
         ylab = 'Concordance (1 - MAE/MAEmax)')
    
    arrows(x0 = x0[trend.list[[i]]], x1 = x1[trend.list[[i]]],
           y0 = y0[trend.list[[i]]], y1 = y1[trend.list[[i]]], 
           lty = 1, lwd = 0.4, length = 0.1, angle = 10, col = alpha('black', 0.3))
    
    points(x = x0[trend.list[[i]]], y = y0[trend.list[[i]]], col = alpha('deepskyblue3', 0.7), pch = 19, cex = 0.5)
    
    points(x = x1[trend.list[[i]]], y = y1[trend.list[[i]]], col = alpha('firebrick3', 0.7), pch = 19, cex = 0.5)
    
  }
  par(mfrow = c(1,1))
  return(list(names(trend1), names(trend2), names(trend3), names(trend4)))
}



#############################################   Technical noise on E-risk - nBead, detP and ICC   #############################################

# Read data
setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/res/pvals/')
stochCpG <- as.vector(read.table(file = 'stochCpG.txt')$V1)
sig = stochCpG


setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/data/qc/erisk/')
nBead <- fread(input = '2019-09-04_nBead.txt')
nBead <- process.beta.fread(nBead); dim(nBead) # 485512    852

detP <- fread(input = '2019-09-04_detP.txt')
detP <- process.beta.fread(detP); dim(detP) # 485512    852

setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/res/pvals/')
pvals1 <- readRDS('y_rTOST_SQN_ComBat_cellcomp_pvalues.rds')
bg = names(pvals1)

setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/res/pvals/')
stochCpG <- as.vector(read.table(file = 'stochCpG.txt')$V1)
sig = stochCpG
nosig = bg[!(bg %in% sig)]

# nBeads
x <- (nBead[nosig,])
y <- (nBead[sig,])
hist(x, col = alpha('red3', 0.5), breaks = max(nBead[nosig,]), freq = F, main = '', xlab = 'nBead')
hist(y, col = alpha('blue3', 0.5), breaks = max(nBead[sig,]), add = T, freq = F)
(a <- rtost(x, y, alpha = 0.05, epsilon = 1))
a$p.value # 0

# detP
x <- rowSums(detP[nosig,] > 0.0001)/ncol(detP)
y <- rowSums(detP[sig,] > 0.0001)/ncol(detP)
boxplot(list(x,y), ylab = 'Proportion of samples with p-val < 0.0001', names = c('Not Significant','Significant'))
hist(x, col = alpha('red3', 0.5), breaks = max(ncol(detP)*x), freq = F)
hist(y, col = alpha('blue3', 0.5), breaks = max(ncol(detP)*y), add = T, freq = F)
(a <- rtost(x, y, alpha = 0.05, epsilon = 0.001))
a$p.value # 4.232592e-54

# ICC
setwd('/media/ultron/2tb_disk2/PROCESSED_DATA/2018/Twin_project/ICC_replicates/') # Bose M et al. BMC Bioinformatics. 2014; 15(1): 312. 
dataset <- fread('ICC_values.csv')
dataset <- dataset[, c(1,4)]
dataset <- as.data.frame(dataset)
colnames(dataset) <- c('ID', 'ICC')
rownames(dataset) <- dataset$ID

tag <- rep('Excluded', times = nrow(dataset))
tag[rownames(dataset) %in% bg] <- 'Discovery Set - Not Significant'
tag[rownames(dataset) %in% sig] <- 'Discovery Set - Significant'
dataset$Group <- tag
dataset <- dataset[,-1]

# Basic violin plot
data_summary <- function(x) 
{
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/res/')
tiff(filename = paste('ICC', 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
ggplot(dataset, aes(x=Group, y=ICC, fill=Group)) + geom_violin() + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  stat_summary(fun.data=data_summary) + labs(title="ICC Distribution",x="Group", y = "ICC") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#############################################   Technical noise on E-risk - per-stage noise assessment   #############################################

setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/data/beta/E_risk/')
matching <- read.table(file = 'matching_MZ.txt', header = T)

raw <- fread('2019-08-22_raw.txt')
raw <- process.beta.fread(raw)
raw <- arrange.beta(raw, matching)
dim(raw) # 346555    852
densityPlot(raw)
delta.raw <- abs(raw[,seq(1,ncol(raw),2)] - raw[,seq(2,ncol(raw),2)])
densityPlot(delta.raw)
dim(delta.raw) # 346555    426
rm(raw, delta.raw); gc()

# SQN
sqn1 <- fread('2019-08-21_SQN_nocombat_nocellcomp.txt')
sqn1 <- process.beta.fread(sqn1)
sqn1 <- arrange.beta(sqn1, matching)
delta.sqn1 <- abs(sqn1[,seq(1,ncol(sqn1),2)] - sqn1[,seq(2,ncol(sqn1),2)])
densityPlot(sqn1)
densityPlot(delta.sqn1, xlim = c(0,1))
rm(sqn1, delta.sqn1); gc()

sqn2 <- fread('2019-08-21_SQN_combat_nocellcomp.txt')
sqn2 <- process.beta.fread(sqn2)
sqn2 <- arrange.beta(sqn2, matching)
delta.sqn2 <- abs(sqn2[,seq(1,ncol(sqn2),2)] - sqn2[,seq(2,ncol(sqn2),2)])
densityPlot(sqn2)
densityPlot(delta.sqn2, xlim = c(0,1))
rm(sqn2, delta.sqn2); gc()

sqn3 <- fread('2019-08-21_SQN_combat_cellcomp.txt')
sqn3 <- process.beta.fread(sqn3)
sqn3 <- arrange.beta(sqn3, matching)
delta.sqn3 <- abs(sqn3[,seq(1,ncol(sqn3),2)] - sqn3[,seq(2,ncol(sqn3),2)])
densityPlot(sqn3)
densityPlot(delta.sqn3, xlim = c(0,1))
rm(sqn3, delta.sqn3); gc()

# dasen
dasen1 <- fread('2019-08-21_dasen_nocombat_nocellcomp.txt')
dasen1 <- process.beta.fread(dasen1)
dasen1 <- arrange.beta(dasen1, matching)
delta.dasen1 <- abs(dasen1[,seq(1,ncol(dasen1),2)] - dasen1[,seq(2,ncol(dasen1),2)])
densityPlot(dasen1)
densityPlot(delta.dasen1, xlim = c(0,1))
rm(dasen1, delta.dasen1); gc()

dasen2 <- fread('2019-08-21_dasen_combat_nocellcomp.txt')
dasen2 <- process.beta.fread(dasen2)
dasen2 <- arrange.beta(dasen2, matching)
delta.dasen2 <- abs(dasen2[,seq(1,ncol(dasen2),2)] - dasen2[,seq(2,ncol(dasen2),2)])
densityPlot(dasen2)
densityPlot(delta.dasen2, xlim = c(0,1))
rm(dasen2, delta.dasen2); gc()

dasen3 <- fread('2019-08-21_dasen_combat_cellcomp.txt')
dasen3 <- process.beta.fread(dasen3)
dasen3 <- arrange.beta(dasen3, matching)
delta.dasen3 <- abs(dasen3[,seq(1,ncol(dasen3),2)] - dasen3[,seq(2,ncol(dasen3),2)])
densityPlot(dasen3)
densityPlot(delta.dasen3, xlim = c(0,1))
rm(dasen3, delta.dasen3); gc()


# oob_Relic_qn_bmiq
bmiq1 <- fread('2019-08-22_oobrelicbmiqqn_nocombat_nocellcomp.txt')
bmiq1 <- process.beta.fread(bmiq1)
bmiq1 <- arrange.beta(bmiq1, matching)
delta.bmiq1 <- abs(bmiq1[,seq(1,ncol(bmiq1),2)] - bmiq1[,seq(2,ncol(bmiq1),2)])
densityPlot(bmiq1)
densityPlot(delta.bmiq1, xlim = c(0,1))
rm(bmiq1, delta.bmiq1); gc()

bmiq2 <- fread('2019-08-22_oobrelicbmiqqn_combat_nocellcomp.txt')
bmiq2 <- process.beta.fread(bmiq2)
bmiq2 <- arrange.beta(bmiq2, matching)
delta.bmiq2 <- abs(bmiq2[,seq(1,ncol(bmiq2),2)] - bmiq2[,seq(2,ncol(bmiq2),2)])
densityPlot(bmiq2)
densityPlot(delta.bmiq2, xlim = c(0,1))
rm(bmiq2, delta.bmiq2); gc()

bmiq3 <- fread('2019-08-22_oobrelicbmiqqn_combat_cellcomp.txt')
bmiq3 <- process.beta.fread(bmiq3)
bmiq3 <- arrange.beta(bmiq3, matching)
delta.bmiq3 <- abs(bmiq3[,seq(1,ncol(bmiq3),2)] - bmiq3[,seq(2,ncol(bmiq3),2)])
densityPlot(bmiq3)
densityPlot(delta.bmiq3, xlim = c(0,1))
rm(bmiq3, delta.bmiq3); gc()


#############################################   Technical noise on Danish Cohort   #############################################


# Danish
setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/data/beta/Danish_cohort/')
danish <- fread('2019-09-09_sqn_combat_cellcomp.txt')
danish <- process.beta.fread(danish)
dim(danish) # 345757    302

matching_danish <- read.table('matching_danish.txt', header = T)
matching_replicates <- data.frame(V1 = c('GSM1506278', 'GSM1506284', 'GSM1506336', 'GSM1506342', 'GSM1506333', 'GSM1506333', 'GSM1506333', 'GSM1506438', 'GSM1506438', 'GSM1506500', 'GSM1506327', 'GSM1506327', 'GSM1506327','GSM1506432', 'GSM1506432', 'GSM1506494'),
                                  V2 = c('GSM1506581', 'GSM1506587', 'GSM1506580', 'GSM1506586', 'GSM1506438', 'GSM1506500', 'GSM1506560', 'GSM1506500', 'GSM1506560', 'GSM1506560', 'GSM1506432', 'GSM1506494', 'GSM1506554', 'GSM1506494', 'GSM1506554', 'GSM1506554'))

rgSamples <- unlist(lapply(strsplit(colnames(danish), split = '_'), function(x) x[1]))
temp <- danish
colnames(temp) <- rgSamples


Replicates <- arrange.beta(temp, matching_replicates); rm(temp)
Twins <- arrange.beta(danish, matching_danish)
dim(Twins) # 345757    292
dim(Replicates) # 345757     32
cross_sig <- sig[sig %in% rownames(Twins)]
length(cross_sig) # 329

# Build delta beta replicate and delta beta twin
delta.rep <- abs(Replicates[cross_sig,seq(1,ncol(Replicates),2)]-
                   Replicates[cross_sig,seq(2,ncol(Replicates),2)])

delta.twin <- abs(Twins[cross_sig,seq(1,ncol(Twins),2)]-
                    Twins[cross_sig,seq(2,ncol(Twins),2)])

dim(delta.rep) # 329  16
dim(delta.twin) #329 146

# Kolmogorov-Smirnov test + ecdf plot
setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/res/')
tiff(filename = paste('technical_KS', 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
pval_KS <- ks.test(as.numeric(delta.rep), as.numeric(delta.twin), alternative = 'greater')$p.value
plot(ecdf(delta.twin), col = 'firebrick3', main = paste('-log10(p-val_KS) =', round(-log10(pval_KS), 4), '\n',  length(cross_sig), 'out of', length(sig), 'stochCpGs', sep = ' '),
     xlab = expression(paste("|", Delta, beta, '|')), ylab = 'Empirical Cumulative Distribution Function',
     lwd = 2, cex = 2)
lines(ecdf(delta.rep), lwd = 2, col = 'deepskyblue3')
legend(x = 'right', col = c('deepskyblue3', 'firebrick3'), legend = c('Technical Replicates', 'MZ twins'), cex = 2, bty = 'n', pch = 19)
dev.off()




# Vector plot
IQR_MAE_rep <- concordance_range(sig, Replicates)
IQR_MAE_twin <- concordance_range(sig, Twins)

setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/res/')
tiff(filename = paste('technical_vector', 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
res <- plot.arrow.field(IQR_MAE_rep, IQR_MAE_twin, c(0,0.7), c(0.4,1))
sapply(res, length) # 50 259   4  16
barplot(sapply(res, length), names.arg = c('A', 'B', 'C', 'D'), ylab = c('Number of CpGs'),
        col = alpha('lightblue3', 0.7))
dev.off()



########################################## Technical noise 2 ###########################

## 
process.beta.fread <- function(beta)
{
  cpgs <- as.vector(beta$rn)
  beta <- beta[,-(1:2)]
  beta <- as.matrix(beta)
  rownames(beta) <- cpgs
  return(beta)
}

# Danish

setwd('/media/ben/DATA/Ben/1_evCpGs/discovery/')
stochCpG <- as.vector(read.table(file = 'stochCpG.txt')$V1)
sig = stochCpG

setwd('/media/ben/DATA/Ben/1_evCpGs/data/mQTL/')
# mQTL <- fread('F7.ALL.M.tab'); dim(mQTL) # 9902081       6
# mQTL <- mQTL[order(mQTL$FDR),]
# fdr_order = unique(mQTL$gene)
# where = na.omit(match(fdr_order, tested))[1:1000]
# control = tested[where]
# control = control[!(control %in% stochCpG)]
# intersect(control, stochCpG) # character(0)
# length(control) # 998
control <- as.vector(read.table(file = 'control.txt')$V1)
length(control) # 998


setwd('/media/ben/DATA/Ben/1_evCpGs/data/beta/Danish_cohort/')
danish <- fread('2019-09-09_sqn_combat_cellcomp.txt')
danish <- process.beta.fread(danish)
dim(danish) # 345757    302
matching_danish <- read.table('matching_danish.txt', header = T)
dim(matching_danish) # 146  2
IDs <- unlist(lapply(X = strsplit(as.vector(colnames(danish)), split = '_'), function(x) x[1]))
colnames(danish) <- IDs

case1 <- danish[, c('GSM1506333', 'GSM1506438', 'GSM1506500', 'GSM1506560',
                    'GSM1506327', 'GSM1506432', 'GSM1506494', 'GSM1506554')]
case2 <- danish[, c('GSM1506342', 'GSM1506586', 
                    'GSM1506336', 'GSM1506580')]
case3 <- danish[, c('GSM1506284', 'GSM1506587', 
                    'GSM1506278', 'GSM1506581')]
colnames(case1) = c(paste('Twin1A', 'rep', 1:4, sep = '_'), paste('Twin1B', 'rep',  1:4, sep = '_'))
colnames(case2) = c(paste('Twin2A', 'rep', 1:2, sep = '_'), paste('Twin2B', 'rep',  1:2, sep = '_'))
colnames(case3) = c(paste('Twin3A', 'rep', 1:2, sep = '_'), paste('Twin3B', 'rep',  1:2, sep = '_'))

cross_stoch = stochCpG[stochCpG %in% rownames(danish)]
length(cross_stoch) # 329
cross_control = control[control %in% rownames(danish)]
length(cross_control) # 992
cross_control = cross_control[1:length(cross_stoch)]
length(cross_control) # 329
# 


heatmap.2(Reduce(cbind, list(case1[cross_stoch,], case2[cross_stoch,], case3[cross_stoch,])),
          tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 1,
          breaks=seq(0, 1, 0.05), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)
heatmap.2(Reduce(cbind, list(case1[cross_control,], case2[cross_control,], case3[cross_control,])),
          tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 1, scale = 'none',
          breaks=seq(0, 1, 0.05), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90)





Names = c(paste('Twin1A', 'rep', 1:4, sep = '_'), paste('Twin1B', 'rep',  1:4, sep = '_'), 
          paste('Twin2A', 'rep', 1:2, sep = '_'), paste('Twin2B', 'rep',  1:2, sep = '_'),
          paste('Twin3A', 'rep', 1:2, sep = '_'), paste('Twin3B', 'rep',  1:2, sep = '_'))
pair = sapply(strsplit(Names, split = "_"), function(x) x[1])
pair = as.factor(pair)
levels(pair) = c("cornflowerblue", "cyan3", "coral1", "brown3", "chartreuse4", "darkolivegreen2")
pair = as.character(pair)


heatmap.2(Reduce(cbind, list(case1[cross_stoch,], case2[cross_stoch,], case3[cross_stoch,])),
          tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 1,
          breaks=seq(0, 1, 0.05), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90,
          ColSideColors = pair)

heatmap.2(Reduce(cbind, list(case1[cross_control,], case2[cross_control,], case3[cross_control,])),
          tracecol = NULL, trace = 'none', mar=c(10.1, 4.1), cexCol = 1, scale = 'none',
          breaks=seq(0, 1, 0.05), dendrogram = 'column', density.info = 'none', 
          offsetRow = -55, keysize = 1, labRow = '', cexRow = 2, srtRow = 90, 
          ColSideColors = pair)


