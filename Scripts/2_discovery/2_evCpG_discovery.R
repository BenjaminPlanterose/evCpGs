############################################################################
############################################################################
###########                                                      ###########
###########            evCpG Discovery (E-risk study)            ###########
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
library(parallel)
library(gplots)
library(minfi)
library(equivalence)
library(DescTools)
library(gplots)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(qqman)
library(VennDiagram)

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
  
  np <- detectCores(logical = FALSE)
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
  
  np <- detectCores(logical = FALSE)
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

# Wrapper that performs Manhattan plots on 45K data based on the qqman library. As this library is 
# originally designded for GWAS data, some preprocessing was required to make it compatible
create.manhattan.plot <- function(pvals, hvCpG, main, alpha_corrected)
{
  CpGs <- names(pvals)
  # Extract annotation. The easiest way to do so, strangely, is to create a directory with 
  # a single sample IDAT (both channels). Read the IDAT and extract the annotation with
  # minfi functions.
  old.dir <- getwd(); # setwd("where")
  RGSET <- read.metharray.exp(getwd()); setwd(old.dir)
  annotation <- getAnnotation(RGSET)
  names <- annotation$Name
  pos <- annotation$pos[na.omit(match(names, CpGs))]
  
  # Chr
  chr <- annotation$chr[match(CpGs, names)]
  chr <- unlist(lapply(strsplit(chr, split = 'chr'), FUN = function(X) X[2]))
  chr <- as.numeric(chr)
  lev <- unique(chr)
  
  # Build data frame in a GWAS-like formatting just for recycling qqman function
  print(length(CpGs))
  print(length(chr))
  print(length(pos))
  print(length(pvals))
  df <- data.frame(SNP = CpGs, CHR = chr, BP = pos, P = pvals)
  
  # In order to place labels side ways so that they are not cropped.
  x_axis <- list()
  for(i in 1:length(lev))
  {
    x_axis[[i]] <- range(df[df$CHR == lev[i], 3])
  }
  ticks <- numeric()
  ticks[1] <- mean(x_axis[[1]])
  acum = 0
  for(i in 2:length(x_axis))
  {
    acum <- acum + x_axis[[i-1]][2]
    ticks[i] <- acum + mean(x_axis[[i]])
  }
  # Perfoms manhattan plot
  qqman::manhattan(df, cex = 0.3, xaxt = 'n', yaxt = 'n',
                   cex.axis = 0.2, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = -log10(alpha_corrected), highlight = hvCpG)
  # Adds labels
  axis(2,at=seq(0, 25, length.out = 6),labels= seq(0, 25, length.out = 6), tck = -0.01, cex.axis = 0.5)
  axis(1,at=ticks,labels= lev, tck = -0.01, las = 2, cex.axis = 0.5)
}

# Employed to visualize all 333 positives
visualize.CpG.testing <- function(hvCpG, beta1, beta2, beta3, ylim, cex)
{
  beta_list <- list(beta1, beta2, beta3)
  title <- c('StrQN_ComBat', 'Dasen_ComBat', 'oob_RELIC_BMIQ_QN')
  epsilon <- c(0.02706526, 0.02739861, 0.03566612)
  
  par(mar=c(5.1, 4.1, 7.1, 2.1), mgp=c(3, 1, 0), las=0)
  par(mfrow = c(3,1))
  
  for(i in 1:3)
  {
    row <- beta_list[[i]][hvCpG,]
    twin <- row[seq(1, length(row), 2)] - row[seq(2, length(row), 2)]
    notwin <- compute.delta.beta.notwin(row)
    
    TOST <- rtost(x = abs(twin), y = abs(notwin), paired = F, epsilon = epsilon[i], 
                  var.equal = FALSE, alpha = 0.05, tr = 0.2)
    res <- TOST$p.value
    
    plot(density(abs(twin)), lwd = 2, xlim = c(-0.1,1), ylim = ylim, main = title[i], xlab = expression(paste('|', Delta, beta, '|')), cex = cex)
    lines(density(abs(notwin)), lwd = 2, col = 'red', lty = 2)
    text(x = c(0.47), y = c(7.5), round(-log10(res), 4), cex = cex)
    text(x = c(0.47), y = c(9.5), '-log10(p-val)', font = 2, cex = cex)
    legend(x = 'right', legend = c(paste('Twin pairs, n =', length(twin)), paste('Unrelated pairs, n =', length(notwin))), fill = c('black', 'red'), bty = 'n', cex = cex)
  }
  
  mtext(hvCpG, side = 3, line = -2, at = c(0.52), outer = TRUE, font = 2)
  
  par(mfrow = c(1,1))
  par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
}



#############################################   Read processed data (see data_prep.R)   #############################################

# Set working directory
# setwd("where")

# Read matching twins
matching <- read.table(file = "matching_MZ.txt", header = T)
head(matching)

# Read raw data (all three normalisations)
# setwd("where")
beta1 <- fread('2019-08-21_SQN_combat_cellcomp.txt', nThread = 4)
beta1 <- process.beta.fread(beta1); dim(beta1) # 346555    852
beta2 <- fread('2019-08-21_dasen_combat_cellcomp.txt', nThread = 4)
beta2 <- process.beta.fread(beta2); dim(beta2) # 346555    852
beta3 <- fread('2019-08-22_oobrelicbmiqqn_combat_cellcomp.txt', nThread = 4)
beta3 <- process.beta.fread(beta3); dim(beta3) # 346555    852


#############################################   Post-normalisation filters   #############################################

# Low ICC probes
# setwd("where")
dataset <- fread('ICC_values.csv')
probes2remove <- as.vector(dataset$`ilmnid(CpG site)`[dataset$`ICC value` < 0.37]); length(probes2remove) # 270527

# Get rid of low ICC probes
beta1 <- beta1[-na.omit(match(probes2remove, rownames(beta1))),]; dim(beta1) # 146833    852
beta2 <- beta2[-na.omit(match(probes2remove, rownames(beta2))),]; dim(beta2) # 146833    852
beta3 <- beta3[-na.omit(match(probes2remove, rownames(beta3))),]; dim(beta3) # 146833    852

# Compute IQRs per CpG
IQRs1 <- sapply(X = 1:nrow(beta1), FUN = function(x) IQR(beta1[x,]))
names(IQRs1) <- rownames(beta1)
IQRs2 <- sapply(X = 1:nrow(beta2), FUN = function(x) IQR(beta2[x,]))
names(IQRs2) <- rownames(beta2)
IQRs3 <- sapply(X = 1:nrow(beta3), FUN = function(x) IQR(beta3[x,]))
names(IQRs3) <- rownames(beta3)

# Justify an IQR threshold
mu = 0.5; sigma = 0.05
alpha = ((1-mu)/(sigma^2) - 1/mu)*mu^2 # estimate alpha parameter from a beta distribution function
beta = alpha*(1/mu - 1) # estimate beta parameter from a beta distribution function
dist <- rbeta(n = 10000000, shape1 = alpha, shape2 = beta)
plot(density(dist))
IQR_threshold <- round(IQR(dist), digits = 2); print(IQR_threshold) # 0.07

#PLOT
par(mfrow = c(1,1))
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(density(IQRs1),ylim = c(0,52), col = 'aquamarine3', lwd = 3, main = 'IQR density plot')
lines(density(IQRs2), col = 'brown2', lwd = 2)
lines(density(IQRs3), col = 'chartreuse4', lwd = 2)
legend(x = 'right', fill = c('aquamarine3', 'brown2', 'chartreuse4'), legend = c('SQN', 'dasen', 'oob_relic_qn_bmiq'),
       bty = 'n', cex = 1.5)
abline(v = IQR_threshold, lty = 2)
#PLOT

# CRITERIA: All norm higher than 0.07
intersection <- Reduce(intersect, list(names(which(IQRs1 > IQR_threshold)), names(which(IQRs2 > IQR_threshold)),names(which(IQRs3 > IQR_threshold))))
length(intersection) # 4652
beta1 <- beta1[intersection,]; dim(beta1) # 4652  852
beta2 <- beta2[intersection,]; dim(beta2) # 4652  852
beta3 <- beta3[intersection,]; dim(beta3) # 4652  852


#############################################   Establishing an epsilon threshold   #############################################

means1 <- hvCpG.epsilon(beta1, matching); gc()
means2 <- hvCpG.epsilon(beta2, matching); gc()
means3 <- hvCpG.epsilon(beta3, matching); gc()

par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
par(mfrow = c(3,1))
SQN_twin <- unlist(lapply(means1, function(x) x[1]))
SQN_notwin <- unlist(lapply(means1, function(x) x[2]))
plot(density(SQN_twin), xlim = c(0,0.4), lwd = 2, col = 'chocolate', main = paste('SQN, mean(delta.beta.twin) =', round(mean(abs(SQN_twin)), 4), ', mean(delta.beta.unrelated) =', round(mean(abs(SQN_notwin)), 4)))
lines(density(SQN_notwin), lwd = 2, col = 'chocolate2')
deltamu_SQN <- SQN_notwin - SQN_twin
abline(v = median(SQN_twin), lty = 2)

dasen_twin <- unlist(lapply(means2, function(x) x[1]))
dasen_notwin <- unlist(lapply(means2, function(x) x[2]))
plot(density(dasen_twin),  lwd = 2, col = 'aquamarine4', xlim = c(0,0.4), main = paste('dasen, mean(delta.beta.twin) =', round(mean(abs(dasen_twin)), 4), ', mean(delta.beta.unrelated) =', round(mean(abs(dasen_notwin)), 4)))
lines(density(dasen_notwin),  lwd = 2, col = 'aquamarine4')
deltamu_dasen <- dasen_notwin - dasen_twin
abline(v = median(dasen_twin), lty = 2)

bmiq_twin <- unlist(lapply(means3, function(x) x[1]))
bmiq_notwin <- unlist(lapply(means3, function(x) x[2]))
plot(density(bmiq_twin), lwd = 2, col = 'chartreuse4', xlim = c(0,0.4), main = paste('oob_relic_qn_bmiq, mean(delta.beta.twin) =', round(mean(abs(bmiq_twin)), 4), ', mean(delta.beta.unrelated) =', round(mean(abs(bmiq_notwin)), 4)))
lines(density(bmiq_notwin), lwd = 2, col = 'chartreuse4')

deltamu_bmiq <- bmiq_notwin - bmiq_twin
abline(v = median(deltamu_bmiq), lty = 2)
abline(v = median(bmiq_twin), lty = 2)
par(mfrow = c(1,3))

plot(density(deltamu_SQN), main = 'SQN'); abline(v = median(deltamu_SQN), lty = 2)
plot(density(deltamu_dasen), main = 'dasen'); abline(v = median(deltamu_dasen), lty = 2)
plot(density(deltamu_bmiq), main = 'BMIQ'); abline(v = median(deltamu_bmiq), lty = 2)

#############################################   evCpG discovery   #############################################

# Criteria: median of delta twin
epsilon = median(SQN_twin); print(epsilon) # 0.02706526
pvals1 <- hvCpG.discovery(beta1, matching, epsilon); gc()

epsilon = median(dasen_twin); print(epsilon) # 0.02739861
pvals2 <- hvCpG.discovery(beta2, matching, epsilon); gc()

epsilon = median(bmiq_twin); print(epsilon) # 0.03566612
pvals3 <- hvCpG.discovery(beta3, matching, epsilon); gc()

# setwd("where")
saveRDS(pvals1, "y_rTOST_SQN_ComBat_cellcomp_pvalues.rds")
saveRDS(pvals2, "y_rTOST_dasen_ComBat_cellcomp_pvalues.rds")
saveRDS(pvals3, "y_rTOST_oobRELICQNBMIQ_ComBat_cellcomp_pvalues.rds")


#############################################   Discovery follow-up   #############################################

# setwd("where")

pvals1 <- readRDS('y_rTOST_SQN_ComBat_cellcomp_pvalues.rds')
pvals2 <- readRDS('y_rTOST_dasen_ComBat_cellcomp_pvalues.rds')
pvals3 <- readRDS('y_rTOST_oobRELICQNBMIQ_ComBat_cellcomp_pvalues.rds')


pvals <- sapply(1:length(pvals1), function(x) max(c(pvals1[[x]], pvals2[[x]], pvals3[[x]])))
names(pvals) <- names(pvals1)

# setwd("where")
write.table(x = pvals[sig], col.names = F, row.names = T, quote = F, sep = '\t', file = 'pvalues_discovery.txt')



# Extract significant in all normalizations
alpha <- 0.05/length(pvals1) # Bonferroni threshold
sig1 <- names(which(unlist(pvals1) < alpha))
sig2 <- names(which(unlist(pvals2) < alpha))
sig3 <- names(which(unlist(pvals3) < alpha))
sig <- Reduce(intersect, list(sig1, sig2, sig3))
length(sig) # 333
# setwd("where")
write.table(x = sig, file = 'stochCpG.txt', quote = F, row.names = F, col.names = F, sep = '\t')


#############################################   Discovery visualizations   #############################################

# Make Venn diagram for different normalisations
a <- venn(list(sig1, sig2, sig3))
a <- attr(a, 'intersections')
# setwd("where")
tiff(filename = paste('venn', 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
grid.newpage()
draw.triple.venn(area1 = length(sig1), area2 = length(sig2), area3 = length(sig3), 
                 n12 = length(intersect(sig1, sig2)), n13 = length(intersect(sig1, sig3)), 
                 n23 = length(intersect(sig2, sig3)), n123 = length(sig), 
                 category = c("StrQN", "Dasen", "oob_RELIC_BMIQ_QN"), lty = 'blank', 
                 fill = c("skyblue", "goldenrod1", "darkorchid1"), cex = 2, cat.cex = 2)
dev.off()


# Make separate manhattan plots
par(mfrow = c(3,1))
# setwd("where")
tiff(filename = paste('SQN_manhattan', 'tiff', sep = '.'), width = 7.5, height = 10, units = 'in', res = 300, compression = 'none')
create.manhattan.plot(unlist(pvals1), sig, 'SQN', alpha)
dev.off()

# setwd("where")
tiff(filename = paste('dasen_manhattan', 'tiff', sep = '.'), width = 7.5, height = 10, units = 'in', res = 300, compression = 'none')
create.manhattan.plot(unlist(pvals2), sig, 'Dasen', alpha)
dev.off()

# setwd("where")
tiff(filename = paste('bmiq_manhattan', 'tiff', sep = '.'), width = 7.5, height = 10, units = 'in', res = 300, compression = 'none')
create.manhattan.plot(unlist(pvals3), sig, 'BMIQ', alpha)
dev.off()


# Make combined manhattan plot - criteria pval-j = maxi(pval-ij); i: normalisation, j: CpG
pvals_list <- list(unlist(pvals1), unlist(pvals2), unlist(pvals3))
pvals_combined <- lapply(1:length(pvals1), function(x) max(c(pvals_list[[1]][x], pvals_list[[2]][x], pvals_list[[3]][x])))
names(pvals_combined) <- names(pvals1)
alpha <- 0.05/length(pvals_combined) # Bonferroni threshold
# setwd("where")
tiff(filename = paste('combined_manhattan', 'tiff', sep = '.'), width = 3.346, height = 5.238, units = 'in', res = 300, compression = 'none')
par(mar=c(2.1, 2.1, 0, 2.1), mgp=c(3, 1, 0), las=0)
create.manhattan.plot(unlist(pvals_combined), sig, 'Combined', alpha)
dev.off()


# Visualize a TP and a FP
TP = 'cg03557725'
TN = 'cg23733394' 
epsilon <- c(0.02706526, 0.02739861, 0.03566612)
# TP
row <- beta1[TP,]
twin <- row[seq(1, length(row), 2)] - row[seq(2, length(row), 2)]
notwin <- compute.delta.beta.notwin(row)
TOST <- rtost(x = abs(twin), y = abs(notwin), paired = F, epsilon = epsilon[1], 
              var.equal = FALSE, alpha = 0.05, tr = 0.2)
res <- TOST$p.value

tiff(filename = paste('TP', 'tiff', sep = '.'), width = 3.346, height = 2.619, units = 'in', res = 300, compression = 'none')
par(mar = c(1.5, 1.5, 1.5, 0.5), oma = c(0, 0, 0, 0))
plot(density(abs(twin)), lwd = 2, xlim = c(-0.05,0.6), ylim = c(0, 26), main = TP, cex.main = 0.7, cex.axis = 0.7, xlab = '',
     xaxt = 'n', yaxt = 'n', ylab = '')
lines(density(abs(notwin)), lwd = 2, col = 'red', lty = 2)
# Axis
axis(1, at = seq(0,0.6,0.1), label = rep('', length(seq(0, 0.6, 0.1))), tck = -0.01)
axis(1, at = seq(0,0.6,0.1), label = as.character(seq(0, 0.6, 0.1)), tck = -0.01, line = -1.2, lwd = 0, cex.axis = 0.5)
axis(2, at = seq(0,25,5), label = rep('', length(seq(0,25,5))), tck = -0.01)
axis(2, at = seq(0,25,5), label = as.character(seq(0,25,5)), tck = -0.01, line = -1, lwd = 0, cex.axis = 0.5)
# axis labels
mtext(expression(paste('|', Delta, beta, '|')), side = 1, outer = F, cex = 0.6, line = 0.6, col = "black")
mtext('Density', side = 2, outer = F, cex = 0.6, line = 0.6, col = "black")
# p-values
text(x = c(0.3), y = c(14), round(-log10(res), 4), cex = 0.5)
text(x = c(0.3), y = c(16), '-log10(p-val)', font = 2, cex = 0.5)
# legend
legend(x = 'top', legend = c(paste('Twin pairs, n =', length(twin)), paste('Unrelated pairs, n =', length(notwin))), 
       fill = c('black', 'red'), bty = 'n', cex = 0.5)
dev.off()

# TN
row <- beta1[TN,]
twin <- row[seq(1, length(row), 2)] - row[seq(2, length(row), 2)]
notwin <- compute.delta.beta.notwin(row)
TOST <- rtost(x = abs(twin), y = abs(notwin), paired = F, epsilon = epsilon[1], 
              var.equal = FALSE, alpha = 0.05, tr = 0.2)
res <- TOST$p.value
tiff(filename = paste('TN', 'tiff', sep = '.'), width = 3.346, height = 2.619, units = 'in', res = 300, compression = 'none')
par(mar = c(1.5, 1.5, 1.5, 0.5), oma = c(0, 0, 0, 0))
plot(density(abs(twin)), lwd = 2, xlim = c(-0.05,0.6), ylim = c(0, 26), main = TP, cex.main = 0.7, cex.axis = 0.7, xlab = '',
     xaxt = 'n', yaxt = 'n', ylab = '')
lines(density(abs(notwin)), lwd = 2, col = 'red', lty = 2)
# Axis
axis(1, at = seq(0,0.6,0.1), label = rep('', length(seq(0, 0.6, 0.1))), tck = -0.01)
axis(1, at = seq(0,0.6,0.1), label = as.character(seq(0, 0.6, 0.1)), tck = -0.01, line = -1.2, lwd = 0, cex.axis = 0.5)
axis(2, at = seq(0,25,5), label = rep('', length(seq(0,25,5))), tck = -0.01)
axis(2, at = seq(0,25,5), label = as.character(seq(0,25,5)), tck = -0.01, line = -1, lwd = 0, cex.axis = 0.5)
# axis labels
mtext(expression(paste('|', Delta, beta, '|')), side = 1, outer = F, cex = 0.6, line = 0.6, col = "black")
mtext('Density', side = 2, outer = F, cex = 0.6, line = 0.6, col = "black")
# p-values
text(x = c(0.3), y = c(14), round(-log10(res), 4), cex = 0.5)
text(x = c(0.3), y = c(16), '-log10(p-val)', font = 2, cex = 0.5)
# legend
legend(x = 'top', legend = c(paste('Twin pairs, n =', length(twin)), paste('Unrelated pairs, n =', length(notwin))), 
       fill = c('black', 'red'), bty = 'n', cex = 0.5)
dev.off()

# Visualize all positives
beta1 <- beta1[stochCpG,]
beta2 <- beta2[stochCpG,]
beta3 <- beta3[stochCpG,]
# setwd("where")
for(i in 1:length(stochCpG))
{
  tiff(filename = paste(stochCpG[i], 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
  visualize.CpG.testing(stochCpG[i], beta1, beta2, beta3, c(0,20), 1.5)
  dev.off()
  print(i)
}