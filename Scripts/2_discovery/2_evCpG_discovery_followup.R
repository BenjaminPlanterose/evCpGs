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
library(reshape2)
library(ggplot2)
library(RColorBrewer)


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

# Shuffle matching
shuffle_matching <- function(matching)
{
  matching <- as.matrix(matching)
  row_order <- sample(x = 1:nrow(matching), size = nrow(matching), replace = F)
  ordered_matching <- matching[row_order,]
  order_list <- lapply(1:nrow(matching), function(x) sample(x = c(1,2), size = 2, replace = F))
  posA <- unlist(lapply(order_list, function(x) x[1]))
  posB <- unlist(lapply(order_list, function(x) x[2]))
  shuffled_ordered_matching <- as.data.frame(matrix(nrow = nrow(matching), ncol = ncol(matching)))
  as_list <- lapply(1:nrow(ordered_matching), function(x) ordered_matching[x,])
  shuffled_ordered_matching$V1 <- as.vector(unlist(lapply(1:length(as_list), function(x) as_list[[x]][posA[x]])))
  shuffled_ordered_matching$V2 <- as.vector(unlist(lapply(1:length(as_list), function(x) as_list[[x]][posB[x]])))
  shuffled_ordered_matching$V2 <- as.vector(shuffled_ordered_matching[c(nrow(shuffled_ordered_matching), 1:(nrow(shuffled_ordered_matching)-1)),2])
  return(shuffled_ordered_matching)
}

# Subfunction employed by hvCpG.epsilon that actually performs evCpG discovery 
hvCpG.discovery2 <- function(beta, matching, epsilon, B)
{
  pvals <- matrix(nrow = nrow(beta), ncol = B)
  beta_H1 <- arrange.beta(beta, matching)
  for(i in 1:B)
  {
    set.seed(i); beta_H0 <- arrange.beta(beta, shuffle_matching(matching))
    pvals[,i] <- sapply(1:nrow(beta_H1), function(x) process.CpG(beta_H0[x,], beta_H1[x,]))
    print(i)
  }
  rownames(pvals) <- rownames(beta)
  colnames(pvals) <- 1:B
  return(pvals)
}


# Subfunction employed by hvCpG.discovery that actually performs equivalence testing based on
# a TOST that employs Yuen´s t-test
process.CpG <- function(row_h0, row_h1)
{
  notwin <- row_h0[seq(1, length(row_h0), 2)] - row_h0[seq(2, length(row_h0), 2)]
  twin <- row_h1[seq(1, length(row_h1), 2)] - row_h1[seq(2, length(row_h1), 2)]
  
  TOST <- rtost(x = abs(twin), y = abs(notwin), paired = F, epsilon = epsilon, 
                var.equal = FALSE, alpha = 0.05, tr = 0.2)
  res <- TOST$p.value
  
  return(res)
}

# Max per position for three matrices with equal dimensions
process.max <- function(pvals1, pvals2, pvals3)
{
  pvals <- matrix(nrow = nrow(pvals1), ncol = ncol(pvals1))
  rownames(pvals) <- rownames(pvals1); colnames(pvals) <- colnames(pvals1)
  for(j in 1:ncol(pvals1))
  {
    for(i in 1:nrow(pvals1))
    {
      pvals[i,j] <- max(c(pvals1[i,j], pvals2[i,j], pvals3[i,j]))
    }
    print(j)
  }
  return(pvals)
}


#############################################   Read processed data (see data_prep.R)   #############################################

# Set working directory
# setwd("where")

# Read matching twins
matching <- read.table(file = "matching_MZ.txt", header = T)

# Read raw data (all three normalisations)
# setwd("where")
beta1 <- fread('2019-08-21_SQN_combat_cellcomp.txt', nThread = 4)
beta1 <- process.beta.fread(beta1); dim(beta1) # 346555    852
beta2 <- fread('2019-08-21_dasen_combat_cellcomp.txt', nThread = 4)
beta2 <- process.beta.fread(beta2); dim(beta2) # 346555    852
beta3 <- fread('2019-08-22_oobrelicbmiqqn_combat_cellcomp.txt', nThread = 4)
beta3 <- process.beta.fread(beta3); dim(beta3) # 346555    852

# Read data
# setwd("where")
stochCpG <- as.vector(read.table(file = 'stochCpG.txt')$V1)
sig = stochCpG

# 
# setwd("where")
pvals1 <- readRDS('y_rTOST_SQN_ComBat_cellcomp_pvalues.rds')
all <- names(pvals1)
nosig = all[!(all %in% sig)]

#set.seed(1994); nosig <- nosig[sample(1:length(nosig), length(sig), replace = F)]
set.seed(1995); nosig <- nosig[sample(1:length(nosig), length(sig), replace = F)]

beta1_cont <- beta1[nosig,]
beta2_cont <- beta2[nosig,]
beta3_cont <- beta3[nosig,]

beta1 <- beta1[stochCpG,]
beta2 <- beta2[stochCpG,]
beta3 <- beta3[stochCpG,]


#############################################   Run inflation control   #############################################

# Criteria: median of delta twin
B = 100

epsilon = 0.02706526
pvals1 <- hvCpG.discovery2(beta1, matching, epsilon, B)

epsilon = 0.02739861
pvals2 <- hvCpG.discovery2(beta2, matching, epsilon, B)

epsilon = 0.03566612
pvals3 <- hvCpG.discovery2(beta3, matching, epsilon, B)

pvals <- process.max(pvals1, pvals2, pvals3)

# Control

# Criteria: median of delta twin
epsilon = 0.02706526
pvals1 <- hvCpG.discovery2(beta1_cont, matching, epsilon, B)

epsilon = 0.02739861
pvals2 <- hvCpG.discovery2(beta2_cont, matching, epsilon, B)

epsilon = 0.03566612
pvals3 <- hvCpG.discovery2(beta3_cont, matching, epsilon, B)
pvals_cont <- process.max(pvals1, pvals2, pvals3)




#############################################   Visualization   #############################################



# Quantitative
longData <- melt(-log10(pvals))
longData$value <- as.numeric(longData$value)
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  labs(x="B", y="CpGs", title="") + 
  scale_x_discrete(limits = 10*(1:10)) + 
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=4),
                     plot.title=element_text(size=11)) +
  scale_fill_gradient2(low="white", high="red3", #colors in the scale
                         breaks=seq(0,19,3), #breaks in the scale bar
                         limits=c(0, 19))

longData <- melt(-log10(pvals_cont))
longData$value <- as.numeric(longData$value)
ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  labs(x="B", y="CpGs", title="") + scale_x_discrete(limits = 10*(1:10)) +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=4),
                     plot.title=element_text(size=11)) +
  scale_fill_gradient2(low="white", high="red3", #colors in the scale
                       breaks=seq(0,19,3), #breaks in the scale bar
                       limits=c(0, 19))


# Qualitative
longData2 <- melt(pvals < 0.05/4652)
longData2$value <- as.factor(as.numeric(longData2$value))

hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')
ggplot(longData2, aes(x = Var2, y = Var1, fill = value)) + geom_tile() +
  scale_fill_brewer() + labs(x="B", y="CpGs", title="") + theme_bw()  + 
  theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
          axis.text.y=element_text(size=4),
          plot.title=element_text(size=11)) + scale_x_discrete(limits = 10*(1:10)) 


longData2 <- melt(pvals_cont < 0.05/4652)
longData2$value <- as.factor(as.numeric(longData2$value))

hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')
ggplot(longData2, aes(x = Var2, y = Var1, fill = value)) + geom_tile() +
  scale_fill_brewer() + labs(x="B", y="CpGs", title="") + theme_bw()  + 
  theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
        axis.text.y=element_text(size=4),
        plot.title=element_text(size=11)) + scale_x_discrete(limits = 10*(1:10)) 

# Summary statistics

# evCpGs
summary <- rowSums(pvals < 0.05/4652)
x <- barplot(sort(summary, decreasing = T), xaxt="n", ylab = 'Number of significant Bs', main = 'evCpGs',
             col = 'cornsilk4', border = 'cornsilk3', ylim = c(0,100))
labs <- names(sort(summary, decreasing = T))
text(x = x , y = -2, labs, xpd=TRUE, srt=90, cex = 0.2)

# control CpGs
summary <- rowSums(pvals_cont < 0.05/(4652-333))
x <- barplot(sort(summary, decreasing = T), xaxt="n", ylab = 'Number of significant Bs', main = 'controlCpGs',
             col = 'cornsilk4', border = 'cornsilk3', ylim = c(0,100))
labs <- names(sort(summary, decreasing = T))
text(x = x , y = -2, labs, xpd=TRUE, srt=90, cex = 0.2)

