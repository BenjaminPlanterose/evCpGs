############################################################################
############################################################################
###########                                                      ###########
###########       Relationship with age (wide-age range)         ###########
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
library(car)
library(lmtest)
library(wateRmelon)
library(data.table)
library(CpGassoc)
library(qqman)
library(minfi)
library(RColorBrewer)
library(ggridges)
library(ggplot2)
library(gplots)

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

# Performs white test correcting for covariate gender in both primary and auxiliary linear models
perform.whitetest <- function(hvCpG, Population, age, gender)
{
  cross_hvCpG <- hvCpG[hvCpG %in% rownames(Population)]
  pvals1 <- numeric()
  for(i in 1:length(cross_hvCpG))
  {
    #plot(age, Population[hvCpG[i],], ylim = c(0,1), xlim = c(10, 100))
    df <- data.frame(CpG = Population[cross_hvCpG[i],], age = age, gender = gender)
    
    lm.fit <- lm(CpG ~ age + gender, df)
    pvals1[i] <- bptest(lm.fit, varformula = ~ age + gender + I(age^2), data = df)$p.value
    print(i)
    
    # Equivalent to:
    #R1 <- resid(lm.fit)  # Extract  the residuals
    #aux <- R1^2   # Square the residuals
    #aux_lm <- lm(aux ~ age + gender + I(age^2)) # Run the auxiliary regression
    #print(summary(aux_lm))       # Print the results
    #obs <- length(aux)
    #data2 <- summary(aux_lm)
    #R_Squared <- data2$r.squared
    #Test_Stat <- obs*R_Squared (n*R^2)
    #pchisq(Test_Stat, degree.freedom, lower.tail = FALSE)
    
  }
  names(pvals1) <- cross_hvCpG
  return(pvals1)
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

# Computes sd(beta) for windows of age defined by an offset and a bin.span
var_across_age <- function(x, age, bin.span, offset)
{
  hist_temp <- hist(age, breaks = seq(min(age), max(age), 1), plot = F)
  breaks <- hist_temp$breaks
  count <- hist_temp$counts
  
  window.size = bin.span + 1
  init = 1
  end = window.size
  meth_vec <- numeric()
  sd_vec <- numeric()
  age_vec <- numeric()
  n_i <- numeric()
  
  l = length(breaks)
  n = (l-window.size)/offset
  
  for(i in 1:n)
  {
    what = breaks[init:end]
    meth_vec[i] <- mean(x[age %in% what])
    sd_vec[i] <- sd(x[age %in% what])
    age_vec[i] <- mean(age[age %in% what])
    n_i[i] <- sum(age %in% what)
    
    init = init + offset
    end = end + offset
    print(what)
  }
  
  
  var_age = list(age = age_vec, sd = sd_vec, p = meth_vec, n = n_i)
  plot(var_age[[1]], var_age[[2]], xlab = 'Age', ylab = 'Var', type = 'o')
  
  return(var_age)
}

#############################################   Read data/phenotypes   #############################################

# Read evCpGs
# setwd("where")
stochCpG <- as.vector(read.table(file = 'evCpGs.txt')$V1)
sig = stochCpG

a = rowMeans(beta1[stochCpG,])
# setwd("where")
names(a) <- stochCpG
write.table(a, file = 'means_erisk.txt', row.names = T, col.names = F, quote = F, sep = '\t')


##### Population datasets ####

# Population
# setwd("where")
Population <- fread('2019-09-10_SQN_combat_cellcomp.txt', nThread = 4)
Population <- process.beta.fread(Population)
dim(Population) # 346469    727
iqrs1 = rowIQRs(Population)
names(iqrs1) = rownames(Population)
rm(Population); gc()

# setwd("where")
phenotype <- getGEO('GSE87571', destdir=".")
pheno <- phenotype[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)
pheno <- pheno[, c(37,39,40)]
IDs <- unlist(lapply(strsplit(colnames(Population), split = '_'), function(x) x[1]))
age <- as.numeric(pheno[IDs, 1])
gender <- as.factor(pheno[IDs, 2])


# Gambia
# setwd("where")
gambia <- fread('2019-09-10_SQN_combat_cellcomp_Gambia.txt')
gambia <- process.beta.fread(gambia)
dim(gambia) # 346099    240
iqrs2 = rowIQRs(gambia)
names(iqrs2) = rownames(gambia)
rm(gambia); gc()


##### Twin datasets ####
# TwinsUK
# setwd("where")
TwinsUK <- fread('2019-09-09_sqn_combat_cellcomp.txt', nThread = 4)
TwinsUK <- process.beta.fread(TwinsUK)
dim(TwinsUK) # 346705    656
matching_UK <- read.table('twin_matching.txt', header = F)
TwinsUK <- arrange.beta(TwinsUK, matching_UK)
iqrs3 = rowIQRs(TwinsUK)
names(iqrs3) = rownames(TwinsUK)
rm(TwinsUK); gc()


# setwd("where")
phenotype <- read.table(file = 'phenotype.txt', header = T)
barcodes <- as.vector(phenotype$Barcode)
age_UK <- phenotype$Age[barcodes %in% colnames(TwinsUK)]
names(age_UK) <- barcodes
TwinsUK <- TwinsUK[,names(age_UK)]

# E-risk
# setwd("where")
beta1 <- fread('2019-08-21_SQN_combat_cellcomp.txt', nThread = 4)
beta1 <- process.beta.fread(beta1); dim(beta1) # 346555    852
matching <- read.table('matching_MZ.txt', header = T)
beta1 <- arrange.beta(beta1, matching)

iqrs4 = rowIQRs(beta1)
names(iqrs4) = rownames(beta1)
rm(beta1); gc()


# Danish
# setwd("where")
danish <- fread('2019-09-09_sqn_combat_cellcomp.txt')
danish <- process.beta.fread(danish)
replicates_excluded <- c('GSM1506278', 'GSM1506587',
                         'GSM1506580', 'GSM1506342',
                         'GSM1506333', 'GSM1506438', 'GSM1506500',
                         'GSM1506327', 'GSM1506494', 'GSM1506554')
rgSamples <- unlist(lapply(strsplit(colnames(danish), split = '_'), function(x) x[1]))
danish <- danish[,!(rgSamples %in% replicates_excluded)]
dim(danish) # 345757    292
iqrs5 = rowIQRs(danish)
names(iqrs5) = rownames(danish)
rm(danish); gc()



matching_danish <- read.table('matching_danish.txt', header = T)
danish <- arrange.beta(danish, matching_danish)
rgSamples <- unlist(lapply(strsplit(colnames(danish), split = '_'), function(x) x[1]))

# setwd("where")
phenotype <- getGEO('GSE61496', destdir=".")
pheno <- phenotype[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)

dim(pheno) # 312 40
pheno <- pheno[,c(36, 38, 39)]
colnames(pheno) <- c('Age', 'Pair_id', 'Sex')

pheno <- pheno[rgSamples,]
dim(pheno) # 292 3
age_danish <- as.numeric(pheno$Age)
plot(age_danish[seq(1,length(age_danish),2)], age_danish[seq(2,length(age_danish),2)])

# Some minor discordance between ages (maybe age at sampling)
table(abs(age_danish[seq(1,length(age_danish),2)] - age_danish[seq(2,length(age_danish),2)]))
#   0   1   2 
# 138   7   1 


# Small population children
# setwd("where")
beta.sqn = fread("2020-04-17_SQN_nocombat_cellcomp.txt")
beta.sqn = process.beta.fread(beta.sqn)
dim(beta.sqn) # 347264     48
iqrs6 = rowIQRs(beta.sqn)
names(iqrs6) = rownames(beta.sqn)
rm(beta.sqn); gc()

list(iqrs1, iqrs2, iqrs3, iqrs4, iqrs5, iqrs6)

# setwd("where")
saveRDS(list(iqrs1, iqrs2, iqrs3, iqrs4, iqrs5, iqrs6), "IQR_list.rds")


#############################################   Age relationship in stochCpGs   #############################################

# Find available in Population dataset
cross_stoch <- stochCpG[stochCpG %in% rownames(Population)]
length(cross_stoch) # 331 out of 333

# EWAS age
model1 <- cpg.assoc(beta.val = Population[cross_stoch,],
                    indep = age,
                    covariates = gender, fdr.method = 'bonferroni')
res <- model1$results
pvals <- res$P.value
names(pvals) <- res$CPG.Labels
sig_age <- names(which(p.adjust(pvals, 'bonferroni') < 0.05))

# setwd("where")
a = pvals[stochCpG]
names(a) <- stochCpG
write.table(a, file = 'pvals_assoc.txt', row.names = T, col.names = F, quote = F, sep = '\t')


# Visualize
summary(model1)
plot(model1, tplot = T, classic = F)
create.manhattan.plot(pvals, sig_age, 'Age dependency', 0.05/length(pvals1))
# setwd("where")
for(i in 1:length(cross_stoch))
{
  tiff(filename = paste(cross_stoch[i], 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
  
  plot(age, Population[cross_stoch[i],], pch = 19, cex = 0.5, ylim = c(0,1), col = gender, xlim = c(0,100),
       xlab = 'Age', ylab = 'Beta-value')
  legend(x = 'top', col = c('red', 'black'), legend = c('Male', 'Female'), bty = 'n', pch =19, ncol = 2, cex = 1.2)
  # F = black
  # R = male
  dev.off()
  print(i)
  
}


#############################################   Age heteroscedasticity in stochCpGs   #############################################


pvals.heter <- perform.whitetest(stochCpG, Population, age, gender)
# setwd("where")
a = pvals.heter[stochCpG]
names(a) <- stochCpG
write.table(a, file = 'pvals_het.txt', row.names = T, col.names = F, quote = F, sep = '\t')


sum(p.adjust(pvals.heter, method = 'bonferroni') < 0.05) # 130
heteros <- names(which(p.adjust(pvals.heter, 'bonferroni') < 0.05))
a = venn(list(heteros = heteros, age = sig_age))

# Visualize
create.manhattan.plot(pvals.heter, names(which(pvals.heter < 0.05/length(pvals.heter))), 
                      'Age heteroscedasticity', 0.05/length(pvals.heter))

par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
colours = rep('', times = length(x))
colours[p.adjust(pvals, 'bonferroni') < 0.05 & p.adjust(pvals.heter, 'bonferroni') < 0.05] <- 'purple'
colours[p.adjust(pvals, 'bonferroni') < 0.05 & p.adjust(pvals.heter, 'bonferroni') >= 0.05] <- 'red3'
colours[p.adjust(pvals, 'bonferroni') >= 0.05 & p.adjust(pvals.heter, 'bonferroni') < 0.05] <- 'blue3'
colours[p.adjust(pvals, 'bonferroni') >= 0.05 & p.adjust(pvals.heter, 'bonferroni') >= 0.05] <- 'black'

table(colours)/length(colours)
# colours
# black     blue3    purple      red3 
# 0.2386707 0.1903323 0.2024169 0.3685801 

thr = -log10(0.05/length(pvals.heter))
x = -log10(pvals.heter)
y = -log10(pvals)
# setwd("where")
tiff(filename = paste('pvals_age_heteroscedasticity', 'tiff', sep = '.'), width = 2.23, height = 2.23, units = 'in', res = 300, compression = 'none')

#############################################   Comparing heteroscedasticity and association   #############################################

par(mar=c(1, 1, 1, 1), mgp=c(3, 1, 0), las=0)
plot(x, y, pch = 19, xlab = '-log10(p-val) (Age Heteroscedasticity)',
     ylab = '-log10(p-val) (Age dependency)', col = alpha(colours, 0.7),
     main = '', cex.main = 0.5, cex = 0.2, xaxt = 'n', yaxt = 'n')
axis(side = 1,at =  seq(0, 14, length.out = 3), labels = seq(0, 14, length.out = 3), tck = -0.01, cex.axis = 0.5)
axis(side = 2,at =  seq(0, 150, length.out = 3), labels = seq(0, 150, length.out = 3), tck = -0.01, cex.axis = 0.5)
abline(h = thr, lty = 2); abline(v = thr, lty = 2)
dev.off()


# 2 most remarkable examples

# Example 1
# setwd("where")
tiff(filename = paste('epigenetic_drift', 'tiff', sep = '.'), width = 2.23, height = 2.23, units = 'in', res = 300, compression = 'none')
i = 'cg00639615'
var_age <- var_across_age(Population[i,], age, 10, 2)
colours <- gender
levels(colours) <- c('deeppink3', 'dodgerblue3')
par(mar=c(2, 2, 2, 1), mgp=c(3, 1, 0), las=0)
plot(age, Population[i,], ylim = c(0.2, 0.8), pch = 19, cex = 0.2, 
     main = '', ylab = 'Beta-value', xlab = 'Age (years)', col = alpha(colours, 0.5),
     xaxt = 'n',yaxt = 'n', lwd = 0.5, bty = 'n')
legend(x ='top', col = c('deeppink3', 'dodgerblue3'), legend = c('', ''), pch = 19,
       cex = 0.2, bty = 'n', ncol = 2)
axis(2,at=seq(0.2, 0.7, length.out = 3),labels= seq(0.2, 0.7, length.out = 3), tck = -0.01, cex.axis = 0.5)
axis(1,at=seq(0,100, length.out = 6),labels= seq(0,100, length.out = 6), tck = -0.01, cex.axis = 0.5)
lines(var_age$age, (var_age$p) + (var_age$sd), col = 'red2', lty = 2, lwd = 1)
lines(var_age$age, (var_age$p) - (var_age$sd), col = 'red2', lty = 2, lwd = 1)
dev.off()

# Example 2
# setwd("where")
tiff(filename = paste('epigenetic_drift', 'tiff', sep = '.'), width = 2.23, height = 2.23, units = 'in', res = 300, compression = 'none')
i = 'cg23479922'
var_age <- var_across_age(Population[i,], age, 10, 2)
colours <- gender
levels(colours) <- c('deeppink3', 'dodgerblue3')
par(mar=c(2, 2, 2, 1), mgp=c(3, 1, 0), las=0)
plot(age, Population[i,], ylim = c(0.15, 0.75), pch = 19, cex = 0.2,
     main = '', ylab = 'Beta-value', xlab = 'Age (years)', col = alpha(colours, 0.5),
     xaxt = 'n',yaxt = 'n', lwd = 0.5, bty = 'n')
legend(x ='top', col = c('deeppink3', 'dodgerblue3'), legend = c('', ''), pch = 19,
       cex = 0.2, bty = 'n', ncol = 2)
axis(2,at=seq(0.2, 0.8, length.out = 4),labels= seq(0.2, 0.8, length.out = 4), tck = -0.01, cex.axis = 0.5)
axis(1,at=seq(0,100, length.out = 6),labels= seq(0,100, length.out = 6), tck = -0.01, cex.axis = 0.5)
lines(var_age$age, (var_age$p) + (var_age$sd), col = 'red2', lty = 2, lwd = 1)
lines(var_age$age, (var_age$p) - (var_age$sd), col = 'red2', lty = 2, lwd = 1)
dev.off()


#############################################   evCpG variation in TwinsUK is greater than in E-risk   #############################################

# evCpGs not discarded by QC in TwinsUK
cross_sig <- Reduce(intersect, list(sig, rownames(TwinsUK), rownames(beta1)))

# Compute |delta beta|
delta.UK <- abs(TwinsUK[cross_sig, seq(1,ncol(TwinsUK),2)]-
                  TwinsUK[cross_sig, seq(2,ncol(TwinsUK),2)])
delta.erisk <- abs(beta1[cross_sig, seq(1,ncol(beta1),2)]-
                   beta1[cross_sig, seq(2,ncol(beta1),2)])

# Perform Kolmogorov-Smirnov test
# setwd("where")
tiff(filename = paste('erisk_uk_deltabeta', 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
pval_KS <- ks.test(as.numeric(delta.UK), as.numeric(delta.erisk), alternative = 'greater')$p.value
plot(ecdf(delta.UK), col = 'firebrick3', main = paste('-log10(p-val_KS) =', round(-log10(pval_KS), 4), '\n',  length(cross_sig), 'out of', length(sig), 'stochCpGs', sep = ' '),
     xlab = expression(paste("|", Delta, beta, '|')), ylab = 'Empirical Cumulative Distribution Function',
     lwd = 2, cex = 2)
lines(ecdf(delta.erisk), lwd = 2, col = 'deepskyblue3')
legend(x = 'right', col = c('deepskyblue3', 'firebrick3'), legend = c('E-risk', 'TwinsUK'), cex = 2, bty = 'n', pch = 19)
dev.off()

#############################################   IQR per dataset   #############################################

# Compute IQRs
iqrs_population <- rowIQRs(Population[stochCpG[stochCpG %in% rownames(Population)],])
names(iqrs_population) <- stochCpG[stochCpG %in% rownames(Population)]

iqrs_gambia <- rowIQRs(gambia[stochCpG[stochCpG %in% rownames(gambia)],])
names(iqrs_gambia) <- stochCpG[stochCpG %in% rownames(gambia)]

iqrs_UK <- rowIQRs(TwinsUK[stochCpG[stochCpG %in% rownames(TwinsUK)],])
names(iqrs_UK) <- stochCpG[stochCpG %in% rownames(TwinsUK)]

iqrs_erisk <- rowIQRs(beta1[stochCpG[stochCpG %in% rownames(beta1)],])
names(iqrs_erisk) <- stochCpG[stochCpG %in% rownames(beta1)]


iqrs_danish <- rowIQRs(danish[stochCpG[stochCpG %in% rownames(danish)],])
names(iqrs_danish) <- stochCpG[stochCpG %in% rownames(danish)]

##### Ridge plot ####

iqrs1 = iqrs1[names(iqrs1) %in% stochCpG]
iqrs2 = iqrs2[names(iqrs2) %in% stochCpG]
iqrs3 = iqrs3[names(iqrs3) %in% stochCpG]
iqrs4 = iqrs4[names(iqrs4) %in% stochCpG]
iqrs5 = iqrs5[names(iqrs5) %in% stochCpG]
iqrs6 = iqrs6[names(iqrs6) %in% stochCpG]



# Prepare data
df <- matrix(NA, ncol = 6, nrow = length(stochCpG))
df <- as.data.frame(df)
colnames(df) <- c('Erisk', 'Danish', 'TwinsUK', 'Population', 'Gambia', "Children")
rownames(df) <- stochCpG
df[names(iqrs4), 1] <- iqrs4 
df[names(iqrs5), 2] <- iqrs5
df[names(iqrs3), 3] <- iqrs3
df[names(iqrs1), 4] <- iqrs1
df[names(iqrs2), 5] <- iqrs2
df[names(iqrs6), 6] <- iqrs6

melted <- melt(df)
status <- as.vector(melted$variable)
status[status %in% c('Erisk', 'Danish', 'TwinsUK')] <- 'Twin study'
status[status %in% c('Population', 'Gambia', "Children")] <- 'Non-twin study'
melted$study_type <- status
melted <- na.omit(melted)
melted$variable <- factor(melted$variable, levels = c("Population","Children", "Gambia", "Danish", "TwinsUK", "Erisk"))

# Perform plot
# setwd("where")
tiff(filename = paste('ridge_age', 'tiff', sep = '.'), width = 2.5, height = 2.5, units = 'in', res = 300, compression = 'none')
ggplot(melted, aes(x = value, y = variable, fill = study_type)) +
  geom_density_ridges(alpha =0.5, size = 0.1) + theme_ridges(grid = T, center_axis_labels = T, font_size = 4, line_size = 0.25) +
  scale_x_continuous(breaks = seq(0, 0.3, length.out = 4), limits = c(0,0.3)) +
  xlab('IQR') + ylab('Dataset') + scale_fill_brewer(palette = "Accent") + theme(legend.title = element_text(size = 1), 
                                                                                legend.text = element_text(size = 1))
  #geom_vline(xintercept = 0.07, lty = 2)
dev.off()


