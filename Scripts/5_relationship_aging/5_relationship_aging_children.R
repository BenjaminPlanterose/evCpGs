############################################################################
############################################################################
###########                                                      ###########
###########             Relationship with age (children)         ###########
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
library(lmtest)
library(CpGassoc)
library(GEOquery)
library(car)
library(wateRmelon)
library(data.table)
library(minfi)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(data.table)

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
length(stochCpG) # 333

##### Population datasets ####
# setwd("where")
beta.sqn = fread("2020-04-17_SQN_nocombat_cellcomp.txt")
beta.sqn = process.beta.fread(beta.sqn)

# Children
# setwd("where")
phenotype <- getGEO('GSE104812', destdir=".")
pheno <- phenotype[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)
pheno <- pheno[, c(33:34)]
colnames(pheno) = c("age", "sex")
pheno$age = as.numeric(pheno$age)
pheno$sex = as.factor(pheno$sex)
rgSamples <- unlist(lapply(X = strsplit(colnames(beta.sqn), split = '_'), function(X) { paste(X[[1]])}))
pheno = pheno[rgSamples,]

age <- as.numeric(pheno$age)
gender <- as.factor(pheno$sex)
table(gender)/length(gender)
mean(age); sd(age)
min(age); max(age)
dim(beta.sqn)


#############################################   Age relationship in stochCpGs   #############################################

# Find available in Population dataset
cross_stoch <- stochCpG[stochCpG %in% rownames(beta.sqn)]
length(cross_stoch) # 333

# EWAS age
model1 <- cpg.assoc(beta.val = beta.sqn[cross_stoch,],
                    indep = age,
                    covariates = gender, fdr.method = 'bonferroni')
# Visualize
summary(model1)
plot(model1, tplot = T, classic = F)

# Res
res <- model1$results
pvals <- res$P.value
names(pvals) <- res$CPG.Labels

hist(-log10(pvals)); abline(v = -log10(0.05/333), lty = 2, col = "red2")
sig_age <- names(which(p.adjust(pvals, 'bonferroni') < 0.05)) # "cg24892069"
hist(-log10(pvals)); abline(v = -log10(0.05/333), lty = 2, col = "red2")


i = 13
plot(age, beta.sqn[cross_stoch[1],], pch = 19, cex = 0.5, ylim = c(0,1), col = gender, xlim = c(0,100),
     xlab = 'Age', ylab = 'Beta-value')
legend(x = 'top', col = c('red', 'black'), legend = c('Male', 'Female'), bty = 'n', pch =19, ncol = 2, cex = 1.2)

plot(age, beta.sqn["cg24892069",], pch = 19, cex = 0.5, ylim = c(0,1), col = gender, xlim = c(0,100),
     xlab = 'Age', ylab = 'Beta-value')
legend(x = 'top', col = c('red', 'black'), legend = c('Male', 'Female'), bty = 'n', pch =19, ncol = 2, cex = 1.2)


#############################################   Age heteroscedasticity in stochCpGs   #############################################


pvals.heter <- perform.whitetest(stochCpG, beta.sqn, age, gender)

sum(p.adjust(pvals.heter, method = 'bonferroni') < 0.05) # 2
heteros <- names(which(p.adjust(pvals.heter, 'bonferroni') < 0.05))
#  "cg20860188" "cg10908869"

thr = -log10(0.05/length(pvals.heter))
x = -log10(pvals.heter)
y = -log10(pvals)


par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
colours = rep('', times = length(x))
colours[p.adjust(pvals, 'bonferroni') < 0.05 & p.adjust(pvals.heter, 'bonferroni') < 0.05] <- 'purple'
colours[p.adjust(pvals, 'bonferroni') < 0.05 & p.adjust(pvals.heter, 'bonferroni') >= 0.05] <- 'red3'
colours[p.adjust(pvals, 'bonferroni') >= 0.05 & p.adjust(pvals.heter, 'bonferroni') < 0.05] <- 'blue3'
colours[p.adjust(pvals, 'bonferroni') >= 0.05 & p.adjust(pvals.heter, 'bonferroni') >= 0.05] <- 'black'

table(colours)/length(colours)
# colours
# black       blue3        red3 
# 0.990990991 0.006006006 0.003003003 

#############################################   Comparing heteroscedasticity and association   #############################################

# setwd("where")
tiff(filename = paste('pvals_age_heteroscedasticity', 'tiff', sep = '.'), width = 2.23, height = 2.23, units = 'in', res = 300, compression = 'none')
par(mar=c(1, 1, 1, 1), mgp=c(1, 0.015, 0), las=0)
plot(x, y, pch = 19, xlab = '-log10(p-val) (Age Heteroscedasticity)',
     ylab = '-log10(p-val) (Age dependency)', col = alpha(colours, 0.7),
     main = '', cex.main = 0.5, cex = 0.2, xaxt = 'n', yaxt = 'n')
axis(side = 1,at =  seq(0, 5, length.out = 3), labels = seq(0, 5, length.out = 3), tck = -0.01, cex.axis = 0.5)
axis(side = 2,at =  seq(0, 4, length.out = 3), labels = seq(0, 4, length.out = 3), tck = -0.01, cex.axis = 0.5)
abline(h = thr, lty = 2); abline(v = thr, lty = 2)
dev.off()

min(age) # 6.4
max(age) # 14.6





