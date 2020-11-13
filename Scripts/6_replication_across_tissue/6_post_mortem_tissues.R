############################################################################
############################################################################
###########                                                      ###########
###########                Post-mortem tissues                   ###########
###########             Author: Benjamin Planterose              ###########
###########                                                      ###########
###########        Erasmus MC University Medical Centre          ###########
###########               Rotterdam, The Netherlands             ###########
###########                                                      ###########
###########             b.planterose@erasmusmc.nl                ###########
###########                                                      ###########
############################################################################
############################################################################


library(GEOquery)
library(data.table)
library(minfi)
library(gplots)
library(prepocessCore)

recode <- function(tissues)
{
  layers <- character(length = length(tissues))
  for(i in 1:length(tissues))
  {
    if(sum(tissues[i] == ectoderm) != 0)
    {
      layers[i] <- 'ectoderm'
    }
    else if(sum(tissues[i] == mesoderm) != 0)
    {
      layers[i] <- 'mesoderm'
    }
    else if(sum(tissues[i] == endoderm) != 0)
    {
      layers[i] <- 'endoderm'
    }
  }
  return(layers)
}

preprocess_pheno <- function(x, class)
{
  blank <- as.list(pheno$title)
  individuals <- character()
  organs <- character()
  for(i in 1:length(class))
  {
    cut <- strsplit(pheno$title, split = class[i])
    for(j in 1:length(cut))
    {
      if(paste(cut[[j]], collapse = "") !=  blank[[j]])
      {
        individuals[j] <- class[i]
        organs[j] <- cut[[j]][1]
      }
    }
    print(i) 
  }
  organs <- unlist(lapply(X = strsplit(organs, split = '_'), function(x) paste(x, collapse = '')))
  organs <- unlist(lapply(X = strsplit(organs, split = ' '), function(x) paste(x, collapse = '')))
  result <- paste(individuals, organs, sep = '_')
  return(result)
}

process.beta.fread <- function(beta)
{
  cpgs <- as.vector(beta$rn)
  beta <- beta[,-(1:2)]
  beta <- as.matrix(beta)
  rownames(beta) <- cpgs
  return(beta)
}



# Read data
# setwd("where")
stochCpG <- as.vector(read.table(file = 'evCpGs.txt')$V1)
length(stochCpG) # 333

# setwd("where")
pvals1 <- readRDS('y_rTOST_SQN_ComBat_cellcomp_pvalues.rds')
tested0 = names(pvals1)
tested = tested0[!(tested0 %in% stochCpG)]
length(tested0) # 4652
length(tested) # 4319

# setwd("where")
control <- as.vector(read.table(file = 'control.txt')$V1)
length(control) # 998

# Tissues - embryonic layers
ectoderm <- c('Ischiaticnerve', 'Medullaoblongata')

mesoderm <- c('Adiposeabdominal', 'Adiposesubcutaenous', 'Aortaabdominal', 'Aortathoracic',
              'Bone', 'Bonemarrowred', 'Bonemarrowyellow', 'Coronaryartery', 'Jointcartilage',
              'Lymphnode', 'Splenicartery', 'Tonsils')

endoderm <- c('Bladder', 'Gallbladder', 'Gastricmucosa')

#######################################################################

## QC
# Data prep
# setwd("where")
raw  = fread("GSE50192_GPL13534_Matrix_signal_intensities_raw_data.txt", sep = "\t")
dim(raw) # 485577    211
CpGs = raw$`Sample name`
raw = raw[,-1]
raw = data.matrix(raw)
detP = raw[, seq(3, ncol(raw), 3)]
dim(detP) # 485577     70
rm(raw); gc()
#
PF_rate = colMeans(detP > 1E-06, na.rm = T)
sum(PF_rate > 0.05) # No samples excluded
#
PF_rate = rowMeans(detP > 1E-06, na.rm = T)
exclude = CpGs[which(PF_rate > 0.05)]
rm(detP); gc()
##


# setwd("where")
phenotype <- getGEO('GSE50192', destdir=".")
pheno <- phenotype[[1]]
pheno <- phenoData(pheno)
pheno <- pData(pheno)
pheno <- pheno[, c(1, 32, 33)]
pheno$title <- as.vector(pheno$title)
pheno$title <- preprocess_pheno(pheno$title, c('BM419/4', 'SJ600-5', 'KT538', 'KA522'))
individuals <- unlist(lapply(strsplit(pheno$title, split = '_'), function(x) x[1]))
tissues <- unlist(lapply(strsplit(pheno$title, split = '_'), function(x) x[2]))
layers <- recode(tissues)

# Extract raw methylation data
eset <- as.matrix(phenotype$GSE50192_series_matrix.txt.gz)
eset <- na.omit(eset)
dim(eset) # 484479     70
eset <- eset[!(rownames(eset) %in% exclude),]
dim(eset) # 484479     70

# Extract SNP data
sum(startsWith(rownames(eset), "rs")) # 65
rs =  rownames(eset)[which(startsWith(rownames(eset), "rs"))]

# Preprocess QN
CpGs = rownames(eset)
samples = colnames(eset)
eset = normalize.quantiles(eset)
rownames(eset) = CpGs
colnames(eset) = samples

#######################################  Produce MDS plots  ################################

# Find available features
cross_stoch = stochCpG[stochCpG %in% rownames(eset)]
length(cross_stoch) # 333
cross_control = control[control %in% rownames(eset)]
length(cross_control) # 995

# 
X = eset[cross_stoch,]
colnames(X) <- pheno$title
Y = eset[cross_control[1:length(cross_stoch)],]
colnames(Y) <- pheno$title
Z = eset[rs,]
colnames(Z) <- pheno$title
#

dim(X) # 333  70
dim(Y) # 333  70
dim(Z) # 65 70


#### Z
eig = cmdscale(dist(t(Z)), k = 2, eig = T)$eig
(comp = 100*eig[1:2]/sum(eig)) # 39.43820 32.53389
mdsPlot(Z, numPositions = nrow(Y), sampGroups = individuals, pch = 19, legendNCol = 1, 
        legendPos = 'bottomleft')


#### Y
eig = cmdscale(dist(t(Y)), k = 2, eig = T)$eig
(comp = 100*eig[1:2]/sum(eig)) # 28.76183 18.91611
mdsPlot(Y, numPositions = nrow(Y), sampGroups = individuals, pch = 19, legendNCol = 1, 
        legendPos = 'bottomleft')
mdsPlot(Y, numPositions = nrow(Y), sampGroups = layers, pch = 19, legendNCol = 2, 
        legendPos = 'bottomright')

##### X
eig = cmdscale(dist(t(X)), k = 2, eig = T)$eig
(comp = 100*eig[1:2]/sum(eig)) # 32.468568  8.897601
mdsPlot(X, numPositions = nrow(X), sampGroups = individuals, pch = 19, legendNCol = 2, 
        legendPos = 'topright')
mdsPlot(X, numPositions = nrow(X), sampGroups = layers, pch = 19, legendNCol = 2, 
        legendPos = 'topright')
