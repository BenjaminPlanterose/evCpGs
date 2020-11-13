############################################################################
############################################################################
###########                                                      ###########
###########               WGBS analysis (adipose)                ###########
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
library(RColorBrewer)
library(scales)
library(matrixStats)
library(gplots)


## Load functions ##

# Reading files
read.files <- function()
{
  files <- list.files() # in order: MZ11, MZ12, MZ21, MZ22, etc
  a <- 1
  for(i in seq(from = 1, to = length(files), by = 2))
  { 
    assign(paste('twin', a, 'L', sep = ''), fread(files[i], nThread = 4),
           envir = parent.frame())
    assign(paste('twin', a, 'R', sep = ''), fread(files[i + 1], nThread = 4),
           envir = parent.frame())
    
    print(a)
    a <- a + 1
    
  }
  gc()
  
  
}
read.files2 <- function()
{
  files <- list.files(pattern = '*.txt') # in order: MZ11, MZ12, MZ21, MZ22, etc
  a <- 1
  for(i in seq(from = 1, to = length(files), by = 2))
  { 
    assign(paste('twin', a, 'L', sep = ''), fread(files[i], nThread = 4),
           envir = parent.frame())
    assign(paste('twin', a, 'R', sep = ''), fread(files[i + 1], nThread = 4),
           envir = parent.frame())
    
    print(a)
    a <- a + 1
    
  }
  gc()
  
  
}

# Processing files
process.files <- function()
{
  files <- list.files()
  
  # Processing: 
  # i) Select positions with less than 20 % methylation difference between strands. 
  # ii) Filter positions with coverage less than two reads per strand
  # iii) add ID
  
  for(i in 1:(length(files)/2))
  {
    temp <- get(paste('twin', i, 'L', sep = ''))
    temp <- temp[abs(temp$meth_fw - temp$meth_rv) < 20,]
    temp <- temp[temp$num_total_fw > 1 & temp$num_total_rv > 1, ]
    temp$ID <- paste(paste(temp$`#chr`, temp$start, sep = ':'), temp$end, sep = '-')
    assign(paste('twin', i, 'L', sep = ''), temp, envir = parent.frame()); gc()
    
    temp <- get(paste('twin', i, 'R', sep = '')); gc()
    temp <- temp[abs(temp$meth_fw - temp$meth_rv) < 20,]
    temp <- temp[temp$num_total_fw > 1 & temp$num_total_rv > 1, ]
    temp$ID <- paste(paste(temp$`#chr`, temp$start, sep = ':'), temp$end, sep = '-')
    assign(paste('twin', i, 'R', sep = ''), temp,  envir = parent.frame()); gc()
    
    print(i)
  }
}
process.files2 <- function()
{
  files <- list.files(pattern = '*.txt')
  
  for(i in 1:(length(files)/2))
  {
    
    # Twin L
    temp <- get(paste('twin', i, 'L', sep = ''))
    chr <- temp$`#chr`
    start <- temp$start
    IDs <- temp$ID
    eliminate <- character()
    
    
    eliminate <- lapply(1:nrow(blacklist), function(j)
    {
      chri <- blacklist[[j, 1]]
      starti <- blacklist[[j, 2]]
      endi <- blacklist[[j, 3]]
      
      ID_j <- IDs[chr == chri]
      start_j <- start[chr == chri]
      ID_j[which(start_j > starti & start_j < endi)]
    })
    
    eliminate <- unique(unlist(eliminate))
    
    temp <- temp[-match(eliminate, IDs),]
    assign(paste('twin', i, 'L', sep = ''), temp, envir = parent.frame())
    
    
    # Twin R
    temp <- get(paste('twin', i, 'R', sep = ''))
    chr <- temp$`#chr`
    start <- temp$start
    IDs <- temp$ID
    eliminate <- character()
    eliminate <- lapply(1:nrow(blacklist), function(j)
    {
      chri <- blacklist[[j, 1]]
      starti <- blacklist[[j, 2]]
      endi <- blacklist[[j, 3]]
      
      ID_j <- IDs[chr == chri]
      start_j <- start[chr == chri]
      ID_j[which(start_j > starti & start_j < endi)]
    })
    
    eliminate <- unique(unlist(eliminate))
    
    temp <- temp[-match(eliminate, IDs),]
    assign(paste('twin', i, 'R', sep = ''), temp, envir = parent.frame())
    
    print(i); gc()
  }
}
process.files3 <- function()
{
  l.files <- 7
  
  low_filt = 10
  
  for(i in 1:l.files)
  {
    temp <- get(paste('twin', i, 'L', sep = ''))
    high_filt <- quantile(temp$total, 0.999)
    temp <- temp[temp$total < high_filt & temp$total > low_filt, ]
    assign(paste('twin', i, 'L', sep = ''), temp, envir = parent.frame())
    
    
    temp <- get(paste('twin', i, 'R', sep = ''))
    high_filt <- quantile(temp$total, 0.999)
    temp <- temp[temp$total < high_filt & temp$total > low_filt, ]
    assign(paste('twin', i, 'R', sep = ''), temp,  envir = parent.frame()); gc()
    
    print(i)
  }
}
process.common.twin <- function(twinL, twinR, common)
{
  whereL <- twinL$ID %chin% common
  whereR <- twinR$ID %chin% common
  
  assign(x = deparse(substitute(twinL)), value = twinL[whereL,], envir = parent.frame())
  assign(x = deparse(substitute(twinR)), value = twinR[whereR,], envir = parent.frame())
}

# Writing files
write.files <- function()
{
  files <- 1:14
  
  for(i in 1:(length(files)/2))
  {
    fwrite(x = get(paste('twin', i, 'L', sep = '')), file = paste('/media/ultron/2tb_disk2/PROCESSED_DATA/2018/Twin_project/WGBS/adipose/', 'twin', i, 'L', '.txt',sep = ''), nThread = 4)
    fwrite(x = get(paste('twin', i, 'R', sep = '')), file = paste('/media/ultron/2tb_disk2/PROCESSED_DATA/2018/Twin_project/WGBS/adipose/', 'twin', i, 'R', '.txt',sep = ''), nThread = 4)
    gc()
    print(i)
  }
}
write.files2 <- function()
{
  l.files <- 7
  
  for(i in 1:l.files)
  {
    fwrite(x = get(paste('twin', i, 'L', sep = '')), file = paste('/media/ultron/2tb_disk2/PROCESSED_DATA/2018/Twin_project/WGBS/proc2/', 'twin', i, 'L', '.txt',sep = ''), nThread = 4)
    fwrite(x = get(paste('twin', i, 'R', sep = '')), file = paste('/media/ultron/2tb_disk2/PROCESSED_DATA/2018/Twin_project/WGBS/proc2/', 'twin', i, 'R', '.txt',sep = ''), nThread = 4)
    gc()
    print(i)
  }
}

# Finding common positions for a given twin pair
getCommon2 <- function(i)
{
  survival <- intersect(get(paste('twin', i, 'L', sep = ''))$ID, get(paste('twin', i, 'R', sep = ''))$ID)
  
  return(survival)
}

# Untargetted positional enrichment analysis
targetted_positional_enrichment_analysis <- function(twinL, twinR, epsilon, window, chr_i)
{
  # Check everything is prepared
  if(nrow(twinL) != nrow(twinR))
  {
    message('different number of row')
    stop()
  }
  if(sum(twinL$ID == twinR$ID) != nrow(twinR))
  {
    message('not ordered')
    stop()
  }
  
  # Prepare IDs/chr/start
  IDs <- twinL$ID
  chr <- twinL$`#chr`
  start <- twinL$start
  
  names(chr) <- IDs
  names(start) <- IDs
  
  # Obtain target from epsilon threshold
  delta_beta <- abs(twinL$total_meth - twinR$total_meth)/100
  names(delta_beta) <- IDs
  target <- names(which(delta_beta >= epsilon))
  bg <- IDs[!(IDs %in% target)] # Bg and target must be mutually exclusive
  
  
  #
  chr_bg <- chr[bg]
  start_bg <- start[bg]
  chr_target <- chr[target]
  start_target <- start[target]
  
  
  
  # Subselect chromosome
  BG <- data.table(chr_ = chr_bg[chr_bg == chr_i], start = start_bg[chr_bg == chr_i], IDs = bg[chr_bg == chr_i])
  TG <- data.table(chr_ = chr_target[chr_target == chr_i], start = start_target[chr_target == chr_i], IDs = target[chr_target == chr_i])
  a = sum(TG$start > window[1] & TG$start < window[2])
  b = sum(BG$start > window[1] & BG$start < window[2])
  
  
  
  contrast <- matrix(c(a, length(target) - a, b, length(bg) - b), nrow = 2);
  print(contrast)
  
  res <- fisher.test(contrast)
  print(res)
  
  return(list(BG = bg, TG = target))
}



# setwd("where")
phenotype <- fread('E-MTAB-3549.sdrf.txt')
phenotype <- phenotype[,c(1,3,4,5,6,8,9,32)]


phenotype <- phenotype[phenotype$`Characteristics[zygosity]` == "monozygotic" & phenotype$`Characteristics[organism part]` == "adipose tissue of abdominal region"]

phenotype
table(phenotype$`Characteristics[zygosity]`) # 32 DZ; 29 MZ
table(phenotype$`Characteristics[sex]`) # 61 Females
table(phenotype$`Characteristics[age]`) # From 39 to 76
table(phenotype$`Characteristics[organism part]`) # 27 blood and 34 adipose







# Read raw files
# setwd("where")
list.files()
read.files()

dim(twin1L) # 25302909  
dim(twin1R) # 24275124
dim(twin2L) # 26026365
dim(twin2R) # 26339470
dim(twin3L) # 26312125
dim(twin3R) # 26319494
dim(twin4L) # 25274262 
dim(twin4R) # 24992386
dim(twin5L) # 25410400 
dim(twin5R) # 25470165 
dim(twin6L) # 26040048
dim(twin6R) # 26047078
dim(twin7L) # 24592446
dim(twin7R) # 23370518


# Processing 1 consists of: 
# i) Selecting positions with less than 20 % methylation difference between strands. 
# ii) Filtering positions with coverage less than two reads per strand
# iii) add position ID column


# setwd("where")
process.files(); gc() # Takes long!

dim(twin1L) # 7081324
dim(twin1R) # 4482956
dim(twin2L) # 11026329 
dim(twin2R) # 13977644
dim(twin3L) # 13198116
dim(twin3R) # 13719458
dim(twin4L) # 6780503
dim(twin4R) # 6670298
dim(twin5L) # 8199686
dim(twin5R) # 8276537
dim(twin6L) # 11010806
dim(twin6R) # 11209299
dim(twin7L) # 4770622
dim(twin7R) # 3580812

# Export files
write.files()


##################################

# Proc 2

# Read blacklisted regions
# setwd("where")
DBR <- fread('consensusBlacklist.bed')
DER <- fread('dukeExcludeRegions.bed')
blacklist <- rbind(DBR, DER)
dim(blacklist) # 2060 6

# Processing 2 consists of: 
# i) Excluding DBR and DER blacklisted regions

# setwd("where")
list.files(pattern = '*.txt')
read.files2()


# setwd("where")
process.files2()

dim(twin1L) # 7060553
dim(twin1R) # 4469168
dim(twin2L) # 10993977
dim(twin2R) # 13935434
dim(twin3L) # 13158240  
dim(twin3R) # 13677467
dim(twin4L) # 6759070  
dim(twin4R) # 6647607
dim(twin5L) # 8172020
dim(twin5R) # 8250357
dim(twin6L) # 10977280
dim(twin6R) # 11177342
dim(twin7L) # 4755042 
dim(twin7R) # 3567144 

# Export files
write.files()


#######################################


# PROC3
# setwd("where")
list.files(pattern = '*.txt') # character(0)
read.files2()

# Processing 4 consists of: 
# i) Applying low end coverage filter: more than 10 reads per position
# ii) Applying high end coverage filter: per sample, filter positions with more than the 99.9 quantile of coverage
process.files3()

dim(twin1L) # 1262361
dim(twin1R) # 273814
dim(twin2L) # 4444770
dim(twin2R) # 8348458
dim(twin3L) # 7504753
dim(twin3R) # 7721147
dim(twin4L) # 936950
dim(twin4R) # 1694686
dim(twin5L) # 2264448
dim(twin5R) # 2130354
dim(twin6L) # 6721038
dim(twin6R) # 6140012
dim(twin7L) # 774091
dim(twin7R) # 490061



# Export files
write.files()



##################################################




# Find common positions for each twin
common1 <- getCommon2(1); length(common1) # 73977
common2 <- getCommon2(2); length(common2) # 2638201
common3 <- getCommon2(3); length(common3) # 3843141
common4 <- getCommon2(4); length(common4) # 269438
common5 <- getCommon2(5); length(common5) # 700185
common6 <- getCommon2(6); length(common6) # 3415115
common7 <- getCommon2(7); length(common7) # 124649

# Select common positions between twins
process.common.twin(twin1L, twin1R, common1)
process.common.twin(twin2L, twin2R, common2)
process.common.twin(twin3L, twin3R, common3)
process.common.twin(twin4L, twin4R, common4)
process.common.twin(twin5L, twin5R, common5)
process.common.twin(twin6L, twin6R, common6)
process.common.twin(twin7L, twin7R, common7)

# Export files
write.files()


#############################################   Data analysis   #############################################


# Read common positions
# setwd("where")
list.files(pattern = '*.txt')
read.files2()

# Visualize data
library(scales)
library(RColorBrewer)

coverage <- c(6.12, 7.42, 28.96, 19.08, 19.22, 18.6, 7.58, 7.92, 7.84, 7.94, 0.7, 6.02, 7.4, 6.56)
nfeat <- c(199573, 199573,17067861, 17067861,
           13786481, 13786481, 464975, 464975,
           610434, 610434, 0, 0,
           350716, 350716)
barplot(nfeat, col = rep(brewer.pal(7, "Dark2"), each = 2))
plot(coverage, nfeat, xlab= 'Sequencing coverage', ylab = 'number of CpGs covered (after filters)',
     pch = 19, col = alpha(rep(brewer.pal(7, "Dark2"), each = 2), 0.8), cex = 2)



coverage <- c(6.66, 4.68, 9.66, 12.92, 12.62, 12.38, 6.22, 7.3, 7.9, 7.66, 5.6, 5.04, 5.86, 5.12)
nfeat <- c(73977, 73977,2638201, 2638201,
           3843141, 3843141, 269438, 269438,
           700185, 700185, 3415115, 3415115,
           124649, 124649)

col =  rep(c(brewer.pal(5, "Dark2"), "darkgreen", "gold"), each = 2)

barplot(nfeat, col = col)
plot(coverage, nfeat, xlab= 'Sequencing coverage', ylab = 'number of CpGs covered (after filters)',
     pch = 19, col = alpha(col, 0.8), cex = 2)



#####################################################################################33


# Targetted positional enrichment analysis
PCDHA = c(140165876, 140391929)
PCDHB = c(140430979, 140627802)
PCDHG = c(140710252, 140892546)
window = c(PCDHA[1], PCDHG[2])



res = targetted_positional_enrichment_analysis(twin1L, twin1R, 0.4, window, 'chr5')
# [,1]  [,2]
# [1,]    0    20
# [2,]  666 73291
# 
# Fisher's Exact Test for Count Data
# 
# data:  contrast
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.00000 22.34777
# sample estimates:
# odds ratio 
#          0 

res = targetted_positional_enrichment_analysis(twin2L, twin2R, 0.4, window, 'chr5')
# [,1]    [,2]
# [1,]    14     680
# [2,] 17699 2619808
# 
# Fisher's Exact Test for Count Data
# 
# data:  contrast
# p-value = 0.0003332
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.656862 5.150141
# sample estimates:
# odds ratio 
#   3.047394 

res = targetted_positional_enrichment_analysis(twin3L, twin3R, 0.4, window, 'chr5')
#    [,1]    [,2]
# [1,]    28    1045
# [2,] 23238 3818830
# 

# Fisher's Exact Test for Count Data
# 
# data:  contrast
# p-value = 2.896e-10
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.910846 6.403084
# sample estimates:
# odds ratio 
#   4.403432 

res = targetted_positional_enrichment_analysis(twin4L, twin4R, 0.4, window, 'chr5')
# [,1]   [,2]
# [1,]    1     63
# [2,] 2130 267244
# 
# Fisher's Exact Test for Count Data
# 
# data:  contrast
# p-value = 0.3985
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.04962026 11.52273199
# sample estimates:
# odds ratio 
#   1.991527 


res = targetted_positional_enrichment_analysis(twin5L, twin5R, 0.4, window, 'chr5')
# [,1]   [,2]
# [1,]    3    162
# [2,] 4573 695447
# 
# Fisher's Exact Test for Count Data
# 
# data:  contrast
# p-value = 0.09464
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.5745811 8.3847956
# sample estimates:
# odds ratio 
#   2.816172 


res = targetted_positional_enrichment_analysis(twin6L, twin6R, 0.4, window, 'chr5')
# [,1]    [,2]
# [1,]    11     813
# [2,] 14769 3399522
# 
# Fisher's Exact Test for Count Data
# 
# data:  contrast
# p-value = 0.001143
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.548103 5.603994
# sample estimates:
# odds ratio 
#   3.114295 


res = targetted_positional_enrichment_analysis(twin7L, twin7R, 0.4, window, 'chr5')
# [,1]   [,2]
# [1,]    0     14
# [2,]  917 123718
# 
# Fisher's Exact Test for Count Data
# 
# data:  contrast
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.00000 40.72197
# sample estimates:
# odds ratio 
#          0 



########




# setwd("where")
tiff(filename = "adipose_blood_wgbs.tiff", width = 6, height = 4, units = "in", res = 300)
par(mar=c(2,2,1,1))
CEX = 1.4
cov <- c(1202, 11200, 17974, 2676, 3557, NA, 1219)
sig <- as.factor(c(F, T, T, T, F, NA, F))
levels(sig) = c("firebrick3", "forestgreen")
sig = as.character(sig)
sig[is.na(sig)] = "black"

set.seed(9)
x = rnorm(n = length(cov), mean = 0, sd = 0.07)
plot(cov, x + rep(1 , times = length(cov)), col = as.character(sig),
     xlim = c(0,25000), ylim = c(0.5, 2.5), yaxt = "n", xlab = "Number of Cs covered in the cPCDH region", ylab = "",
     pch = as.character(c(1:5, 8, 9)), cex = CEX, axes = F)
axis(side = 2, at = c(1, 2), labels = c("Blood", "Adipose"))
axis(side = 1, at = seq(0, 25000, 5000))

#coverage <- c(6.66, 4.68, 9.66, 12.92, 12.62, 12.38, 6.22, 7.3, 7.9, 7.66, 5.6, 5.04, 5.86, 5.12)
#cov <- sapply(seq(1, length(coverage), 2), function(x) min(coverage[x:(x+1)]))
cov <- c(666, 17713, 23266, 2131, 4576, 14780, 917)
sig <- as.factor(c(F, T, T, F, F, T, F))
#levels(sig) = c("red3", "green3")
levels(sig) = c("firebrick3", "forestgreen")
set.seed(2)
x = rnorm(n = length(cov), mean = 0, sd = 0.07)
points(cov, x + rep(2, times = length(cov)), pch = as.character(1:7), col = as.character(sig),
     xlim = c(0,20), cex = CEX)
dev.off()

#

coverage <- c(6.12, 7.42, 28.96, 19.08, 19.22, 18.6, 7.58, 7.92, 7.84, 7.94, 0.7, 6.02, 7.4, 6.56,
              6.66, 4.68, 9.66, 12.92, 12.62, 12.38, 6.22, 7.3, 7.9, 7.66, 5.6, 5.04, 5.86, 5.12)
cov <- sapply(seq(1, length(coverage), 2), function(x) min(coverage[x:(x+1)]))
cov2 <- c(1202, 11200, 17974, 2676, 3557, NA, 1219,
          666, 17713, 23266, 2131, 4576, 14780, 917)
plot(cov, cov2)





library(data.table)
# setwd("where")
data = fread("E-MTAB-3549.sdrf.txt")
head(data)
adipose = data[data$`Characteristics[zygosity]` == "monozygotic" & data$`Characteristics[organism part]` == "adipose tissue of abdominal region",]
blood = data[data$`Characteristics[zygosity]` == "monozygotic" & data$`Characteristics[organism part]` == "blood",]
adipose = adipose[!(adipose$`Characteristics[co-twin]` == "not applicable"),]
blood = blood[!(blood$`Characteristics[co-twin]` == "not applicable"),]
adipose = adipose[order(adipose$`Characteristics[co-twin]`),]
blood = blood[order(blood$`Characteristics[co-twin]`),]

adipose[, c(3:4, 6)]
blood[, c(3:4, 6)]

















