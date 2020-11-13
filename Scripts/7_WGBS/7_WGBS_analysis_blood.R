############################################################################
############################################################################
###########                                                      ###########
###########               WGBS analysis blood                    ###########
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
    fwrite(x = get(paste('twin', i, 'L', sep = '')), file = paste('/media/ultron/2tb_disk2/PROCESSED_DATA/2018/Twin_project/WGBS/', 'twin', i, 'L', '.txt',sep = ''), nThread = 4)
    fwrite(x = get(paste('twin', i, 'R', sep = '')), file = paste('/media/ultron/2tb_disk2/PROCESSED_DATA/2018/Twin_project/WGBS/', 'twin', i, 'R', '.txt',sep = ''), nThread = 4)
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


#############################################   Prepare phenotypes   #############################################

# setwd("where")
phenotype <- fread('E-MTAB-3549.sdrf.txt')
phenotype <- phenotype[,c(1,3,4,5,6,8,9,32)]

table(phenotype$`Characteristics[zygosity]`) # 32 DZ; 29 MZ
table(phenotype$`Characteristics[sex]`) # 61 Females
table(phenotype$`Characteristics[age]`) # From 39 to 76
table(phenotype$`Characteristics[organism part]`) # 27 blood and 34 adipose

# Extract Blood from MZ
phenotype <- phenotype[phenotype$`Characteristics[zygosity]` == 'monozygotic' &
                         phenotype$`Characteristics[organism part]` == 'blood',]
dim(phenotype) # 7 MZ twin pairs are available for whole-blood

# Write bash script to move samples of interest to a folder called MZ_twins in current directory
write.table(paste('mv', phenotype$`Derived Array Data File`, './MZ_twins'), quote = F, row.names = F, col.names = F)

#############################################   Pre-processing 1   #############################################

# Read raw files
# setwd("where")
read.files()

dim(twin1L) # 25,505,737 positions
dim(twin1R) # 26,013,357 positions
dim(twin2L) # 26,804,319 positions
dim(twin2R) # 26,686,723 positions
dim(twin3L) # 26,705,981 positions
dim(twin3R) # 26,673,411 positions
dim(twin4L) # 26,163,271 positions
dim(twin4R) # 25,665,732 positions
dim(twin5L) # 26,160,875 positions
dim(twin5R) # 25,868,966 positions
dim(twin6L) # 12,142,848 positions *** Lower coverage
dim(twin6R) # 25,686,220 positions
dim(twin7L) # 26,164,359 positions
dim(twin7R) # 26,167,333 positions

# Processing 1 consists of: 
# i) Selecting positions with less than 20 % methylation difference between strands. 
# ii) Filtering positions with coverage less than two reads per strand
# iii) add position ID column

# setwd("where")
process.files(); gc() # Takes long!

dim(twin1L) # 7,176,329 positions
dim(twin1R) # 9,822,360 positions
dim(twin2L) # 22,198,454 positions
dim(twin2R) # 20,478,731 positions
dim(twin3L) # 20,194,175 positions
dim(twin3R) # 19,252,913 positions
dim(twin4L) # 10,394,267 positions
dim(twin4R) # 8,282,004 positions
dim(twin5L) # 10,435,814 positions
dim(twin5R) # 9,337,371 positions
dim(twin6L) # 40,182 positions *** Lower coverage
dim(twin6R) # 6,930,329 positions
dim(twin7L) # 10,503,618 positions
dim(twin7R) # 9,897,931 positions

# Export files
write.files()

#############################################   Pre-processing 2   #############################################

# Read process1 files
# setwd("where")
read.files2()

# Read blacklisted regions
# setwd("where")
DBR <- fread('consensusBlacklist.bed')
DER <- fread('dukeExcludeRegions.bed')
blacklist <- rbind(DBR, DER)
dim(blacklist) # 2060 6

# Processing 2 consists of: 
# i) Excluding DBR and DER blacklisted regions

# setwd("where")
process.files2()

dim(twin1L) # 7155734 positions
dim(twin1R) # 9794859 positions
dim(twin2L) # 22131827 positions
dim(twin2R) # 20418513 positions
dim(twin3L) # 20135549 positions
dim(twin3R) # 19196071 positions
dim(twin4L) # 10363754 positions
dim(twin4R) # 8255881 positions
dim(twin5L) # 10403335 positions
dim(twin5R) # 9308105 positions
dim(twin6L) # 39127 *** Lower coverage
dim(twin6R) # 6908656 positions
dim(twin7L) # 10472334 positions
dim(twin7R) # 9868506 positions

# Export files
# setwd("where")
list.files(pattern = '*.txt') # character(0)
write.files()

#############################################   Pre-processing 3   #############################################

# Read process2 files
# setwd("where")
list.files(pattern = '*.txt')
read.files2()

# Processing 4 consists of: 
# i) Applying low end coverage filter: more than 10 reads per position
# ii) Applying high end coverage filter: per sample, filter positions with more than the 99.9 quantile of coverage
process.files3()

dim(twin1L) # 860698 positions
dim(twin1R) # 2123108 positions
dim(twin2L) # 21605001 positions
dim(twin2R) # 19247617 positions
dim(twin3L) # 18435192 positions
dim(twin3R) # 17447081 positions
dim(twin4L) # 2516881 positions
dim(twin4R) # 2160317 positions
dim(twin5L) # 2715423 positions
dim(twin5R) # 2418088 positions
dim(twin6L) # 0 *** Lower coverage
dim(twin6R) # 756554 positions
dim(twin7L) # 2363754 positions
dim(twin7R) # 2025865 positions

# Export files
# setwd("where")
list.files(pattern = '*.txt') # character(0)
write.files()

#############################################   Common for each twin pair   #############################################

# Find common positions for each twin
common1 <- getCommon2(1); length(common1) # 199,573
common2 <- getCommon2(2); length(common2) # 17,067,861
common3 <- getCommon2(3); length(common3) # 13,786,481
common4 <- getCommon2(4); length(common4) # 464,975
common5 <- getCommon2(5); length(common5) # 610,434
common6 <- getCommon2(6); length(common6) # 0
common7 <- getCommon2(7); length(common7) # 350,716

# Select common positions between twins
process.common.twin(twin1L, twin1R, common1)
process.common.twin(twin2L, twin2R, common2)
process.common.twin(twin3L, twin3R, common3)
process.common.twin(twin4L, twin4R, common4)
process.common.twin(twin5L, twin5R, common5)
process.common.twin(twin6L, twin6R, common6)
process.common.twin(twin7L, twin7R, common7)

# Export files
# setwd("where")
list.files(pattern = '*.txt') # character(0)
write.files()


#############################################   Simulation to decided delta_beta threshold   #############################################

# At least 10X means that for CpG with beta = 0.5, the expected sd will be:
sqrt(0.5*(1-0.5)/10) # 0.1581139

# This can be simulated:

subset <- lapply(1:10000, function(x) sample(c(0,1), 10, replace = T))
#subset <- lapply(1:10000, function(x) rbinom(n = 10, size = 1, prob = 0.5))

sampled.betas <- sapply(subset, mean)
sd(sampled.betas) # 0.1551832

# As well, we can simulate the expected difference distribution
sampled.diff <- unlist(lapply(1:(length(sampled.betas)-1), function(x) 
{
  sampled.betas[(x+1):length(sampled.betas)] - sampled.betas[x]
}))
length(sampled.diff) # 49995000
choose(length(sampled.betas), 2) # 49995000

# Visualize distribution delta beta
hist(sampled.diff)

# Visualize distribution |delta beta|
hist(abs(sampled.diff)); median(abs(sampled.diff)); abline(v = median(abs(sampled.diff)), col = 'blue3') # 0.1
quantile(abs(sampled.diff), 0.95) # 0.4
# Selecting the target set with a difference threshold of 0.4 will help filter out
# most of the differences that arose via sampling bias.

#############################################   Data analysis   #############################################

# Read common positions
# setwd("where")
list.files(pattern = '*.txt')
read.files2()

# Visualize data
coverage <- c(6.12, 7.42, 28.96, 19.08, 19.22, 18.6, 7.58, 7.92, 7.84, 7.94, 0.7, 6.02, 7.4, 6.56)
nfeat <- c(nrow(twin1L), nrow(twin1R), nrow(twin2L), nrow(twin2R),
           nrow(twin3L), nrow(twin3R), nrow(twin4L), nrow(twin4R),
           nrow(twin5L), nrow(twin5R), nrow(twin6L), nrow(twin6R),
           nrow(twin7L), nrow(twin7R))

plot(coverage, nfeat, xlab= 'Sequencing coverage', ylab = 'number of CpGs covered (after filters)',
     pch = 19, col = alpha(rep(brewer.pal(7, "Dark2"), each = 2), 0.8), cex = 2)
barplot(nfeat, col = rep(brewer.pal(7, "Dark2"), each = 2))

hist(twin2L$total, breaks = 100, col = brewer.pal(7, "Dark2")[2], right = F)
hist(twin2L$total_meth/100, breaks = 100, col = brewer.pal(7, "Dark2")[2], right = F)
hist(twin2R$total, breaks = 100, col = brewer.pal(7, "Dark2")[2], right = F)
hist(twin2R$total_meth/100, breaks = 100, col = brewer.pal(7, "Dark2")[2], right = F)

hist(twin3L$total, breaks = 100, col = brewer.pal(7, "Dark2")[3], right = F)
hist(twin3L$total_meth/100, breaks = 100, col = brewer.pal(7, "Dark2")[3], right = F)
hist(twin3R$total, breaks = 100, col = brewer.pal(7, "Dark2")[3], right = F)
hist(twin3R$total_meth/100, breaks = 100, col = brewer.pal(7, "Dark2")[3], right = F)

# Targetted positional enrichment analysis
PCDHA = c(140165876, 140391929)
PCDHB = c(140430979, 140627802)
PCDHG = c(140710252, 140892546)
window = c(PCDHA[1], PCDHG[2])

### High coverage twins ###
# Twin pair 2

targetted_positional_enrichment_analysis(twin2L, twin2R, 0.4, window, 'chr5')
#       [,1]     [,2]
#[1,]    28     6687
#[2,] 11172 17049974
#Fisher's Exact Test for Count Data
#data:  contrast
#p-value = 4.785e-14
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 4.241085 9.249179
#sample estimates:
#odds ratio 
#  6.390182 
100*28/(28+11172) # 0.25 %
mat <- matrix(c(28, 6687, 11172, 17049974), byrow = T, nrow = 2)
fisher.test(mat)
mat[1,] <- mat[1,]/sum(mat[1,])
mat[2,] <- mat[2,]/sum(mat[2,])
rownames(mat) <- c('Bg', 'Target')
rownames(mat) <- c('within PCDH', 'outside')
balloonplot(as.table(mat), main = '', ylab = 'Island Status', xlab = 'Set')

# Twin pair 3

targetted_positional_enrichment_analysis(twin3L, twin3R, 0.4, window, 'chr5')
#      [,1]     [,2]
#[1,]    46     4628
#[2,] 17928 13763879
#Fisher's Exact Test for Count Data
#data:  contrast
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  5.576757 10.200848
#sample estimates:
#odds ratio 
#  7.630906 
100*46/(46+17928) # 0.25 %
mat <- matrix(c(46, 4628, 17928, 13763879), byrow = T, nrow = 2)
fisher.test(mat)
mat[1,] <- mat[1,]/sum(mat[1,])
mat[2,] <- mat[2,]/sum(mat[2,])
rownames(mat) <- c('Bg', 'Target')
rownames(mat) <- c('within PCDH', 'outside')
balloonplot(as.table(mat), main = '', ylab = 'Island Status', xlab = 'Set')

### Lower coverage twins ###

# Twin pair 1
targetted_positional_enrichment_analysis(twin1L, twin1R, 0.4, window, 'chr5')
#     [,1]   [,2]
#[1,]    1     36
#[2,] 1201 198335
# Fisher's Exact Test for Count Data
# data:  contrast
# p-value = 0.2003
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.1129464 27.3088075
# sample estimates:
# odds ratio 
#   4.587163 

# Twin pair 4
targetted_positional_enrichment_analysis(twin4L, twin4R, 0.4, window, 'chr5')
#     [,1]   [,2]
#[1,]    7    107
#[2,] 2669 462192
#Fisher's Exact Test for Count Data
# data:  contrast
# p-value = 4.993e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   4.443767 24.195065
# sample estimates:
# odds ratio 
#   11.32856 

# Twin pair 5
targetted_positional_enrichment_analysis(twin5L, twin5R, 0.4, window, 'chr5')
#      [,1]   [,2]
# [1,]    0    153
# [2,] 3557 606724

# Fisher's Exact Test for Count Data

# data:  contrast
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.000000 4.164648
# sample estimates:
# odds ratio 
#         0 

# Twin pair 7
targetted_positional_enrichment_analysis(twin7L, twin7R, 0.4, window, 'chr5')
#      [,1]   [,2]
# [1,]    0     73
# [2,] 1219 349424
# Fisher's Exact Test for Count Data
# data:  contrast
# p-value = 1
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.00000 14.87632
# sample estimates:
# odds ratio 
#          0 

#############################################   Export files for visualizing in IGV   #############################################

# Perform targetted positional enrichment analysis, now exporting features
list_1 <- targetted_positional_enrichment_analysis(twin2L, twin2R, 0.4, window, 'chr5')
list_2 <- targetted_positional_enrichment_analysis(twin3L, twin3R, 0.4, window, 'chr5')


a = list(BG = c("chr1:122-143", "chr1:122-143", "chr1:122-143"), TG = c("chr1:122-143","chr1:122-143","chr1:122-143"))

prepare.bed <- function(list_tg_bg)
{
  bg = list_tg_bg[[1]]
  tg = list_tg_bg[[2]]
  
  b_chr = unlist(lapply(strsplit(bg, split = ":"), function(x) x[1]))
  b_res = unlist(lapply(strsplit(bg, split = ":"), function(x) x[2]))
  t_chr = unlist(lapply(strsplit(tg, split = ":"), function(x) x[1]))
  t_res = unlist(lapply(strsplit(tg, split = ":"), function(x) x[2]))
  
  b_start = as.integer(unlist(lapply(strsplit(b_res, split = "-"), function(x) x[1])))
  b_end = as.integer(unlist(lapply(strsplit(b_res, split = "-"), function(x) x[2])))
  t_start = unlist(lapply(strsplit(t_res, split = "-"), function(x) x[1]))
  t_end = unlist(lapply(strsplit(t_res, split = "-"), function(x) x[2]))
  
  bed1 = data.frame(chr = b_chr, start = b_start, end = b_end, ID = bg, value = 1000)
  bed2 = data.frame(chr = t_chr, start = t_start, end = t_end, ID = bg, value = 1000)
  
  return(list(bed1, bed2))
  
}

prepare.bed(a)

# Prepare bed-like file
bed_1 <- prepare.bed(list_1)
bed_2 <- prepare.bed(list_2)

# Initialize files with header
# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1', file = 'wgbs_target_twin2.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1', file = 'wgbs_target_twin3.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)

# Write bed
fwrite(x = bed_1[[2]], file = 'wgbs_target_twin2.bed', nThread = 4, sep = '\t', col.names = F, append = T)
fwrite(x = bed_2[[2]], file = 'wgbs_target_twin3.bed', nThread = 4, sep = '\t', col.names = F, append = T)

# As Bgs are too big for IGV, we subselected sites at Chr5
bg1 <- bed_1[[1]]
bg1 <- bg1[bg1$seqname == 'chr5',]
bg2 <- bed_2[[1]]
bg2 <- bg2[bg2$seqname == 'chr5',]

# As blocks are not very informative, line plot of normalized coverage were prefered
bins <- seq(from = 1, to = 180915260, by = 1000) # http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
length(bins) # 180916

# Twin 2
counts <- binCounts(x = bg1$start, bx = bins)
sum(counts)
coverage_track <- data.frame(seqname = 'chr5',
                         start = as.integer(bins[-length(bins)]), end = as.integer(bins[2:length(bins)] - 1),
                         name = paste('chr5', paste(bins[-length(bins)], bins[2:length(bins)] - 1, sep = '-'), sep = ':'),
                         score = 1000*counts/max(counts))
# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1 graphType=line viewLimits=0:1000', file = 'coverage_bg_withouttarget_twin2_chr5.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = coverage_track, file = 'coverage_bg_withouttarget_twin2_chr5.bed', nThread = 4, sep = '\t', col.names = F, append = T)

# Twin 3
counts <- binCounts(x = bg2$start, bx = bins)
coverage_track <- data.frame(seqname = 'chr5',
                             start = as.integer(bins[-length(bins)]), end = as.integer(bins[2:length(bins)] - 1),
                             name = paste('chr5', paste(bins[-length(bins)], bins[2:length(bins)] - 1, sep = '-'), sep = ':'),
                             score = 1000*counts/max(counts))
# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1 graphType=line viewLimits=0:1000', file = 'coverage_bg_withouttarget_twin3_chr5.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = coverage_track, file = 'coverage_bg_withouttarget_twin3_chr5.bed', nThread = 4, sep = '\t', col.names = F, append = T)


# Bg and target included together
# Twin 2
counts <- binCounts(x = c(bg1$start, tg1$start), bx = bins)
coverage_track <- data.frame(seqname = 'chr5',
                             start = as.integer(bins[-length(bins)]), end = as.integer(bins[2:length(bins)] - 1),
                             name = paste('chr5', paste(bins[-length(bins)], bins[2:length(bins)] - 1, sep = '-'), sep = ':'),
                             score = 1000*counts/max(counts))
# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1 graphType=line viewLimits=0:1000', file = 'coverage_bg_twin2_chr5.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = coverage_track, file = 'coverage_bg_twin2_chr5.bed', nThread = 4, sep = '\t', col.names = F, append = T)

# Twin 3
counts <- binCounts(x = c(bg2$start, tg2$start), bx = bins)
coverage_track <- data.frame(seqname = 'chr5',
                             start = as.integer(bins[-length(bins)]), end = as.integer(bins[2:length(bins)] - 1),
                             name = paste('chr5', paste(bins[-length(bins)], bins[2:length(bins)] - 1, sep = '-'), sep = ':'),
                             score = 1000*counts/max(counts))
# setwd("where")
write.table(x = 'track name=pairedReads description=Clone Paired Reads useScore=1 graphType=line viewLimits=0:1000', file = 'coverage_bg_twin3_chr5.bed',
            quote = F, sep = '\t', row.names = F, col.names = F)
fwrite(x = coverage_track, file = 'coverage_bg_twin3_chr5.bed', nThread = 4, sep = '\t', col.names = F, append = T)


