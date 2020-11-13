############################################################################
############################################################################
###########                                                      ###########
###########             DNA motif enrichment analysis            ###########
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

library(minfi)
library(data.table)
library(seqinr)

## Load functions ##

# Help parse coordinates for input into samtools
parse.coordinates.out <- function(annotation, distance, name)
{
  seqDA <- as.integer(annotation$pos - distance)
  seqDB <- as.integer(annotation$pos + distance)
  seqD <- paste(seqDA, seqDB, sep = '-')
  seqD <- paste(as.vector(annotation$chr), seqD, sep = ':')
  fwrite(x = as.data.frame(seqD), file = paste(name, distance, sep = ''), quote = F, row.names = F, col.names = F)
}

#############################################   Parse coordinates to extract   #############################################

old.dir <- getwd()
# setwd("where")
RGSET <- read.metharray.exp(getwd())
full_annot <- getAnnotation(RGSET)
annotation <- functional.annotation.hg19(stochCpG)

# Read raw data
# setwd("where")
beta1 <- fread('2019-08-21_SQN_combat_cellcomp.txt', nThread = 4)
beta1 <- process.beta.fread(beta1); dim(beta1) # 346555    852

full_annot <- full_annot[rownames(beta1),]

# setwd("where")

# Target
parse.coordinates.out(annotation, 500, 'evCpGs')

# Bg
parse.coordinates.out(full_annot, 500, 'bg')

#############################################   Motif enrichment analysis   #############################################

## To extract sequences for target/Bg we employed samtools:

# samtools faidx hg19.fasta.gz -r evCpGs500 > evCpGs500.fa
# samtools faidx hg19.fasta.gz -r bg500 > bg500.fa

## To validate the the right number of entries exist:

# tr -cd '>' < evCpGs500.fa | wc -c
# 333
# tr -cd '>' < bg500.fa | wc -c
# 485512

## We then performed motif enrichment analysis with Homer by calling:
# findMotifs.pl evCpGs500.fa fasta . -fasta bg500.fa -p 4 -humanGO > log.txt


#############################################   GC content   #############################################


# setwd("where")

target <- read.fasta('stochCpGs500.fa', seqtype = 'DNA')
length(target) # 333

bg <- read.fasta('450K_500.fa', seqtype = 'DNA')
length(bg) # 346555

gc_target <- sapply(1:length(target), function(x) mean(target[[x]] %in% c('g', 'c')))
gc_bg <- sapply(1:length(bg), function(x) mean(bg[[x]] %in% c('g', 'c')))

# setwd("where")
tiff(filename = paste('GC', 'tiff', sep = '.'), width = 10, height = 10, units = 'in', res = 300, compression = 'none')
plot(density(gc_bg), xlim = c(0,1), ylim = c(0,4.5), lty = 2, col = 'green4', lwd = 3,
     xlab = 'GC content', main = 'CpG site ± 500 bp')
lines(density(gc_target), lty = 2, col= 'blue3', lwd = 3)
legend(x = 'right', bty = 'n', legend = c('450K (filtered)', 'stochCpGs'), col = c('green4', 'blue3'), pch = 19, cex = 1.2)
abline(v = 0.461, lty = 3)
a <- legend(x = 0.461, y = 4.5, 'Average human [G+C]', box.col = "black", bg = "white", adj = 0, plot = F)
legend(x = 0.461-(a$rect[[1]])/2, y = 4.5, 'Average human [G+C]', box.col = "black", bg = "white", adj = 2*(a$text[[1]]-a$rect[[3]]), plot = T)
dev.off()

wilcox.test(gc_target, gc_bg, alternative = 'less') # p-value = 1.063e-08