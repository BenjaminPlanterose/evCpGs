############################################################################
############################################################################
###########                                                      ###########
###########                   Rep-evCpG Annotation               ###########
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
library(VennDiagram)
library(data.table)
library(genefilter)
library(matrixStats)
library(minfi)
library(lattice)
library(RColorBrewer)
library(scales)
library(missMethyl)
library(PAMES)
library(gplots)
library(ggplot2)
library(hrbrthemes)

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

# Annotate CpGs based on Illumina annotation (hg19)
functional.annotation.hg19 <- function(hvCpGs)
{
  old.dir <- getwd()
  setwd('/media/ben/DATA/Ben/1_evCpGs/data/test/')
  RGSET <- read.metharray.exp(getwd())
  annotation <- getAnnotation(RGSET)
  names <- rownames(annotation)
  manifest <- getManifest(RGSET)
  
  # Chr
  chr <- annotation$chr[match(hvCpGs, names)]
  pos <- annotation$pos[match(hvCpGs, names)]
  
  # strand
  strand <- annotation$strand[match(hvCpGs, as.vector(annotation$Name))]
  
  # Type
  I <- getProbeInfo(manifest, type = c('I'))
  I_annot <- as.vector(I$nCpG)
  names(I_annot) <- I$Name
  
  II <- getProbeInfo(manifest, type = c('II'))
  II_annot <- as.vector(II$nCpG)
  names(II_annot) <- II$Name
  
  # Types
  types <- character(length = length(hvCpGs))
  types[na.omit(match(names(I_annot), hvCpGs ))] <- 'I'
  types[na.omit(match(names(II_annot), hvCpGs ))] <- 'II'
  
  # Number of CpGs assessed by probes
  n.CpGs <- numeric(length = length(hvCpGs))
  n.CpGs[na.omit(match(names(I_annot), hvCpGs ))] <- as.vector(na.omit(I_annot[hvCpGs]))
  n.CpGs[na.omit(match(names(II_annot), hvCpGs ))] <- as.vector(na.omit(II_annot[hvCpGs]))
  
  # Island status
  island <- annotation$Relation_to_Island[match(hvCpGs, as.vector(annotation$Name))]
  
  # Feature group
  feature <- annotation$Regulatory_Feature_Group[match(hvCpGs, as.vector(annotation$Name))]
  feature[which(feature == '')] <- NA
  
  # Enhancer
  enhancer <- annotation$Enhancer[match(hvCpGs, as.vector(annotation$Name))]
  enhancer[which(enhancer == '')] <- 'FALSE'
  enhancer <- as.logical(enhancer)
  
  # Annot
  UCSC_name <- annotation$UCSC_RefGene_Name[match(hvCpGs, as.vector(annotation$Name))]
  #UCSC_accession <- annotation$UCSC_RefGene_Accession[match(hvCpGs, as.vector(annotation$Name))]
  UCSC_group <- annotation$UCSC_RefGene_Group[match(hvCpGs, as.vector(annotation$Name))]
  
  # DMR
  dmr <- annotation$DMR[match(hvCpGs, as.vector(annotation$Name))]
  
  # Probe MAF
  MAF <- annotation$Probe_maf[match(hvCpGs, as.vector(annotation$Name))]
  
  # DNAse-I hypersensitivity
  DHS <- annotation$DHS[match(hvCpGs, as.vector(annotation$Name))]
  DHS <- as.logical(DHS)
  DHS[which(DHS == '')] <- NA
  
  final.annot <- data.frame(chr = chr, pos= pos, strand = strand, type = types, nCpGs = n.CpGs,
                            island.status = island, feature.group = feature,
                            enhancer = enhancer, UCSC.name = UCSC_name, 
                            UCSC.group = UCSC_group, DMR = dmr, MAF = MAF,
                            DHS = DHS)
  
  rownames(final.annot) <- hvCpGs
  setwd(old.dir)
  return(final.annot)                    
  
  
}

# Annotate CpGs based on Illumina annotation (liftover-ed to hg38)
functional.annotation.hg38 <- function(hvCpGs)
{
  old.dir <- getwd()
  setwd('/media/ben/DATA/Ben/1_evCpGs/data/test/')
  RGSET <- read.metharray.exp(getwd())
  annotation <- getAnnotation(RGSET)
  annotation2 <- PAMES::illumina450k_hg38
  names <- annotation2$`Composite Element REF`
  manifest <- getManifest(RGSET)
  
  # Chr
  chr <- annotation2$Chromosome[match(hvCpGs, names)]
  pos <- annotation2$Start[match(hvCpGs, names)]
  
  # strand
  strand <- annotation$strand[match(hvCpGs, as.vector(annotation$Name))]
  
  # Type
  I <- getProbeInfo(manifest, type = c('I'))
  I_annot <- as.vector(I$nCpG)
  names(I_annot) <- I$Name
  II <- getProbeInfo(manifest, type = c('II'))
  II_annot <- as.vector(II$nCpG)
  names(II_annot) <- II$Name
  
  # Types
  types <- character(length = length(hvCpGs))
  types[na.omit(match(names(I_annot), hvCpGs ))] <- 'I'
  types[na.omit(match(names(II_annot), hvCpGs ))] <- 'II'
  
  # Number of CpGs assessed by probes
  n.CpGs <- numeric(length = length(hvCpGs))
  n.CpGs[na.omit(match(names(I_annot), hvCpGs ))] <- as.vector(na.omit(I_annot[hvCpGs]))
  n.CpGs[na.omit(match(names(II_annot), hvCpGs ))] <- as.vector(na.omit(II_annot[hvCpGs]))
  
  # Island status
  island <- annotation$Relation_to_Island[match(hvCpGs, as.vector(annotation$Name))]
  
  # Feature group
  feature <- annotation$Regulatory_Feature_Group[match(hvCpGs, as.vector(annotation$Name))]
  feature[which(feature == '')] <- NA
  
  # Enhancer
  enhancer <- annotation$Enhancer[match(hvCpGs, as.vector(annotation$Name))]
  enhancer[which(enhancer == '')] <- 'FALSE'
  enhancer <- as.logical(enhancer)
  
  # Annot
  UCSC_name <- annotation$UCSC_RefGene_Name[match(hvCpGs, as.vector(annotation$Name))]
  #UCSC_accession <- annotation$UCSC_RefGene_Accession[match(hvCpGs, as.vector(annotation$Name))]
  UCSC_group <- annotation$UCSC_RefGene_Group[match(hvCpGs, as.vector(annotation$Name))]
  
  # DMR
  dmr <- annotation$DMR[match(hvCpGs, as.vector(annotation$Name))]
  
  # Probe MAF
  MAF <- annotation$Probe_maf[match(hvCpGs, as.vector(annotation$Name))]
  
  # DNAse-I hypersensitivity
  DHS <- annotation$DHS[match(hvCpGs, as.vector(annotation$Name))]
  DHS <- as.logical(DHS)
  DHS[which(DHS == '')] <- NA
  
  Gene_type <- annotation2$Gene_Type[match(hvCpGs, names)]
  
  final.annot <- data.frame(chr = chr, pos= pos, strand = strand, type = types, nCpGs = n.CpGs,
                            island.status = island, feature.group = feature,
                            enhancer = enhancer, UCSC.name = UCSC_name, 
                            UCSC.group = UCSC_group, DMR = dmr, MAF = MAF,
                            DHS = DHS, Gene_type = Gene_type)
  rownames(final.annot) <- hvCpGs
  setwd(old.dir)
  return(final.annot)                    
}

# Custom deep annotation based on UCSC refGene (hg19). Include intron/exon information
functional.annotation.deep <- function(hvCpG, d)
{
  setwd('/media/ben/DATA/Ben/1_evCpGs/data/test/')
  RGSET <- read.metharray.exp(getwd())
  annotation <- getAnnotation(RGSET)
  annotation <- as.data.frame(annotation[hvCpG, 1:3])
  
  setwd('/media/ben/DATA/Ben/1_evCpGs/annotation/')
  annot <- fread('refGene.txt')
  
  annotated1 <- annot.within(annotation, annot)
  annotated2 <- annot.nearby(annotation, annot, d)
  
  return(cbind(annotated1, annotated2))
}

# Subfunction for functional.annotation.deep
annot.nearby <- function(annotation, annot, d)
{
  annotated2 <- annotation
  head(annotated2)
  colnames(annotated2) <- c("chr",    "pos CpG",    "strand")
  annotated2$gene <- NA
  annotated2$pos <- NA
  annotated2$strandgene <- NA
  annotated2$refseq <- NA
  
  for(i in 1L:(nrow(annotation)))
  {
    chr <- annotation[i,1]
    pos <- annotation[i,2]
    subset <- annot[annot$V3 == chr,]
    subset1 <- subset[subset$V4 == '+',]
    subset2 <- subset[subset$V4 == '-',]
    subset1 <- subset1[pos > subset1$V5 - d & pos < subset1$V5,]
    subset2 <- subset2[pos > subset2$V6 & pos < subset2$V6 + d]
    
    acum_strand <- character()
    acum_pos <- character()
    acum_gene <- character()
    acum_refseq <- character()
    if(nrow(subset1) >= 1)
    {
      acum_strand <- c(acum_strand, subset1$V4)
      acum_pos <- c(acum_pos, paste(subset1$V5, subset1$V6, sep = '-'))
      acum_gene <- c(acum_gene, subset1$V13)
      acum_refseq <- c(acum_refseq, subset1$V2)
    }
    if(nrow(subset2) >= 1)
    {
      acum_strand <- c(acum_strand, subset2$V4)
      acum_pos <- c(acum_pos, paste(subset2$V5, subset2$V6, sep = '-'))
      acum_gene <- c(acum_gene, subset2$V13)
      acum_refseq <- c(acum_refseq, subset2$V2)
    }
    if(length(acum_strand) >= 1)
    {
      annotated2$gene[i] <- paste(acum_gene, collapse = ' ')
      annotated2$pos[i] <- paste(acum_pos, collapse = ' ')
      annotated2$strandgene[i] <- paste(acum_strand, collapse = ' ')
      annotated2$refseq[i] <- paste(acum_refseq, collapse = ' ')
    }
    print(i)
  }
  return(annotated2)
}

# Subfunction for functional.annotation.deep
annot.within <- function(annotation, annot)
{
  annotated <<- annotation
  colnames(annotated) <<- c("chr",    "pos CpG",    "strand")
  annotated$gene <<- NA
  annotated$pos <<- NA
  annotated$region <<- NA
  annotated$strandgene <<- NA
  annotated$refseq <<- NA
  
  for(i in 1L:(nrow(annotation)))
  {
    chr <- annotation[i,1]
    pos <- annotation[i,2]
    subset <- annot[annot$V3 == chr,]
    subset <- subset[subset$V5 <= pos & subset$V6 >= pos,]
    
    if(nrow(subset) == 1)
    {
      single.entry(subset, annot, i, pos)
    }
    
    else if(nrow(subset) > 1)
    {
      multi.entry(subset, annot, i, pos)
    }
    print(i)
  }
  return(annotated)
}

# Subfunction for functional.annotation.deep
single.entry <- function(subset, annot, i, pos)
{
  annotated$gene[i] <<- subset$V13
  annotated$pos[i] <<- paste(subset$V5, subset$V6, sep = '-')
  annotated$strandgene[i] <<- unique(subset$V4)
  annotated$refseq[i] <<- subset$V2
  
  A <- as.numeric(unlist(strsplit(subset$V10, split = ',')))
  B <- as.numeric(unlist(strsplit(subset$V11, split = ',')))
  
  if(annotated$strandgene[i] == '+')
  {
    single.entry.plus(subset, annot, i, A, B, pos)
  }
  
  else
  {
    single.entry.minus(subset, annot, i, A, B, pos)
  }
}

# Subfunction for functional.annotation.deep
multi.entry <- function(subset, annot, i, pos)
{
  annotated$gene[i] <<- paste(subset$V13, collapse = ' ')
  annotated$strandgene[i] <<- paste(unique(subset$V4), collapse = ' ')
  annotated$refseq[i] <<- paste(subset$V2, collapse = ' ')
  
  if(annotated$strandgene[i] == '+')
  {
    multi.entry.plus(subset, annot, i, pos)
  }
  else
  {
    multi.entry.minus(subset, annot, i, pos)
  }
  
}

# Subfunction for functional.annotation.deep
single.entry.plus <- function(subset, annot, i, A, B, pos)
{
  index <- which(sort(c(pos, A, B)) == pos)
  if(index %/% 2 == index/2) # It is in an exon or 5´UTR
  {
    if(index == 2 | index == length(c(A,B))) # If in the first or last exon
    {
      res <- UTR_test_plus(subset, pos) # Check if it is in UTR
      if(res[[1]] == TRUE)
      {
        annotated$region[i] <<- res[[2]]
      }
      else
      {
        annotated$region[i] <<- paste('exon', index/2, 'of', subset$V9)
      }
    }
    
    else
    {
      annotated$region[i] <<- paste('exon', index/2, 'of', subset$V9)
    }
  }
  else # It is in an intron
  {
    annotated$region[i] <<- paste('intron', floor(index/2), 'of', subset$V9 -1)
  }
}

# Subfunction for functional.annotation.deep
single.entry.minus <- function(subset, annot, i, A, B, pos)
{
  index <- which(sort(c(pos, A, B), decreasing = T) == pos)
  if(index %/% 2 == index/2) # It is in an exon
  {
    annotated$region[i] <<- paste('exon', index/2, 'of', subset$V9)
    if(index == 2 | index == length(c(A,B))) # If in the first or last exon
    {
      res <- UTR_test_minus(subset, pos) # Check if it is in UTR
      if(res[[1]] == TRUE)
      {
        annotated$region[i] <<- res[[2]] # Change your mind
      }
    }
  }
  else # It is in an intron
  {
    annotated$region[i] <<- paste('intron', floor(index/2), 'of', subset$V9 -1)
  }
}

# Subfunction for functional.annotation.deep
multi.entry.plus <- function(subset, annot, i, pos)
{
  acum <- ''
  acum2 <- ''
  for(j in 1:nrow(subset))
  {
    #j <- 1
    A <- as.numeric(unlist(strsplit(subset$V10[j], split = ',')))
    B <- as.numeric(unlist(strsplit(subset$V11[j], split = ',')))
    index <- which(sort(c(pos, A, B)) == pos)
    temp <- paste(subset$V5[j], subset$V6[j], sep = '-')
    acum2 <- paste(acum2, temp)
    if(index %/% 2 == index/2) # It is in an exon or UTR
    {
      if(index == 2 | index == length(c(A,B))) # If in the first or last exon
      {
        res <- UTR_test_plus(subset[j,], pos) # Check if it is in UTR
        if(res[[1]] == TRUE)
        {
          acum <- paste(acum, ',', res[[2]]) # Change your mind
        }
        else
        {
          acum <- paste(acum, ',', 'exon', index/2, 'of', subset$V9[j])
        }
      }
      else
      {
        acum <- paste(acum, ',', 'exon', index/2, 'of', subset$V9[j])
      }
    }
    else # It is in an intron
    {
      acum <- paste(acum, ',', 'intron', floor(index/2), 'of', subset$V9[j])
    }
  }
  annotated$region[i] <<- substr(acum, start = 3, stop = nchar(acum))
  annotated$pos[i] <<- acum2
}

# Subfunction for functional.annotation.deep
multi.entry.minus <- function(subset, annot, i, pos)
{
  acum <- ''
  acum2 <- ''
  for(j in 1:nrow(subset))
  {
    A <- as.numeric(unlist(strsplit(subset$V10[j], split = ',')))
    B <- as.numeric(unlist(strsplit(subset$V11[j], split = ',')))
    index <- which(sort(c(pos, A, B), decreasing = T) == pos)
    temp <- paste(subset$V5[j], subset$V6[j], sep = '-')
    acum2 <- paste(acum2, temp)
    if(index %/% 2 == index/2) # It is in an exon
    {
      if(index == 2 | index == length(c(A,B))) # If in the first or last exon
      {
        res <- UTR_test_plus(subset[j,], pos) # Check if it is in UTR
        if(res[[1]] == TRUE)
        {
          acum <- paste(acum, ',', res[[2]]) # Change your mind
        }
        else
        {
          acum <- paste(acum, ',', 'exon', index/2, 'of', subset$V9[j])
        }
      }
      else
      {
        acum <- paste(acum, ',', 'exon', index/2, 'of', subset$V9[j])
      }
    }
    else # It is in an intron
    {
      acum <- paste(acum, ',', 'intron', floor(index/2), 'of', subset$V9[j])
    }
  }
  annotated$region[i] <<- substr(acum, start = 3, stop = nchar(acum))
  annotated$pos[i] <<- acum2
}

# Subfunction for functional.annotation.deep
UTR_test_plus <- function(subset, pos)
{
  cds.init <- subset$V7
  cds.end <- subset$V8
  if(pos < cds.init && cds.init != cds.end)
  {
    res <- list(TRUE, "5'-UTR")
  }
  else if(pos > cds.end && cds.init != cds.end)
  {
    res <- list(TRUE, "3'-UTR")
  }
  else
  {
    res <- list(FALSE, "Exon")
  }
}

# Subfunction for functional.annotation.deep
UTR_test_minus <- function(subset, pos)
{
  cds.init <- subset$V7
  cds.end <- subset$V8
  if(pos < cds.init && cds.init != cds.end)
  {
    res <- list(TRUE, "3'-UTR")
  }
  else if(pos > cds.end && cds.init != cds.end)
  {
    res <- list(TRUE, "5'-UTR")
  }
  else
  {
    res <- list(FALSE, "Exon")
  }
}

# Performs enrichment per category for CpG island status/functional status
enrichment.per.category <- function(mat)
{
  pvals <- numeric(length = ncol(mat))
  odd_ratio <- numeric(length = ncol(mat))
  names(pvals) <- colnames(mat)
  names(odd_ratio) <- colnames(mat)
  
  for(i in 1:ncol(mat))
  {
    category <- colnames(mat)[i]
    data <- cbind(mat[,i], rowSums(mat[,-i]))
    rownames(data) <- c('Bg', 'target')
    colnames(data) <- c(category, paste('Not', category))
    message(category)
    print(data)
    res <- fisher.test(data)
    print(res)
    pvals[i] <- res$p.value
    odd_ratio[i] <-res$estimate
  }
  output <- list(pvals, odd_ratio)
  names(output) <- c('pvals', 'Odd_ratios')
  return(output)
}
enrichment.per.category <- function(mat)
{
  pvals <- numeric(length = ncol(mat))
  odd_ratio <- numeric(length = ncol(mat))
  names(pvals) <- colnames(mat)
  names(odd_ratio) <- colnames(mat)
  
  for(i in 1:ncol(mat))
  {
    category <- colnames(mat)[i]
    data <- cbind(mat[,i], rowSums(mat[,-i]))
    rownames(data) <- rownames(mat)
    colnames(data) <- c(category, paste('Not', category))
    message(category)
    print(data)
    res <- fisher.test(data)
    print(res)
    pvals[i] <- res$p.value
    odd_ratio[i] <-res$estimate
  }
  output <- list(pvals, odd_ratio)
  names(output) <- c('pvals', 'Odd_ratios')
  return(output)
}

# Chooses highest priority function status term
process.priority <- function(func, priority)
{
  split_class <- lapply(strsplit(as.vector(func), split = ';'), unique)
  split_class[unlist(lapply(split_class, length)) == 0] <- 'Not gene associated'
  where <- split_class[unlist(lapply(split_class, length)) > 1]
  where_priority <- choose_priority(where, priority)
  split_class[unlist(lapply(split_class, length)) > 1] <- where_priority
  split_class <- table(unlist(split_class))
  return(split_class[priority])
}

# Subfunction for process.priority
choose_priority <- function(choice_list, priority)
{
  for(i in 1:length(choice_list))
  {
    choice_list[[i]] <- priority[min(match(choice_list[[i]], priority))]
  }
  return(choice_list)
}

# Performs positional enrichment analysis
positional_enrichment_analysis <- function(target, bg, window)
{
  # Extract annotation
  old.dir <- getwd()
  setwd('/media/ben/DATA/Ben/1_evCpGs/data/test/')
  RGSET <- read.metharray.exp(getwd())
  annotation <- getAnnotation(RGSET)
  setwd(old.dir)
  
  # Extract annotation for target and bg
  annotation_target <- annotation[target,]
  annotation_bg <- annotation[bg,]
  
  # Segregate annotation per chromosome
  levels <- unique(annotation_target$chr)
  annotation_target <- lapply(1:length(levels), function(x) annotation_target[annotation_target$chr == levels[x],1:2])
  names(annotation_target) <- levels
  
  levels <- unique(annotation_bg$chr)
  annotation_bg <- lapply(1:length(levels), function(x) annotation_bg[annotation_bg$chr == levels[x],1:2])
  names(annotation_bg) <- levels
  
  p_val <- numeric()
  for(i in 1:length(target))
  {
    chr <- annotation[target[i], 1]
    pos <- annotation[target[i], 2]
    annot_chr_bg <- annotation_bg[[chr]]
    annot_chr_target <- annotation_target[[chr]]
    
    count_bg <- sum(annot_chr_bg$pos < pos + window/2  & annot_chr_bg$pos > pos - window/2)
    count_target <- sum(annot_chr_target$pos < pos + window/2  & annot_chr_target$pos > pos - window/2)
    
    contrast <- matrix(c(count_target, nrow(annot_chr_target) - count_target,
                         count_bg, nrow(annot_chr_bg) - count_bg), nrow = 2)
    
    p_val[i] <- (fisher.test(contrast))$p.value
    if(i %/% 100 == i/100)
    {
      print(i)
    }
  }
  names(p_val) <- target
  return(p_val)
}

# Target enrichment analysis
enrichment_pergene <- function(target, bg, key)
{
  a <- sum(startsWith(target, prefix = key))
  b <- sum(startsWith(bg, prefix = key))
  c <- length(target) - a
  d <- length(bg) - b
  contrast <- matrix(c(a,c,b,d), nrow = 2)
  colnames(contrast) <- c('target', 'bg')
  rownames(contrast) <- c(key, paste('Not', key, sep = ' '))
  res <- fisher.test(contrast)
  
  print(contrast)
  message('p-val enrichment')
  print(res$p.value)
  message('MLE of odds ratio')
  print(res$estimate)
}


#############################################   Annotation   #############################################

# Read data
setwd("/media/ben/DATA/Ben/1_evCpGs/data/adipose/")
stochCpG <- as.vector(read.table(file = 'replicated.txt', header = F)$V1); length(stochCpG) # 154
sig = stochCpG

# Illumina annotation (hg19)
annotation <- functional.annotation.hg19(stochCpG)
annotation <- annotation[order(annotation$chr, annotation$pos),]
setwd('~/Desktop')
write.table(file = 'annot.txt', annotation, row.names = T, col.names = T, sep = '\t', quote = F)

# PAMES annotation (hg38)
annotation2 <- functional.annotation.hg38(stochCpG)
annotation2 <- annotation2[order(annotation2$chr, annotation2$pos),]
setwd('~/Desktop')
write.table(file = 'annot2.txt', annotation2, row.names = T, col.names = T, sep = '\t', quote = F)

# Deep annotation (hg19)
annotation3 <- functional.annotation.deep(stochCpG, 1500)
annotation3 <- annotation3[order(annotation3$chr, annotation3$`pos CpG`),]
setwd('~/Desktop')
write.table(file = 'annot3.txt', annotation3, row.names = T, col.names = T, sep = '\t', quote = F)

#############################################   Enrichment   #############################################

# Read whole annotation
setwd('/media/ben/DATA/Ben/1_evCpGs/data/test/')
RGSET <- read.metharray.exp(getwd())
complete_annot <- getAnnotation(RGSET)
annotation <- functional.annotation.hg19(stochCpG)


# Define bg0
Y_probes = rownames(complete_annot)[complete_annot$chr == "chrY"]
X_probes = rownames(complete_annot)[complete_annot$chr == "chrX"]
data(SNPs.147CommonSingle)
f.SNP <- c(rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$Probe_maf >= 0.01],
           rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$CpG_maf > 0],
           rownames(SNPs.147CommonSingle)[SNPs.147CommonSingle$SBE_maf > 0])
SNP_probes <- na.omit(unique(f.SNP)); length(SNP_probes) # 99337

setwd("/media/ben/DATA/Ben/1_evCpGs/annotation/")
CR_1 <- as.vector(read.table('crossreactive_Chen.txt', header = F)$TargetID) # Chen YA et al. Epigenetics. 2013 Feb;8(2):203-9. doi: 10.4161/epi.23470. Epub 2013 Jan 11.
kobor <- fread('GPL16304-47833.txt') # Price ME et al. Epigenetics Chromatin. 2013 Mar 3;6(1):4. doi: 10.1186/1756-8935-6-4.
CR_2 <- unique(c(kobor$ID[kobor$Autosomal_Hits == 'A_YES'], kobor$ID[kobor$XY_Hits == 'XY_YES']))
CR_probes <- unique(c(CR_1, CR_2)); length(CR_probes) # 41937
probes2remove <- unique(c(Y_probes, X_probes, SNP_probes, CR_probes))
length(probes2remove) # 137884

bg0 = rownames(complete_annot)[!(rownames(complete_annot) %in% probes2remove)]
length(bg0) # 347628
bg = bg0[!(bg0 %in% sig)] # Bg and target must be mutually exclusive
length(bg) #  347474
bg_complete = rownames(complete_annot)

#############################################   Imprinted genes   #############################################


# Read imprinted genes, from http://www.geneimprint.com/site/genes-by-species
setwd('/media/ben/DATA/Ben/1_evCpGs/annotation/')
imprinted <- fread('IMPRINTED')
# Any alias will do. Combine all into a single vector
imprinted_genes <- unique(c(imprinted$Gene, unlist(strsplit(imprinted$Aliases, split = ', '))))
length(imprinted_genes) # 1004 aliases

# Extract associated genes in target
associated_genes <- unlist(lapply(strsplit(as.character(annotation$UCSC.name), ';'), function(x) x[1]))
associated_genes[associated_genes %in% imprinted_genes]
# BEFORE: "FUCA1" "HOXA5" "IGF2"  "HOXC4" "HTR2A" "DLK1"  "GNAS" 
# AFTER: "HOXA5" "HOXA5" "HOXC4" "DLK1"  "GNAS"

# Extract associated genes in bg
bg_associated_genes <- unlist(lapply(strsplit(as.character(complete_annot[bg,]$UCSC_RefGene_Name), ';'), function(x) x[1]))

# Count number of genes that are imprinted
a = sum(associated_genes %in% imprinted_genes) # 7
b = sum(bg_associated_genes %in% imprinted_genes) # 6940

# Build matrix
mat <- matrix(c(a, length(associated_genes) - a, b, length(bg_associated_genes)-b), ncol = 2, byrow = F)
#      [,1]   [,2]
# [1,]    5   6971
# [2,]  149 340503

# Perform fisher´s exact test
fisher.test(mat) 
# Fisher's Exact Test for Count Data
# 
# data:  mat
# p-value = 0.2425
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.5243717 3.9130298
# sample estimates:
# odds ratio 
#   1.639098 

#############################################   Metastable epialleles   #############################################


# Metastable epialleles
setwd("/media/ben/DATA/Ben/1_evCpGs/annotation/")
metastable <- as.vector(read.table('metastable', header = F)$V1) # Harris, RA et al (2013) Epigenetics. 2013 Feb 1; 8(2): 157–163. 

venn(list(sig, metastable)) # 4 CpGs
functional.annotation.hg19(sig[sig %in% metastable])[,c(1:2, 9)]
# cg16181718 chr4 190731655                                                
# cg01575930 chr5 140733385 PCDHGA2;PCDHGA4;PCDHGA1;PCDHGA4;PCDHGB1;PCDHGA3
# cg27341636 chr6  27236324                                                
# cg01696605 chr7 150871692 

# Count number of metastables in background and target
a = sum(stochCpG %in% metastable) # 4
b = sum(bg %in% metastable) # 1218

# Build matrix
mat <- matrix(c(a, length(associated_genes) - a, b, length(bg_associated_genes)-b), ncol = 2, byrow = F)
#      [,1]   [,2]
# [1,]    4   1218
# [2,]  150 346256

# Perform Fisher´s exact test
fisher.test(mat) 
# 	Fisher's Exact Test for Count Data
# data:  mat
# p-value = 0.002256
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   2.035738 19.865667
# sample estimates:
#   odds ratio 
# 7.580866 


#############################################   EWAS trait association   #############################################

setwd("/media/ben/DATA/Ben/1_evCpGs/annotation/") # https://bigd.big.ac.cn/ewas/index
#assoc <- fread('association.tsv', sep = '\t'); dim(assoc) # 416331      8
assoc <- fread('association.tsv', sep = '\t'); dim(assoc) # 534437      11
head(assoc)

assoc <- as.data.frame(assoc)
dim(assoc) # 416331     11

# Terms to test for enrichment
terms <- unique(assoc[assoc$`Probe id` %in% sig, 2]); length(terms) # 81

# Extract terms in target and background
target.t <- assoc[assoc$`Probe id` %in% sig, 2]
length(target.t) # 502
bg.t <- assoc[assoc$`Probe id` %in% bg, 2]
length(bg.t) # 293341

# Show most common terms
head(sort(table(target.t), decreasing = T), 12)


# Perform enrichment for all terms
pvali <- numeric(length = length(terms))
names(pvali) <- terms
OR <- numeric(length = length(terms))
names(OR) <- terms
for(i in 1:length(terms))
{
  a = sum(target.t %in% terms[i])
  b = sum(bg.t %in% terms[i])
  
  mat <- matrix(c(a, length(target.t) - a, b, length(bg.t)-b), ncol = 2, byrow = T)
  test = fisher.test(mat, alternative = 'greater')
  pvali[i] <- test$p.value
  OR[i] <- test$estimate
  if(pvali[i] < 0.05/length(terms))
  {
    print(terms[i])
    print(mat)
    print(test)
  }
  
}

# Adjust for multiple testing correction
pval.adj <- p.adjust(pvali, 'bonferroni')


# Visualize
par(mar=c(5.1, 20, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
barplot(sort(-log10(pval.adj[pval.adj < 0.05]), decreasing = F), horiz = T, srt = 90, las = 2, col = 'dodgerblue')
barplot((OR[pval.adj < 0.05])[order(-log10(pval.adj[pval.adj < 0.05]))], horiz = T, srt = 90, las = 2, col = 'dodgerblue')

#############################################   RNA expression   #############################################


setwd('/media/ben/DATA/Ben/1_evCpGs/annotation/')
genes <- unique(unlist(strsplit(as.vector(annotation$UCSC.name), split = ';')))
length(genes) # 121

gene_expression <- fread('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct')
gene_expression <- gene_expression[gene_expression$Description %in% genes,]
dim(gene_expression) # 110 55

X <- melt(gene_expression)
X$value <- log2(X$value + 1)

setwd('~/Desktop/')
tiff(filename = paste('RNA', 'tiff', sep = '.'), width = 8, height = 10, units = 'in', res = 300, compression = 'none')
ggplot(X, aes(variable, Description, fill= value)) + 
  geom_tile() + scale_fill_gradient(low="white", high="blue") +
  theme_ipsum(base_size = 3, axis_title_size = 0) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



#############################################   Island status enrichment   #############################################

# Extract bg/target Island status
island_bg <- complete_annot[bg, ]$Relation_to_Island
island_target <- complete_annot[sig, ]$Relation_to_Island

# Build contigency table
mat1 <- rbind(table(island_bg), table(island_target))
rownames(mat1) <- c('Bg', 'target')
#        Island N_Shelf N_Shore OpenSea S_Shelf S_Shore
# Bg     111727   16868   46593  121003   15005   36278
# target     45       6      28      41       4      30
# Global enrichment
set.seed(1); fisher.test(mat1, simulate.p.value = T, B = 100000) # p-value = 0.00041

# Per category
res1 <- enrichment.per.category(mat1)
p.adjust(res1$pvals, method = 'bonferroni') < 0.05
# Island N_Shelf N_Shore OpenSea S_Shelf S_Shore 
# FALSE   FALSE   FALSE   FALSE   FALSE    TRUE 

y <- barplot(1/res1$Odd_ratios-1, horiz = T, xaxt = 'n', xlab = 'Target to Background Odds ratio', las = 1,
             col = c('brown3', 'brown3','chartreuse4', 'brown3', 'brown3', 'chartreuse4'), xlim = c(1/2-1, 1/0.5-1))
# axis(side = 1, labels = c(0.5, 1, 1.5), at = c(-0.5, 0, 0.5))
# abline(v = 0, lty = 1)
# text(x = c(0.1, 0.1, -0.1, 0.1, 0.1, -0.1), y = y, c('ns', 'ns', 'ns', 'ns', 'ns', '*'))


mat1 <- rbind(table(island_bg), table(island_target))
mat1[1,] <- mat1[1,]/sum(mat1[1,])
mat1[2,] <- mat1[2,]/sum(mat1[2,])
rownames(mat1) <- c('Bg', 'Target')
mat1 <- mat1[,ncol(mat1):1]
balloonplot(as.table(mat1), main = '', ylab = 'Island Status', xlab = 'Set')


#############################################   Functional status enrichment   #############################################

# Extract bg/target functional status
func_bg1 <- complete_annot[bg, ]$UCSC_RefGene_Group
func_target <- complete_annot[sig, ]$UCSC_RefGene_Group

# Define priority
priority <- c("5'UTR","TSS200", "TSS1500", "1stExon", "Body", "3'UTR", "Not gene associated")

# Decide which term remains on target/bg
count_bg1 <- process.priority(func = func_bg1, priority = priority)
count_target <- process.priority(func = func_target, priority = priority)

# Build contingency table
mat1 <- rbind(count_bg1, count_target)
rownames(mat1) <- c('Bg', 'target')
#        5'UTR TSS200 TSS1500 1stExon   Body 3'UTR Not gene associated
# Bg     49226  41022   50467    8163 110195 11454               76947
# target     9     10      39      10     37     3                  46


# Global enrichment
set.seed(2)
fisher.test(mat1, simulate.p.value = T, B = 100000) # p-value = 1e-05

# Perform enrichment per category
res1 <- enrichment.per.category(mat1)
p.adjust(res1$pvals, method = 'bonferroni') < 0.05
# 5'UTR              TSS200             TSS1500             1stExon 
#                TRUE               FALSE                TRUE                TRUE 
#                Body               3'UTR Not gene associated 
# FALSE               FALSE               FALSE 

# Plot
par(mar=c(5.1, 10.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
y <- barplot(1/res1$Odd_ratios-1, horiz = T, xaxt = 'n', xlab = 'Target to Background Odds ratio', las = 1,
             col = c('brown3', 'brown3','chartreuse4', 'chartreuse4','brown3', 'brown3', 'chartreuse4'),
             xlim = c(-0.5, 1.5))
axis(side = 1, labels = c(0.5, 1, 1.5, 2, 2.5), at = c(-0.5, 0, 0.5, 1, 1.5))
abline(v = 0, lty = 1)
# 
# p.adjust(res1$pvals, method = 'bonferroni') < 0.001
# p.adjust(res1$pvals, method = 'bonferroni') < 0.01
# p.adjust(res1$pvals, method = 'bonferroni') < 0.05
# text(x = c(0.1, 0.1, -0.1, -0.1, 0.1, 0.1, -0.1), y = y, c('*', '**', '*', '*', 'ns', 'ns', '*'))

mat1[1,] <- mat1[1,]/sum(mat1[1,])
mat1[2,] <- mat1[2,]/sum(mat1[2,])
rownames(mat1) <- c('Bg', 'Target')
mat1 <- mat1[,ncol(mat1):1]
balloonplot(as.table(mat1), main = '', ylab = 'Functional Status', xlab = 'Set')



#############################################   ChromHMM enrichment analysis   #############################################

setwd('/media/ultron/2tb_disk2/PROCESSED_DATA/2019/HMM_annot/PBMC/')
full = fread('PBMC_15stat_annot.txt')
full = as.data.frame(full)
CpGs = full$CpG
full = full[,-1]
rownames(full) <- CpGs

full$HMM_450K <- as.factor(x = full$HMM_450K)

target_hmm = full[rownames(full) %in% stochCpG, 3]
bg_hmm = full[rownames(full) %in% bg, 3]

counts = rbind(table(bg_hmm), table(target_hmm))
rownames(counts) <- c('Bg', 'Target')
rowSums(counts)
#     Bg Target 
# 345963    333 (259 unmapped elements were substracted from the bg)


counts <- counts[,-ncol(counts)]
counts
#         10_TssBiv 11_BivFlnk 12_EnhBiv 13_ReprPC 14_ReprPCWk 15_Quies 1_TssA 2_TssAFlnk 3_TxFlnk  4_Tx 5_TxWk
# Bg          9795      10454      5581     33219       37902    74233  26541      62878      969 33635  32102
# Target         3          7         1        66          69       70     20         29        0    11     22
#        6_EnhG 7_Enh 8_ZNF/Rpts 9_Het
# Bg       3364 12411        617  2262
# Target      0    15          3    17

set.seed(1994)
fisher.test(counts, simulate.p.value = T, B = 100000)
# 	Fisher's Exact Test for Count Data with simulated p-value (based on 1e+05 replicates)
# data:  counts
# p-value = 1e-05
# alternative hypothesis: two.sided

# Performs enrichment per category for CpG island status/functional status

res <- enrichment.per.category(counts)
mat1 <- counts
mat1[1,] <- mat1[1,]/sum(mat1[1,])
mat1[2,] <- mat1[2,]/sum(mat1[2,])

order_lev <- as.numeric(lapply(strsplit(colnames(mat1), split = '_'), function(x) x[1]))
mat1 <- mat1[,order(order_lev)]
balloonplot(as.table(mat1), main = '', ylab = 'Island Status', xlab = 'Set')



par(mar=c(5.1, 10.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
order_lev <- as.numeric(lapply(strsplit(names(res$Odd_ratios), split = '_'), function(x) x[1]))
x = 1/res$Odd_ratios[order(order_lev, decreasing = T)]
x[x == 0] = NA

x = log2(x)

y <- barplot(x, horiz = T, xlab = 'log2(Target to Background Odds ratio)', las = 1,
             col = c('brown3', 'chartreuse4','chartreuse4', 
                     'brown3','brown3', 'brown3',
                     'chartreuse4', 'chartreuse4', 'chartreuse4',
                     'brown3', 'brown3', 'brown3', 'brown3','brown3', 'brown3'),
             xlim = c(-4, 4))
abline(v = 0, lty = 1)

p.adjust(res$pvals[order(order_lev, decreasing = T)], method = 'bonferroni') < 0.001
p.adjust(res$pvals[order(order_lev, decreasing = T)], method = 'bonferroni') < 0.01
p.adjust(res$pvals[order(order_lev, decreasing = T)], method = 'bonferroni') < 0.05
text(x = c(0.5, -0.5, -0.5, 
           0.5, 0.5, 0.5,
           -0.5, -0.5, -0.5,
           rep(0.5, 6)), 
     y = y, 
     c('ns', '***', '***', 'ns', 'ns', 
       'ns', '***', 'ns', 'ns', 'ns',
       'ns', '***', 'ns', '***', 'ns'))


#############################################   Positional enrichment analysis   #############################################



# 500 bp up-/down-stream
pvals_enrich <- positional_enrichment_analysis(sig, bg, 1000)
pos_enrich_res <- p.adjust(pvals_enrich, 'bonferroni')
complete_annot[names(which(pos_enrich_res < 0.05)),]$UCSC_RefGene_Name
# [1] "CACNA1E"                                        
# [2] "LIMS1"                                          
# [3] ""                                               
# [4] ""                                               
# [5] "ZIC4;ZIC4;ZIC4"                                 
# [6] "PCDHGA2;PCDHGA4;PCDHGA1;PCDHGA4;PCDHGB1;PCDHGA3"
# [7] "PCDHGA2;PCDHGA4;PCDHGA1;PCDHGA4;PCDHGB1;PCDHGA3"
# [8] ""                                               
# [9] ""                                               
# [10] "LOC389458"                                      
# [11] "LOC389458"                                      
# [12] ""                                               
# [13] ""                                               
# [14] "OR4A16"                                         
# [15] "TMEM126A"                                       
# [16] ""                                               
# [17] ""                                               
# [18] "TNFAIP8L3"  

a <- as.data.frame(annotation[names(which(pos_enrich_res < 0.05)),])
a$pval <- pvals_enrich[rownames(a)]
# setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/res/supp/')
# write.table(x = a, file = 'regions_1kb_positional_enrichment.txt', quote = F, sep = '\t')


# 1 Mb
pvals_enrich <- positional_enrichment_analysis(sig, bg, 10^6)
pos_enrich_res <- p.adjust(pvals_enrich, 'bonferroni')
setwd('/media/ultron/2tb_disk2/Papers/Genome_independent_interindividual_variation_in_dna_methylation_2/res/supp/')
complete_annot[names(which(pos_enrich_res < 0.05)),]$UCSC_RefGene_Name
# [1] "VTRNA1-3"                                                                                                    
# [2] "PCDHA2;PCDHA1;PCDHA1;PCDHA2;PCDHA2"                                                                          
# [3] "PCDHA6;PCDHA2;PCDHA1;PCDHA9;PCDHA7;PCDHA1;PCDHA11;PCDHA6;PCDHA5;PCDHA11;PCDHA10;PCDHA3;PCDHA4;PCDHA10;PCDHA8"
# [4] "PCDHB4"                                                                                                      
# [5] "PCDHB5;PCDHB5"                                                                                               
# [6] "PCDHB6"                                                                                                      
# [7] "PCDHB17"                                                                                                     
# [8] "PCDHB7"                                                                                                      
# [9] "PCDHB16"                                                                                                     
# [10] "PCDHB14;PCDHB14"                                                                                             
# [11] "PCDHGA2;PCDHGA4;PCDHGA1;PCDHGA4;PCDHGB1;PCDHGA3"                                                             
# [12] "PCDHGA2;PCDHGA4;PCDHGA1;PCDHGA4;PCDHGB1;PCDHGA3"                                                             
# [13] "PCDHGA4;PCDHGA2;PCDHGA5;PCDHGB2;PCDHGA1;PCDHGB1;PCDHGA3;PCDHGA5"                                             
# [14] ""                                                                                                            
# [15] ""                                                                                                            
# [16] ""                                                                                                            
# [17] ""                              

# a <- as.data.frame(annotation[names(which(pos_enrich_res < 0.05)),])
# a$pval <- pvals_enrich[rownames(a)]
# write.table(x = a, file = 'regions_1Mb_positional_enrichment.txt', quote = F, sep = '\t')

