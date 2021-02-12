library(vcfR)
library(tidyverse) 
library("DESeq2")

setwd("/Users/pamelacamejo/Documents/IBIO/Rodrigo_Gutierrez/Projects/Javier_Paz-Ares/Analyses/ASEReadCounter/")

# Import vcf file
vcf <- read.vcfR("../Freebayes/Ler/fb.LeLe.LlLl.filt.final.genome.vcf")

# Remove indels from vcf
vcf_filt <- extract.indels(vcf)
vcf_fix <- data.frame(vcf_filt@fix)
colnames(vcf_fix) <- c("contig","position","ID","refAllele","altAllele",colnames(vcf_fix)[6:8])

# Merge counts file to create a table with all samples
vcf_fix_counts_merge <- vcf_fix 
file.names <- dir("Counts_Hisat_Ler_genome", pattern =".out",full.names = TRUE)
for(i in 1:length(file.names)){
  sample <- unlist(strsplit(unlist(strsplit(file.names[i],"/"))[2],"_"))[1]
  replicate <- unlist(strsplit(unlist(strsplit(file.names[i],"/"))[2],".",fixed=TRUE))[1]
  counts <- read.table(file.names[i],header = TRUE)
  counts_selected <- counts[,c(1:2,4:8,11)]
  colnames(counts_selected) <- c("contig","position","refAllele","altAllele",paste0("refCount_",replicate),paste0("altCount_",replicate),paste0("totalCount_",replicate), paste0("rawDepth_",replicate))
  
  vcf_fix_counts_merge <- merge(vcf_fix_counts_merge, counts_selected, by=c("contig","position","refAllele","altAllele"),all.x = TRUE) 
}

# SNPs with no counts is zero.
vcf_fix_counts_merge[is.na(vcf_fix_counts_merge)] <- 0

# Export file with counts
write.table(vcf_fix_counts_merge,"Results/all_discriminants_snps_counts_refLegenome.tsv",quote = FALSE, row.names = FALSE, sep = "\t")

# Filter low-confidence SNPs
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

vcf_fix_counts_merge_filt <- vcf_fix_counts_merge[,1:8]
i=9
while(i< ncol(vcf_fix_counts_merge)){
  headers <- vcf_fix_counts_merge[,i:(i+11)]
  total_counts <- headers[,grepl("totalCount",colnames(headers))]
  #Only kept SNPs present in all three replicates
  headers[rowSums(total_counts > 0) < 3,] <- 0
  #Remove SNPs where the sum of mapped reference and alternate allele counts was <90% of the total read depth in in at least one replicate.
  rawDepth<- headers[,grepl("rawDepth",colnames(headers))]
  ratio <- total_counts/rawDepth
  ratio[is.nan(ratio)] <- 0
  headers[rowSums(ratio < 0.9) > 0,] <- 0
  vcf_fix_counts_merge_filt <- cbind(vcf_fix_counts_merge_filt,headers)
  i = i+12
}

# Remove variants with no counts
vcf_fix_counts_merge_filt2 <- vcf_fix_counts_merge_filt[rowSums(vcf_fix_counts_merge_filt[,9:ncol(vcf_fix_counts_merge_filt)]) > 0,]

# Remove SNPs in the hybrid ASE dataset where the proportion of reference alleles in pure Landsberg erecta was <95%.
LeLe_ref_counts <- rowSums(vcf_fix_counts_merge_filt2[,grepl("refCount_LeLe",colnames(vcf_fix_counts_merge_filt2))])
LeLe_rawDepth <- rowSums(vcf_fix_counts_merge_filt2[,grepl("rawDepth_LeLe",colnames(vcf_fix_counts_merge_filt2))])
ratio_LeLe <- LeLe_ref_counts/LeLe_rawDepth
ratio_LeLe[is.na(ratio_LeLe)] <- 0

# Remove SNPs in the hybrid ASE dataset where the proportion of alternate alleles in pure Llagostera was <95%.
LlLl_ref_counts <- rowSums(vcf_fix_counts_merge_filt2[,grepl("altCount_LlLl",colnames(vcf_fix_counts_merge_filt2))])
LlLl_rawDepth <- rowSums(vcf_fix_counts_merge_filt2[,grepl("rawDepth_LlLl",colnames(vcf_fix_counts_merge_filt2))])
ratio_LlLl <- LlLl_ref_counts/LlLl_rawDepth
ratio_LlLl[is.na(ratio_LlLl)] <- 0

vcf_fix_counts_merge_filt3 <- vcf_fix_counts_merge_filt2[ratio_LeLe > 0.95,]
vcf_fix_counts_merge_filt3 <- vcf_fix_counts_merge_filt2[ratio_LlLl > 0.95,]

#Keep only SNPs with hybrid counts.
ASE_hybrid <- vcf_fix_counts_merge_filt3 %>% select(contains(c("LeLl","LlLe")))
vcf_fix_counts_merge_filt4 <- vcf_fix_counts_merge_filt3[rowSums(ASE_hybrid) > 0,]

#Export file with filtered counts
write.table(vcf_fix_counts_merge_filt4,"Results/filtered_discriminants_snps_counts_refLe_genome.tsv",quote = FALSE, row.names = FALSE, sep = "\t")


# Test DE in all samples

x <- data.frame(Sample=character(),
                DE_Genes=integer(),
                DE_Genes_UP_REF=integer(),
                DE_Genes_UP_ALT=integer())

for (sample in c("LeLlRH","LeLlRL","LeLlSH","LeLlSL","LlLeRH","LlLeRL","LlLeSH","LlLeSL")) {
  
  sample_data <- ASE_hybrid %>% select(contains(sample)) %>% select(contains(c("ref","alt")))
  rownames(sample_data) <- do.call(paste, c(vcf_fix_counts_merge_filt3[c("contig","position")], sep="_"))
  
  condition <- c("Le","Le","Le","Ll","Ll","Ll")
  coldata <- data.frame(row.names=colnames(sample_data), condition)
  
  # DE analysis
  dds <- DESeqDataSetFromMatrix(countData=sample_data,  colData=coldata, design=~condition)
  dds <- estimateSizeFactors(dds)
  
  dds_cond <- DESeq(dds, test="LRT", reduced= ~ 1)
  
  resLRT <- results(dds_cond,contrast=c("condition","Le","Ll"))
  res <- na.exclude(as.data.frame(resLRT))
  filter <- res[(res$padj<0.05 & abs(res$log2FoldChange) > 0.6),]
  
  write.table(filter,paste("../DeSeq/",sample,"_DESeq.out.tsv",sep="")
  
  #Number of DE variants
  nde <- nrow(filter)
  
  #Number of genes up in REF
  ref_up <- sum(filter$log2FoldChange > 0.6)
  
  #Number of genes up in REF
  alt_up <-sum(filter$log2FoldChange < -0.6)
  
  x[sample,] <- c(sample,nde,ref_up,alt_up)
  
}

write.table(x,"ASE_statistics_refLer_genome.tsv",quote=FALSE,sep="\t",row.names = TRUE)
