# Manisha Munasinghe - Last Updated: 6/27/19
# Script used for:
#   - Differential Expression Analysis
# For More Details See:
#   - deseq2_personal_func.R
#   - 


.libPaths('/workdir/mam737/Rlibs')
#Set this to wherever you store your Rlibs#

required_packages <- c("geneplotter","ggplot2","plyr","LSD","DESeq2","gplots",
	"RColorBrewer","stringr","topGO","genefilter","biomaRt","dplyr","EDASeq",
	"fdrtool","reshape2","pheatmap","TxDb.Dmelanogaster.UCSC.dm6.ensGene","gridExtra","goseq")
lapply(required_packages,library,character.only=TRUE)

#Seed set for reproducibility
set.seed(1994)



## Relevant Directories ##
setwd('/workdir/arh223/brc_data/')
source('./deseq2_personal_func.R')
fastqDir <- file.path('./raw_files')
htseqDir <- file.path('./htseq_count_bam')
## Relevant Directories ##

# Load & Format Count Data 

fastq <- list.files(fastqDir, pattern = '*.fastq')
sampleID <- unlist(lapply(strsplit(fastq,"_"), function(x) x[1]))
condition <- unlist(lapply(strsplit(sampleID,'-'), function(x) x[1]))
replicate <- rep(c('1','2'),15)[-c(8,18)]
sex <- rep(c('F','F','M','M'),8)[-c(8,18,19,20)]
phenotype <- c(rep(c('SE'),4),rep(c('SE'),3),rep(c('LE'),8),rep(c('SE'),5),rep(c('LE'),8))
libraryName <- paste(condition,'_',replicate,'_',sex,sep="")
countFile <- list.files(htseqDir, pattern = '*.counts')

metadata <- data.frame(sampleID=sampleID, condition=condition, replicate=replicate, sex=sex,phenotype=phenotype,libraryName=libraryName, fastq=fastq, countFile=countFile)
sampleTable <- data.frame(sampleName = metadata$libraryName, fileName= metadata$countFile, sampleID = metadata$sampleID, condition = metadata$condition, sex=metadata$sex,replicate = metadata$replicate,phenotype=metadata$phenotype,fastq=metadata$fastq)

DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = htseqDir, design = ~ phenotype + sex + phenotype:sex)
print("Number of Genes from HTSEQ Count")
rowData(DESeq2Table)

## QC & Normalization of the Count Data ##

raw_counts <- counts(DESeq2Table) ##Raw Counts Matrix

# Remove rows where all values are zero
idx.nz <- apply(raw_counts[,-1],1,function(x) !all(x==0))
print("Number of genes with non-zero counts across samples")
sum(idx.nz)
DESeq2Table<- DESeq2Table[idx.nz,]


# Normalization Method

DESeq2Table <- estimateSizeFactors(DESeq2Table)
#sizeFactors(DESeq2Table)
norm_counts <- counts(DESeq2Table, normalized = TRUE) #<- normalized data

# Filter out genes with low counts
norm_counts2<-norm_counts[rowSums(norm_counts > 10) >=14,]
DESeq2Table<-DESeq2Table[rowSums(norm_counts > 10) >=14,]
dim(norm_counts2)
rld_blind <- rlog(DESeq2Table, blind=TRUE)

#If you want to visualize the impact of
# normalizing run the following chunk of code below

norm_plots <- generate_normalization_plots(raw_counts,norm_counts,norm_counts2)
pdf("./Rplots/final/normalization_plots.pdf")
for (i in seq(1, length(norm_plots), 1)) {
	grid.arrange(grobs=norm_plots[i],ncol=1,nrow=1)
}
dev.off()


#If you want to visualize the PCA
#  run the following chunk of code below

pdf('./Rplots/final/pca_plot.pdf')
generate_pca_plots(rld_blind)
dev.off()


#### RUN DESEQ2's Differential Expression Analysis ####

# Test for differences between the two phenotypes
# i.e differences between LE and SE
# design = ~ phenotype + sex vs ~ sex

dds_phenotype <- DESeq(DESeq2Table,test='LRT',reduced=~sex)
filtered_res_phenotype <- results(dds_phenotype,contrast=c('phenotype', 'SE', 'LE'), alpha=0.05)
sum(filtered_res_phenotype$padj < 0.05, na.rm=TRUE)
ordered_filtered_res_phenotype <- filtered_res_phenotype[order(filtered_res_phenotype$padj),]

#If you want to visual the top 10 hits
#  either via scatterplots or a boxplot
#  run the following chunk of code below 

pdf('./Rplots/final/phenotype_top10hits_counts.pdf')
generate_count_plots(dds_phenotype, ordered_filtered_res_phenotype,'scatter')
dev.off()

pdf('./Rplots/final/phenotype_top10hits_boxplots.pdf')
generate_count_plots(dds_phenotype, ordered_filtered_res_phenotype,'boxplot')
dev.off()

#If you want to generate a QQplot for the LRT of DE between phenotypes
#  run the following chunk of code below

pdf('./Rplots/final/phenotypeLRT_qq.pdf')
generate_qq_plots(filtered_res_phenotype)
dev.off()


#If you want to generate a heatmap of specific genes
#  follow the chunk of code below
#  (the specific chunk below is to generate it for specific Hsp genes)

specified_genes <- c("FBgn0013278","FBgn0001225","FBgn0001230","FBgn0001226","FBgn0013279")
specified_gene_names <- c('Hsp70Bc','Hsp70Bb','Hsp27','Hsp26','Hsp68')

pdf('./Rplots/final/hsp_heatmap.pdf')
generate_heatmap(specified_genes,specified_gene_names,rld_blind)
dev.off()


#If you want to overlay scatterplot w/ boxplot for a specific gene
#  follow the chunk of code below
#  
pdf('./Rplots/final/Hsp70Bc_counts.pdf')
generate_joint_count_plots(dds_phenotype,gene_id='FBgn0013279',gene_name="Hsp70Bc")
dev.off()

