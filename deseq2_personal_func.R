# Manisha Munasinghe - Last Updated: 6/27/19
# Script containing functions used for:
#   - Visualizing Differential Expression QC & Analysis
# For More Details See:
#   - deseq2_analysis.R
#   - 


# The following function generates either
#   a boxplot or a scatterplot visualizing
#   normalized counts + pseudocount of 0.5 for a list of genes
generate_count_plots <- function(dds_object,ordered_results,plot_type,var) {
	plotlist <- list()
	boxplotlist <- list()
	layout <- rbind(c(1,2),
                    c(1,2),
                    c(1,2),
                    c(1,2),
                    c(1,2),
                    c(1,2),
                    c(1,2),
                    c(1,2),
                    c(1,2),
                    c(1,2))
	for (i in head(rownames(ordered_results),n=10)) {
		d <- plotCounts(dds_object,gene=i,intgroup=c('phenotype','condition'),returnData=TRUE)
		count_plot <- ggplot(d,aes(x=phenotype,y=count,color=condition))+scale_y_log10(breaks=c(25,100,400))+geom_point(position=position_jitter(w=0.1,h=0),size=3) + labs(title=i,xlab='Phenotype',ylab='Counts')+
  theme(plot.title = element_text(hjust = 0.5))
		plotlist[[i]] <- count_plot
		box_plot <- ggplot(d,aes(x=phenotype,y=count)) + geom_boxplot()+scale_y_log10(breaks=c(25,100,400)) + labs(title=i,xlab='Phenotype',ylab='Counts')+theme(plot.title = element_text(hjust = 0.5))
		boxplotlist[[i]] <- box_plot
	}
	if (plot_type =='scatter') {
		marrangeGrob(grobs=plotlist,nrow=1,ncol=2,layout_matrix=layout)
	} else if (plot_type == 'boxplot') {
		marrangeGrob(grobs=boxplotlist,nrow=1,ncol=2,layout_matrix=layout)
	}
	do.call('grid.arrange',c(plotlist,ncol=2))
}


# The following function generates a
#   a boxplot overlaid with a scatterplot visualizing
#   normalized counts + pseudocount of 0.5 for a specific gene
generate_joint_count_plots <- function(dds_obj,gene_id,gene_name) {
	d <- plotCounts(dds_obj,gene=gene_id,intgroup=c('phenotype','condition'),returnData=TRUE)
	d$condition <- factor(d$condition,levels=c('26','45','83','85','595','639','850','859'))
	d$phenotype <- factor(d$phenotype,levels=c('SE','LE'))

	y_title <- expression(paste("Read Counts - ", italic(gene_name)))

	print(ggplot(d,aes(x=phenotype,y=count))+geom_boxplot(width=0.4)+scale_y_log10(breaks=c(25,100,400))+ 
		geom_point(aes(color=factor(condition)),position=position_jitter(w=0.1,h=0),size=3)+
		labs(x='Phenotype',y=y_title)+
		scale_x_discrete(labels=c("SE" = "Small Eyes", "LE" = "Large Eyes"))+
		scale_color_manual(name='Line',values=c("#F6E8C3","#DFC27D","#BF812D","#8C510A","#E6F5D0"
			,"#B8E186","#7FBC41","#4D9221")))

}


# The following function generates
#   2 different ways of visualizing the raw data
generate_density_functions <- function(DESeq2Table_obj) {
	jpeg('./Rplots/ecdf.jpg')
	multiecdf( counts(DESeq2Table_obj, normalized = T),xlab='Mean Counts', xlim = c(0,2000))
	dev.off()#

	jpeg('./Rplots/densities.jpg')
	multidensity( counts(DESeq2Table_obj, normalized = T),xlab='Mean Counts', xlim = c(0,2000))
	dev.off()	
}

# The following function generates
#   a series of histograms showing raw counts
#   and various normalization & filtering steps
generate_normalization_plots <- function(raw_counts_matrix,norm_counts_matrix,norm_counts2_matrix) {
	plotList <- list()

	pseudoRawCounts <- log2(raw_counts_matrix +1 )
	melt_raw_counts <- melt(raw_counts_matrix)
	melt_pseudoRawCounts <- melt(pseudoRawCounts)

	pseudoNormCounts <- log2(norm_counts_matrix +1 )
	melt_norm_counts <- melt(norm_counts_matrix)
	melt_pseudoNormCounts <- melt(pseudoNormCounts)

	pseudoNormCounts2 <- log2(norm_counts2_matrix + 1)
	melt_norm_counts2 <- melt(norm_counts2_matrix)
	melt_pseudoNormCounts2 <- melt(pseudoNormCounts2)	

	plotList[['a']] <- ggplot(melt_raw_counts, aes(x=value)) + geom_histogram(breaks=seq(0,400000,by=5000)) + ggtitle("Histogram of Raw HTSEQ Counts") + labs(x='HTSEQ Counts',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['b']] <- ggplot(melt_pseudoRawCounts, aes(x=value)) + geom_histogram(breaks=seq(0,20,by=1)) + ggtitle("Histogram of -log2(Raw HTSEQ Counts)") + labs(x='HTSEQ Counts',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['c']] <- ggplot(melt_norm_counts, aes(x=value)) + geom_histogram(breaks=seq(0,400000,by=10000)) + ggtitle("Histogram of Normalized HTSEQ Counts") + labs(x='Normalized HTSEQ Counts',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['d']] <- ggplot(melt_pseudoNormCounts, aes(x=value)) + geom_histogram(breaks=seq(0,20,by=1)) + ggtitle("Histogram of -log2(Normalized HTSEQ Counts)") + labs(x='Normalzied HTSEQ Counts',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['e']] <- ggplot(melt_pseudoRawCounts, aes(x=Var2,y=value,fill=Var2))+geom_boxplot()+ggtitle('Boxplot Distribution of Normalized Counts in Each Sample')+labs(x='Sample',y='log2(Normalized Count+1)')+scale_color_brewer(palette='Spectral')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position='none')
	plotList[['f']] <- ggplot(melt_norm_counts2, aes(x=value)) + geom_histogram(breaks=seq(0,400000,by=10000)) + ggtitle("Histogram of Normalized HTSEQ Counts: At least 14 samples have norm counts > 10") + labs(x='Normalized HTSEQ Counts ',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['g']] <- ggplot(melt_pseudoNormCounts2, aes(x=value)) + geom_histogram(breaks=seq(0,20,by=1)) + ggtitle("Histogram of -log2(Normalized HTSEQ Counts): At least 14 samples have norm counts >10") + labs(x='Normalzied HTSEQ Counts',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['h']] <- ggplot(melt_pseudoNormCounts2, aes(x=Var2,y=value,fill=Var2))+geom_boxplot()+ggtitle('Boxplot Distribution of Normalized Counts in Each Sample: At least 10 samples have norm counts >2 ')+labs(x='Sample',y='log2(Normalized Count+1)')+scale_color_brewer(palette='Spectral')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position='none')

	return(plotList)
}


# The following function generates
#   a PCA plot of top 500 most variable genes
generate_pca_plots <- function(rld_object) {
	
	pcaData <- plotPCA(rld_object, intgroup=c('phenotype','condition','replicate','sex'), returnData=TRUE)
	pcaData$condition <- factor(pcaData$condition,levels=c('26','45','83','85','595','639','850','859'))
	percentVar <- round(100*attr(pcaData,"percentVar"))
	print(ggplot(pcaData, aes(x=PC1,y=PC2,shape=sex,color=condition))+geom_point(size=3)+xlab(paste0("PC1: ",percentVar[1],"% Variance"))+ylab(paste0("PC2: ", percentVar[2],"% Variance"))+coord_fixed()+scale_shape_manual(name='Sex',values=c(16,17))+
  scale_color_manual(name='Line',values=c("#F6E8C3","#DFC27D","#BF812D","#8C510A","#E6F5D0","#B8E186","#7FBC41","#4D9221")))

}

# The following function generates
#   a QQ plot for LRT used for DE analysis
generate_qq_plots <- function(filtered_res_obj) {
	observed_pvals <- sort(-log10(filtered_res_obj$pvalue),decreasing = FALSE)
	expected_pvals <- (sort(runif(length(na.omit(filtered_res_obj$pvalue))), decreasing = FALSE))
	expected_pvals <- sort(-log10(expected_pvals),decreasing=FALSE)
	plot(expected_pvals, observed_pvals, main = "Phenotype QQplot for RNASeq ",xlab='Expected Pvals', ylab='Observed Pvals')
	abline( a= 0 , b = 1, col = "red")

}

# The following function generates
#   a heatmap for a given set of genes 
generate_heatmap <- function(vec_of_genes, vec_of_gene_names, rld_obj) {
	specified_rows <- match(vec_of_genes, row.names(rld_obj))
	specified_phenotype_mat <- assay(rld_obj)[specified_rows,]
	specified_phenotype_mat <- specified_phenotype_mat - rowMeans(specified_phenotype_mat)
	specified_phenotype_mat_col.order <- c("26_1_F","26_2_F" ,"26_1_M" ,"26_2_M",
		"45_1_F" ,"45_2_F" ,"45_1_M","83_1_F"  ,"85_1_F"  ,"85_2_F","85_1_M"  ,"85_2_M",
		"595_1_F","595_2_F" ,"595_1_M" ,"595_2_M","639_1_F" ,"639_2_F","639_1_M","639_2_M",
		"850_1_F","850_2_F","850_1_M" ,"850_2_M" ,"859_1_F","859_2_F" ,"859_1_M" ,"859_2_M")

	df <- as.data.frame(colData(rld_blind)[,c("condition","phenotype")])

	cond_color_vec <- c("#F6E8C3","#DFC27D","#BF812D","#8C510A","#E6F5D0","#B8E186","#7FBC41","#4D9221")
	names(cond_color_vec) <- c('26','45','83','85','595','639','850','859')
	colScale <- scale_colour_manual(name = "condition",values = cond_color_vec)

	color_list = list(phenotype = c(SE="#A65628",LE="#276419"), condition = cond_color_vec)
	breaksList = seq(-1.5, 1.5, length.out = 100)

	pheatmap(specified_phenotype_mat[,specified_phenotype_mat_col.order],breaks=breaksList, annotation_col=df,fontsize=9,fontsize_number=0.4*fontsize,cluster_cols=FALSE,cellheight=50,annotation_colors=color_list,
	labels_row = vec_of_gene_names,
	labels_col = c('Line 26','Line 26','Line 26','Line 26','Line 45','Line 45','Line 45',
		'Line 83','Line 85','Line 85','Line 85','Line 85','Line 595','Line 595','Line 595','Line 595',
		'Line 639','Line 639','Line 639','Line 639','Line 850','Line 850','Line 850','Line 850',
		'Line 859','Line 859','Line 859','Line 859'))

}


