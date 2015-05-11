library(lattice); # If you do not have lattice library installed, simply run install.packages("lattice")
library(ggplot2);   # If you do not have ggplot2 library installed, simply run install.packages("ggplot2")
library(reshape2);


###
# To install phyloseq package, the following dependent packages need to be installed before head.
# https://github.com/joey711/phyloseq/issues/329
# 
# And may require an update-to-date BioGeneric (version >=0.14) 
# wget http://www.bioconductor.org/packages/release/bioc/src/contrib/BiocGenerics_0.14.0.tar.gz
# install.packages("BiocGenerics_0.14.0.tar.gz", repos = NULL)
#
# library("devtools")
# install_github("phyloseq", "joey711")
# 
# If the following error occurs: 
# >Error in system(full, intern = quiet, ignore.stderr = quiet, ...) : 
# >  error in running command
#
# Setting options(unzip="internal") can solve the error
###



###
if(FALSE) {
 df <- read.table("test.dat", sep="\t", header=T)
 #x_label="Sample";y_label="Code";value_label="Coverage";excluding_zero=TRUE; xlabels=character(0); ylabels=character(0); xtitle=character(0); ytitle=character(0);pdf_fn=character(0)
 plot_bubble(df, "Sample", "Code", "Coverage")
}

###
plot_bubble <- function(df, x_label, y_label, value_label, excluding_zero=TRUE, xlabels=character(0), ylabels=character(0), xtitle=character(0), ytitle=character(0), pdf_fn=character(0))
{
	require(reshape2)
	require(ggplot2)
	
	plot_df <- melt(df, ids=c(y_label), value.name=value_label)

	#plot_df <- plot_df[which(plot_df$variable == value_label), c(x_label, y_label, value_label)]
	colnames(plot_df)[(colnames(plot_df) == "variable")] <- x_label
	
	plot_df <- plot_df[, c(y_label, x_label, value_label)]
	
	print(paste(ncol(plot_df), ": ", nrow(plot_df), sep=""))
	

	if(excluding_zero)
	{
		plot_df <- plot_df[which(plot_df[,value_label] > 0), ]
	}

	
	#gg <- ggplot(plot_df, aes(x=plot_df[,x_label], y=plot_df[,y_label], size=plot_df[,value_label]),guide=FALSE) +
	gg <- ggplot(plot_df, aes_string(x=x_label, y=y_label, size=value_label),guide=FALSE) +
		geom_point(colour="white", fill="black", shape=21) 


	if(length(xtitle) > 0)
	{
		gg <- gg + scale_x_discrete(name=xtitle)	
	} else {
		gg <- gg + scale_x_discrete(name="")
	}

	
	if(length(ytitle) > 0)
	{
		gg <- gg + scale_y_discrete(name=ytitle)	
	} else {
		gg <- gg + scale_y_discrete(name="")
	}

	gg <- gg + theme_bw() + theme(axis.text.x = element_text(angle = -45, hjust = 0)) 
	print(gg)
}


### 
calculate_bray_curtis <- function(vegan_biom, normalization=T, dist_method="bray") {
	require(vegan)
	require(phyloseq)
	
	if(normalization)
	{
		normalized_vegan_biom <- decostand(vegan_biom, "total")
	} else {
		normalized_vegan_biom <- vegan_biom
	}
	
	vegan_biom.dist <- vegdist(normalized_vegan_biom, method=dist_method)
	
	return(vegan_biom.dist)
}



calculate_pcoa <- function(vegan_biom, dist_method="bray") {
	require(vegan)
	require(phyloseq)
	require(ggplot2)
	
	biom.dist <- calculate_bray_curtis(vegan_biom, dist_method=dist_method)
	biom.pcoa <- capscale(biom.dist ~ 1, distance=dist_method)
	
	# Variance explained: MDS1, MDS2
	mds1_ve_1 <- summary(biom.pcoa)$cont$importance[2,1] * 100
	mds1_ve_2 <- summary(biom.pcoa)$cont$importance[2,2] * 100
	
	# Coordinates of the samples in PCoA
	sample_coords <- as.data.frame(scores(biom.pcoa)$sites)
	
	# Plot the PCoA
	ggplot(sample_coords, aes_string(x="MDS1", y="MDS2")) + geom_point() + theme_bw() + xlab(paste("PCoA_1 (variance explained=",mds1_ve_1,"%",")",sep="" )) + ylab(paste("PCoA_2 (variance explained=",mds1_ve_2,"%",")",sep="" )) + geom_text(aes(label=rownames(sample_coords)),hjust=-0.1, vjust=0.5, cex=3) + geom_vline(xintercept = 0, lty = "dotted") + geom_hline(yintercept = 0, lty = "dotted") 
	ggplot(sample_coords, aes_string(x="MDS1", y="MDS2")) + geom_point() + theme_bw() + xlab(paste("PCoA_1 (variance explained=",mds1_ve_1,"%",")",sep="" )) + ylab(paste("PCoA_2 (variance explained=",mds1_ve_2,"%",")",sep="" )) + geom_vline(xintercept = 0, lty = "dotted") + geom_hline(yintercept = 0, lty = "dotted") 
	
	
	return (biom.pcoa)
}


#
# 
#
filter_phyloseq_by_read_count <- function(phyloseq_biom)
{
	
}


#
#
#


# Obtained from http://joey711.github.io/phyloseq-demo/phyloseq-demo.html
# This function converts Phyloseq object into Vegan compatible data.
#
#
veganotu = function(phyloseq_biom) {
	require("vegan")
	OTU = otu_table(phyloseq_biom)
	if (taxa_are_rows(OTU)) {
		OTU = t(OTU)
	}
	return(as(OTU, "matrix"))p
}



initial <- function()
{
	require(phyloseq)
	require(rgl)
	
	#sample_id <- "S1-1B"
	sample_id <- "SWH-Cell_Y2"
	min_len <- 5000
	
	#cov_fn <- paste(sample_id, ".coverage.summary", sep="")
	cov_fn <- paste(sample_id, ".sorted.bam.coverage.summary", sep="")
	cov <- read.table(cov_fn, sep="\t", header=T, row.names=1, stringsAsFactors=F)
	#colnames(cov) <- c("coverage")
	
	#tetra_fn <- paste(sample_id, ".scaffold.fa.tetra_freq", sep="")
	tetra_fn <- paste(sample_id, ".fa.tetra_freq", sep="")
	tetra <- read.table(tetra_fn, sep="\t", header=T, row.names=1, stringsAsFactors=F)
	if(length(which(colnames(tetra) == "X")) == 1)
	{
		tetra <- tetra[, -which(colnames(tetra) == "X")]
	}
	
	lineage_fn <- paste(sample_id, "+nr.m8.lineage", sep="")
	lineages <- read.table(lineage_fn, sep="\t", header=T, row.names=1, stringsAsFactors=F)
	
	tetra.pcoa <- calculate_pcoa(tetra, dist_method="bray")
	tetra.coord <- as.data.frame(scores(tetra.pcoa)$sites)
	tetra.cov <- merge(tetra.coord, cov, by="row.names", all=TRUE)
	row.names(tetra.cov) <- tetra.cov[,1]
	tetra.cov <- tetra.cov[,-which(colnames(tetra.cov) == "Row.names")]
	
	tetra.cov_lineage_pre <- merge(tetra.cov, lineages, by="row.names", all=TRUE)
	
	# Update for min_len
	tetra.cov_lineage <- tetra.cov_lineage_pre[which(tetra.cov_lineage_pre$length > min_len),]
	tetra.cov_lineage$log_length <- log(tetra.cov_lineage$length)
	
	
	labels <- unique(tetra.cov_lineage$lineage)
	labels <- labels[-which(labels == "")]
	cols = rainbow(length(labels))
	
	tetra.cov_lineage$col <- rep("#00000000")
	for(i in 1 : length(labels))
	{
		label = labels[i]
		col <- cols[i]
		tetra.cov_lineage$col[which(tetra.cov_lineage$lineage == label)] <- col
	}
	
	col_alpha <- 1
	# Plot
	with(tetra.cov_lineage, plot3d(MDS1, MDS2, log(coverage), size=0, col=col, alpha=col_alpha, colkey=F, main=sample_id))

	# 
	#http://stackoverflow.com/questions/10341963/3d-scatterplot-in-r-using-rgl-plot3d-different-size-for-each-data-point
	for(i in 1:nrow(tetra.cov_lineage)) {
		points3d(tetra.cov_lineage$MDS1[i], tetra.cov_lineage$MDS2[i], log(tetra.cov_lineage$coverage[i]), size=(tetra.cov_lineage$length[i])/10000, col=tetra.cov_lineage$col[i], alpha=col_alpha)
	}

	text_threshold = 10
	#tetra.cov_lineage.text = tetra.cov_lineage[tetra.cov_lineage$coverage >= text_threshold,]
	
	tetra.cov_lineage.text = tetra.cov_lineage[nchar(tetra.cov_lineage$lineage) > 0,]
	text3d(tetra.cov_lineage.text$MDS1, tetra.cov_lineage.text$MDS2, log(tetra.cov_lineage.text$coverage), paste(tetra.cov_lineage.text$Row.names,";", tetra.cov_lineage.text$lineage, sep=""), cex=0.6, adj=c(1.3,1.0))
	
	
	# add legend
	legend("topright", legend = labels, pch = 16, col = cols, cex=1, inset=c(0.02))
	
	
#	# Plot
#	with(tetra.cov, plot3d(MDS1, MDS2, log(coverage), size=0))
#	
#	#http://stackoverflow.com/questions/10341963/3d-scatterplot-in-r-using-rgl-plot3d-different-size-for-each-data-point
#	for(i in seq_along(tetra.cov)) {
#		points3d(tetra.cov$MDS1[i], tetra.cov$MDS2[i], log(tetra.cov$coverage[i]), size=log(tetra.cov$length[i]))
#	}
#	
	
}
