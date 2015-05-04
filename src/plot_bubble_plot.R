library(lattice); # If you do not have lattice library installed, simply run install.packages("lattice")
library(ggplot2);   # If you do not have ggplot2 library installed, simply run install.packages("ggplot2")
library(reshape2);

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
calculate_bray_curtis <- function(vegan_biom, normalization=T) {
	require(vegan)
	
	if(normalization)
	{
		normal_vegan_biom <- decostand(vegan_biom, "total")
	} else {
		normal_vegan_biom <- vegan_biom
	}
	
	bc_dist <- vegdist(normal_vegan_biom, "bray")
	
	return(bc_dist)
}



###
filter_phyloseq_by_read_count <- function(phyloseq_biom)
{
	
}



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
	return(as(OTU, "matrix"))
}



initial <- function()
{
	require(phyloseq)
	
}
