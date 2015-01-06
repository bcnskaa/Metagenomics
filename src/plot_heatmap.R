library(lattice); # If you do not have lattice library installed, simply run install.packages("lattice")
library(ggplot2);   # If you do not have ggplot2 library installed, simply run install.packages("ggplot2")
library(reshape2);

# Process combined
if(FALSE)
{
	source("plot_heatmap.R")
	cols <- c(rep('character', 5), 'numeric')
	mtx_fn = "bla.combined.mtx"

	bla_mtx <- read.table(mtx_fn, sep="\t", header=F, stringsAsFactors=F)
	
	bla_mtx$V4 <- round(bla_mtx$V4, 3)
	
	#ids = sort(unique(paste(bla_mtx$V2, bla_mtx$V3, sep=".")))
	ids = sort(unique(c(bla_mtx$V2, bla_mtx$V3)))
	
	mtx <- matrix(0.0, length(ids), length(ids))
	
	colnames(mtx) <- ids
	rownames(mtx) <- ids
	
	for(i in 1 : nrow(bla_mtx))
	{
		cid <- bla_mtx[i, 2]
		rid <- bla_mtx[i, 3]
		val <- bla_mtx[i, 4];
		
		val[which(val > 1)] <- 1.0 
		mtx[cid,rid] <- val;
		mtx[rid,cid] <- val;
	}
	
	plot_heatmap_mtx(mtx, out_fn="bla_mtx.combined.mtx.pdf",height=5, width=5, pdf_output=T, midpoint=0.5, colorbar_scheme=c("steelblue", "yellow", "darkred"), plot_label=T, label_size=3)

	
	library(ape)
	out_fn="bla_mtx.combined.tree.pdf"
	pdf(out_fn, height=5, width=5)
	plot.phylo(as.phylo(hclust(dist(mtx))), type="fan", cex=0.8, use.edge.length = TRUE)
	dev.off()
	
	pie_chart_data <- data.frame(bla_mtx$V1, )

}


if(FALSE)
{
	pdf("contig_info.overview.pdf", width=10, height=6)
	#ggplot(contig_info2, aes(x=V2.y, y=V3.y, colour=V3.x)) + geom_point(aes(size=(V2.x), cex=2)) + geom_text(aes(label=V8), cex=3, hjust=-0.1,vjust=0) + xlim(10, 200) + ylim(20000,240000) + theme(panel.grid = element_blank(), panel.background = element_blank()) + scale_size_area()
	ggplot(contig_info2, aes(x=V2.y, y=V3.y)) + geom_point(aes(size=(V2.x), cex=2, alpha=0.8)) + xlim(10, 200) + ylim(20000,240000) + theme(panel.grid = element_blank(), panel.background = element_blank()) + scale_size_area() + xlab("Number of Contig") + ylab("Mean Contig Length (bp)") # + geom_text(aes(label=V8), cex=3, hjust=-0.1,vjust=0)
	
	dev.off()
	
	pdf("contig_info.zoom.pdf", width=10, height=6)
	ggplot(contig_info2, aes(x=V2.y, y=V3.y, colour=V3.x)) + geom_point(aes(size=(V2.x), cex=2)) + geom_text(aes(label=V8), cex=3, hjust=-0.1,vjust=0) + xlim(10, 75) + ylim(60000,240000) + theme(panel.grid = element_blank(), panel.background = element_blank()) + scale_size_area()
	dev.off()
}


if(FALSE)
{
	source("workspace/Metagenomics/src/plot_heatmap.R")
	mtx_fn <- "score_mtx.2000"
	bla_mtx <- read.table(mtx_fn, sep="\t", header=T, stringsAsFactors=F)
	labels <- bla_mtx$Label
	plot_heatmap_mtx(bla_mtx, x_axis_labels=labels, y_axis_labels=labels, title="Unweighted Score of Genomic Similarity", out_fn=paste(mtx_fn, ".mtx", sep=""), height=5, width=5, pdf_output=T, midpoint=0.5, colorbar_scheme=c("steelblue", "yellow", "red"), plot_label=T, label_size=3, value_decimal_len=4)
	
	mtx_fn <- "weighted_score_mtx.2000"
	bla_mtx <- read.table(mtx_fn, sep="\t", header=T, stringsAsFactors=F)
	labels <- bla_mtx$Label
	plot_heatmap_mtx(bla_mtx, x_axis_labels=labels, y_axis_labels=labels, title="Weighted Score of Genomic Similarity", out_fn=paste(mtx_fn, ".mtx", sep=""), height=5, width=5, pdf_output=T, midpoint=0.5, colorbar_scheme=c("steelblue", "yellow", "red"), plot_label=T, label_size=3, value_decimal_len=4)
	
	
}


if(FALSE)
{
source("~/Desktop/workspace/Metagenomics/src/plot_heatmap.R")
cols <- c(rep('character', 5), 'numeric')
mtx_fn = "bla.mtx"
bla_mtx <- read.table(mtx_fn, sep="\t", header=F, stringsAsFactors=F, colClasses = cols)

bin_val_threshold <- 8
bin_1_vals <- as.integer(bla_mtx$V3);
bin_2_vals <- as.integer(bla_mtx$V5);
selected_idx <- which(bin_1_vals < bin_val_threshold & bin_2_vals < bin_val_threshold)

bla_mtx <- bla_mtx[selected_idx,]

ids = sort(unique(paste(bla_mtx$V2, bla_mtx$V3, sep=".")))

mtx <- matrix(0.0, length(ids), length(ids))
colnames(mtx) <- ids
rownames(mtx) <- ids

for(i in 1 : nrow(bla_mtx))
{
	
	cid <- paste(bla_mtx[i, 2], bla_mtx[i, 3], sep=".");
	rid <- paste(bla_mtx[i, 4], bla_mtx[i, 5], sep=".");
	val <- bla_mtx[i, 6];
	if(val > 1){
		val <- 1;
	}
	
	mtx[cid,rid] <- 1 - val;
	mtx[rid,cid] <- 1 - val;
}
plot_heatmap_mtx(mtx, out_fn="All",height=25, width=25, pdf_output=T, midpoint=0.5, colorbar_scheme=c("blue", "orange", "red"))

library(ape)
plot(as.phylo(hclust(dist(mtx))), type="fan", cex=0.5, use.edge.length = TRUE)

selected_group_id <- "SWH-Cell"
smtx <- mtx[grep(selected_group_id, colnames(mtx)), grep(selected_group_id, colnames(mtx))]
plot_heatmap_mtx(smtx, out_fn=paste(selected_group_id, ".mtx", sep=""),height=10, width=10, pdf_output=T, midpoint=0.5, colorbar_scheme=c("blue", "orange", "red"))

library(ape)
pdf(paste(selected_group_id, ".cluster.pdf", sep=""))
plot(hclust(dist(smtx)))
dev.off()


# Specify the row and column
selected_row_id <- "GZ-Cell_Y1"
selected_col_id <- "GZ-Cell_Y0"
smtx <- mtx[grep(selected_row_id, colnames(mtx)), grep(selected_col_id, colnames(mtx))]
plot_heatmap_mtx(smtx, out_fn=paste(selected_row_id, ".", selected_col_id, ".mtx", sep=""),height=10, width=10, pdf_output=T, midpoint=0.5, colorbar_scheme=c("blue", "orange", "red"))



}


# To use it,
plot_heatmap_mtx <- function(mtx, x_axis_labels=character(0), y_axis_labels=character(0), out_fn="output", export_to_file=T, height=2.8, width=3.5, colorbar_scheme=c("red", "yellow", "green"), midpoint=0, colorbar_witdh = 8.3, title=character(0), xtitle=character(0),  ytitle=character(0), legend_title=character(0), font="Courier", delim="\t", pdf_output=F, range_limit=character(0), plot_label=F, label_size=1, value_decimal_len=2)
{
	library(ggplot2)
	library(reshape2)
	#legend_position = c(1.0, 0.0);
	barwitdh = 10;
	image_scale = 100;
	pdf_scale = 1.4;
	
	
	#mtx <- read.table(mtx_fn, header=T, sep=delim, quote="", na.strings="", stringsAsFactors=F, strip.white=TRUE);
	
	plot_mtx <- melt(mtx);
	#x_label <- "variable";
	#y_label <- "label";
	
	colnames(plot_mtx) <- c("x_label", "y_label", "value");
	plot_mtx$value_str <- format(round(plot_mtx$value, value_decimal_len), nsmall = value_decimal_len)

	if(length(x_axis_labels) > 0)
	{
		plot_mtx$x_label <- factor(plot_mtx$x_label, levels=x_axis_labels, labels=x_axis_labels)
	}
	
	if(length(y_axis_labels) > 0)
	{
		plot_mtx$y_label <- factor(plot_mtx$y_label, labels=y_axis_labels)
	}
	
#	if(length(range_limit) == 0)
#	{
#		max_guide_range <- max(abs(min(plot_mtx$value)), plot_mtx$value);
#		range_limit = c(-1 * max_guide_range, max_guide_range);
#	}

	#g <- ggplot(plot_mtx, aes(x=Var1, y=Var2, fill=value)) + geom_tile(aes(height=0.97, width=0.97)) +
	g <- ggplot(plot_mtx, aes(x=y_label, y=x_label, fill=value)) + geom_tile(aes(height=0.97, width=0.97)) +
			theme(panel.background=element_blank(), axis.ticks=element_blank()) +
			theme(legend.position="bottom", axis.text=element_text(family="Courier")) +
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			#scale_fill_gradient2(name="", low=colorbar_scheme[1], mid=colorbar_scheme[2], high=colorbar_scheme[3], limits=range_limit) +
			scale_fill_gradient2(name="", low=colorbar_scheme[1], mid=colorbar_scheme[2], high=colorbar_scheme[3], midpoint=midpoint) +
			guides(fill = guide_colorbar(barwidth=colorbar_witdh, title.position = "bottom", direction = "horizontal")) 
	
#	if(length(x_axis_labels) > 0)
#	{
#		g <- g + xlim(levels(x_axis_labels))
#	}
#	if(length(y_axis_labels) > 0)
#	{
#		g <- g + ylim(levels(y_axis_labels))
#	}
	
	if(plot_label)
	{
		g <- g + geom_text(aes(label=value_str), size=label_size)
	}
	
	
	if(length(title) > 0)
	{
		g <- g + ggtitle(title);
	}
		
		
	if(length(xtitle) > 0)
	{
		g <- g + labs(x=xtitle);
	} else {
		g <- g + theme(axis.title.x=element_blank())
	}
	
	if(length(ytitle) > 0)
	{
		g <- g + labs(y=ytitle);
	} else {
		g <- g + theme(axis.title.y=element_blank())
	}
	
#	if(length(legend_title) > 0)
#	{
#		g <- g + labs(fill=legend_title);
#		g <- g + scale_fill_gradient(name=legend_title, colours=color_palette)
#	} else {
#		g <- g + theme(legend.title=element_blank())
#	}
	
	
	if(export_to_file)
	{
		if(pdf_output) {
			pdf(paste(out_fn, ".pdf", sep=""), height * pdf_scale, width * pdf_scale);
		} else {
			png(paste(out_fn, ".png", sep=""), height * image_scale, width * image_scale);
		}
	}

	if(export_to_file)
	{
		print(g, yscale.components=NULL);
		dev.off();
	}
	
	return(g);
}	
	
	

# To use it,
plot_heatmap_2 <- function(mtx_fn, export_to_file=T, height=2.8, width=3.5, colorbar_scheme=c("red", "yellow", "green"), midpoint=0, colorbar_witdh = 8.3, title=character(0), xtitle=character(0),  ytitle=character(0), legend_title=character(0), font="Courier", delim="\t", pdf_output=F, range_limit=character(0))	
{
	#legend_position = c(1.0, 0.0);
	barwitdh = 10;
	image_scale = 100;
	pdf_scale = 1.4;
	
	mtx <- read.table(mtx_fn, header=T, sep=delim, quote="", na.strings="", stringsAsFactors=F, strip.white=TRUE);
	
	plot_mtx <- melt(mtx);
	colnames(plot_mtx) <- c("label", "variable", "value");
	
#	if(length(range_limit) == 0)
#	{
#		max_guide_range <- max(abs(min(plot_mtx$value)), plot_mtx$value);
#		range_limit = c(-1 * max_guide_range, max_guide_range);
#	}
	
	g <- ggplot(plot_mtx, aes(x=variable, y=label, fill=value)) + geom_tile(aes(height=0.97, width=0.97)) +
			theme(panel.background=element_blank(), axis.ticks=element_blank()) +
			theme(legend.position="bottom", axis.text=element_text(family="Courier")) +
			#scale_fill_gradient2(name="", low=colorbar_scheme[1], mid=colorbar_scheme[2], high=colorbar_scheme[3], limits=range_limit) +
			scale_fill_gradient2(name="", low=colorbar_scheme[1], mid=colorbar_scheme[2], high=colorbar_scheme[3], midpoint=midpoint) +
			guides(fill = guide_colorbar(barwidth=colorbar_witdh, title.position = "bottom", direction = "horizontal")) 
	
	if(length(xtitle) > 0)
	{
		g <- g + labs(x=xtitle);
	} else {
		g <- g + theme(axis.title.x=element_blank())
	}
	
	if(length(ytitle) > 0)
	{
		g <- g + labs(y=ytitle);
	} else {
		g <- g + theme(axis.title.y=element_blank())
	}
	
#	if(length(legend_title) > 0)
#	{
#		g <- g + labs(fill=legend_title);
#		g <- g + scale_fill_gradient(name=legend_title, colours=color_palette)
#	} else {
#		g <- g + theme(legend.title=element_blank())
#	}
	
	
	if(export_to_file)
	{
		if(pdf_output) {
			pdf(paste(mtx_fn, ".pdf", sep=""), height * pdf_scale, width * pdf_scale);
		} else {
			png(paste(mtx_fn, ".png", sep=""), height * image_scale, width * image_scale);
		}
	}
	
	if(export_to_file)
	{
		print(g, yscale.components=NULL);
		dev.off();
	}
	
	return(g);
}	




# To use it,
plot_heatmap <- function(mtx_fn, export_to_file=T, height=500, weight=350, font="Courier", delim=",", pdf_output=F)
{
	mtx <- read.table(mtx_fn, header=T, sep=delim, quote="", na.strings="", stringsAsFactors=F);
	
	mtx$selected[which(nchar(mtx$selected) == 0)] <- " "
	mtx[is.na(mtx$selected),1] <- " ";
	
	max_id_len <- max(nchar(mtx$id))
	ids <- mtx$id;
	for(i in 1 : length(ids))
	{
		id <- ids[i]
		if(nchar(id) != max_id_len)
			ids[i] <- paste(id, rep(" ", max_id_len - nchar(id)), sep="");
	}
	
	ylabels <- do.call("expression", lapply(1:nrow(mtx), function(i) substitute(X ~ italic(Y), list(X=paste(mtx$selected[i], mtx$name[i], " "), Y=ids[i]))))
	plot_mtx <- mtx[,4:6];
	rownames(plot_mtx) <- paste(mtx$selected, " ", mtx$name, mtx$id, sep="");
	plot_mtx <- as.matrix(plot_mtx);
	
	max_range <- max(abs(min(plot_mtx)), plot_mtx);
	
	rgb.palette <- colorRampPalette(c("red", "white", "blue"), space = "rgb");
	
	at = seq(-max_range, max_range, length = 100)

	if(export_to_file)
	{
		if(pdf_output) {
			pdf(paste(mtx_fn, ".pdf", sep=""), 10, 5);
		} else {
			png(paste(mtx_fn, ".png", sep=""), height, weight);
		}
	}
	
	g <- levelplot(plot_mtx, regions=F, col.regions=rgb.palette, border="white", outer=FALSE, scale=list(x=list(labels=ylabels, rot=90, cex=1), y=list(rot=90,alternating=1, cex=1.1)), xlab=NULL, ylab=NULL, at=seq(-1 * max_range, max_range, length=100), par.settings=list(axis.text=list(fontfamily=font)), colorkey=list(labels=list(rot=90, las=0)));
	
	if(export_to_file)
	{
		print(g, yscale.components=NULL);
		dev.off();
	}
	return(g);
}



list_font <- function()
{
	print(names(pdfFonts()));
}


test <- data.frame(x=c(1,1,2,2), y=c(1,2,1,2), w=c(0.2,0.1,0.4,0.3), h=c(0.2,0.1,0.4,0.3))
test <- data.frame(x=c("A","A","B","B"), y=c("C","D","C","D"), w=c(0.2,0.1,0.4,0.3), h=c(0.2,0.1,0.4,0.3))

library(ggplot2)
ggplot(test, aes(x=x, y=y), xlim=c(0,2), ylim=c(0,2)) + geom_tile(aes(width=w, height=h))


