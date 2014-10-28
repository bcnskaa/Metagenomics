library(lattice); # If you do not have lattice library installed, simply run install.packages("lattice")
library(ggplot2);   # If you do not have ggplot2 library installed, simply run install.packages("ggplot2")
library(reshape2);



# To use it,
plot_heatmap_2 <- function(mtx_fn, export_to_file=T, height=2.8, width=3.5, colorbar_scheme=c("red", "yellow", "green"), midpoint=0, colorbar_witdh = 8.3, xtitle=character(0),  ytitle=character(0), legend_title=character(0), font="Courier", delim="\t", pdf_output=F, range_limit=character(0))
{
	#legend_position = c(1.0, 0.0);
	barwitdh = 10;
	image_scale = 100;
	pdf_scale = 1.4;
	
	mtx <- read.table(mtx_fn, header=T, sep=delim, quote="", na.strings="", stringsAsFactors=F, strip.white=TRUE);
	
	plot_mtx <- melt(mtx);
	
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

