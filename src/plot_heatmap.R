library(lattice); # If you do not have lattice library installed, simply run install.packages("lattice")
library(ggplot2);   # If you do not have ggplot2 library installed, simply run install.packages("ggplot2")
library(reshape2);



# To use it,
plot_heatmap <- function(mtx_fn, export_to_file=T, height=500, weight=350, tile_width=0.97, tile_height=0.97, xtitle=character(0),  ytitle=character(0), legend_title=character(0), font="Courier", delim="\t", pdf_output=F, color_palette=0, symmetric_guide=T)
{
	mtx <- read.table(mtx_fn, header=T, sep=delim, quote="", na.strings="", stringsAsFactors=F, strip.white=TRUE);
	
	plot_mtx <- melt(mtx);
	
	tileHeight <- 0.97;
	tileWidth <- 0.97;
	
	if(symmetric_guide)
	{	
		max_guide_range <- max(abs(min(plot_mtx$value)), plot_mtx$value);
	}
	
	color_palette <- colorRampPalette(c("red", "white", "blue"))(n = 20);

	g <- ggplot(plot_mtx, aes(x=variable, y=label, fill=value)) + geom_tile(aes(height=0.97, width=0.97)) +
			theme(panel.background=element_blank(), axis.ticks=element_blank()) +
			scale_fill_gradientn(colours=color_palette);
	
	if(length(xtitle) > 0)
	{
		g <- g + labs(x=xtitle) + theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust=1, family="Courier"));
	} else {
		g <- g + theme(axis.title.x=element_blank())
	}
	
	
	if(length(ytitle) > 0)
	{
		g <- g + labs(y=ytitle);
	} else {
		g <- g + theme(axis.title.y=element_blank())
	}
	
	if(length(legend_title) > 0)
	{
		g <- g + labs(fill=legend_title);
		g <- g + scale_fill_gradient(name=legend_title, colours=color_palette)
	} else {
		g <- g + theme(legend.title=element_blank())
	}
	
	
	if(export_to_file)
	{
		if(pdf_output) {
			pdf(paste(mtx_fn, ".pdf", sep=""), 10, 5);
		} else {
			png(paste(mtx_fn, ".png", sep=""), height, weight);
		}
	}
	
	#g <- levelplot(plot_mtx, regions=F, col.regions=rgb.palette, border="white", outer=FALSE, scale=list(x=list(labels=ylabels, rot=90, cex=1), y=list(rot=90,alternating=1, cex=1.1)), xlab=NULL, ylab=NULL, at=seq(-1 * max_range, max_range, length=100), par.settings=list(axis.text=list(fontfamily=font)), colorkey=list(labels=list(rot=90, las=0)));
	
	if(export_to_file)
	{
		print(g, yscale.components=NULL);
		dev.off();
	}
	
	return(g);
}	
	
	


#
## To use it,
#plot_heatmap <- function(mtx_fn, export_to_file=T, height=500, weight=350, font="Courier", delim=",", pdf_output=F)
#{
#	
#	mtx <- read.table(mtx_fn, header=T, sep=delim, quote="", na.strings="", stringsAsFactors=F);
#	
#	mtx$selected[which(nchar(mtx$selected) == 0)] <- " "
#	mtx[is.na(mtx$selected),1] <- " ";
#	
#	
#	max_id_len <- max(nchar(mtx$id))
#	ids <- mtx$id;
#	for(i in 1 : length(ids))
#	{
#		id <- ids[i]
#		if(nchar(id) != max_id_len)
#			ids[i] <- paste(id, rep(" ", max_id_len - nchar(id)), sep="");
#	}
#	
#	ylabels <- do.call("expression", lapply(1:nrow(mtx), function(i) substitute(X ~ italic(Y), list(X=paste(mtx$selected[i], mtx$name[i], " "), Y=ids[i]))))
#	plot_mtx <- mtx[,4:6];
#	rownames(plot_mtx) <- paste(mtx$selected, " ", mtx$name, mtx$id, sep="");
#	plot_mtx <- as.matrix(plot_mtx);
#	
#	max_range <- max(abs(min(plot_mtx)), plot_mtx);
#	
#	rgb.palette <- colorRampPalette(c("red", "white", "blue"), space = "rgb");
#	
#	at = seq(-max_range, max_range, length = 100)
#
#	if(export_to_file)
#	{
#		if(pdf_output) {
#			pdf(paste(mtx_fn, ".pdf", sep=""), 10, 5);
#		} else {
#			png(paste(mtx_fn, ".png", sep=""), height, weight);
#		}
#	}
#	
#	#g <- levelplot(plot_mtx, regions=F, col.regions=rgb.palette, border="white", outer=FALSE, scale=list(x=list(labels=ylabels, rot=90, cex=1), y=list(rot=90,alternating=1, cex=1.1)), xlab=NULL, ylab=NULL, at=seq(-1 * max_range, max_range, length=100), par.settings=list(axis.text=list(fontfamily=font)), colorkey=list(labels=list(rot=90, las=0)));
#	
#	if(export_to_file)
#	{
#		print(g, yscale.components=NULL);
#		dev.off();
#	}
#	return(g);
#}

list_font <- function()
{
	print(names(pdfFonts()));
}

