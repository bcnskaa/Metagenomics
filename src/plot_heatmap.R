library(lattice);

# To use it, simply 
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

