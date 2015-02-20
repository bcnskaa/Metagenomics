
#
# selected_sample_ids <- vector(mode="list", length=3)
# names(selected_sample_ids) <- c("Xylan", "Seed", "Cell")
# selected_sample_ids$Xylan <- c("GZ-Xyl_Y1", "GZ-Xyl_Y2")
# selected_sample_ids$Cell <- c("GZ-Cell_Y1", "GZ-Cell_Y2")
# selected_sample_ids$Seed <- c("GZ-See_Y0")
# grouping 
#
import_merged_tbl <- function(tbl_merged_fn, normalization=F, selected_sample_ids=character(0), selected_tax_ids=character(0), selected_gene_id=character(0), colorbar_scheme=c("steelblue", "yellow", "darkred"), pdf_fn=character(0))
{
	library(ggplot2)
	
	tbl <- read.table(tbl_merged_fn, sep="\t", header=F, stringsAsFactors=F);
	colnames(tbl) <- c("sample_id", "tax_id", "gene_id", "abundance", "count")
	group_ids <- unique(tbl$sample_id)
	
	if(length(selected_sample_ids) > 0)
	{
		patterns = paste(selected_sample_ids, collapse="|", sep="")
		tbl <- tbl[which(grepl(patterns, tbl$sample_id)), ]
	}
	tbl$group_id <- paste(tbl$sample_id, tbl$tax_id, sep="@")	
	
	
#	if(length(selected_sample_ids) > 0)	
#	{
#		tbl$group_id <- rep("", nrow(tbl))
#		for(group_id in names(selected_sample_ids))
#		{
#			print(paste("Processing ", group_id, sep=""))
#			sample_id_to_be_grouped <- selected_sample_ids[group_id]
#			tbl$group_id[which(tbl$sample_id %in% as.character(unlist(sample_id_to_be_grouped)))] <- group_id
#		}
#		tbl <- tbl[which(tbl$group_id != ""), ]
#		
#		library(plyr)
#		merge_tbl <- ddply(tbl, .(group_id, tax_id), abundance=sum(abundance), count=sum(count))
#		
#	} else {
#		tbl$group_id <- paste(tbl$sample_id, tbl$tax_id, sep="@")	
#	}

	
	tbl <- tbl[which(tbl$gene_id %in% tbl$gene_id[grepl("GH", tbl$gene_id)]), ]
	
	
	
	
	
	nr_gene_ids <- as.character(unique(tbl$gene_id))
	excluded_gene_ids <- rep(F, length(nr_gene_ids))
	count_cutoff <- 3
	#count_cutoff <- 10
	
	for(i in 1 : length(nr_gene_ids))
	{
		nr_gene_id <- nr_gene_ids[i]
		if(sum(tbl$count[which(tbl$gene_id == nr_gene_id)]) < count_cutoff)
		{
			excluded_gene_ids[i] <- T
		}
	}
	
	tbl <- tbl[which(tbl$gene_id %in% nr_gene_ids[!excluded_gene_ids]), ]

	g <- ggplot(data=tbl, aes(x=gene_id, y=group_id, fill=abundance)) + geom_tile() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
	scale_fill_gradient2(name="", low=colorbar_scheme[1], mid=colorbar_scheme[2], high=colorbar_scheme[3], midpoint=5000)
	# (max(tbl$count) - min(tbl$count) / 2)
	if(length(pdf_fn) > 0)
	{
		pdf(pdf_fn, width=15, height=8)
		print(g);
		dev.off();
	}
	g;
	
}

#
# tbl_merged_fn = "dbCAN.bin.tbl.consolidated.melt.abund_merged"
# import_merged_tbl(tbl_merged_fn, pdf_fn="all.pdf", selected_sample_ids=c("SWH-Seed_Y0", "SWH-Xyl_Y1", "SWH-Xyl_Y2", "SWH-Cell_Y1", "SWH-Cell_Y2"))
# import_merged_tbl(tbl_merged_fn, pdf_fn="seed.pdf", selected_sample_ids=c("SWH-Seed_Y0", "GZ-Seed_Y0"))
# import_merged_tbl(tbl_merged_fn, pdf_fn="xyl.pdf", selected_sample_ids=c("SWH-Xyl_Y2", "GZ-Xyl_Y2"))
# import_merged_tbl(tbl_merged_fn, pdf_fn="cell.pdf", selected_sample_ids=c("SWH-Cell_Y2", "GZ-Cell_Y2"))
plot_gradient_tile <- function(gradient_mtx, weight_mtx, x_labels=character(0), y_labels=character(0), out_fn=character(0))
{
	
}