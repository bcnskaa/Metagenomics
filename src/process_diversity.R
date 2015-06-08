library(fBasics)
library(reshape2)

script_home = "/home/bcnskaa/projects/Metagenomics/trunk/src/"
source(paste(script_home, "plot_heatmap.R", sep=""))
source(paste(script_home, "plot_bubble_plot.R", sep=""))



pfam_fn <- "/home/bcnskaa/projects/Metagenomics_WD/FunctionalDiversity/all_samples/Pfam/all_samples.renamed+Pfam.dom.tbl.summary"
trait_fn <- "samples-clustered-tax_id+go.traits"
tx_fn_s <- "samples-tax_id+go.samples"

tx_fn_g <- "samples-tax_id+go.samples.g.clustered"
tx_fn_g_meta <- "samples-tax_id+go.samples.g.clustered.meta"

tx_fn <- "samples-clustered-tax_id+go.samples"
#tx_fn <- "all_samples.tax.species.txt"
fd_fn <- "all_samples.seed.functions.txt"

log_scale <- FALSE
selected_level <- 30

wk_fn <- tx_fn_g

# Read counts assigned to SEED subsystems
fd <- read.table(wk_fn, sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)
if(wk_fn == tx_fn_g)
{
	fd_meta <- read.table(tx_fn_g_meta, sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)
}
fd_bak <- fd




if(wk_fn == "samples-clustered-tax_id+go.samples" ) {
	###################
	# Min-count
	min_count <- 5000
	species_counts <- colSums(fd)
	
	excluded_species <- names(which(species_counts < min_count))
	excluded_species <- c(excluded_species, "o__")
	
	fd <- fd[, which(! (colnames(fd) %in% excluded_species))]
	
	excluded_samples = c("all_samples+nr.renamed.m8")
	
	fd <- fd[which(! (rownames(fd) %in% excluded_samples)), ]
	
	rownames(fd) <- gsub("+nr.renamed.m8", "", rownames(fd), fixed=T)
	fd <- t(fd)
}



if(grep(".clustered", wk_fn)) {
	###################
	# Min-count
	min_count <- 5000
	species_counts <- colSums(fd)
	
	excluded_species <- names(which(species_counts < min_count))
	# Remove species without kingdom assignment
	excluded_species <- rbind(excluded_species, names(fd[grep("k__.p_", rownames(fd))]))
	excluded_species <- rbind(excluded_species, "X.1.k__.p__.c__.o__.f__.g__")
	
	excluded_fd <- fd[, which((colnames(fd) %in% excluded_species))]
	fd <- fd[, which(! (colnames(fd) %in% excluded_species))]
	
	
	excluded_samples = c("all_samples+nr.renamed.m8")
	
	fd <- fd[which(! (rownames(fd) %in% excluded_samples)), ]
	
	rownames(fd) <- gsub("+nr.m8", "", rownames(fd), fixed=T)
	fd <- t(fd)
}



# Remove NA
if(wk_fn == "samples-clustered-tax_id+go.traits") {
	###################
	# Min-count
	min_global_trait_count <- 1000
	min_species_trait_count <- 300
	
	trait_counts <- colSums(fd)
	length(which(trait_counts < min_global_trait_count))
	excluded_traits <- colnames(fd)[which(trait_counts < min_global_trait_count)]
	
	
	# Exclude taxa whose triat count is too small
	trait_counts <- rowSums(fd)
	excluded_species <- rownames(fd)[which(trait_counts < min_species_trait_count)]


	fd <- fd[which(! (rownames(fd) %in% excluded_species)), ! (colnames(fd) %in% excluded_traits)]
	fd <- fd[,- which(colnames(fd) == "NA.")]
	fd <- fd[- which(rowSums(fd) == 0),]	
	
	fd <- t(fd)
}


# Processing the pfam data
if(wk_fn == pfam_fn)
{
	min_global_trait_count <- 20
	min_species_trait_count <- 100
	
	desc <- fd[, "DESC"]
	fd <- fd[, which(!colnames(fd) %in% "DESC")]
	
	min_count <- 2
	
	fd <- t(fd[which(rowSums(fd) > min_count), ])
}



#
#sorted_subsystems <- sort(rowSums(fd_normalized), decreasing=T)
#
#selected_level <- min(selected_level, nrow(fd_normalized))
#selected_fd <- fd_normalized[names(sorted_subsystems)[1:selected_level],]
#

sorted_subsystems <- sort(rowSums(fd), decreasing=T)

selected_level <- min(selected_level, nrow(fd))
fd <- fd[names(sorted_subsystems)[1:selected_level],]




if(wk_fn == "all_samples.renamed+Pfam.dom.tbl.summary")
{
	fd_normalized <- fd
	
} else {
	# Normalize the data
	fd_totals <- colSums(fd)
	
	fd_normalized <- (as.data.frame(t(t(fd) / fd_totals))) * 100
	
	# http://www.statsblogs.com/2014/07/14/a-log-transformation-of-positive-and-negative-values/
	if(log_scale)
	{
		#fd_normalized <- log10(fd_normalized)
		fd_normalized <- log10(fd_normalized + 1)
	}
	
} 

selected_fd <- fd_normalized





#selected_fd <- fd_normalized[which(fd_normalized[, "GZ.Cell_Y1.nr"] > 0.5),]
#selected_fd <- fd_normalized[which(fd_normalized[, "GZ.Cell_Y1.nr"] > 0.5),]


# Do hierachical clustering, and extract the labels
library(vegan)
if(wk_fn == "samples-clustered-tax_id+go.traits")
{
	fd_dist <- vegdist(t(selected_fd), "bray")
} else{ 
	fd_dist <- vegdist(t(fd_normalized), "bray")
}

hc <- hclust(fd_dist)


pdf(paste(wk_fn, "hclust", "pdf", sep="."), w=14, h=5)
plot(hc)
dev.off()

sample_ids <- hc$labels[hc$order]
selected_fd <- selected_fd[sample_ids]

#selected_fd <- selected_fd[c("GZ.Seed_Y0.nr", "SWH.Seed_Y0.nr", "SWH.Cell55_Y2.nr", "GZ.Cell_Y1.nr", "GZ.Cell_Y2.nr","GZ.Xyl_Y1.nr", "GZ.Xyl_Y2.nr", "SWH.Cell_Y1.nr", "SWH.Cell_Y2.nr","SWH.Xyl_Y1.nr", "SWH.Xyl_Y2.nr")


pdf(paste(wk_fn,"mtx","pdf",sep="."), w=14, h=16)
plot_mtx <- as.matrix(selected_fd)
plot_heatmap_mtx(plot_mtx, plot_label=T, midpoint=(max(plot_mtx) - min(plot_mtx)) / 2, x_axis_labels=rev(rownames(selected_fd)), xtitle="Sample", ytitle="SEED Subsystem", colorbar_scheme=c("steelblue", "yellow", "darkred"), label_size=3, value_decimal_len=2)
dev.off()


dev.off()

selected_fd$subsystem <- rownames(selected_fd)
plot_df <- melt(selected_fd)
colnames(plot_df) <- c("system", "sample", "read_count_perc")

pdf(paste(wk_fn,"bubble","pdf",sep="."), w=14, h=16)
df <- plot_bubble(plot_df, "sample", "system", "read_count_perc", xtitle="Sample", ytitle="SEED System", ylabels=rev(rownames(selected_fd)))
dev.off()



# Group analysis
if(grep("sample", wk_fn) & grep("clustered", wk_fn))
{
	group_cell <- colnames(plot_mtx)[grep("Cell_", colnames(plot_mtx))]
	group_xyl <- colnames(plot_mtx)[grep("Xyl", colnames(plot_mtx))]
	group_seed <- colnames(plot_mtx)[grep("Seed", colnames(plot_mtx))]
	
	group_cell_avgs <- rowAvgs(plot_mtx[,group_cell])
	group_xyl_avgs <- rowAvgs(plot_mtx[,group_xyl])
	group_seed_avgs <- rowAvgs(plot_mtx[,group_seed])
	
	# Normal to zero
	group_cell_avgs <- group_cell_avgs - group_seed_avgs
	group_xyl_avgs <- group_xyl_avgs - group_seed_avgs
	group_seed_avgs <- group_seed_avgs - group_seed_avgs
			
	group_avgs <- data.frame(rownames(plot_mtx), group_cell_avgs, group_xyl_avgs, group_seed_avgs)
	colnames(group_avgs) <- c("name", "Cell", "Xyl", "Seed")
	
	melt_group_avgs <- melt(group_avgs)
	
	
	#melt_group_avgs$value <- sign(melt_group_avgs$value) * log2(abs(1 + melt_group_avgs$value))
	
	
	pdf(paste(wk_fn,"group_compare","pdf",sep="."), w=14, h=16)
	xyplot(name ~ value, data=melt_group_avgs, groups=melt_group_avgs$variable, col=c("steelblue", "red", "black"), cex=1, pch=19)
	dev.off()
}




# Read counts assigned to taxonomy
#tx <- read.table("all_samples.tax.genus.txt", sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)




# Do PCoA
library(rgl)

plot_mtx.pca <- princomp(plot_mtx)
pdf(paste(wk_fn,"pca_biplot","pdf",sep="."), w=8, h=16)
#biplot(plot_mtx.pca, scale=T, type = c("text", "points"))
biplot(plot_mtx.pca, scale=T)

dev.off()


plot3d(plot_mtx.pca$scores[,1:3])
text3d(plot_mtx.pca$scores[,1:3], texts=rownames(plot_mtx))
text3d(plot_mtx.pca$loading[,1:3], texts=rownames(plot_mtx.pca$loadings), col="red")
coords <- NULL
for(i in 1 : nrow(plot_mtx.pca$loadings))
{
	coords <- rbind(coords, rbind(c(0,0,0), plot_mtx.pca$loadings[i, 1:3]))
}
lines3d(coords, col="red", lwd=4)



library(FD)




# http://stats.stackexchange.com/questions/79688/using-anosim-or-permanova-to-abiotics-data
# http://www.researchgate.net/post/What_is_the_best_method_to_analyze_change_in_species_composition_over_three_years_in_different_sites
# 1. Calculate the Bray-Curtis dissimilarity matrix using species abundance
# 2. PERMANOVA on the matrices to test for differences and 
calculate_adonis <- function(fd, fd_meta)
{
	require(vegan)
	adonis(fd ~ inoculum * substrate * time_point * temperature, data=fd_meta, permuations=100)
	adonis(fd ~ temperature * substrate, data=fd_meta, permuations=100)
	adonis(fd ~ substrate, data=fd_meta, permuations=100)
	adonis(fd ~ temperature, data=fd_meta, permuations=100)
}


process_functional_diversity <- function(tx_fn="FD.traits.txt", fd_fn="FD.example.sample.txt", normalized=FALSE)
{
	require(fBasics)
	require(vegan)
	
	digits=2
	#options(digits=digits)
	
	samples <- read.table(fd_fn, header=T, row.names=1)
	traits <- read.table(tx_fn, header=T, row.names=1)
	
	# Initialize a new list for storing the results
	sample_traits <- list()
	
	
	sample_traits[["normalized"]] <- normalized
	if(normalized) 
	{
		samples <- decostand(samples, method="total", MARGIN=1)
		# Generate a weight matrix
		sample_traits[['weight']] <- decostand(samples, method="total", MARGIN=1)
		
		# Community weighted mean
		sample_traits[['cwm']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))		
	} else {
		sample_traits[['mean']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))
		sample_traits[['sd']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))
		sample_traits[['kurtosis']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))
		sample_traits[['skewness']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))	
		
		#sample_traits[['weight']] <- matrix(1.0, nrow(samples), ncol(samples), dimnames=list(rownames(samples), colnames(samples)))
	}

	
	#sample_traits[['print']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))
	#rownames(sample_traits[['print']]) <- rownames(samples)
	#colnames(sample_traits[['print']]) <- colnames(traits)

	#sample_traits[['mean']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))
	#sample_traits[['sd']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))
	#sample_traits[['kurtosis']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))
	#sample_traits[['skewness']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))
	


	
	sample_n <- nrow(samples)
	for(i in 1 : sample_n)
	{
		sample_id <- rownames(samples)[i]
		
		if(normalized)
		{
			# Check the orders of species and traits
			#### !!!Potential bug!!! ###
		
			for(j in 1 : ncol(selected_traits))
			{
				trait_name <- colnames(selected_traits)[j]
				sample_traits[['cwm']][i,trait_name] <- weighted.mean(traits[colnames(samples), trait_name], sample_traits[['weight']][sample_id, colnames(samples)])
			}

		} else {
			sample_traits[['print']] <- matrix(0.0, nrow(samples), ncol(traits), dimnames=list(rownames(samples), colnames(traits)))
			
			
			# Select the species observed in the current sample
			selected_species <- names(samples[sample_id, samples[sample_id,] > 0])	
			selected_traits <- traits[selected_species, ]	
			
			sample_traits[['kurtosis']][i,] <- as.numeric(format(round(colKurtosis(selected_traits, method="moment"), digits), nsmall = digits))
			sample_traits[['sd']][i,] <- as.numeric(format(round(colSds(selected_traits), digits), nsmall = digits))
			sample_traits[['mean']][i,] <- as.numeric(format(round(colAvgs(selected_traits), digits), nsmall = digits))
			sample_traits[['skewness']][i,] <- as.numeric(format(round(colSkewness(selected_traits, methods="moment"), digits), nsmall = digits))
			
			sample_traits[['print']][i,] <- paste(sample_traits[['mean']][i,], " (sd=", sample_traits[['sd']][i,], ", kurt=", sample_traits[['kurtosis']][i,], ", skw=", sample_traits[['skewness']][i,], ")", sep="")
			
			sample_traits[['kurtosis']][i,] <- colKurtosis(selected_traits, method="moment")
			sample_traits[['sd']][i,] <- colSds(selected_traits)
			sample_traits[['mean']][i,] <- colAvgs(selected_traits)
			sample_traits[['skewness']][i,] <- colSkewness(selected_traits, methods="moment")
		}
		
		#k_vals <- format(round(colKurtosis(selected_traits, method="moment"), digits), nsmall = digits)
		#s_vals <- format(round(colSds(selected_traits), digits), nsmall = digits)
		#m_vals <- format(round(colAvgs(selected_traits), digits), nsmall = digits)
		#w_vals <- format(round(colSkewness(selected_traits, methods="moment"), digits), nsmall = digits)
		
		#sample_traits[i,] <- paste(m_vals, " (sd=", s_vals, ", kurt=", k_vals, ", skw=", w_vals, ")", sep="")
		#selected_traits <- traits[names(sample[sample_id, sample[sample_id, ] > 0]), ]
		#k_val <- kurtosis(as.numeric(unlist(selected_traits)), method="moment")
		#s_val <- sd(as.numeric(unlist(selected_traits)))
		#m_val <- mean(as.numeric(unlist(selected_traits)))
		
		
		#print(paste(m_vals, " (sd=", s_vals, ", kurt=", k_vals, ", skw=", w_vals, ")", sep=""), sep=": ", collapse="\t")
	}
	
	return(sample_traits)
}




ss <- process_functional_diversity()
nss <- process_functional_diversity(normalized=T)




do_coocurrence_analysis("samples-tax_id+go.samples.g.clustered.picked")


do_coocurrence_analysis <- function(picked_fn)
{
	#picked_fn = "samples-tax_id+go.samples.g.clustered.picked"
	
	
	wk_fn <- picked_fn
	species_min_hit <- 500000
	normalized <- T
	
	fd <- read.table(wk_fn, sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)
	
	
# Function for calculating the column max
	colMaxs <- function(fd) { return(sapply(fd, max))} 
	
	C55_group <- rownames(fd)[grep("Cell55_", rownames(fd))]
	seed_group <- rownames(fd)[grep("Seed_", rownames(fd))]
	C35_group <- rownames(fd)[grep("Cell_", rownames(fd))]
	X35_group <- rownames(fd)[grep("Xyl_", rownames(fd))]
	
	
#filtered_fd <- fd[which(!rownames(fd) %in% c(seed_group, C55_group)),]
	filtered_fd <- fd[which(rownames(fd) %in% X35_group),]
	
	
	excluding_list <- which(colMaxs(filtered_fd) < species_min_hit)
	
# Subset fd
	filtered_fd <- filtered_fd[, - excluding_list]
	
	if(normalized)
	{
		sums <- rowSums(filtered_fd)
		filtered_fd <- filtered_fd / sums
	}
	
	
# Compute the correlation
	filtered_fd.cor <- cor(filtered_fd)
	
# Plot the matrix
	library(lattice)
	
	rgb.palette <- colorRampPalette(c("red", "white", "steelblue"), space = "rgb")
#levelplot(filtered_fd.cor, col.regions=rgb.palette(120), scales=list(x=list(rot=90)), at=seq(-1,1,0.05))
	levelplot(filtered_fd.cor, col.regions=rgb.palette(120),  title=wk_fn, scales=list(x=list(rot=90)))
	
	
	levelplot(filtered_fd, col.regions=rgb.palette(120), title=wk_fn, scales=list(x=list(rot=90)))
	

}



###

fn = "samples-tax_id+go.samples.g.clustered.transposed.filtered_norm_0.05"
plot_diversity_mtx(fn)
 
###
plot_diversity_mtx <- function(fn)
{
	require(reshape2)
	require(ggplot2)
	require(plyr)
	
	wk_fn = fn
	level = "g__"
	
	#df <- read.table(wk_fn, sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)
	df <- read.table(wk_fn, sep="\t", header=T, stringsAsFactors=F, comment.char="@")
		
	colnames(df)[1] <- "tax_name"
	tax_names <- unique(df$tax_name)
	#levels(df$tax_name) <- tax_names
	#df$tax_name <- reorder(df$tax_name, rowSums(df[,2:ncol(df)]))
	
	
	plot_df <- melt(df)	
	colnames(plot_df) <- c("tax_name", "sample", "relative_abundance")
	
	plot_df <- ddply(plot_df, .(sample), transform, pos = cumsum(relative_abundance) - (0.5 * relative_abundance))
	
	#plot_df <- ddply(plot_df, .(sample), transform, label= strsplit(strsplit(tax_name, level)[[1]][2], ";")[[1]][1])
	plot_df <- ddply(plot_df, .(sample), transform, label= extract_taxid(tax_name, paste(format(round(relative_abundance * 100, 2), nsmall = 2), "%", sep=""), level))
	
	
	#levels(plot_df$tax_name) <- unique(plot_df$tax_name)
	pdf(paste(fn, ".diversity_barplot.pdf" , sep=""), w=30, h=10)
	g <- ggplot(plot_df, aes_string(x="sample", y="relative_abundance", fill="tax_name")) + geom_bar(stat='identity') +
			#scale_y_continuous(limits = c(0, 1)) +
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
			geom_text(aes_string(label = "label", y = "pos"), size = 3)
	print(g)
	dev.off()
}





#
# extract_taxid(plot_df$tax_name, "g__")

#
extract_taxid <- function(tx, val, level) 
{
	taxids = rep("NA", length(tx))
	#print(length(tx))
	for(i in 1 : length(tx))
	{
		t <- tx[i]
		v <- val[i]
		ids <- unlist(strsplit(t, level))
		if (length(ids) > 1)
		{
			#print(ids[2])
			taxids[i] <- paste(ids[2], v, sep=": ")
		} else {
			taxids[i] <- paste(t, v, sep=": ")
		}
	}
	
	return(taxids)
}





# do_metagenome_tax_diversity("samples-tax_id+go.samples.g.clustered", "samples-tax_id+go.samples.g.clustered.meta")
do_metagenome_tax_diversity <- function(tx_fn_g, tx_fn_g_meta) 
{
	library(fBasics)
	library(reshape2)
	
	script_home = "/home/bcnskaa/projects/Metagenomics/trunk/src/"
	source(paste(script_home, "plot_heatmap.R", sep=""))
	source(paste(script_home, "plot_bubble_plot.R", sep=""))
	
	#tx_fn_g <- "samples-tax_id+go.samples.g.clustered"
	#tx_fn_g_meta <- "samples-tax_id+go.samples.g.clustered.meta"

	discard_others <- FALSE
	log_scale <- FALSE
	selected_level <- 30
	
	wk_fn <- tx_fn_g
	
	
# Read counts assigned to SEED subsystems
	fd <- read.table(wk_fn, sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)
	fd_meta <- read.table(tx_fn_g_meta, sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)
	
	tax_level = "Genus"
	fd_bak <- fd
	
	
	# Filtering the dataset by minimal abundance
	if(grep(".clustered", wk_fn))
	{
		###################
		# Min-count
		min_count <- 5000
		species_counts <- colSums(fd)
		
		excluded_species <- names(which(species_counts < min_count))
		# Remove species without kingdom assignment
		excluded_species <- rbind(excluded_species, names(fd[grep("k__.p_", rownames(fd))]))
		excluded_species <- rbind(excluded_species, "X.1.k__.p__.c__.o__.f__.g__")
		
		# Backup the excluded_fd
		excluded_fd <- fd[, which((colnames(fd) %in% excluded_species))]
		
		
		fd <- fd[, which(! (colnames(fd) %in% excluded_species))]
		
		
		excluded_samples = c("all_samples+nr.renamed.m8")
		
		fd <- fd[which(! (rownames(fd) %in% excluded_samples)), ]
		
		rownames(fd) <- gsub("+nr.m8", "", rownames(fd), fixed=T)
		colnames(fd) <- extract_tax_name(colnames(fd))
		fd <- t(fd)
	}
	
	
	# Sort the 
	sorted_subsystems <- sort(rowSums(fd), decreasing=T)
	
	selected_level <- min(selected_level, nrow(fd))
	selected_subsystem_ids <- names(sorted_subsystems)[1:selected_level]
	fd <- fd[selected_subsystem_ids,]
	

	# Normalize the data
	fd_totals <- colSums(fd)
		
	fd_normalized <- (as.data.frame(t(t(fd) / fd_totals))) * 100



	# http://www.statsblogs.com/2014/07/14/a-log-transformation-of-positive-and-negative-values/
	if(log_scale)
	{
		#fd_normalized <- log10(fd_normalized)
		fd_normalized <- log10(fd_normalized + 1)
	}

	selected_fd <- fd_normalized

	# Do hierachical clustering, and extract the labels
	library(vegan)

	fd_dist <- vegdist(t(fd_normalized), "bray")
	hc <- hclust(fd_dist)
	
	
	cell_height <- 0.3
	cell_width <- 0.4
	pdf_width <- 3 + (cell_width * ncol(selected_fd))
	pdf_height <- cell_height * nrow(selected_fd)
	
	pdf(paste(wk_fn, "hclust", "pdf", sep="."), w=pdf_width, h=4)
	plot(hc)
	dev.off()
	
	sample_ids <- hc$labels[hc$order]
	selected_fd <- selected_fd[sample_ids]

#selected_fd <- selected_fd[c("GZ.Seed_Y0.nr", "SWH.Seed_Y0.nr", "SWH.Cell55_Y2.nr", "GZ.Cell_Y1.nr", "GZ.Cell_Y2.nr","GZ.Xyl_Y1.nr", "GZ.Xyl_Y2.nr", "SWH.Cell_Y1.nr", "SWH.Cell_Y2.nr","SWH.Xyl_Y1.nr", "SWH.Xyl_Y2.nr")
	
	
	
	pdf(paste(wk_fn,"mtx","pdf",sep="."), w=pdf_width, h=pdf_height)
	
	plot_mtx <- as.matrix(selected_fd)
	plot_heatmap_mtx(plot_mtx, plot_label=T, midpoint=(max(plot_mtx) - min(plot_mtx)) / 2, x_axis_labels=rev(rownames(selected_fd)), xtitle="Sample", ytitle=tax_level, colorbar_scheme=c("steelblue", "yellow", "darkred"), label_size=3, value_decimal_len=2)
	dev.off()
	

	selected_fd$class <- rownames(selected_fd)
	
	#selected_fd$class <- factor(selected_fd$x_label, levels=x_axis_labels)
	
	plot_df <- melt(selected_fd)
	colnames(plot_df) <- c("class", "sample", "read_count_perc")
	
	pdf(paste(wk_fn,"bubble","pdf",sep="."), w=pdf_width, h=pdf_height)
	df <- plot_bubble(plot_df, "sample", "class", "read_count_perc", xtitle="Sample", ytitle=tax_level, ylabels=rev(rownames(selected_fd)))
	dev.off()
	
	
# Group analysis
	

	group_cell <- colnames(plot_mtx)[grep("Cell_", colnames(plot_mtx))]
	group_xyl <- colnames(plot_mtx)[grep("Xyl", colnames(plot_mtx))]
	group_seed <- colnames(plot_mtx)[grep("Seed", colnames(plot_mtx))]
	group_cell55 <- colnames(plot_mtx)[grep("Cell55", colnames(plot_mtx))]
	
	
	group_cell_avgs <- rowAvgs(plot_mtx[,group_cell])
	group_xyl_avgs <- rowAvgs(plot_mtx[,group_xyl])
	group_seed_avgs <- rowAvgs(plot_mtx[,group_seed])
	group_cell55_avgs <- plot_mtx[,group_cell55]
	
	fold_change <- TRUE
	
	# Normal to zero
	if(fold_change)
	{
		require(gtools)
		group_cell_avgs <- log(foldchange(group_cell_avgs, group_seed_avgs))
		group_xyl_avgs <- log(foldchange(group_xyl_avgs, group_seed_avgs))
		group_seed_avgs <- log(foldchange(group_seed_avgs, group_seed_avgs))
		group_cell55_avgs <- log(foldchange(group_cell55_avgs, group_seed_avgs))
	} else {
		group_cell_avgs <- group_cell_avgs - group_seed_avgs
		group_xyl_avgs <- group_xyl_avgs - group_seed_avgs
		group_seed_avgs <- group_seed_avgs - group_seed_avgs
		group_cell55_avgs <- group_cell55_avgs - group_seed_avgs		
	}
	
	
	
	#group_avgs <- data.frame(rownames(plot_mtx), group_cell_avgs, group_xyl_avgs, group_cell55_avgs, group_seed_avgs)
	#colnames(group_avgs) <- c("name", "Cell", "Xyl", "Cell55", "Seed")
	group_avgs <- data.frame(rownames(plot_mtx), group_cell_avgs, group_xyl_avgs, group_cell55_avgs)
	colnames(group_avgs) <- c("name", "Cell", "Xyl", "Cell55")

	melt_group_avgs <- melt(group_avgs)
	melt_group_avgs$name <- factor(melt_group_avgs$name, levels=rev(selected_subsystem_ids))
	
	#melt_group_avgs$value <- sign(melt_group_avgs$value) * log2(abs(1 + melt_group_avgs$value))
	
	if(fold_change)
	{	
		pdf(paste(wk_fn,"group_compare_fold","pdf",sep="."), w=5, h=pdf_height)	
	} else {
		pdf(paste(wk_fn,"group_compare","pdf",sep="."), w=5, h=pdf_height)
	}
	#xyplot(name ~ value, data=melt_group_avgs, groups=melt_group_avgs$variable, col=c("red", "steelblue", "green"), cex=1, pch=19)
	xyplot(name ~ value, data=melt_group_avgs, groups=melt_group_avgs$variable, alpha=0.6, col=c("red", "steelblue", "green"), cex=1, pch=19, panel = function(...) {
				panel.abline(v=0.5, lty = "dotted", col = "grey")
				panel.xyplot(...)
			})
	dev.off()

	
	
	
}



extract_tax_name <- function(lineages, missing="Others") 
{
	lineage_names <- strsplit(lineages, ".", fixed=T)
	tax_names <- rep("", length(lineage_names))
	tax_level="g__"

	for(i in 1 : length(tax_names))
	{
		lineage = lineage_names[[i]]
		tax_names[i] = lineage[length(lineage)]
		tax_names[i] = gsub(tax_level, "", tax_names[i], fixed=T)
		
		if(nchar(tax_names[i]) == 0)
		{
			tax_names[i] = missing
		}
		print(tax_names[i])
	}
	#tax_names <- gsub(tax_level, "", tax_names, fixed=T)
	
	return(tax_names)
}
