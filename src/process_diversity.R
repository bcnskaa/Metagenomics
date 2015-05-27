library(fBasics)
library(reshape2)

script_home = "/home/bcnskaa/projects/Metagenomics/trunk/src/"
source(paste(script_home, "plot_heatmap.R", sep=""))
source(paste(script_home, "plot_bubble_plot.R", sep=""))



pfam_fn <- "/home/bcnskaa/projects/Metagenomics_WD/FunctionalDiversity/all_samples/Pfam/all_samples.renamed+Pfam.dom.tbl.summary"
trait_fn <- "samples-clustered-tax_id+go.traits"
tx_fn_s <- "samples-tax_id+go.samples"
tx_fn_g <- "samples-tax_id+go.samples.g.clustered"
#tx_fn <- "samples-clustered-tax_id+go.samples"
#tx_fn <- "all_samples.tax.species.txt"
#fd_fn <- "all_samples.seed.functions.txt"

log_scale <- FALSE
selected_level <- 30

wk_fn <- tx_fn_g

# Read counts assigned to SEED subsystems
fd <- read.table(wk_fn, sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)

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
	min_global_trait_count <- 20
	min_species_trait_count <- 100
	
	trait_counts <- colSums(fd)
	length(which(trait_counts < min_global_trait_count))
	excluded_traits <- colnames(fd)[which(trait_counts < min_global_trait_count)]
	
	
	# Exclude taxa whose triat count is too small
	trait_counts <- rowSums(fd)
	excluded_species <- rownames(fd)[which(trait_counts < min_species_trait_count)]


	fd <- fd[which(! (rownames(fd) %in% excluded_species)), ! (colnames(fd) %in% excluded_traits)]
	fd <- fd[,- which(colnames(fd) == "NA.")]
	
	fd <- t(fd)
}



if(wk_fn == pfam_fn)
{
	min_global_trait_count <- 20
	min_species_trait_count <- 100
}


# Processing the pfam data
if(wk_fn == "samples-clustered-tax_id+go.traits")
{
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


pdf(paste(wk_fn,"hclust","pdf",sep="."), w=14, h=5)
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

pdf(paste(wk_fn,"bubble","pdf",sep="."), w=8, h=16)
df <- plot_bubble(plot_df, "sample", "system", "read_count_perc", xtitle="Sample", ytitle="SEED System", ylabels=rev(rownames(selected_fd)))
dev.off()





# Read counts assigned to taxonomy
tx <- read.table("all_samples.tax.genus.txt", sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)



# Do PCoA
library(rgl)

plot_mtx.pca <- princomp(plot_mtx)
pdf(paste(wk_fn,"pca_biplot","pdf",sep="."), w=8, h=16)
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
