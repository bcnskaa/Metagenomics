library(fBasics)
library(reshape2)

script_home = "/home/bcnskaa/projects/Metagenomics/trunk/src/"
source(paste(script_home, "plot_heatmap.R", sep=""))
source(paste(script_home, "plot_bubble_plot.R", sep=""))


tx_fn <- "all_samples.tax.genus.txt"
fd_fn <- "all_samples.seed.functions.txt"

log_scale <- FALSE




# Read counts assigned to SEED subsystems
fd <- read.table(fd_fn, sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)

fd_totals <- colSums(fd)

fd_normalized <- (as.data.frame(t(t(fd) / fd_totals))) * 100

# http://www.statsblogs.com/2014/07/14/a-log-transformation-of-positive-and-negative-values/
if(log_scale)
{
	#fd_normalized <- log10(fd_normalized)
	fd_normalized <- log10(fd_normalized + 1)
}

sorted_subsystems <- sort(rowSums(fd_normalized), decreasing=T)

selected_fd <- fd_normalized[names(sorted_subsystems)[1:60],]

#selected_fd <- fd_normalized[which(fd_normalized[, "GZ.Cell_Y1.nr"] > 0.5),]
#selected_fd <- fd_normalized[which(fd_normalized[, "GZ.Cell_Y1.nr"] > 0.5),]


# Do hierachical clustering, and extract the labels
library(vegan)
fd_dist <- vegdist(t(fd_normalized), "bray")
hc <- hclust(fd_dist)

pdf(paste(tx_fn,"hclust","pdf",sep="."), w=14, h=5)
plot(hc)
dev.off()
sample_ids <- hc$labels[hc$order]

selected_fd <- selected_fd[sample_ids]

#selected_fd <- selected_fd[c("GZ.Seed_Y0.nr", "SWH.Seed_Y0.nr", "SWH.Cell55_Y2.nr", "GZ.Cell_Y1.nr", "GZ.Cell_Y2.nr","GZ.Xyl_Y1.nr", "GZ.Xyl_Y2.nr", "SWH.Cell_Y1.nr", "SWH.Cell_Y2.nr","SWH.Xyl_Y1.nr", "SWH.Xyl_Y2.nr")


pdf(paste(tx_fn,"mtx","pdf",sep="."), w=14, h=16)
plot_mtx <- as.matrix(selected_fd)
plot_heatmap_mtx(plot_mtx, plot_label=T, midpoint=(max(plot_mtx) - min(plot_mtx)) / 2, x_axis_labels=rev(rownames(selected_fd)), xtitle="Sample", ytitle="SEED Subsystem", colorbar_scheme=c("steelblue", "yellow", "darkred"), label_size=3, value_decimal_len=2)
dev.off()



selected_fd$subsystem <- rownames(selected_fd)
plot_df <- melt(selected_fd)
colnames(plot_df) <- c("system", "sample", "read_count")

df <- plot_bubble(plot_df, "sample", "system", "read_count", xtitle="Sample", ytitle="SEED System", ylabels=rev(rownames(selected_fd)))






# Read counts assigned to taxonomy
tx <- read.table("all_samples.tax.genus.txt", sep="\t", header=T, stringsAsFactors=F, comment.char="@", row.names=1)









library(FD)

