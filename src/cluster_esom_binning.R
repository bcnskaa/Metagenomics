# 2014 SKWoolf bcnskaa AT gmail DOT com


#
#
process_esom <- function(path_to_esom_dir, output_prefix=character(0), output_dir=character(0), selected_ids=character(0)) 
{	
	# We are about to import class definitions
	cls_fn <- list.files(path=path_to_esom_dir, pattern="cls");	

	if(length(cls_fn) != 1)
	{
		print(paste("Unable to read class definition from ", path_to_esom_dir, sep=""));
		#return(NULL);
	}
	#cls_fn <- cls_fn[0];
	
	print_msg(paste("Reading class definition from", cls_fn));
	fh<-file(cls_fn);
	open(fh, "r");
	cls_lines <- readLines(fh);
	close(fh);
	print_msg(paste("Number of lines read:", length(cls_lines)));
	
	
	# Import coordinates of bestmatches
	bm_fn <- list.files(path=path_to_esom_dir, pattern="bm");	
	
	if(length(bm_fn) != 1)
	{
		print(paste("Unable to read coordinates of bestmatches from ", path_to_esom_dir, sep=""));
		return(NULL);
	}

	print_msg(paste("Reading coordinates of bestmatches from", bm_fn));
	fh<-file(bm_fn);
	open(fh, "r");
	bm_lines <- readLines(fh);
	close(fh);
	print_msg(paste("Number of coordinates read:", length(bm_lines)));
	
	# Discard the first two lines
	bm_lines <- bm_lines[c(-1, -2)];
	
	
#	# Import class names
#	names_fn <- list.files(path=path_to_esom_dir, pattern="names");	
#	
#	if(length(names_fn) != 1)
#	{
#		print(paste("Unable to read class names from ", path_to_esom_dir, sep=""));
#		return(NULL);
#	}
#	#names_fn <- names_fn[0];
#	
#	print_msg(paste("Reading class names from", names_fn));
#	fh<-file(names_fn);
#	open(fh, "r");
#	names_lines <- readLines(fh);
#	close(fh);
#	print_msg(paste("Number of lines read:", length(names_lines)));
#	
	
	# create a dataframe for holding class definition
	# Get a hint on the number of bestmatch to be processed
	bm_n <- as.integer(gsub("% ", "", cls_lines[1]));
	# Remove the first line from cls_lines
	cls_lines <- cls_lines[-1];

	
	# Classes to which all bestmatches belong
	esom_cls <- data.frame(class_id=character(0), name=character(0), stringsAsFactors=F);
	# Bestmatch and their coordinates
	#esom_bm <- data.frame(bm_id=seq(1, bm_n), class_id=character(bm_n), x=integer(bm_n), y=integer(bm_n), stringsAsFactors=F);
	esom_bm <- data.frame(bm_id=seq(1, bm_n), class_id=character(bm_n), class_name=character(bm_n), x=rep(-1, bm_n), y=rep(-1, bm_n), stringsAsFactors=F);
	
	
	line_n <- 0;
	for(line in cls_lines)
	{
		if(nchar(line) < 3)
			next;
		
		if(length(grep("^%[0-9]", line)) > 0) # Class definition
		{
			#print(line);
			items = unlist(strsplit(gsub("%", "", line), "\t"));
			print_msg(paste(items[1], items[2], sep=":"));	
			esom_cls <- rbind(esom_cls, data.frame(class_id=items[1], name=items[2], stringsAsFactors=F));
		} else {
			ids <- unlist(strsplit(line, "\t"));
			esom_bm$class_id[as.integer(ids[1])] <- as.character(ids[2]);
			
			if(line_n %% 1000 == 0)
			{
				print_msg(paste("Line: ", line_n, " / ", bm_n, sep=""))
			}
		}
		line_n <- line_n + 1
	}
	print_msg(paste("Line: ", line_n, " / ", bm_n, " Done.", sep=""))
	
	
	
	# Parse coordinates of bestmatches
	line_n <- 0;
	for(line in bm_lines)
	{
		if(nchar(line) < 3)
			next;
		
		ids <- unlist(strsplit(line, "\t"));
		
		bm_idx <- as.integer(ids[1]);
		
		esom_bm$x[bm_idx] <- as.integer(ids[2]);
		esom_bm$y[bm_idx] <- as.integer(ids[3]);
		
		
		if(line_n %% 1000 == 0)
		{
			print_msg(paste("Line: ", line_n, " / ", bm_n, sep=""))
		}
		line_n <- line_n + 1
	}
	print_msg(paste("Line: ", line_n, " / ", bm_n, " Done.", sep=""))
	
	
	# Assign class name
	#apply(esom_cls, 1, function(cls) { esom_bm[which(esom_bm$class_id==cls["class_id"]), "class_name"] <- cls["name"] })
	
	for(i in 1 : nrow(esom_cls))
	{
		print(length(esom_bm$class_name[which(esom_bm$class_id==esom_cls$class_id[i])]))
		
		esom_bm$class_name[which(esom_bm$class_id==esom_cls$class_id[i])] <- as.character(esom_cls$name[i]);
	}

	
	# Release resources
	rm(cls_lines, bm_lines);
	
	
	pdf_fn <- character(0);
	outdir <- character(0);
	if(length(output_dir) > 0)
	{
		outdir <- paste(output_dir, "/", sep="");
	}
	
	if(length(output_prefix) > 0)
	{
		pdf_fn <- paste(outdir, output_prefix, ".mtx.pdf" ,sep="")
	}
	
	if(length(selected_ids) > 0)
	{
		row_ids <- selected_ids;
	}
	
	# Calculate the distribution differernces
	distance_mtx <- group_diff(esom, row_ids=row_ids, mtx_pdf_fn=pdf_fn, pdf_title=output_prefix);

	# Estimate grouping
	hist(distance_mtx[distance_mtx >= 0])
	fitdist(abs((distance_mtx[distance_mtx >= 0] / max(distance_mtx)) - 1) , "gamma")
	for(id in row_ids) {
		cat(paste(id, "\n", sep="")); 
		cat(rownames(distance_mtx)[which(distance_mtx[,id] < 0.175 & distance_mtx[,id] >= 0)]);
		cat("\n")
	}
	

	
	return(fns);
}



# 
plot_esom_binning <- function(distance_mtx, abundances, completeness, genome_sizes, gc_contents, pdf_fn=character(0), pdf_title=character(0))
{
	require(ggplot2);
	require(reshape2);
	

	rgb.palette <- colorRampPalette(c("yellow", "blue"), space = "rgb");
	#levelplot(distance_mtx, main="", xlab="", ylab="", col.regions=rgb.palette(120), at=seq(0,1,0.01), scales=list(x=list(rot=90)))
	#plot_mtx <- levelplot(distance_mtx, main="Distance", xlab="", ylab="", scales=list(x=list(cex=.8, rot=90),y=list(cex=.8)), at=seq(0,max(distance_mtx),0.01), col.regions=heat.colors)
	
	melt_distance_mtx <- melt(distance_mtx);
	colnames(melt_distance_mtx) <- c("x", "y", "distance");
	#plot_mtx <- ggplot(melt_distance_mtx, aes(x, y, fill=distance)) + geom_tile();
	plot_mtx <- ggplot(melt_distance_mtx, aes(x, y, fill=distance)) + 
			theme_bw() + 
			geom_tile() + scale_fill_gradientn(guide = FALSE, colours = c("red", "orange", "yellow", "white"), values=c(0,0.3, 0.6,1), limits=c(0,max(melt_distance_mtx$distance))) + 
			theme(plot.margin = unit(c(0,0.1,-0.25,-0.25), "cm")) + 
			#theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
			theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=8), axis.text.y = element_text(size=8), axis.title.x = element_blank(), axis.title.y = element_blank());
	
			#theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8), axis.text.y = element_text(size=8), axis.title.x = element_blank(), axis.title.y = element_blank());
	
	print_msg(paste("Ploting distance matrix", sep=""));
	
	g1_xlab <- get_grob_xlab(plot_mtx);
	g1_ylab <- get_grob_ylab(plot_mtx);
	plot_mtx <- plot_mtx + labs(x=NULL, y=NULL) + theme(axis.text = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank());
	g1 <- ggplotGrob(plot_mtx);
	
	
	plot_abundances <- generate_binning_hist(data.frame(name=names(abundances), value=as.double(abundances)));
	#plot_abundances <- ggplot(data=data.frame(name=names(abundances), value=as.double(abundances)), aes(y=value, x=name)) + geom_bar(stat="identity") + coord_flip() + #ggtitle("Abundance") +
	#theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank());
			
	print_msg(paste("Ploting bin groups' abundances", sep=""));
	g2 <- ggplotGrob(plot_abundances);
	
	plot_completeness <- generate_binning_hist(data.frame(name=names(completeness), value=as.double(completeness)));
	#plot_completeness <- ggplot(data=data.frame(name=names(completeness), value=as.double(completeness)), aes(y=value, x=name)) + geom_bar(stat="identity") + coord_flip() + #ggtitle("Completeness") +
	#theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank());
	
	print_msg(paste("Ploting bin groups' completeness", sep=""));
	g3 <- ggplotGrob(plot_completeness);
	
	plot_genome_sizes <- generate_binning_hist(data.frame(name=names(genome_sizes), value=as.double(genome_sizes)));
	#plot_genome_sizes <- ggplot(data=data.frame(name=names(genome_sizes), value=as.double(genome_sizes)), aes(y=value, x=name)) + geom_bar(stat="identity") + coord_flip() + #ggtitle("Genome Size") +
	#theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank());
	
	print_msg(paste("Ploting bin groups' genome_sizes", sep=""));
	g4 <- ggplotGrob(plot_genome_sizes);
	
	plot_gc_contents <- generate_binning_hist(data.frame(name=names(gc_contents), value=as.double(gc_contents)));
	# Slight leave some space to the right
	#plot_gc_contents <- plot_gc_contents + theme(plot.margin = unit(c(0.1,0.1,0.1,-0.3), "cm"));
	plot_gc_contents <- plot_gc_contents + theme(plot.margin = unit(c(0.0,0.0,-0.25,-0.3), "cm"));


	#plot_gc_contents <- ggplot(data=data.frame(name=names(gc_contents), value=as.double(gc_contents)), aes(y=value, x=name)) + geom_bar(stat="identity") + coord_flip() + #ggtitle("GC%") + 
	#theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank());
	print_msg(paste("Ploting bin groups' gc_contents", sep=""));
	g5 <- ggplotGrob(plot_gc_contents);
	

	#plots <- list(plot_mtx, plot_abundances, plot_completeness, plot_genome_sizes, plot_gc_contents);
	gplots <- list(g1_ylab, g1,g2,g3,g4,g5,g1_xlab);
	
	print_msg(paste("Number of plots=", length(gplots), sep=""));
	
	# For setting up a plot layout system
	require(gtable);
	
	print_msg(paste("Preparing Layout...", sep=""));
	
	#glayout <- lapply(1:length(plots), function(i) ggplotGrob(plots[i]));
	
	#print_msg(paste("Number of glayout=", length(glayout), sep=""));
	
	
	gtable <- gtable(widths=unit(rep(1,12), "null"), heights=unit(rep(1,7), "null"));
	gtable <- gtable_add_grob(gtable, gplots, l=c(1,2,9,10,11,12,2), r=c(1,8,9,10,11,12,8), t=c(1,1,1,1,1,1,6), b=c(5,5,5,5,5,5,6));
	
	
	if(length(pdf_fn) > 0)
	{
		pdf_w <- (0.2 * nrow(distance_mtx)) + (2 * 4);
		pdf_h <- 0.3 * ncol(distance_mtx);
			
		print_msg(paste("Plots will be exported to ", pdf_fn, sep=""));
		pdf(pdf_fn, width=pdf_w, height=pdf_h);
		
		grid.newpage();
		grid.draw(gtable);
		#print(pp);
		
		dev.off();
	} else {
		print_msg(paste("Plots will be exported to ", pdf_fn, sep=""));
		grid.newpage();
		grid.draw(gtable);
	}
}


generate_binning_hist <- function(df)
{
	gp <- ggplot(data=df, aes(y=value, x=name)) + 
			theme_bw() + 
			geom_bar(stat="identity") + coord_flip() + #ggtitle("GC%") + 
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) + 
			theme(plot.margin = unit(c(0,0.1,-0.25,-0.2), "cm")) + 
			labs(x=NULL, y=NULL) + 
			theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank());
	return (gp);
}


#
# 
#
get_grob_object <- function(ggplot_object, grob_label="guide-box")
{
	#require(gridExtra)
	tmp <- ggplot_gtable(ggplot_build(ggplot_object))
	selected_idx <- grep(grob_label, sapply(tmp$grobs, function(x) x$name))
	print(selected_idx)
	if(length(selected_idx) > 0)
	{
		obj <- tmp$grobs[[selected_idx]]
		return(obj);
	} else {
		return(NULL);
	}
}


get_grob_xlab <- function(ggplot_object)
{
	#require(gridExtra)
	tmp <- ggplot_gtable(ggplot_build(ggplot_object))
	selected_idx <- grep("axis.title.x", sapply(tmp$grobs, function(x) x$name))
	print(selected_idx)
	if(length(selected_idx) > 0)
	{
		obj <- tmp$grobs[[selected_idx - 1]]
		return(obj);
	} else {
		return(NULL);
	}
	
	#return(get_grob_object(ggplot_object, "axis.title.x"));
}

get_grob_ylab <- function(ggplot_object)
{
	#require(gridExtra)
	tmp <- ggplot_gtable(ggplot_build(ggplot_object))
	selected_idx <- grep("axis.title.y", sapply(tmp$grobs, function(x) x$name))
	print(selected_idx)
	if(length(selected_idx) > 0)
	{
		obj <- tmp$grobs[[2]]
		return(obj);
	} else {
		return(NULL);
	}
}


get_grob_legend_guide <- function(ggplot_object)
{
	return(get_grob_object(ggplot_object, "guide-box"));
}


estimate_centroid <- function(esom_bm)
{
	class_names <- as.character(unique(esom_bm$class_name));
	
	#centroids <- data.frame(class_name=class_names, x=integer(length(class_names)), y=integer(length(class_names)));
	centroids <- data.frame("class_name"=character(0), x=integer(0), y=integer(0), stringsAsFactors=F);
	
	for(i in 1 : length(class_names))
	{
		class_name <- class_names[i];
		selected_bm <- esom_bm[esom_bm$class_name==class_name,];
		print_msg(paste(class_name, "=",  summary(selected_bm$x)[4], ", ", summary(selected_bm$y)[4], sep=""));
		centroids <- rbind(centroids, data.frame(class_name=class_name, x=as.integer(summary(selected_bm$x)[4]), y=as.integer(summary(selected_bm$y)[4]), stringsAsFactors=F));
	}
	
	return(centroids);
}

# Calculate the differences of bestmatches
#
# Cartesian coordinates of bestmatches are first separated into x and y components. Distribution of each components
# are then constructed
# 
# http://stackoverflow.com/questions/5453336/plot-correlation-matrix-into-a-graph
# Use: group_diff(esom_bm, row_ids=unique(esom_bm$class_name)[grep("contig", unique(esom_bm$class_name))])
#group_diff <- function(esom_bm, row_ids=character(0), distance_thres=0.1, mtx_pdf_fn=character(0), pdf_title=character(0), step_size=5, compare_myself=F, na.value=-1)
group_diff <- function(esom_bm, row_ids=character(0), distance_thres=0.1, step_size=5, compare_myself=F, na.value=-1)		
{
	require(lattice);
	
	class_names <- as.character(unique(esom_bm$class_name));
	class_n <- length(class_names);
	
	xrange <- range(esom_bm$x);
	yrange <- range(esom_bm$y);
	
	#break_n_x <- 100;  
	#break_n_y <- 100;
	breaks_x <- c(seq(range(esom_bm$x)[1], range(esom_bm$x)[2], step_size), max(esom_bm$x));
	breaks_y <- c(seq(range(esom_bm$y)[1], range(esom_bm$y)[2], step_size), max(esom_bm$y));
	
	# Creat an empty distance matrix
	#distance_mtx <- data.frame(matrix(-1, class_n, class_n));
	#colnames(distance_mtx) <- class_names;
	#rownames(distance_mtx) <- class_names;
	distance_mtx <- matrix(rep(na.value, class_n*class_n), nrow=class_n, ncol=class_n, dimnames = list(class_names, class_names))
	
	print_msg("Calculating distance matrix");
	
	# Go through all the matrix
	for(i in 1 : (class_n - 1))
	{
		query_name <- class_names[i];
		query_subset <- esom_bm[esom_bm$class_name == query_name,];
		
		query_hist_x <- hist(query_subset$x, breaks=breaks_x, plot=F);
		query_hist_y <- hist(query_subset$y, breaks=breaks_y, plot=F);
		
		# Go through all the classes
		for(j in i: class_n)
		{
			if(!compare_myself){
				if(i == j)
					next;
			}
			
			subject_name <- class_names[j];
			subject_subset <- esom_bm[esom_bm$class_name == subject_name,];
	
			subject_hist_x <- hist(subject_subset$x, breaks=breaks_x, plot=F);
			subject_hist_y <- hist(subject_subset$y, breaks=breaks_y, plot=F);
			
			x_diff <- sum(abs(query_hist_x$density - subject_hist_x$density));
			y_diff <- sum(abs(query_hist_y$density - subject_hist_y$density));
			
			#print_msg(paste(query_name, " vs ", subject_name, ": x_diff=", x_diff, ", y_diff=", y_diff, sep=""));
		
			#distance_mtx[query_name, subject_name] <- (x_diff + y_diff) / 2;
			#distance_mtx[subject_name, query_name] <- (x_diff + y_diff) / 2;
			distance_mtx[query_name, subject_name] <- max(c(x_diff, y_diff));
			distance_mtx[subject_name, query_name] <- max(c(x_diff, y_diff));	
		}
	}
	
	
	print_msg(paste("Selected row: ", length(row_ids), sep=""));
	if(length(row_ids) > 0)
	{
		print_msg(paste("Selecting rows", sep=""));
		distance_mtx <- distance_mtx[, colnames(distance_mtx) %in% row_ids];
	}
	
	#pdf_w <- 0.2 * nrow(distance_mtx);
	#pdf_h <- 0.3 * ncol(distance_mtx);
	
	print_msg("Generating distance matrix");
#	if(length(mtx_pdf_fn) > 0)
#	{
#		print_msg(paste("Result will be exported to ", mtx_pdf_fn, sep=""));
#		pdf(mtx_pdf_fn, width=pdf_w, height=pdf_h);
#	}
	
	# Create a color scale
	rgb.palette <- colorRampPalette(c("yellow", "blue"), space = "rgb");
	#levelplot(distance_mtx, main="", xlab="", ylab="", col.regions=rgb.palette(120), at=seq(0,1,0.01), scales=list(x=list(rot=90)))
	pp <- levelplot(distance_mtx, main=pdf_title, xlab="", ylab="", scales=list(x=list(cex=.8, rot=90),y=list(cex=.8)), at=seq(0,max(distance_mtx),0.01), col.regions=heat.colors)
#	if(length(mtx_pdf_fn) > 0)
#	{
#		print(pp);
#		dev.off();
#	}
	pp;
	
	
	return(distance_mtx);
}



print_msg <- function(msg)
{
	cat(paste(msg, "\n", sep=""));
}