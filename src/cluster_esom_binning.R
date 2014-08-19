# 2014 SKWoolf bcnskaa AT gmail DOT com


#
#
prepare_esom_entity <- function(path_to_esom_dir) 
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
	
	
	
	
	return(fns);
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


print_msg <- function(msg)
{
	cat(paste(msg, "\n", sep=""));
}