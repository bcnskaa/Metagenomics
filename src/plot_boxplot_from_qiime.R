# 2014 SKWoolf bcnskaa AT gmail DOT com
library(ggplot2)
library(reshape)

# The group info file should be in tabular format and contains the following headers
group_info_field_labels <- c("Label", "Place", "Treatment", "Year", "Temperature", "Replicate", "Excluded");
# Name of alpha diversities
alpha_diversities <- c("PD_whole_tree", "chao1", "fisher_alpha", "shannon", "observed_species");

# 
import_group_info <- function(info_fn)
{
	group_info <- read.table(info_fn, header=T, sep="\t", stringsAsFactors=F);

	cflag = TRUE;
	for(field_label in group_info_field_labels)
	{
		cflag = cflag & field_label %in% colnames(group_info);
	}
	
	label_cflag = FALSE;
	for(i in 1:length(group_info$Label))
	{
		label <- group_info$Label[i];
		if(grepl("-", label))
		{
			group_info$Label[i] <- gsub("-", ".", label);
			label_cflag = TRUE;
		}
	}
	
	if(label_cflag)
		cat(paste("Warning: label(s) contain illegal character.\n"));	
	
	# 
	if(!cflag){
		cat(paste("Group info ", info_fn, " is not valid, please check the names of column headers.\n", sep=""));
		return(NULL);
	}else{
		cat(paste(nrow(group_info), " items read from ", info_fn, ".\n", sep=""));
		return(group_info);
	}
}


# Collapse and melt a list of dataframes into a single dataframe
melt_df <- function(df_list, new_class_label="class", ignore_na=TRUE) {
	library(reshape2);
	
	class_labels <- names(df_list);
	
	melted_df <- data.frame();
	
	# Go through all labels
	
	for(label in class_labels) {
		print(paste("Processing ", label, "...", sep=""));
		
		cur_df <- df_list[[label]];
		cur_df <- cbind(cur_df, iteration=seq(1, nrow(cur_df)))
		
		# Melt cur_df by holding iteration as an id
		mdf <- melt(cur_df, id=c("iteration"));
		
		# 
		mdf <- cbind(mdf, X______X_____X=rep(label, nrow(mdf)));
		
		if(dim(melted_df)[1] == 0)
		{
			melted_df <- mdf;
		} else {
			melted_df <- rbind(melted_df, mdf);
		}
	}
	
	colnames(melted_df)[colnames(melted_df) == "X______X_____X"] <- new_class_label;
	
	return(melted_df);
}


main2 <- function(group_info_fn, path, mean_replicate=FALSE)
{
	alpha_diversities <- c("PD_whole_tree", "chao1", "fisher_alpha", "shannon", "observed_species");
	# Discard some columns from the datasets
	col_to_be_removed <- c("X", "sequences.per.sample", "iteration", "H.GZ", "H.SHX");

	group_info <- import_group_info(group_info_fn);

	# Import the data
	dfs<- list()
	for(alpha_diversity in alpha_diversities)
	{
		cat(paste("Reading from ", alpha_diversity, ".txt\n", sep=""));
		
		# Id to the current datasets
		id <- alpha_diversity;
		
		# Import dataset
		data <- read.table(paste(path, "/", alpha_diversity, ".txt", sep=""), header=T, sep="\t", stringsAsFactors=F);
		
		# select the rows with a highest sequences per sample value
		data <- subset(data, data$sequences.per.sample == max(data$sequences.per.sample));
		
		# Discard columns
		data <- data[,-which(colnames(data) %in% col_to_be_removed)]
			
		# Insert data into list
		dfs[[id]] <- data;
	}
	
	cat("The following labels are set to be excluded from further analyses.\n");
	excluded_list <- group_info$Label[which(group_info$Excluded == "1")];
	cat("\t");cat(paste(excluded_list,sep=",", join="  "));cat("\n");
	
	
	
	cat(paste("Checking values...\n", sep=""));
	# Raise warning message if the dataset contains N/A
	for(alpha_diversity in alpha_diversities)
	{
		d <- dfs[[alpha_diversity]];	
		# Discard excluded columns
		d <- d[,-which(colnames(d) %in% excluded_list)];
		
		labels <- colnames(d);
		for(label in labels) 
		{
			d[[label]] <- as.double(d[[label]]);
			
			if(length(which(is.na(d[[label]]))) > 0)
			{
				cat(paste("Warning: ", alpha_diversity, "/", label, " contains N/A values\n", sep=""));
			}
		}
		
		
		dfs[[alpha_diversity]] <- d;
		
	}
	
	# Reshape dataframe
	melt_dfs <- melt_df(dfs);
	
	# Discard rows with n/a values
	melt_dfs <- melt_dfs[-which(is.na(melt_dfs$value)), ];

	# Go through all alpha diversities
	for(alpha_diversity in alpha_diversities)
	{
		# Select data on basis of selected alpha diversity
		ss <- melt_dfs[melt_dfs$alpha_diversity == alpha_diversity, ];
		
		
		
		
		
		#melt_ss <- melt(ss, id=c("sample", "experiment", "temperature"))
		ss <- cbind(ss, id=paste(ss$sample, ss$experiment, ss$temperature, sep=":"))
		
		# Take a mean value from replicates
		if(flag_mean)
		{
			require(plyr)
			ss <- ddply(ss, colnames(ss)[which(colnames(ss) != "value")], summarise, value=mean(value))
		}
		
		pdf(paste("combined_replicates.", alpha_diversity, ".", kingdom,".boxplot.pdf", sep=""), width=4, height=6);
		print(qplot(factor(id), value, data = ss, geom = "boxplot") + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
		# qplot(factor(variable), value, data = ss, geom = "boxplot", fill=factor(alpha_diversity)) + geom_jitter(position=position_jitter(w=0.1, h=0.1)) 
		#p <- ggplot(ss, aes(factor(variable), value)) + geom_jitter();
		#p + geom_boxplot(aes(fill=factor(alpha_diversity))) + theme(axis.text.x = element_text(angle = 90, hjust = 1));
		dev.off();
		
		# Combined place
		pdf(paste("combined_place.", alpha_diversity, ".", kingdom,".boxplot.pdf", sep=""), width=4, height=6);
		print(qplot(factor(sample), value, data = ss, geom = "boxplot") + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
		dev.off();
		
		# Combined experiment
		pdf(paste("combined_experiment.", alpha_diversity, ".", kingdom,".boxplot.pdf", sep=""), width=4, height=6);
		print(qplot(factor(experiment), value, data = ss, geom = "boxplot") + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
		dev.off();
		
		
		# Combined temperature
		pdf(paste("combined_temperature.", alpha_diversity, ".", kingdom,".boxplot.pdf", sep=""), width=4, height=6);
		print(qplot(factor(temperature), value, data = ss, geom = "boxplot") + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
		dev.off();
		rm(ss);			
	}
}



# Usage: 
# 
main <- function(path, mean_replicate=False)
{
	# Prepare for importing data
	kingdoms <- c("bacteria", "archaea");
	alpha_diversities <- c("PD_whole_tree", "chao1", "fisher_alpha", "shannon", "observed_species");
	
	# Discard some columns from the datasets
	col_to_be_removed <- c("X", "sequences.per.sample", "iteration", "H.GZ", "H.SHX")
	sample_to_be_ignored <- "gDNA"
	
	# Characters used to delimit id
	ids_delim_chars <- "[.|_]";

	# flag: take a mean value from replicates
	flag_mean = T
	
	# Data split into two kingdoms, import them separately
	dfs<- list()
	for(kingdom in kingdoms)
	{
		# import each alpha diversity
		for(alpha_diversity in alpha_diversities)
		{
			print(paste("Reading from ", kingdom, "/", alpha_diversity, ".txt", sep=""))
			
			# Id to the current datasets
			id <- paste(kingdom, ".", alpha_diversity, sep="");
			
			# Import dataset
			data <- read.table(paste(kingdom, "/", alpha_diversity, ".txt", sep=""), header=T, sep="\t", stringsAsFactors=F);
			
			# select the rows with a highest sequences per sample value
			data <- subset(data, data$sequences.per.sample == max(data$sequences.per.sample));
			
			# Discard columns
			data <- data[,-which(colnames(data) %in% col_to_be_removed)]
			
			
			# Insert data into list
			dfs[[id]] <- data;
		}
	}
	
	# Reshape dataframe
	melt_dfs <- melt_df(dfs);
	
	# Discard rows with n/a values
	melt_dfs <- melt_dfs[-which(melt_dfs$value == "n/a"), ]
	
	# Discard rows containing patterns enlisted in sample_to_be_ignored
	melt_dfs <- melt_dfs[-which(grepl("gDNA", melt_dfs$variable) == TRUE), ]
	
	
	# split variable column into something meaningful
	# As strsplit looking for regexp pattern, we have to put the period into a bracket or to include a fixed flag
	library(stringr)
	meta_info <- str_split_fixed(paste(melt_dfs$class, melt_dfs$variable, sep="."), "[.]", 6);
	
	# Rename the colnames into something meaningful
	colnames(meta_info) <- c("kingdom", "alpha_diversity", "sample", "experiment", "temperature", "replicate_id");
	
	# Merge the new columns into data.frame
	melt_dfs <- cbind(melt_dfs, meta_info);
	
	# Tidy up ids by removing redundant columns
	melt_dfs <- melt_dfs[, -which(colnames(melt_dfs) %in% c("iteration", "class"))]
	
	# Cast values into correct datatype 
	melt_dfs$value <- as.double(melt_dfs$value);
	
	# Factorize values 
	melt_dfs$variable <- factor(melt_dfs$variable);
	melt_dfs$alpha_diversity <- factor(melt_dfs$alpha_diversity);
	
	
	# Plot a histogram of observed values
	pdf("all_values.histogram.pdf")
	hist(as.double(melt_dfs$value), breaks=100, col="steelblue", xlab="Alpha Diversity (Combined all)", main="Spectrum of Alpha Diversity")
	dev.off();
	
	
	# Prepare important variables
	kingdoms <- as.character(unique(melt_dfs$kingdom))
	alpha_diversities <- as.character(unique(melt_dfs$alpha_diversity));
	experiments <- as.character(unique(melt_dfs$experiment));
	samples <- as.character(unique(melt_dfs$sample));
	temperatures <- as.character(unique(melt_dfs$temperature))
	
	
	# Combined all and fixed in default variables
	for(kingdom in kingdoms)
	{
		# Filter data from selected kingdom
		ss <- melt_dfs[melt_dfs$kingdom == kingdom, ];
		
		pdf(paste("all_values.", kingdom,".boxplot.pdf", sep=""), width=20, height=5);
		print(qplot(factor(variable), value, data = ss, geom = "boxplot", fill=factor(alpha_diversity)))
		# qplot(factor(variable), value, data = ss, geom = "boxplot", fill=factor(alpha_diversity)) + geom_jitter(position=position_jitter(w=0.1, h=0.1)) 
		#p <- ggplot(ss, aes(factor(variable), value)) + geom_jitter();
		#p + geom_boxplot(aes(fill=factor(alpha_diversity))) + theme(axis.text.x = element_text(angle = 90, hjust = 1));
		dev.off();
		rm(ss);
	}
	
	
	# Combined replicates 
	for(kingdom in kingdoms)
	{
		# Go through all alpha diversities
		for(alpha_diversity in alpha_diversities)
		{
			# Select data on basis of selected kingdom and alpha diversity
			ss <- melt_dfs[melt_dfs$kingdom == kingdom & melt_dfs$alpha_diversity == alpha_diversity, ];
			
			#melt_ss <- melt(ss, id=c("sample", "experiment", "temperature"))
			ss <- cbind(ss, id=paste(ss$sample, ss$experiment, ss$temperature, sep=":"))
			
			# Take a mean value from replicates
			if(flag_mean)
			{
				require(plyr)
				ss <- ddply(ss, colnames(ss)[which(colnames(ss) != "value")], summarise, value=mean(value))
			}
			
			pdf(paste("combined_replicates.", alpha_diversity, ".", kingdom,".boxplot.pdf", sep=""), width=4, height=6);
			print(qplot(factor(id), value, data = ss, geom = "boxplot") + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
			# qplot(factor(variable), value, data = ss, geom = "boxplot", fill=factor(alpha_diversity)) + geom_jitter(position=position_jitter(w=0.1, h=0.1)) 
			#p <- ggplot(ss, aes(factor(variable), value)) + geom_jitter();
			#p + geom_boxplot(aes(fill=factor(alpha_diversity))) + theme(axis.text.x = element_text(angle = 90, hjust = 1));
			dev.off();
			
			# Combined place
			pdf(paste("combined_place.", alpha_diversity, ".", kingdom,".boxplot.pdf", sep=""), width=4, height=6);
			print(qplot(factor(sample), value, data = ss, geom = "boxplot") + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
			dev.off();
			
			# Combined experiment
			pdf(paste("combined_experiment.", alpha_diversity, ".", kingdom,".boxplot.pdf", sep=""), width=4, height=6);
			print(qplot(factor(experiment), value, data = ss, geom = "boxplot") + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
			dev.off();
			
			
			# Combined temperature
			pdf(paste("combined_temperature.", alpha_diversity, ".", kingdom,".boxplot.pdf", sep=""), width=4, height=6);
			print(qplot(factor(temperature), value, data = ss, geom = "boxplot") + geom_jitter(position=position_jitter(w=0.1, h=0.1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
			dev.off();
			rm(ss);			
		}
	}
}


