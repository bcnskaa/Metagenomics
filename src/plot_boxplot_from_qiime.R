library(ggplot2)
library(reshape)



# Collapse and melt a list of dataframes into a single dataframe
melt_df <- function(df_list, new_class_label="class", ignore_na=TRUE) {
	library(reshape);
	
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



# Prepare for impoting data
kingdoms <- c("bacteria", "archaea");
alpha_diversities <- c("PD_whole_tree", "chao1", "fisher_alpha", "shannon", "observed_species");

# Discard some columns from the datasets
col_to_be_removed <- c("X", "sequences.per.sample", "iteration", "H.GZ", "H.SHX")
sample_to_be_ignored <- "gDNA"
# Characters used to delimit id
ids_delim_chars <- "[.|_]";


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
		
		# select the rows with highest sequences per sample
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

# Discard rows with ids containing patterns enlisted in sample_to_be_ignored
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

# Cast values into correct datatype 
melt_dfs$variable <- factor(melt_dfs$variable);
melt_dfs$alpha_diversity <- factor(melt_dfs$alpha_diversity);


# Plot value histogram
pdf("all_values.histogram.pdf")
hist(as.double(melt_dfs$value), breaks=100, col="steelblue", xlab="Alpha Diversity (Combined all)", main="Spectrum of Alpha Diversity")
dev.off();


# 
kingdoms <- as.character(unique(melt_dfs$kingdom))
alpha_diversities <- as.character(unique(melt_dfs$alpha_diversity));
experiments <- as.character(unique(melt_dfs$experiment));
samples <- as.character(unique(melt_dfs$sample));
temperatures <- as.character(unique(melt_dfs$temperature))


# Combined all and fixed in default variables
for(kingdom in kingdoms)
{
	ss <- melt_dfs[melt_dfs$kingdom == kingdom, ];
	
	#ss[which(ss$alpha_diversity == "chao1"), c(1,2,4)]
	
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

	for(alpha_diversity in alpha_diversities)
	{
		ss <- melt_dfs[melt_dfs$kingdom == kingdom & melt_dfs$alpha_diversity == alpha_diversity, ];
		
		#melt_ss <- melt(ss, id=c("sample", "experiment", "temperature"))
		ss <- cbind(ss, id=paste(ss$sample, ss$experiment, ss$temperature, sep=":"))
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


