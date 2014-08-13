# 2014 SKWoolf bcnskaa AT gmail DOT com


#
#
prepare_esom_entity <- function(path_to_esom_dir) 
{	
	# We are about to import class definitions
	cls_fn <- list.files(path=path_to_esom_dir, pattern=".cls");	

	if(length(cls_fn) != 1)
	{
		print(paste("Unable to read class definition from ", path_to_esom_dir, sep=""));
		return(NULL);
	}
	cls_fns <- cls_fns[0];
	
	print_msg(paste("Reading class definition from", cls_fns);
	cls_lines <- readLines(cls_fns);
	print_msg(paste("Number of lines read:", length(cls_lines)));
	
	
	# Import class names
	names_fn <- list.files(path=path_to_esom_dir, pattern=".names");	
	
	if(length(names_fn) != 1)
	{
		print(paste("Unable to read class names from ", path_to_esom_dir, sep=""));
		return(NULL);
	}
	names_fn <- names_fn[0];
	
	print_msg(paste("Reading class names from", names_fn);
	names_lines <- readLines(names_fn);
	print_msg(paste("Number of lines read:", length(names_lines)));
	
	
	
	# create a dataframe for holding class definition
	class_n <- integer(gsub("% ", "", lines[0]));
	esom_cls <- data.frame(character(class_n), character(class_n));
	for(line in cls_lines)
	{
		#if(length(line) == 0)
		#	next;
		
		if(length(grep("^%[0-9]", line)) > 0) # Class definition
		{
			#print(line);
			items = unlist(strsplit(gsub("%", "", line), "\t"));
			print_msg(paste(items[1], items[2], sep=":"));
		} else {
			
		}
	}
	
	
	return(fns);
}


print_msg <- function(msg)
{
	print(msg);
}