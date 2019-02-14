#sample_id<-"all_to_B" # Ok
default_sample_id<-"";
# default_kmer_fn<-paste(sample_id, ".kmer.tsv", sep="")
# #default_cov_fn<-paste(sample_id, ".cov.tsv", sep="")
# default_gc_fn<-paste(sample_id, ".gc.tsv", sep="")
# default_seqlen_fn<-paste(sample_id, ".seq_len.tsv", sep="")
# default_tax_fn<-paste(sample_id, ".tax.tsv", sep="")
default_min_contig_len<-"5000"
default_dimension<-2;



# main_menu
#   +-- Data (D)
#       +-- import (I)
#       +-- export (E)
#   +-- cluster (C)
#       +-- pick (P): zoom(), ahull(), 
#       +-- explore (X): list_cluster(), select_cluster()
#       +-- modify (M): append_to(), remove_from(), delete_cluster(), rename_cluster()
#   +-- Setting (S)
#       +-- Filtering
#           +-- GC
#           +-- sequence length
#           +-- taxonomic lineage
#       +-- Mode
#           +-- 2-d (1)
#               +-- X and Y (1)
#               +-- X and Z (2)
#               +-- Y and Z (3)
#           +-- 3-d (2) [Enabled only if 3 or more coverages are available]
#               +-- mds (3)
#           +-- ca (4)
#       +-- Plot
#           +-- cluster_size
#           +-- panels: 1,2,3,4
#           +-- font_size
#           +-- log
#           +-- show_cluster
#           +-- show_taxonomy
#   +-- exit (Q)
#


#SEPARATOR<-"===================================="
#SEPARATOR2<-"------------------------------------"
SEP_WIDTH<-60;
SEPARATOR<-paste(rep("=", SEP_WIDTH), sep="", collapse="");
SEPARATOR2<-paste(rep("-", SEP_WIDTH), sep="", collapse="");

STATUS_MAIN<-0;
STATUS_DATA<-1;
STATUS_CLUST<-2;
STATUS_SETTINGS<-3;
STATUS_EXIT<-4;
STATUS_ON_EXIT<-100;
DATA_MODE_KMER <- 1;
DATA_MODE_COV <- 2;

VERSION <- "0.2a";

TAG_GC<-"GC";
TAG_SEQ_LEN<-"SEQ_LEN";
TAG_TAX_LINEAGE<-"TAX_LINEAGE";
TAG_WORKING_CTG<-"WOKRING_CTG";

PDF_HEIGHT<-8;
PDF_WIDTH<-8;


RETRY_N <- 3;
DEFAULT_MIN_CONTIG_LEN<-800;
DEFAULT_TOP_N_MARKER<-10;

MIN_CONTIG_LEN<-DEFAULT_MIN_CONTIG_LEN;
PRINT_DEBUG <- F;
LOG_ENABLED<-TRUE;
SHOW_CLUSTER_BOUND_BOX<-FALSE;
SHOW_MARKER_LEGEND<-TRUE;
NO_TRANSPARENT<-FALSE;
TOP_N_MARKER<-DEFAULT_TOP_N_MARKER;
# 1 = PCA, 2 = Correspondence Analysis
TRANSFORM_METHOD_KMER<-F;
TRANSFORM_METHOD_MDS<-F;
PRINT_LARGE_PDF <- F;

cl_cols <- {};
col_palettes <- {};
col_palettes_map <- {};

check_sys <- function()
{
  if("rgl" %in% rownames(installed.packages()) == FALSE) {install.packages("rgl", repos="http://R-Forge.R-project.org")};
  if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer")};
  if("scatterplot3d" %in% rownames(installed.packages()) == FALSE) {install.packages("scatterplot3d", repos="http://R-Forge.R-project.org")}
  if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")};
  if("plyr" %in% rownames(installed.packages()) == FALSE) {install.packages("plyr")};
  if("vegan" %in% rownames(installed.packages()) == FALSE) {install.packages("vegan")};
  if("zoom" %in% rownames(installed.packages()) == FALSE) {install.packages("zoom")};
  if("alphahull" %in% rownames(installed.packages()) == FALSE) {install.packages("alphahull")};
  if("fpc" %in% rownames(installed.packages()) == FALSE) {install.packages("fpc")};
  #if("md5" %in% rownames(installed.packages()) == FALSE) {install.packages("md5")};
  
  #if("rgl" %in% rownames(installed.packages()) == FALSE) {install.packages("rgl", repos="http://R-Forge.R-project.org")}
  library(rgl)
  #if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer")}
  library(RColorBrewer)
  #if("scatterplot3d" %in% rownames(installed.packages()) == FALSE) {install.packages("scatterplot3d", repos="http://R-Forge.R-project.org")}
  library(scatterplot3d)
  #if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
  library(ggplot2)
  #if("plyr" %in% rownames(installed.packages()) == FALSE) {install.packages("plyr")}
  library(plyr)
  #if("vegan" %in% rownames(installed.packages()) == FALSE) {install.packages("vegan")}
  library(vegan)
  #if("zoom" %in% rownames(installed.packages()) == FALSE) {install.packages("zoom")}
  #if("alphahull" %in% rownames(installed.packages()) == FALSE) {install.packages("alphahull")}
  library(alphahull)
  #if("fpc" %in% rownames(installed.packages()) == FALSE) {install.packages("fpc")}
  library(fpc)
  #library(md5)
}


init <- function()
{
  # colors for 100 genome bins
  color_ramp_step <- 20;
  colorPalette <- colorRampPalette(c("darkred", "orange", "green", "steelblue"));
  cl_cols <<- colorPalette(color_ramp_step)[sample(seq(1,color_ramp_step))];
  
  
  col_palettes$default <<- colorRampPalette(c("darkred", "orange", "green", "steelblue"));
  col_palettes$spectrum <- colorRampPalette(c("purple", "red", "yellow", "green", "blue"));
  col_palettes$thermal <<- colorRampPalette(c("red", "blue"));
  col_palettes$extreme <<- colorRampPalette(c("red", "green", "purple", "yellow", "magenta", "blue", "orange", "darkgrey"));
  #col_palettes$thermal2 <<- colorRampPalette(c("yellow", "steelblue"));
  #col_palettes$thermal3 <<- colorRampPalette(c("red", "darkblue"));
  
  col_palettes_map$default <<- colorRampPalette(c("darkred", "yellow", "steelblue"));
  col_palettes_map$thermal <<- colorRampPalette(c("red", "blue"));
  col_palettes_map$thermal2 <<- colorRampPalette(c("yellow", "steelblue"));
  col_palettes_map$thermal3 <<- colorRampPalette(c("darkblue", "yellow", "darkred"));
  col_palettes_map$highlight_h <<- colorRampPalette(c("darkgray", "darkgray", "darkgray", "darkgray", "red"));
  col_palettes_map$highlight_m <<- colorRampPalette(c("darkgray", "darkgray", "red", "darkgray", "darkgray"));
  col_palettes_map$highlight_l <<- colorRampPalette(c("red", "darkgray", "darkgray", "darkgray", "darkgray"));
  
}


main_ctr <- function()
{
  check_sys();
  init();
  
  data_src <- {};
  status <- STATUS_MAIN;
  if(! interactive()) { display_msg("", mode=-1); display_msg("Please use the command: source(\"plot_differential_bin2.R\") under R environment."); status <- STATUS_ON_EXIT; }
  
  while(status != STATUS_ON_EXIT) {
    display_msg(character(0), mode=-1);
    
    if(status == STATUS_MAIN) {
      status <- mn_main(data_src);
    } else if(status == STATUS_DATA) {
      data_src <- mn_data(data_src);
      if(length(data_src) > 0) { print_debug_msg("main_ctr().STATUS_DATA"); plot_cov(data_src); }
      status <- STATUS_MAIN;
    } else if(status == STATUS_CLUST) {
      data_src <- mn_clustering(data_src);
      status <- STATUS_MAIN;
    } else if(status == STATUS_SETTINGS) {
      if(length(data_src) > 0) {
        update_pd <- mn_settings(data_src);
        if(update_pd$updating) {
          data_src$cov$sel_cov_col_idx <- update_pd$sel_cov_col_idx;
          data_src$pd <- update_pd(data_src);
        }
        plot_cov(data_src); 
      }
      status <- STATUS_MAIN;
    } else if(status == STATUS_EXIT) {
      if(mn_exit(data_src) == 0) { status <- STATUS_ON_EXIT; } else { status <- STATUS_MAIN; }
    } else if(status == STATUS_ON_EXIT) {
      break;
    } else {
      status <- STATUS_MAIN;
    }
  }
  
  quit("no");
}


#
mn_main <- function(ds)
{
  
  if(length(ds) == 0) {   
    display_msg("", mode=-1);
    display_msg("No data imported.");
    return(STATUS_DATA);
  } else {
    opt <- -1;      
    ops <- c("D", "C", "S", "Q");
    ops_status <- c(STATUS_DATA, STATUS_CLUST, STATUS_SETTINGS, STATUS_EXIT);
    while(opt < 0 || opt > length(ops)) {
      display_msg("", mode=-1);
      display_msg("", mode=2);
      display_msg("Main Menu", mode=1)
      display_msg("", mode=2);
      display_msg("", mode=4);
      display_msg("[D] Data", mode=3);
      display_msg("[C] Clustering", mode=3);
      display_msg("[S] Settings", mode=3);
      display_msg("", mode=4);
      display_msg("[Q] Quit", mode=3);
      display_msg("", mode=4);
      display_msg("", mode=4);
      display_msg("", mode=7);
      display_msg("", mode=4);
      display_msg("", mode=6);
      
      choice <- toupper(readline("Please enter an option: "));
      opt <- ifelse(!is.na(match(choice, ops)), match(choice,ops), -1);
    }
    return(ops_status[opt]);
  }
  
}


#   +-- Data (D)
#       +-- import (I)
#       +-- export (E)
mn_data <- function(ds)
{
  # Existing data
  if(length(ds) == 0) {
    display_msg("", mode=-1);
    ds <- mn_data_import();
  } else {
    opt <- -1;      
    ops <- c("N", "I", "E", "Q")
    while(opt < 0 || opt > length(ops)) {
      display_msg("", mode=-1);
      display_msg("", mode=2);
      display_msg("Data Menu", mode=1)
      display_msg("", mode=2);
      display_msg("", mode=4);
      
      display_msg(paste("Data currently in memory:"), mode=3);
      display_msg("", mode=4);
      display_msg(paste("> imported coverage data: ", ds$cov_src), mode=5);
      display_msg(paste("> imported kmer data: ", ds$kmer_src), mode=5);
      display_msg("", mode=4);
      
      display_msg("", mode=4);
      display_msg("[I] Import new data", mode=3);
      #display_msg("Import clusters [I]", mode=3);
      display_msg("[E] Export clusters to files", mode=3);
      display_msg("", mode=4);
      display_msg("[Q] Return [Q]", mode=3);
      display_msg("", mode=4);
      
      display_msg(paste("> imported coverage data: ", ds$cov_src), mode=5);
      display_msg(paste("> imported kmer data: ", ds$kmer_src), mode=5);
      display_msg("", mode=4);
      display_msg("", mode=4);
      display_msg("", mode=6);
      
      choice <- toupper(readline("Please enter an option: "));
      opt <- ifelse(!is.na(match(choice, ops)), match(choice,ops), -1);
    }
    
    # Import new data set
    if(opt == which(ops == "N")) {    
      new_ds <- mn_data_import();
      eval.parent(substitute(ds<-new_ds));
    }
    
    # Import clusters
    if(opt == which(ops == "I")) {
      #new_ds <- mn_data_import();
      eval.parent(substitute(ds<-new_ds));
    }
    
    # Export clusters
    if(opt == which(ops == "E")) { export_clusters_to_file(ds); }
    
    #if(opt == which(ops == "R")) { return(STATUS_MAIN); }
  }
  
  return(ds);
}


#       +-- pick (P): zoom(), ahull(), 
#       +-- explore (X): list_cluster(), select_cluster()
#       +-- modify (M): append_to(), remove_from(), delete_cluster(), rename_cluster()
mn_clustering <- function(ds)
{
  require(zoom);
  paxy <- par("usr"); 
  
  # Existing data
  if(length(ds) == 0)
  {
    display_msg("", mode=-1);
    display_msg("Data is empty", mode=1);
  } else {
    is_zoomed <- F;
    show_gc <- F;
    opt <- -1;      
    ops <- c("E", "P", "M", "L", "G", "Z", "V", "S", "Q");
    while(opt < 0 || opt > length(ops)) {
      display_msg("", mode=-1);
      display_msg("", mode=2);
      display_msg("Clustering Menu", mode=1);
      display_msg("", mode=2);
      display_msg("", mode=4);
      display_msg("[P] Create a new cluster", mode=3);
      
      if(length(ds$clusters) > 0) { 
        #display_msg("[D] Delete existing clusters", mode=3);
        #display_msg("[R] Rename an existing cluster", mode=3);
        display_msg("[M] Modify existing clusters", mode=3);
        display_msg("[E] Export clusters to list files", mode=3);
        
        display_msg("", mode=4);
        display_msg(paste("> Number of clusters: ", length(ds$clusters), sep=""), mode=5);
      }
      
      display_msg("", mode=4);
      display_msg("[L] Toogle marker", mode=3);      
      display_msg("[G] Toogle GC", mode=3);
      display_msg("[Z] Zoom", mode=3);
      if(is_zoomed) { display_msg("[V] Siwtch back to the global view", mode=3); }
      
      if(length(ds$clusters) > 0) { 
        display_msg("[S] Quick save", mode=3);
      }
      
      display_msg("", mode=4);  
      display_msg("[Q] Return to Main Menu", mode=3);
      display_msg("", mode=4);
      
      
      display_msg("", mode=4);
      display_msg("", mode=4);
      display_msg("", mode=6);
      
      choice <- toupper(readline("Please enter an option: "));
      opt <- ifelse(!is.na(match(choice, ops)), match(choice,ops), -1);
      
      # Export clusters
      if(opt == which(ops == "E") && length(ds$clusters) > 0) { export_clusters_to_file(ds); }
      if(opt == which(ops == "P")) { ds <- mn_clustering_pick(ds); }
      if(opt == which(ops == "L") && length(ds$clusters) > 0) { if(length(ds$clusters) > 0) { mn_clustering_list(ds); } else { opt <- -1; } }
      #if(opt == which(ops == "D") && length(ds$clusters) > 0) { ds <- mn_clustering_delete(ds); }
      #if(opt == which(ops == "R") && length(ds$clusters) > 0) { ds <- mn_clustering_rename(ds); }
      if(opt == which(ops == "M") && length(ds$clusters) > 0) { ds <- mn_clustering_modify(ds); }
      if(opt == which(ops == "G")) { show_gc <- !show_gc; }
      if(opt == which(ops == "V")) { is_zoomed <- !is_zoomed; }
      if(opt == which(ops == "S")) { save_clusters_to_file(ds, character(0)); }
      if(opt == which(ops == "Z")) { zm(); paxy <- par("usr"); is_zoomed <- T;}  
      if(opt == which(ops == "Q")) { break; }
      
      if(is_zoomed) { pax<-paxy[1:2]; pay<-paxy[3:4]; } else { pax<-character(0);pay<-character(0); }
      
      if(show_gc) { col_mode <- 1; } else { col_mode <- 0; };
      
      print_debug_msg("mn_clustering().while");
      plot_cov(ds, x_lim=pax, y_lim=pay, col_mode=col_mode);
      opt <- -1;
    }
    #if(opt == which(ops == "R")) { return(STATUS_MAIN); }
  }
  #eval.parent(substitute(ds<-ds));
  return(ds);
}


#
#
mn_clustering_list <- function(ds)
{
  display_msg("", mode=-1);
  #if(dev.list() > 0) { dev.off() };
  
  display_msg("", mode=-1);
  display_msg(SEPARATOR);
  display_msg("List of Clusters", mode=1);
  display_msg(SEPARATOR);  
  display_msg("", mode=4);
  
  if(length(ds$clusters) > 0) {  
    for(i in 1:length(ds$clusters)) {
      cl_id <- names(ds$clusters)[i];
      cl_name <- ds$cluster_names[[cl_id]];
      cl_ctgs <- ds$clusters[[cl_id]];
      display_msg(paste(i, ".) ", cl_name, ": contig#=", length(cl_ctgs), ", length=", sum(ds$seq_len[match(cl_ctgs, rownames(ds$cov$dat))]), " bp", sep=""), mode=5);
    }
    display_msg("", mode=4);
    display_msg("", mode=4);
    
  } else {
    display_msg("No cluster in the list", mode=5);
  }
  
  
  print_debug_msg("mn_clustering_list()");
  plot_cov(ds, col_mode=0, hl_mode=1);
  
  display_msg("", mode=6);
  
  choice <- toupper(readline("Return to Main Menu [Enter]: "));
  
}


#
# Menu for modifying the existing clusters
#
mn_clustering_modify <- function(ds)
{
  opt <- -1;
  ops <- c("D", "R", "M", "Q");
  while(opt == -1) {
    display_msg("", mode=-1);
    #if(dev.list() > 0) { dev.off() };
    
    display_msg("", mode=-1);
    display_msg(SEPARATOR);
    display_msg("Modifying Clusters", mode=1);
    display_msg(SEPARATOR);  
    display_msg("", mode=4);
    display_msg("", mode=4);
    display_msg("List of Clusters:", mode=3);
    display_msg("", mode=4);
    if(length(ds$clusters) > 0) { 
      for(i in 1:length(ds$clusters)) {
        cl_id <- names(ds$clusters)[i];
        cl_name <- ds$cluster_names[[cl_id]];
        cl_ctgs <- ds$clusters[[cl_id]];
        display_msg(paste("[", i, "] ", cl_name, ": contig#=", length(cl_ctgs), ", length=", sum(ds$seq_len[match(cl_ctgs, rownames(ds$cov$dat))]), " bp", sep=""), mode=5);
      }
      display_msg("", mode=4);
    } else {
      display_msg("No cluster in the list", mode=5);
    }
    
    display_msg("", mode=4);  
    #display_msg("List clusters [L]", mode=3);
    display_msg("[R] Rename a cluster", mode=3);
    display_msg("[D] Delete a cluster", mode=3);
    display_msg("", mode=4);
    display_msg("[Q] Return to previous menu", mode=3);
    display_msg("", mode=4);
    display_msg("", mode=4);
    display_msg("", mode=6);
    
    choice <- toupper(readline("Please select an option: "));
    opt <- ifelse(!is.na(match(choice, ops)), match(choice,ops), -1);
    
    pass <- F;
    if(opt == which(ops == "D")) {
      cl_idx <- -1;
      while(cl_idx == -1) {
        cl_idx <- as.numeric(readline(paste("Please select a cluster [1-", length(ds$clusters),"] to delete: ", sep="")));
        if (!is.na(cl_idx)) {
          if(cl_idx > 0 && cl_idx <= length(ds$clusters)) {
            cl_id <- names(ds$clusters)[cl_idx];
            cl_name <- ds$cluster_names[cl_id];
            choice <- toupper(readline(paste("Are you sure to delete ", cl_name, " (warning: the operation cannot be undone)? [Y/N] ", sep="")));
            
            if(choice == "Y") {
              ds$clusters <- ds$clusters[-which(names(ds$clusters) == cl_id)];
              ds$cluster_names <- ds$cluster_names[-which(names(ds$cluster_names) == cl_id)];
              pass <- T;
              break;
            }
            if(choice == "N") { break; }
          }
        }
        if(!pass) {
          display_msg(paste(cl_idx, "<-- Invalid cluster id selected.", sep=""), mode=3);
          cl_idx <- -1;  
        }
      }
      opt <- -1;
      
    }
    if(opt == which(ops == "R")) { 
      cl_idx <- -1;
      while(cl_idx == -1) {
        cl_idx <- as.numeric(readline(paste("Please select a cluster [1-", length(ds$clusters),"] to rename: ", sep="")));
        if (!is.na(cl_idx)) {
          if(cl_idx > 0 && cl_idx <= length(ds$clusters)) {
            cl_id <- names(ds$clusters)[cl_idx];
            cl_name <- ds$cluster_names[cl_id];
            
            while(T) {
              new_name <- readline(paste("Please enter a new name for the cluster, ", cl_name, "(hash_code=", cl_id,"): ", sep=""));
              
              if(nchar(new_name) == 0) { break; }
              if(check_name(new_name)) {
                if(! new_name %in% unlist(ds$cluster_names))  {
                  ds$cluster_names[which(names(ds$cluster_names) == cl_id)] <- new_name;
                  pass <- T;
                  break;                  
                } else {
                  display_msg(paste(new_name, "<-- Name is in used, please enter another name.", sep=""), mode=3);
                }
              }
            }
          }
        }
        if(!pass) {
          display_msg(paste(cl_idx, "<-- Invalid cluster id selected.", sep=""), mode=3);
          cl_idx <- -1;  
        }
      }
      opt <- -1;
    }
    
    if(opt == which(ops == "Q")) { break; }
  }
  
  return(ds);
}



mn_clustering_delete <- function(ds)
{
  #display_msg("List clusters [L]", mode=3);
  display_msg("Rename cluster [R]", mode=3);
  display_msg("Remove points from a cluster [R]", mode=3);
  display_msg("Delete cluster [D]", mode=3);
}

#11 
#
mn_clustering_pick <- function(ds)
{
  global_cov_paxy <- par("usr");
  
  #   display_msg(paste("xlim=", paxy[1:2], ", ylim=", paxy[3:4], collapse=",", sep=""), mode=3);
  #   
  #   display_msg("Please select 6 points that form a region covering points of interest: Manual\n");
  #   bin_def <- NULL
  #   bin_def <- ahull(locator(6, type="p", pch=20), alpha=100);
  #   plot(bin_def, add=T, col="black");
  #   # Extract contigs belonging to selected genome bin
  #   #pd <- ds$cov$dat[, c(ds$cov$sel_cov_col_idx[1], ds$cov$sel_cov_col_idx[2])];
  #   pd <- ds$pd$dat;
  #   picked_pts <- {}
  #   for(i in 1:nrow(pd)) { if(inahull(bin_def, c(pd$x[i], pd$y[i]))) picked_pts <- rbind(picked_pts, pd[i,])};
  #  
  
  #picked_pts <- pick_points(ds$pd$dat);
  if(length(ds$clusters_working) > 0) { picked_pts <- ds$clusters_working; } else { picked_pts <- character(0); }
  
  
  marker_col_palette <- colorRampPalette(c("purple", "red", "yellow", "green", "blue"));
  
  opt <- -1;
  
  v_params <- {};
  v_params$show_kmer <- F;
  v_params$show_picked <- T;
  v_params$show_zoomed <- T;
  v_params$show_cluster <- F;
  v_params$show_gc_col <- F;
  v_params$show_markers <- F; 
  # Current paxy
  v_params$paxy <- global_cov_paxy;
  v_params$zoomed_paxy <- NULL;
  v_params$col_mode <- 0; 
  v_params$hl_mode <- 0;
  v_params$saved_dat_mode <- 0;
  v_params$dat_mode <- 0;
  v_params$marker_idx <- 0;
  v_params$marker_col_palette_idx <- 1;
  v_params$map_idx <- 0;
  v_params$map_col_palette_idx <- 1;
  v_params$sel_cov_pair_idx <- 1;
  v_params$graph_style <- F;
  v_params$cov_pairs <- generate_cov_pairs(ds, length(ds$cov$sel_cov_col_idx));
  saved_v_params <- v_params;
  
  pdf_outfile_name <- character(0);
  
  #v_params$saved_cur_paxy <- paxy;
  #v_params$saved_zoomed_paxy <- zoomed_paxy;
  
  #   show_kmer <- F;
  #   show_picked <- T;
  #   show_zoomed <- T;
  #   show_cluster <- F;
  #   show_gc_col <- F;
  # 
  #   # Current paxy
  #   paxy <- global_cov_paxy;
  #   zoomed_paxy <- global_cov_paxy;
  #   saved_cur_paxy <- paxy;
  #   saved_zoomed_paxy <- zoomed_paxy;
  # col_mode <- 0; hl_mode <- 0; saved_dat_mode <- 0; dat_mode <- 0;
  ops <- c("A", "D", "P", "X", "G", "T", "O", "K", "V", "Z", "C", "Q", "E", "S", "N", "I", "J", "U", "F");
  while(opt < 0 || opt > length(ops)) {
    picked_pt_names <- character(0);
    
    display_msg("", mode=-1);
    display_msg("", mode=2);
    display_msg("Cluster Picking Menu", mode=1);
    display_msg("", mode=2);
    display_msg("", mode=4);
    
    display_msg("[A] Pick new points to the list", mode=3);
    if(length(picked_pts) > 0) { 
      display_msg("[X] Remove points from the list", mode=3); 
      display_msg("[E] Clear the list", mode=3);
      display_msg("", mode=4);
      display_msg(paste("> Number of points in the list: ", length(picked_pts), sep=""), mode=5);
      display_msg(paste("> Total size (bp): ", sum(ds$seq_len[match(picked_pts, rownames(ds$cov$dat))]), sep=""), mode=5);
      display_msg("", mode=4);
      display_msg("[N] Create a new cluster with picked points", mode=3);
      display_msg("", mode=4);
    }
    
    display_msg("", mode=4);
    
    
    display_msg(paste("[G] Toogle mapping: ", ifelse(v_params$map_idx > 0, names(ds$maps)[v_params$map_idx], "None"), sep=""), mode=3);
    display_msg(paste("[T] Toogle marker: ", ifelse(v_params$marker_idx > 0, names(ds$markers)[v_params$marker_idx], "None"), sep=""), mode=3);
    
    display_msg("[P] Toogle picked points", mode=3);
    if(length(ds$clusters) > 0) {
      display_msg("[C] Toogle clusters", mode=3);
      display_msg("", mode=4);
      display_msg(paste("> Number of clusters: ", length(ds$clusters), sep=""), mode=5);
      
      if(v_params$show_cluster) {
        for(i in 1:length(ds$clusters)) {
          cl_id <- names(ds$clusters)[i];
          cl_name <- ds$cluster_names[[cl_id]];
          cl_ctgs <- ds$clusters[[cl_id]];
          display_msg(paste(i, ".) ", cl_name, ": contig#=", length(cl_ctgs), ", length=", sum(ds$seq_len[match(cl_ctgs, rownames(ds$cov$dat))]), " bp", sep=""), mode=5);
        }
      }
    }
    display_msg("", mode=4);
    
    # Check if coverage pairs switching is available
    if(nrow(v_params$cov_pairs) > 1) {
      # Print cov paris
      for(i in 1 : nrow(v_params$cov_pairs)) {
        print_debug_msg(paste(i, ") ", paste(v_params$cov_pairs[i, ], collapse=", "), sep=""));
      }
      
      # Switch between coverage pairs
      display_msg("[D] Switch dataset", mode=3); 
      display_msg("", mode=4);
      display_msg(paste("> Number of coverage pairs: ", nrow(v_params$cov_pairs), sep=""), mode=5);
      display_msg(paste("> Current: ", paste(colnames(ds$cov$dat)[v_params$cov_pairs[v_params$sel_cov_pair_idx, ]], collapse=", "), sep=""), mode=5);
    }
    display_msg("", mode=4);
    
    
    if(v_params$show_kmer) { display_msg("[K] Switch to the coverage view", mode=3); } else {  display_msg("[K] Switch to the Kmer view", mode=3); }
    if(!is.null(v_params$zoomed_paxy)) {
      if(v_params$show_zoomed) {
        display_msg("[V]Show the global view", mode=3); 
      } else {
        display_msg("[V] Show previous zoomed view", mode=3);
      }
    }
    #if(!show_kmer) { display_msg("Zoom to a new location [Z]", mode=3); }
    display_msg("[Z] Zoom to a new location", mode=3);
    
    
    display_msg("", mode=4);
    display_msg(paste("Color scheme: ", sep=""), mode=3);
    display_msg(paste("  [I] Mapping: ", names(col_palettes_map)[v_params$map_col_palette_idx], sep=""), mode=3);
    display_msg(paste("  [J] Marker: ", names(col_palettes)[v_params$marker_col_palette_idx], sep=""), mode=3);
    display_msg(paste("  [U] Graph style: ", ifelse(v_params$graph_style == T, "Centroid", "Default"), sep=""), mode=3);
    
    display_msg("", mode=4);
    
    
    display_msg("", mode=4); 
    #if(v_params$graph_style) { display_msg("[F] Export the current plot to PDF", mode=3); }
    display_msg("[F] Export the current plot to PDF", mode=3);
    
    display_msg("[S] Quick save", mode=3);
    display_msg("", mode=4);
    display_msg("[Q] Return to Main Menu", mode=3);
    
    display_msg("", mode=4);
    display_msg("", mode=4);
    display_msg("", mode=6);
    
    choice <- toupper(readline("Please enter an option: "));
    opt <- ifelse(!is.na(match(choice, ops)), match(choice,ops), -1);
    
    if(opt == which(ops == "A")) {
      newly_picked_pts <- pick_points(ds$pd$dat, msg="Please select new points to be added to the current set.");
      if(length(newly_picked_pts) > 0) {
        if(length(picked_pts) == 0) { 
          picked_pts <- newly_picked_pts;
        } else {  
          old_nrow <- length(picked_pts);
          #newly_picked_pts <- subset(newly_picked_pts, ! rownames(newly_picked_pts) %in% rownames(picked_pts));
          #picked_pts <- rbind(picked_pts, newly_picked_pts);
          picked_pts <- unique(c(newly_picked_pts, picked_pts));
        }
      }
      v_params$show_picked <- T;
    }
    if(opt == which(ops == "X")) {
      to_be_removed_pts <- pick_points(ds$pd$dat, msg="Please select selected points to be removed from the current set.");
      if(length(to_be_removed_pts) > 0) { 
        old_nrow <- length(picked_pts);
        #picked_pts <- subset(picked_pts, ! rownames(picked_pts) %in% rownames(to_be_removed_pts));
        picked_pts <- unique(picked_pts[!picked_pts %in% to_be_removed_pts]);
        display_msg(paste("> Number of points has been removed: ", old_nrow - length(picked_pts), sep=""), mode=5);
      }
      v_params$show_picked <- TRUE;
    }
    
    if(opt == which(ops == "K")) { 
      #v_params$show_kmer <- !v_params$show_kmer;
      if(!v_params$show_kmer)
      {
        saved_v_params <- v_params;
        saved_v_params$paxy <- par("usr");
        v_params$show_kmer <- T;
        v_params$show_zoomed <- F;
        v_params$dat_mode <- 1;
        v_params$zoomed_paxy <- NULL;
        v_params$paxy <- par("usr");
        if(length(picked_pts) > 0 && v_params$show_picked) { 
          sel_ctgs <- picked_pts;
        } else {
          #cat("Select DS")
          x1<-v_params$paxy[1];x2<-v_params$paxy[2];y1<-v_params$paxy[3];y2<-v_params$paxy[4];
          #cat(paste("x1:", x1, " x2:", x2, " y1:", y1, " y2:", y2, sep=""));
          #cat(paste(colnames(ds$pd$dat)));
          sel_ds <- ds$pd$dat[which(x1 <= ds$pd$dat$x & x2 >= ds$pd$dat$x),];
          #cat(paste("sel_ds: ", nrow(sel_ds), sep=""));
          sel_ds <- sel_ds[which(y1 <= sel_ds$y & y2 >= sel_ds$y),];
          #cat(paste("sel_ds: ", nrow(sel_ds), sep=""));
          sel_ctgs <- rownames(sel_ds);
          #sel_ctgs <- rownames(ds$pd$dat)[which(x1 >= ds$pd$dat$x && x2 <= ds$pd$dat$x && y1 <= ds$pd$dat$y && y2 >= ds$pd$dat$y)];
        }
      } else {
        v_params <- saved_v_params;
        v_params$dat_mode <- 0;
      }
      ds$pd <- update_pd(ds, sel_ctgs=sel_ctgs, dat_mode=v_params$dat_mode);
    }
    
    #if(opt == which(ops == "Z") && !show_kmer) { zm(); paxy <- par("usr"); show_zoomed <- T; }
    if(opt == which(ops == "Z")) { v_params$paxy <- par("usr"); zm(); v_params$zoomed_paxy <- par("usr"); v_params$show_zoomed <- T; }
    if(opt == which(ops == "G")) {
      v_params$map_idx <- v_params$map_idx + 1;
      if(v_params$map_idx > length(ds$maps)) {
        v_params$map_idx <- 0;
      }
      
      #      cat(paste("v_params$map_idx=", v_params$map_idx, sep=""))
      
      #      if(v_params$map_idx > 0) { v_params$show_markers <- T; } else { v_params$show_markers <- F; }
      
      #       v_params$col_mode <- v_params$col_mode + 1;
      #       
      #       if(v_params$col_mode > 1) {
      #         if(length(picked_pts) == 0) { v_params$show_gc_col <- F; v_params$col_mode <- 0; }
      #         if(v_params$col_mode > 2 && length(picked_pts) > 0) { v_params$show_gc_col <- F; v_params$col_mode <- 0; }
      #       }
      #       
      #       #v_params$show_gc_col <- !v_params$show_gc_col;
      #       if(v_params$col_mode > 0) { 
      #         v_params$show_gc_col <- T;
      #         
      #         if(length(picked_pts) == 0 && v_params$col_mode != 1) { v_params$show_gc_col <- F; v_params$col_mode <- 0; }
      #         if(length(picked_pts) > 0 && v_params$col_mode > 1 && !v_params$show_picked) { v_params$show_gc_col <- F; v_params$col_mode <- 0; }
      #       } else { 
      #         v_params$col_mode <- 0; v_params$show_gc_col <- F;
      #       };
      
    }
    
    if(opt == which(ops == "T")) {
      v_params$marker_idx <- v_params$marker_idx + 1;
      if(v_params$marker_idx > length(ds$markers)) {
        v_params$marker_idx <- 0;
      } 
      
      if(v_params$marker_idx > 0) { v_params$show_markers <- T; } else { v_params$show_markers <- F; }
      #v_params$show_markers <- !v_params$show_markers;
    }
    if(opt == which(ops == "E")) { picked_pts <- {}; sel_ctgs <- {}; }
    if(length(ds$clusters) > 0 && opt == which(ops == "U")) { dev.off(); v_params$graph_style <- !v_params$graph_style; }
    
    if(opt == which(ops == "C")) { 
      if(length(ds$clusters) > 0){ v_params$show_cluster <- !v_params$show_cluster; } else { v_params$show_cluster <- FALSE; }
    }
    
    if(opt == which(ops == "S")) { save_clusters_to_file(ds, picked_pts); }
    if(opt == which(ops == "P")) { v_params$show_picked <- !v_params$show_picked; }
    if(opt == which(ops == "V") && !is.null(v_params$zoomed_paxy)) { v_params$show_zoomed <- !v_params$show_zoomed; }
    if(opt == which(ops == "N")) {
      if(length(picked_pts) > 0) {
        cluster_list <- create_new_cluster(ds, picked_pts);
        
        ds$clusters <- cluster_list$clusters;
        ds$cluster_names <- cluster_list$cluster_names;
        
        v_params$show_picked <- F;
        #show_cluster <- ds$cluster_names[length(ds$clusters)];
        picked_pts <- {};
      }
      v_params$show_cluster <- T;
    }
    if(opt == which(ops == "Q")) { break; }
    
    # Color schemes
    if(opt == which(ops == "J")) {
      v_params$marker_col_palette_idx <- v_params$marker_col_palette_idx + 1;
      if(v_params$marker_col_palette_idx > length(col_palettes)) {
        v_params$marker_col_palette_idx <- 1;
      } 
    }
    if(opt == which(ops == "I")) {
      v_params$map_col_palette_idx <- v_params$map_col_palette_idx + 1;
      if(v_params$map_col_palette_idx > length(col_palettes_map)) {
        v_params$map_col_palette_idx <- 1;
      } 
    }
    
    if(opt == which(ops == "F")) {
      while(T) {
        pdf_outfile_name <- readline("Please enter a name for PDF file: ");
        if(check_name(pdf_outfile_name)) { break; }
        display_msg(paste(pdf_outfile_name, " is not valid.", sep=""));
      }
      if(!grepl(".pdf", pdf_outfile_name)) { pdf_outfile_name <- paste(pdf_outfile_name, ".pdf", sep=""); }
    }
    
    #update_pd(ds, sel_cov_col_idx=character(0), sel_ctgs=character(0), pt_weighted=T, log=T, vis_mode=0, dat_mode=0, transform_method=1, min_seq_len=MIN_CONTIG_LEN);     v_params$sel_cov_pair_idx <- 1;
    #v_params$cov_pairs <- generate_cov_pairs(ds, length(ds$cov$sel_cov_col_idx));
    
    # Switch dataset
    if(nrow(v_params$cov_pairs) > 1 && opt == which(ops == "D")) {
      v_params$sel_cov_pair_idx <- v_params$sel_cov_pair_idx + 1;
      if(v_params$sel_cov_pair_idx > nrow(v_params$cov_pairs)) { v_params$sel_cov_pair_idx <- 1; }
      # Update selected coverages
      ds$cov$sel_cov_col_idx <- v_params$cov_pairs[v_params$sel_cov_pair_idx, ];
      ds$pd <- update_pd(ds, sel_ctgs=sel_ctgs, dat_mode=v_params$dat_mode);
    }
    
    
    if(v_params$show_picked) { picked_pt_names <- picked_pts; } else { picked_pt_names<-character(0); }
    #if(v_params$show_gc_col) { v_params$col_mode <- 1; } else { v_params$col_mode <- 0; }
    if(v_params$show_cluster) { v_params$hl_mode <- 1; } else { v_params$hl_mode <- 0; }
    if(v_params$show_zoomed) { pax<-v_params$zoomed_paxy[1:2];pay<-v_params$zoomed_paxy[3:4];} else { pax<-character(0);pay<-character(0); }
    #if(v_params$show_zoomed) { pax<-v_params$zoomed_paxy[1:2];pay<-v_params$zoomed_paxy[3:4]; } else { pax<-v_params$paxy[1:2];pay<-v_params$paxy[3:4]; }
    
    #if(v_params$show_markers) { v_params$marker_idx <- c(TAG_TAX_LINEAGE); } else { v_params$marker_idx <- character(0); } 
    
    print_debug_msg("mn_clustering_pick()");
    
    if(v_params$graph_style) {
      if(length(pdf_outfile_name) > 0) { 
        plot_cluster_diagram(ds, x_lim=v_params$paxy[1:2], y_lim=v_params$paxy[3:4], pdf_ofn=pdf_outfile_name); 
        pdf_outfile_name <- character(0);
      }
      plot_cluster_diagram(ds, x_lim=v_params$paxy[1:2], y_lim=v_params$paxy[3:4]);
      
    } else {
      print_debug_msg(paste("print_debug_msg=", pdf_outfile_name, sep=""));
      if(length(pdf_outfile_name) > 0) { 
        plot_cov(ds, sel_ctgs=picked_pt_names, map_idx=v_params$map_idx, marker_idx=v_params$marker_idx, map_col_palette=col_palettes_map[[v_params$map_col_palette_idx]], top_n_marker=TOP_N_MARKER, marker_col_palette=col_palettes[[v_params$marker_col_palette_idx]], x_lim=pax, y_lim=pay, col_mode=v_params$col_mode, hl_mode=v_params$hl_mode, pdf_ofn=pdf_outfile_name);
        pdf_outfile_name <- character(0);
      }
      #plot_cov(ds, sel_ctgs=picked_pt_names, marker_idx=v_params$marker_idx, marker_col_palette=marker_col_palette, x_lim=pax, y_lim=pay, col_mode=v_params$col_mode, hl_mode=v_params$hl_mode);
      plot_cov(ds, sel_ctgs=picked_pt_names, map_idx=v_params$map_idx, marker_idx=v_params$marker_idx, map_col_palette=col_palettes_map[[v_params$map_col_palette_idx]], top_n_marker=TOP_N_MARKER, marker_col_palette=col_palettes[[v_params$marker_col_palette_idx]], x_lim=pax, y_lim=pay, col_mode=v_params$col_mode, hl_mode=v_params$hl_mode);
    }
    
    opt <- -1;
    cat(length(ds$clusters))
  }
  
  
  if(v_params$dat_mode != 0) { ds$pd <- update_pd(ds, dat_mode=0); }
  
  ds$clusters_working <- picked_pts;
  
  #eval.parent(substitute(ds<-ds));
  return(ds);
}




# dat_src
#   +$cov
#     +$dat
#     +$sel_col_idx
#   +$kmer
#   +$gc
#   +$seq_len
#   +$markers (1..n)
#     +$scmg
#     +$rbp
#     +$16s
#     +$tax
#   +$pd
#   +$clusters
#   +$clusters_working
#   +$cluster_names
#   +src
mn_data_import <- function()
{    
  display_msg("", mode=-1);
  display_msg(SEPARATOR);
  display_msg("Import data", mode=1);
  display_msg(SEPARATOR);
  separator = ","
  
  dat_src <- {}
  dat_src$cov <- {}
  dat_src$markers <- {}
  dat_src$maps <- {}
  dat_src$marker_src <- {}
  dat_src$pd <- {}
  
  wkdir <- getwd();
  display_msg(paste("Current working directory: ", wkdir, sep=""));
  
  sample_id <- ""
  if(nchar(default_sample_id) > 0) {
    sample_id <- default_sample_id;
  } else {
    sample_id <- readline("Please enter a sample id (>SAMPLE_ID<.cov.tsv): ");
  }
  
  # Assign sample id
  dat_src$sample_id <- sample_id;
  
  combined_infile <- paste(sample_id, "combined.tsv+++++", sep=".");
  if(file.exists(combined_infile)) {
    #     dat_src$cov_src <- combined_infile;
    #     combined_dat <- import_data(dat_src$cov_src);
    #     print(paste(names(combined_dat)));
  } else {
    cov_infile <- paste(sample_id, "cov.tsv", sep=".");
    
    fail_n <- 0;
    while(! file.exists(cov_infile) && ! file.exists(paste(cov_infile, "cov.tsv", sep="."))) {
      display_msg(paste("File, \"", cov_infile, "\", does not exist.", sep=""));
      if(fail_n > RETRY_N) { display_msg(paste("Unable to import any *.cov.tsv file, abort now.", sep="")); exit(0)}
      cov_infile <- readline("Please enter a coverage file (tab delimited): ");
      fail_n <- fail_n + 1;
    }
    
    dat_src$cov_src <- cov_infile;
    dat_src$cov$dat <- import_data(dat_src$cov_src);
    print_debug_msg(paste("mn_data_import().dat_src$cov$dat=", nrow(dat_src$cov$dat), "x", ncol(dat_src$cov$dat), sep=""));
    print_debug_msg(paste("mn_data_import().dat_src$cov$dat: ", paste(dat_src$cov$dat[1:3,1], collapse=", "), sep=""));
    print_debug_msg(paste("mn_data_import().dat_src$cov$dat.rownames: ", paste(rownames(dat_src$cov$dat)[1:3], collapse=", "), sep=""));
    
    cluster_fn <- paste(cov_infile, ".clusters", sep="");
    if(file.exists(cluster_fn)) {
      opt <- -1;
      while(opt != "Y" && opt != "N") {
        opt <- toupper(readline("Do you want to import previous cluster data? (Y/N) "));
      }
      
      if(opt == "Y") {
        clusters <- import_cluster_file(cluster_fn);
        dat_src$clusters <- clusters$clusters;
        dat_src$cluster_names <- clusters$cluster_names;
        dat_src$clusters_working <- clusters$working_clusters;
      }
    } else {
      dat_src$clusters <- {};
      dat_src$cluster_names <- {};
    }
    
    cov_n <- 0;
    while(cov_n < 2) {
      #cov_n <- readline("Please enter the number of coverage columns in the data: ");
      cov_n <- default_dimension;
    }
    dat_src$cov$sel_cov_col_idx <- seq(1, as.numeric(cov_n));
    
    display_msg(paste(length(dat_src$cov$sel_cov_col_idx), " columns will be used as coverages: ", paste(colnames(dat_src$cov$dat)[dat_src$cov$sel_cov_col_idx], collapse=", "), sep=""));
    display_msg(paste("Default columns for plotting ", length(dat_src$cov$sel_cov_col_idx), "-D coverage plot: ", paste(colnames(dat_src$cov$dat)[dat_src$cov$sel_cov_col_idx], collapse=", "), sep=""));
    if(length(dat_src$cov$sel_cov_col_idx) > 3) { display_msg(paste("Multi-dimensional scaling will be used to reduce the ", length(dat_src$cov$sel_cov_col_idx)," dimension.")); }
    display_msg(paste("Column for sequence length: ", TAG_SEQ_LEN, " [n=", length(dat_src$maps$gc), "]", sep=""));
    #display_msg(paste("Column for GC: ", TAG_GC, " [n=", length(dat_src$maps$gc), "]", sep=""));
    
    
    #dat_src$cov_col_idx <- seq(1, as.numeric(cov_n));
    
    #dat_src$cov$dat <- dat_src$data[, dat_src$cov$sel_cov_col_idx];
    #dat_src$cov$dat <- dat_src$data;
    
    #   if(nchar(default_gc_fn) > 0) {
    #     gc_infile <- default_gc_fn;
    #   } else {
    #     gc_infile <- "";
    #     gc_infile <- readline("Please enter a file with GC content (tab delimited): ");
    #   }
    
    gc_infile <- paste(sample_id, "gc.tsv", sep=".");
    while(! file.exists(gc_infile)) {
      display_msg(paste("File, \"", gc_infile, "\", does not exist.", sep=""));
      gc_infile <- readline("Please enter a file with GC content (tab delimited): ");
    }
    #dat_src$gc <- import_data(gc_infile);
    gc <- import_data(gc_infile);
    gc <- gc[match(names(gc), rownames(dat_src$cov$dat))];
    dat_src$maps[[TAG_GC]] <- gc;
    
    #gc <- import_data(gc_infile);
    #dat_src$gc <- gc[match(rownames(gc), rownames(dat_src$cov$dat)), 1];
    
    seqlen_infile <- paste(sample_id, "seq_len.tsv", sep=".");
    #   if(nchar(default_seqlen_fn) > 0) {
    #     seqlen_infile <- default_seqlen_fn;
    #   } else {
    #     seqlen_infile <- "";
    #     seqlen_infile <- readline("Please enter a file with sequence length (tab delimited): ");
    #   }
    while(! file.exists(seqlen_infile)) {
      display_msg(paste("File, \"", seqlen_infile, "\", does not exist.", sep=""));
      seqlen_infile <- readline("Please enter a file with sequence length (tab delimited): ");
    }
    seq_len <- import_data(seqlen_infile);
    dat_src$seq_len <- seq_len[match(names(seq_len), rownames(dat_src$cov$dat))];
    #dat_src$seq_len <- import_data(seqlen_infile);
    print_debug_msg(paste("mn_data_import().seq_len=", length(dat_src$seq_len), sep=""));
    #print_debug_msg(paste("mn_data_import().seq_len: ", paste(seq_len[1:3,1], collapse=", "), sep=""));
    #print_debug_msg(paste("mn_data_import().seq_len.rownames: ", paste(rownames(seq_len)[1:3], collapse=", "), sep=""));
    print_debug_msg(paste("mn_data_import().dat_src$seq_len: ", paste(names(dat_src$seq_len)[1:3], collapse=", "), sep=""));
    
    
    #dat_src$seq_len <- as.numeric(dat_src$data[, which(colnames(dat_src$data) == TAG_SEQ_LEN)]);
    #dat_src$gc <- dat_src$data[, TAG_GC];
    
    ##
    #   #dat_src$markers$tax_lineage <- dat_src$data[, which(colnames(dat_src$data) == TAG_TAX_LINEAGE)];
    #   if(nchar(default_tax_fn) > 0) {
    #     tax_infile <- default_tax_fn;
    #   } else {
    #     tax_infile <- "";
    #     tax_infile <- readline("Please enter a file with taxonomic lineage (tab delimited): ");
    #   }
    tax_infile <- paste(sample_id, "tax.tsv", sep=".");
    while(! file.exists(tax_infile)) {
      display_msg(paste("File, \"", tax_infile, "\", does not exist.", sep=""));
      tax_infile <- readline("Please enter a file with taxonomic lineage (tab delimited): ");
    }
    tax <- import_data(tax_infile);
    if(length(tax) > 0)
      dat_src$markers[[TAG_TAX_LINEAGE]] <- tax[match(names(tax), rownames(dat_src$cov$dat))];
    #dat_src$markers[[TAG_TAX_LINEAGE]] <- import_data(tax_infile);
    
    
    #  tax <- import_data(tax_infile);
    #  dat_src$markers[[TAG_TAX_LINEAGE]] <- tax[match(rownames(tax), rownames(dat_src$cov$dat)),1];
    #rownames(dat_src$markers[[TAG_TAX_LINEAGE]]) <- rownames(tax);
    #   display_msg(paste("Column for taxonomic lineage: ", TAG_TAX_LINEAGE, " [n=", nrow(dat_src$markers[[TAG_TAX_LINEAGE]]), "]", sep=""));  
    #   print_debug_msg(paste("First few markers: ", paste(rownames(tax)[1:5], collapse=", "), sep=""));
    #   tax2 <- dat_src$markers[[TAG_TAX_LINEAGE]]; 
    #   print_debug_msg(paste("First few markers: ", paste(rownames(tax2)[1:5], collapse=", "), sep=""));
    #   print_debug_msg(paste("First few markers: ", paste(rownames(dat_src$markers[[TAG_TAX_LINEAGE]])[1:5], collapse=", "), sep=""));
    
  }
  
  ##############################
  # Import auxiliary data files
  ##############################
  
  kmer_infile <- paste(sample_id, "kmer.tsv", sep=".");
  #   if(nchar(default_cov_fn) > 0) {
  #     kmer_infile <- default_kmer_fn;
  #   } else {
  #     kmer_infile <- "";
  #     kmer_infile <- readline("Please enter a kmer file (tab delimited): ");
  #   }
  while(! file.exists(kmer_infile)) {  
    display_msg(paste("File, \"", kmer_infile, "\", does not exist.", sep=""));
    kmer_infile <- readline("Please enter a kmer file (tab delimited): ");
  }
  dat_src$kmer_src <- kmer_infile; 
  kmer <- import_data(dat_src$kmer_src);
  dat_src$kmer <- kmer[match(rownames(kmer), rownames(dat_src$cov$dat)), ];
  #dat_src$kmer <- import_data(dat_src$kmer_src);
  
  #   min_contig_len <- -1
  #   min_contig_len <- DEFAULT_MIN_CONTIG_LEN;
  #   if(nchar(default_min_contig_len) > 0) {
  #     min_contig_len <- as.integer(default_min_contig_len);
  #   } else {
  #     min_contig_len <- readline(paste("Please enter minimum contig size to analyze [default=", min_contig_len,"]: ", sep=""));
  #   }
  
  if(nchar(default_min_contig_len) > 0) { min_contig_len <- as.integer(default_min_contig_len); } else { min_contig_len <- -1; }
  while(min_contig_len < DEFAULT_MIN_CONTIG_LEN) {
    min_contig_len <- as.integer(readline(paste("Please enter a new minimum contig length [integer >", DEFAULT_MIN_CONTIG_LEN,"]: ", sep="")));
    if (!is.na(min_contig_len)) {
      if(min_contig_len >= DEFAULT_MIN_CONTIG_LEN) { break; }
    }
    display_msg(paste(min_contig_len, " is not a valid value.", sep=""));
    min_contig_len <- -1;
  }
  MIN_CONTIG_LEN <<- min_contig_len;
  dat_src$min_ctg_len <- MIN_CONTIG_LEN; 
  
  
  #dat_src$sel_cov_col_idx <- c(dat_src$cov$sel_cov_col_idx[1], dat_src$cov$sel_cov_col_idx[2]);
  #dat_src$pd <- {}
  
  #dat_src$pd <- update_pd(dat_src, c(dat_src$cov$sel_cov_col_idx[1], dat_src$cov$sel_cov_col_idx[2]));
  dat_src$pd <- update_pd(dat_src);
  
  marker_fns <- list.files(wkdir, "*.markers.tsv$")[grep(paste(sample_id, ".", sep=""), list.files(wkdir, "*.markers.tsv$"))];
  # Import markers: 16S, RBP, ESOM, MaxBin
  if(length(marker_fns) > 0) {
    for(k in 1 : length(marker_fns)) {
      marker_fn <- marker_fns[k];
      display_msg(paste("Importing markers from ", marker_fn, sep=""));
      if(file.exists(marker_fn)) {  
        dat_src$marker_src <- c(dat_src$marker_src, marker_fn);
        markers <- import_data(marker_fn, convert_to_list = F);
        for(i in 1 : ncol(markers)) {
          marker_name <- colnames(markers)[i];
          display_msg(paste(i, "=", marker_name, sep=""));
          print_debug_msg(paste("Markers=", paste(rownames(markers)[1:5], collapse=", "), sep=""));
          print_debug_msg(paste("Markers=", paste(markers[1:5,i], collapse=", "), sep=""));
          dat_src$markers[[marker_name]] <- markers[, i];
          names(dat_src$markers[[marker_name]]) <- rownames(markers);
          print_debug_msg(paste("marker_name=", marker_name, sep=""));
          print_debug_msg(paste("Markers=", paste(names(dat_src$markers[[marker_name]])[1:5], collapse=", "), sep=""));
        }
      } 
    }    
  }
  
  display_msg(paste("Number of markers imported: ", length(dat_src$markers), sep=""));
  display_msg(paste("Markers: ", paste(names(dat_src$markers), collapse=", "), sep=""));
  
  # Import color mapping
  # GC, expression level, 
  #   col_map_fn <- gsub(".cov.", ".col_maps.", cov_infile)
  #   if(file.exists(marker_fn)) {  
  #     markers <- import_data(dat_src$kmer_src, convert_to_list = F); 
  #     for(i in 1 : ncol(markers)) {
  #       marker_name <- colnames(markers)[i];
  #       dat_src$maps[[marker_name]] <- markers[, i];
  #     }
  #   }
  
  maps_fns <- list.files(wkdir, "*.maps.tsv$")[grep(paste(sample_id, ".", sep=""), list.files(wkdir, "*.maps.tsv$"))];
  print_debug_msg(paste("maps_fns=", paste(maps_fns, collapse=", "), sep=""));
  print_debug_msg(paste("maps_fns=", length(maps_fns), sep=""));
  
  
  if(length(maps_fns) > 0) {
    # Import maps
    for(k in 1 : length(maps_fns)) {
      maps_fn <- maps_fns[k];
      display_msg(paste("Importing map from ", maps_fn, sep=""));
      if(file.exists(maps_fn)) {  
        dat_src$maps_src <- c(dat_src$maps_src, maps_fn);
        col_maps <- import_data(maps_fn, convert_to_list = F);
        for(i in 1 : ncol(col_maps)) {
          maps_name <- colnames(col_maps)[i];
          display_msg(paste(i, "=", maps_name, sep=""));
          print_debug_msg(paste("Markers=", paste(rownames(col_maps)[1:5], collapse=", "), sep=""));
          print_debug_msg(paste("Markers=", paste(col_maps[1:5,i], collapse=", "), sep=""));
          dat_src$maps[[maps_name]] <- col_maps[, i];
          names(dat_src$maps[[maps_name]]) <- rownames(col_maps);
          print_debug_msg(paste("maps_name=", maps_name, sep=""));
          print_debug_msg(paste("Map=", paste(names(dat_src$maps[[maps_name]])[1:5], collapse=", "), sep=""));
        }
      } 
    }
  }
  display_msg(paste("Number of maps imported: ", length(dat_src$maps), sep=""));
  display_msg(paste("Maps: ", paste(names(dat_src$maps), collapse=", "), sep="")); 
  
  
  display_msg(paste("Complete importing data", sep=""));
  
  return(dat_src);
}



mn_settings <- function(ds)
{
  update_pd <- {};
  update_pd$updating <- F;
  update_pd$MIN_CONTIG_LEN <- MIN_CONTIG_LEN;
  update_pd$SHOW_CLUSTER_BOUND_BOX <- SHOW_CLUSTER_BOUND_BOX;
  update_pd$sel_cov_col_idx <- ds$cov$sel_cov_col_idx;
  
  opt <- -1;
  ops <- c("C", "B", "L", "T", "M", "D", "N", "K", "V", "P", "Q");
  while(opt < 0 || opt > length(ops)) {
    picked_pt_names <- character(0);
    
    display_msg("", mode=-1);
    display_msg("", mode=2);
    display_msg("Settings Menu", mode=1);
    display_msg("", mode=2);
    display_msg("", mode=4);
    
    
    display_msg(paste("[C] Minimun contig length: ", MIN_CONTIG_LEN, sep=""), mode=3);
    display_msg(paste("[D] Set coverages to be plotted: ", sep=""), mode=3);
    display_msg("", mode=4); 
    display_msg(paste("> Current coverage: ", paste(colnames(ds$cov$dat)[ds$cov$sel_cov_col_idx], collapse=", "), sep=""), mode=3);    
    #display_msg(paste("Minimun contig length [C]: ", MIN_CONTIG_LEN, sep=""), mode=3);
    display_msg("", mode=4);  
    
    display_msg(paste("[B] Show cluster boundary: ", as.character(SHOW_CLUSTER_BOUND_BOX), sep=""), mode=3);
    display_msg(paste("[L] Show marker legend: ", as.character(SHOW_MARKER_LEGEND), sep=""), mode=3);
    display_msg(paste("[T] Enable transparent points: ", as.character(!NO_TRANSPARENT), sep=""), mode=3);
    display_msg(paste("[M] Show debug message: ", as.character(PRINT_DEBUG), sep=""), mode=3);
    display_msg("", mode=4);
    
    display_msg(paste("[N] Show marker from top: ", TOP_N_MARKER, sep=""), mode=3);
    display_msg("", mode=4);
    display_msg(paste("[K] Transform method for K-mer profile: ", ifelse(TRANSFORM_METHOD_KMER, "Correspondence Analysis (CA)", "Principle Components Analysis (PCA)"), sep=""), mode=3);
    display_msg(paste("[V] Transform method for multiple coverages (n>2): ", ifelse(TRANSFORM_METHOD_MDS, "Multidimensional Scaling (MDS)", "Principle Components Analysis (PCA)"), sep=""), mode=3);
    
    display_msg("", mode=4);
    
    display_msg(paste("[P] PDF size (W x H) ", PDF_WIDTH, " x ", PDF_HEIGHT, sep=""), mode=3);
    display_msg("", mode=4); 
    display_msg("[Q] Return to Main Menu", mode=3);
    
    display_msg("", mode=4);
    display_msg("", mode=4);
    display_msg("", mode=6);
    choice <- toupper(readline("Please enter an option: "));
    opt <- ifelse(!is.na(match(choice, ops)), match(choice,ops), -1);
    
    
    if(opt == which(ops == "C")) {
      new_ncl <- -1;
      repeat{
        display_msg("", mode=-1);
        new_ncl <- as.numeric(readline(paste("Please enter a new minimum contig length [integer >", DEFAULT_MIN_CONTIG_LEN,"]: ", sep="")));
        if (!is.na(new_ncl)) {
          if(new_ncl > DEFAULT_MIN_CONTIG_LEN) { break; }
        }
        display_msg(paste(new_ncl, " is not a valid value.", sep=""));
      }
      MIN_CONTIG_LEN <<- new_ncl;
      update_pd$updating <- T;
    }
    
    
    if(opt == which(ops == "N")) {
      new_top_n <- -1;
      repeat{
        display_msg("", mode=-1);
        new_top_n <- as.numeric(readline(paste("Please enter a value for top n marker to show [current: top ", TOP_N_MARKER,"]: ", sep="")));
        if (!is.na(new_top_n)) {
          if(new_top_n > 0) { break; }
        }
        display_msg(paste(new_top_n, " is not a valid value.", sep=""));
      }
      TOP_N_MARKER <<- new_top_n;
      update_pd$updating <- T;
    }
    
    
    if(opt == which(ops == "D")) {
      new_sel_cov_col_idx <- rep(F, ncol(ds$cov$dat));
      print_debug_msg(paste("Number of new_sel_cov_col_idx: ", new_sel_cov_col_idx, sep=""));
      for(i in 1 : ncol(ds$cov$dat)) {
        if(i %in% ds$cov$sel_cov_col_idx) { new_sel_cov_col_idx[i] <- T; }
      }
      print_debug_msg(paste("Number of new_sel_cov_col_idx==TRUE: ", length(which(new_sel_cov_col_idx == T)), sep=""));
      
      while(T) {
        display_msg("", mode=-1);
        display_msg(paste("Avaiable coverages:", sep=""));
        display_msg("", mode=6);
        display_msg("", mode=4);
        for(i in 1 : ncol(ds$cov$dat)) {
          cov_name <- colnames(ds$cov$dat)[i];
          if(new_sel_cov_col_idx[i] == T) {
            display_msg(paste(i, ") ", cov_name, " [Selected]", sep=""), mode=3);
          } else {
            display_msg(paste(i, ") ", cov_name, sep=""), mode=3);
          }
        }
        display_msg("", mode=4);
        if(length(which(new_sel_cov_col_idx == T)) > 1) { display_msg(paste("Q) Quit [Please select two or more coverages]", sep=""), mode=3); }
        display_msg("", mode=4);
        display_msg("", mode=6);
        
        cov_idx <- readline(paste("Please select a coverage [1-", ncol(ds$cov$dat), "]: ", sep=""));
        
        if (!is.na(as.numeric(cov_idx))) {
          cov_idx <- as.numeric(cov_idx);
          if(cov_idx > 0 && cov_idx <= ncol(ds$cov$dat)) {
            new_sel_cov_col_idx[cov_idx] <- !new_sel_cov_col_idx[cov_idx];
          } else {
            display_msg(paste("Invalid coverage index", sep=""), mode=3);
          }
        } else {
          if(toupper(cov_idx) == "Q" && length(which(new_sel_cov_col_idx == T)) > 1) { 
            if(length(which(new_sel_cov_col_idx == T)) > 1) { break; }
            display_msg(paste("Error: at least two coverages selected is required.", sep=""), mode=3);
          }
        }
      }
      update_pd$updating <- T;
      update_pd$sel_cov_col_idx <- which(new_sel_cov_col_idx == T);
    }
    
    if(opt == which(ops == "P")) { PRINT_LARGE_PDF <<- !PRINT_LARGE_PDF; }
    if(opt == which(ops == "V")) { TRANSFORM_METHOD_MDS <<- !TRANSFORM_METHOD_MDS; update_pd$updating <- T;}
    if(opt == which(ops == "K")) { TRANSFORM_METHOD_KMER <<- !TRANSFORM_METHOD_KMER; update_pd$updating <- T;}
    if(opt == which(ops == "B")) { SHOW_CLUSTER_BOUND_BOX <<- !SHOW_CLUSTER_BOUND_BOX; }
    if(opt == which(ops == "L")) { SHOW_MARKER_LEGEND <<- !SHOW_MARKER_LEGEND; }
    if(opt == which(ops == "T")) { NO_TRANSPARENT <<- !NO_TRANSPARENT; update_pd$updating <- T;}
    if(opt == which(ops == "M")) { PRINT_DEBUG <<- !PRINT_DEBUG; }    
    if(opt == which(ops == "Q")) { break; }
    
    if(PRINT_LARGE_PDF) { PDF_WIDTH <<- 16; PDF_HEIGHT <<- 16; } else { PDF_WIDTH <<- 8; PDF_HEIGHT <<- 8; }
    
    opt <- -1;
  }
  
  return(update_pd);
}


get_cl_color <- function(i)
{
  return(cl_cols[i]);
}


generate_cov_pairs <- function(ds, dimension=2)
{
  display_msg(paste("dimension=", dimension, sep="")); 
  display_msg(paste("Number of col: ", ncol(ds$cov$dat), sep="")); 
  cov_id_pairs <- t(combn(seq(1, ncol(ds$cov$dat)), dimension));
  for(i in nrow(cov_id_pairs))
  {   display_msg(paste(i, ") ", paste(cov_id_pairs[i, ], collapse=", "), sep=""));  }
  #cov_id_idx <- seq(1, ncol(ds$cov$dat))
  #cov_ids <- colnames(ds$cov$dat);
  #cov_id_pairs <- matrix(c(cov_ids[-length(cov_ids)], cov_ids[-1]), ncol=dimension);
  return(cov_id_pairs);
}

#
# Show exit menu
#
mn_exit <- function(ds)
{   
  display_msg("", mode=-1);
  #dev.off();
  #if(dev.list() > 0) { dev.off(); }
  
  # Existing data
  if(length(ds) == 0) { return(0); } 
  else {
    choice <- "S";
    while(choice != "C" && choice != "Y" && choice != "N" && choice != "Q") {  choice <- toupper(readline("Want to save cluster data? [Y]es/[N]o/[C]ancel ")); }
    if(choice == "N" || choice == "Q") { return(0); }
    if(choice == "C") { return(-1); }
    if(choice == "Y") { 
      display_msg("Saving..."); 
      save_clusters_to_file(ds);
      return(0); 
    }
  }
}

#
# 
# dat_mode:
# 0 - cov
# 1 - kmer
#
# vis_mode:
# 0 - 2d
# 1 - 3d (dat_mode=0 only)
#
# transform_method:
# 0 - PCA
update_pd <- function(ds, sel_ctgs=character(0), pt_weighted=T, log=T, vis_mode=0, dat_mode=0, transform_method_mds_cov=TRANSFORM_METHOD_MDS, transform_method_kmer=TRANSFORM_METHOD_KMER, min_seq_len=MIN_CONTIG_LEN)
  #update_pd <- function(ds, sel_cov_col_idx=character(0), sel_ctgs=character(0), pt_weighted=T, log=T, vis_mode=0, dat_mode=0, transform_method=1, min_seq_len=MIN_CONTIG_LEN)
{
  pd <- {};
  pd$mode <- DATA_MODE_COV;
  if(dat_mode == 1) { # CA on Kmer 
    ctgs_names <- rownames(ds$cov$dat);
    if(length(sel_ctgs) > 0) { ctgs_names <- sel_ctgs; }
    kmer <- ds$kmer[match(ctgs_names, rownames(ds$kmer)), ];
    seq_len <- ds$seq_len[match(ctgs_names, rownames(ds$cov$dat))];
    
    if(transform_method_kmer) {   
      print_debug_msg("CA on kmer");
      ca <- scores(cca(kmer), choices=1:5)$sites;
      pd$dat <- data.frame(x=ca[,"CA1"], y=ca[,"CA2"]); 
      pd$main_lab <- paste("Correspondence Analysis (PCA) on k-mers", sep="");
      pd$x_lab <- "CA1";
      pd$y_lab <- "CA2";
    } else {
      print_debug_msg("PCA on kmer");
      kmer_pca <- prcomp(kmer);
      pd$dat <- data.frame(x=kmer_pca$x[,1], y=kmer_pca$x[,2]);
      pd$main_lab <- paste("Principle Compoenents Analysis (PCA) on k-mers", sep="");
      
      vs <- apply(kmer_pca$x, 2, var);
      ps <- cumsum(vs / sum(vs));
      
      pd$x_lab <- paste("PCA1 (", format(round(ps[1] * 100, 1), nsmall = 2), "%)", sep="");
      pd$y_lab <- paste("PCA2 (", format(round(diff(ps)[1] * 100, 1), nsmall = 2), "%)" , sep="");
    }
    
    rownames(pd$dat) <- ctgs_names;
    #pd$sel_cov_col_idx <- "KMER";
    pd$sel_cov_col_idx <- ds$cov$sel_cov_col_idx;
    pd$mode <- DATA_MODE_KMER;
    
    if(pt_weighted) {
      if(length(seq_len) > 0) {
        pd$pt_size <- sqrt(seq_len)/100;
      }else{
        pd$pt_size <- 1.0;
      }
    }
    
    
  } else { # Normal coverage
    cat("Preparing plot...");
    
    if(length(ds$cov$sel_cov_col_idx) > 2) { # MDS: N->2
      cat("\nPerforming MDS...");
      
      pd$sel_cov_col_idx <- ds$cov$sel_cov_col_idx;
      
      # Filter by 
      coords <- ds$cov$dat[which(ds$seq_len >= MIN_CONTIG_LEN), pd$sel_cov_col_idx];
      
      #if(log) coords <- log(coords + 1);
      if(log) coords <- log_coord(coords);
      
      if(transform_method_mds_cov) {   
        print_debug_msg(paste("\nMDS on ", nrow(coords), " coords...", sep=""));
        coord_dist <- dist(coords, method="euclidean");
        cmd_fit <- cmdscale(coord_dist, k=2, eig=T);
        # Offset the coordinates from negatives
        pd$dat <- data.frame(x=(cmd_fit$points[,1] - min(cmd_fit$points[,1]) + 1), y=(cmd_fit$points[,2] - min(cmd_fit$points[,2]) + 1));
        rownames(pd$dat) <- rownames(cmd_fit$points);
        
        pd$main_lab <- paste("MDS on ", length(pd$sel_cov_col_idx), " dimensional data [", paste(colnames(ds$cov$dat)[ds$cov$sel_cov_col_idx], collapse=", "), "]", sep="");
        
        pd$x_lab <- "MDS1";
        pd$y_lab <- "MDS2";
      } else {
        coords_pca <- prcomp(coords, scaling=T);
        pd$dat <- data.frame(x=coords_pca$x[,1], y=coords_pca$x[,2]);
        
        pd$main_lab <- paste("PCA on ", length(pd$sel_cov_col_idx), " dimensional data [", paste(colnames(ds$cov$dat)[ds$cov$sel_cov_col_idx], collapse=", "), "]", sep="");
        
        vs <- apply(coords_pca$x, 2, var)  ;
        ps <- cumsum(vs / sum(vs));
        
        pd$x_lab <- paste("PCA1 (", format(round(ps[1]*100, 1), nsmall = 2), "%)", sep="");
        pd$y_lab <- paste("PCA2 (", format(round(diff(ps)[1]*100, 1), nsmall = 2), "%)" , sep="");
      }
      
      print_debug_msg(paste("pd$sel_cov_col_idx:", length(pd$sel_cov_col_idx), sep=""))
      
      #pd$dat <- data.frame(x=as.numeric(ds$cov$dat[, ds$cov$sel_cov_col_idx[sel_cov_col_idx[1]]]), y=as.numeric(ds$cov$dat[, ds$cov$sel_cov_col_idx[sel_cov_col_idx[2]]]), stringsAsFactors = F);
      #pd$dat <- data.frame(x=as.numeric(ds$cov$dat[, pd$sel_cov_col_idx[1]]), y=as.numeric(ds$cov$dat[, pd$sel_cov_col_idx[2]]), stringsAsFactors = F);
      
    } else { # Normal
      #if(length(sel_cov_col_idx) == 0) { sel_cov_col_idx <- ds$cov$sel_cov_col_idx; }
      pd$sel_cov_col_idx <- ds$cov$sel_cov_col_idx;
      pd$main_lab <- paste("Differential Coverages between ", paste(colnames(ds$cov$dat)[pd$sel_cov_col_idx], collapse=" and "), sep="");
      
      print_debug_msg(paste("pd$sel_cov_col_idx:", length(pd$sel_cov_col_idx), sep=""))
      
      
      #pd$dat <- data.frame(x=as.numeric(ds$cov$dat[, ds$cov$sel_cov_col_idx[sel_cov_col_idx[1]]]), y=as.numeric(ds$cov$dat[, ds$cov$sel_cov_col_idx[sel_cov_col_idx[2]]]), stringsAsFactors = F);
      pd$dat <- data.frame(x=as.numeric(ds$cov$dat[, pd$sel_cov_col_idx[1]]), y=as.numeric(ds$cov$dat[, pd$sel_cov_col_idx[2]]), stringsAsFactors = F);
      rownames(pd$dat) <- rownames(ds$cov$dat);
      pd$dat <- pd$dat[which(ds$seq_len >= MIN_CONTIG_LEN), ];
      
      # Convert d into plot_data
      #if(log) pd$dat <- log(pd$dat + 1);
      if(log) pd$dat <- log_coord(pd$dat);
      
      pd$x_lab <- colnames(ds$cov$dat)[ds$cov$sel_cov_col_idx[1]];
      pd$y_lab <- colnames(ds$cov$dat)[ds$cov$sel_cov_col_idx[2]];
    }
    
    #colnames(pd$dat) <- c("x", "y");
    #pd$y <- ds$cov$dat[, ds$cov$sel_cov_col_idx[sel_cov_col_idx[2]]];
    
    #if(length(sel_ctgs) > 0) { pd$dat <- pd$dat[which(rownames(pd$dat) %in% sel_ctgs),]; }
    
    if(pt_weighted) {
      if(length(ds$seq_len) > 0) {
        pd$pt_size <- sqrt(ds$seq_len[match(rownames(pd$dat), rownames(ds$cov$dat))])/100;
      }else{
        pd$pt_size <- 1.0;
      }
    }
  }
  
  
  # Filter by seq-len
  print_debug_msg(paste("\nDebug update_pd():"))
  print_debug_msg(paste("Before filtering: ", "update_pd().pd$dat=", nrow(pd$dat), ",", ncol(pd$dat), sep=""));
  print_debug_msg(paste("update_pd().min_seq_len=", MIN_CONTIG_LEN , sep=""));
  print_debug_msg(paste("update_pd().pd$dat$x=", paste(pd$dat$x[1:10], collapse=", "), sep=""));
  print_debug_msg(paste("update_pd().pd$dat$y=", paste(pd$dat$y[1:10], collapse=", "), sep=""));
  #pd_seq_len <- unlist(ds$seq_len[match(rownames(pd$dat), rownames(ds$cov$dat))]);
  pd_seq_len <- unlist(ds$seq_len[match(rownames(pd$dat), rownames(ds$cov$dat))]);
  
  print_debug_msg(paste("update_pd().pd_seq_len.length=", length(pd_seq_len), sep=""));
  print_debug_msg(paste("update_pd().pd_seq_len=", paste(pd_seq_len[1:10], collapse=", "), sep=""));
  pd$dat <- pd$dat[which(pd_seq_len >= MIN_CONTIG_LEN ), c("x","y")];
  print_debug_msg(paste("update_pd().pd$dat$x=", paste(pd$dat$x[1:10], collapse=", "), sep=""));
  print_debug_msg(paste("update_pd().pd$dat$y=", paste(pd$dat$y[1:10], collapse=", "), sep=""));
  print_debug_msg(paste("update_pd().pd$dat.colnames=", paste(colnames(pd$dat), collapse=", "), sep=""));
  print_debug_msg(paste("After filtering: ", "update_pd().pd$dat=", nrow(pd$dat), ",", ncol(pd$dat), sep=""));
  
  #pd$seq_len <- pd$seq_len[pd_seq_len >= min_seq_len];
  #print_debug_msg(paste("pd$seq_len=", length(pd$seq_len), sep=""));
  
  #pd$pt_size <- pd$pt_size[pd_seq_len >= MIN_CONTIG_LEN ];
  
  print_debug_msg(paste("update_pd().pd$pt_size=", length(pd$pt_size), sep=""));
  
  
  return(pd)
}


log_coord <- function(pd_dat) { return(log(pd_dat + 1)); }



# 
# update_pd <- function(ds, sel_cov_col_idx=character(0), pt_weighted=T, log=T, vis_mode=0, dat_mode=0)
# {
#   pd <- {};
#   pd$sel_cov_col_idx <- sel_cov_col_idx;
#   pd$dat <- data.frame(x=as.numeric(ds$cov$dat[, ds$cov_col_idx[sel_cov_col_idx[1]]]), y=as.numeric(ds$cov$dat[, ds$cov_col_idx[sel_cov_col_idx[2]]]), stringsAsFactors = F);
#   rownames(pd$dat) <- rownames(ds$cov$dat);
#   
#   #colnames(pd$dat) <- c("x", "y");
#   #pd$y <- ds$cov$dat[, ds$cov_col_idx[sel_cov_col_idx[2]]];
#   pd$x_lab <- colnames(ds$cov$dat)[ds$cov_col_idx[sel_cov_col_idx[1]]];
#   pd$y_lab <- colnames(ds$cov$dat)[ds$cov_col_idx[sel_cov_col_idx[2]]];
#   
#   
#   if(dat_mode == 1)
#   {
#     # kmer
#     pd$main_lab <- paste("Correspondance Aanlysis (CA) on k-mers", sep="");
#   } else {
#     # cov
#     pd$main_lab <- paste(x_lab, " vs ", y_lab, sep="");
#   }
#   
#   # Convert d into plot_data
#   if(log) pd$dat <- log(pd$dat + 1);
#   
#   if(pt_weighted)
#   {
#     if(length(ds$seq_len) > 0) {
#       pd$pt_size <- sqrt(ds$seq_len)/100;
#     }else{
#       pd$pt_size <- 1.0;
#     }
#   }
#   
#   return(pd)
# }

#
# Generate coverage plot in 2D
# hl_mode: 
# 0 - sel_ctgs
# 1 - clusters
#
# col_mode:
# 0 - No coloring
# 1 - GC
# 
#plot_cov <- function(d, main=character(0), cl=character(0), show_cl_id=FALSE, sel_ctgs=character(0), x_lab=character(0), y_lab=character(0), x_lim=character(0), y_lim=character(0), zm=FALSE, log=LOG_ENABLED, mode=0)
plot_cov <- function(d, main=character(0), sel_cls=character(0), sel_ctgs=character(0), map_idx=0, marker_idx=0, top_n_marker=character(0), map_col_palette=character(0), marker_col_palette=character(0), x_lab=character(0), y_lab=character(0), x_lim=character(0), y_lim=character(0), zm=FALSE, log=LOG_ENABLED, col_mode=0, hl_mode=0, dat_mode=0, pdf_ofn=character(0), pdf_w=PDF_WIDTH, pdf_h=PDF_HEIGHT)
{
  par("plt")
  #par(oma = c(4, 1, 1, 1));
  
  if(length(pdf_ofn) > 0) {
    display_msg(paste("Exporting plot to PDF file", pdf_ofn, "...", sep=""));
    pdf(pdf_ofn, width=pdf_w, height=pdf_h);
  }
  
  
  color_ramp_step <- 10;
  
  # Convert d into plot_data
  #pd <- d$cov$dat[,1:2];
  #pd <- data.frame(x=d$cov$dat[,x_idx], y=d$cov$dat[,y_idx], stringsAsFactors=F);
  pd <- d$pd;
  pd_dat <- pd$dat;
  #gc <- d$maps$gc[match(rownames(pd_dat), rownames(d$cov$dat))];
  
  
  print_debug_msg(paste("plot_cov().pd_dat=", nrow(pd_dat), ", ", ncol(pd_dat), sep=""));
  print_debug_msg(paste("plot_cov().pd_dat.colnames=", paste(colnames(pd_dat), collapse=", "), sep=""));
  
  y_lab <- pd$y_lab;
  x_lab <- pd$x_lab;
  
  # Main title
  main_lab <- pd$main_lab;
  
  print_debug_msg(paste("plot_cov().main_lab=", main_lab, sep=""));
  print_debug_msg(paste("plot_cov().x_lim=", length(y_lim), sep=""));
  print_debug_msg(paste("plot_cov().y_lim=", length(x_lim), sep=""));
  
  step <- 100;
  if(length(y_lim) == 0) { y_lim<-c(min(pd_dat$y) - (diff(range(pd_dat$y)) / step), max(pd_dat$y) + (diff(range(pd_dat$y)) / step)); }
  if(length(x_lim) == 0) { x_lim<-c(min(pd_dat$x) - (diff(range(pd_dat$x)) / step), max(pd_dat$x) + (diff(range(pd_dat$x)) / step)); }
  print_debug_msg(paste("plot_cov().x=", paste(pd_dat$x[1:10], collapse=", "), sep=""));
  print_debug_msg(paste("plot_cov().y=", paste(pd_dat$y[1:10], collapse=", "), sep=""));
  print_debug_msg(paste("plot_cov().x_range=", paste(range(pd_dat$x), collapse=", "), sep=""));
  print_debug_msg(paste("plot_cov().y_range=", paste(range(pd_dat$y), collapse=", "), sep=""));
  print_debug_msg(paste("plot_cov().x_lim=", paste(x_lim, collapse=", "), sep=""));
  print_debug_msg(paste("plot_cov().y_lim=", paste(y_lim, collapse=", "), sep=""));
  
  
  default_pt_col <- "lightgray";
  default_sel_ctgs_col <- "red";
  
  sel_ctgs_col <- default_sel_ctgs_col;
  
  print_debug_msg(paste("plot_cov().marker_idx: ", marker_idx, sep=""));
  
  # Color with GC content
  if(map_idx > 0) { 
    map_id <- names(d$maps)[map_idx];
    
    #colorPalette <- colorRampPalette(c("darkred", "yellow", "steelblue"));
    if(length(map_col_palette) != 0) {
      colorPalette <- map_col_palette;
    } else {
      colorPalette <- col_palettes[[7]];
    }
    
    # Generate colors for all points
    pt_cols <- rep(default_pt_col, nrow(d$cov$dat));
    
    if(length(d$maps[[map_id]]) > 0) {
      map_pt_cols <- colorPalette(color_ramp_step)[cut(d$maps[[map_id]], breaks = color_ramp_step)];
      #pt_cols[which(rownames(d$cov$dat) %in% names(d$maps[[map_id]]))] <- map_pt_cols[which(rownames(d$cov$dat) %in% names(d$maps[[map_id]]))];
      print_debug_msg(paste("map_pt_cols=", length(map_pt_cols), sep=""));
      print_debug_msg(paste("nrow(d$cov$dat)=", nrow(d$cov$dat), sep=""));
      print_debug_msg(length(which(!is.na(match(rownames(d$cov$dat), names(d$maps[[map_id]]))))));
      matched_ids <- match(rownames(d$cov$dat), names(d$maps[[map_id]]));
      pt_cols[which(! is.na(matched_ids))] <- map_pt_cols[matched_ids[which(! is.na(matched_ids))]];
      pt_cols <- pt_cols[match(rownames(pd_dat), rownames(d$cov$dat))];
    }
    
    
    # Belowing lines are working
    #pt_cols <- colorPalette(color_ramp_step)[cut(d$maps[[map_id]], breaks = color_ramp_step)];
    # Have to apply this line to fix kmer GC coloring
    # Subset the colors to fit pd_dat
    #pt_cols <- pt_cols[match(rownames(pd_dat), rownames(d$cov$dat))];
    
    
    if(col_mode == 2) {
      if(length(sel_ctgs) > 0) {  
        pt_cols[! rownames(pd_dat) %in% sel_ctgs] <- default_pt_col;
      }
    }
  } else {
    pt_cols <- "black";
  }
  
  
  #   # Color with GC content
  #   if(col_mode == 1 || col_mode == 2) {
  #     if(length(marker_id) > 0) {
  #       marker_name <- marker_id;
  #     } else {
  #       marker_name <- TAG_GC;
  #       marker_id <- TAG_GC;
  #     }
  #     
  #     colorPalette <- colorRampPalette(c("darkred", "yellow", "steelblue"));
  #     # Generate GC colors for all points
  #     #pt_cols <- colorPalette(color_ramp_step)[cut(d$gc, breaks = color_ramp_step)];
  #     pt_cols <- colorPalette(color_ramp_step)[cut(d$maps[[marker_id]], breaks = color_ramp_step)];
  #     
  #     # Have to apply this line to fix kmer GC coloring
  #     # Subset the colors to fit pd_dat
  #     pt_cols <- pt_cols[match(rownames(pd_dat), rownames(d$cov$dat))];
  #  
  #     # Only color selected points
  #     if(col_mode == 2) {
  #       if(length(sel_ctgs) > 0) {  
  #         pt_cols[! rownames(pd_dat) %in% sel_ctgs] <- default_pt_col;
  #       }
  #     }
  #     
  # #     if(SHOW_MARKER_LEGEND) {
  # #       marker_lbls <- paste(top_n_markers, " (n=", top_n_marker_ns,")", sep="");
  # #       legend("topright", marker_lbls, col=colorPalette, pch=marker_pch, cex=0.6, box.col=adjustcolor("black", 0.1),bg=adjustcolor("white", 0.8), title=legend_marker_title);
  # #     }
  #   } else {
  #     pt_cols <- "black";
  #   }
  
  # Fix for some devices without transparent supports
  if(! NO_TRANSPARENT) {
    if(col_mode == 0 && map_idx == 0) { pt_cols <- adjustcolor(pt_cols, 0.2); } else { pt_cols <- adjustcolor(pt_cols, 0.6); }
    sel_ctgs_col <- adjustcolor(sel_ctgs_col, 0.4);
    #pt_col <- rgb(0.2,0.2,0.2,0.1);
    #sel_ctgs_col <- rgb(1,0,0,0.4);
  }
  
  print_debug_msg("Rendering...")
  
  # Generate the plot
  if(length(d$cov$sel_cov_col_idx) == 2 && pd$mode == DATA_MODE_COV) {
    # Take care of the axis labels
    print(plot(x=pd_dat$x, y=pd_dat$y, cex=pd$pt_size, xlim=x_lim, ylim=y_lim, ylab=y_lab, xlab=x_lab, xaxt="n", yaxt="n", xaxs="i", yaxs="i", main=main_lab, col=pt_cols, pch=19));
    
    tick_n <- 10;
    ticks <- seq(0, tick_n-1);
    if(log) { axis_labels <- (10^seq(0, tick_n, 0.5)[1:tick_n]); axis_labels[seq(2, tick_n, 2)] <- ""; } else { axis_labels <- ticks; }
    axis(1, at=ticks, labels=axis_labels);
    axis(2, at=ticks, labels=axis_labels);
  } else {
    print(plot(x=pd_dat$x, y=pd_dat$y, cex=pd$pt_size, xlim=x_lim, ylim=y_lim, ylab=y_lab, xlab=x_lab, xaxs="i", yaxs="i", main=main_lab, col=pt_cols, pch=19));
    
    if(pd$mode == DATA_MODE_KMER) {
      abline(h=0, col="lightgrey", lty=2);
      abline(v =0, col="lightgrey", lty=2);  
    }
  }
  
  
  # Highlight clusters
  if(hl_mode == 1) {
    print_debug_msg("Showing clusters: ")
    print_debug_msg(paste("Number of clusters: ", length(d$clusters), "\n", sep=""))
    
    if(length(d$clusters) > 0) {
      #cl_color_palette <- colors()[ceiling(runif(10, 1, length(colors()) -1))];
      for(i in 1 : length(d$clusters)) {
        cl_id <- names(d$clusters)[i];
        cl_name <- d$cluster_names[cl_id];
        
        # Exclude unselected cluster
        if(length(sel_cls)>0) { if(! cl_id %in% sel_cls) { next; } }
        
        cl_ctgs <- unlist(d$clusters[cl_id]);
        cl_pts <- d$pd$dat[match(cl_ctgs, rownames(pd_dat)), ];
        print_debug_msg(paste("cl_pts: ", nrow(cl_pts), "\n", sep=""));
        
        cl_pt_size <- pd$pt_size[match(cl_ctgs, rownames(pd_dat))];
        #cl_col <- cl_color_palette[i];
        cl_col <- get_cl_color(i);
        
        #points(x=cl_pts$x, y=cl_pts$y, cex=cl_pt_size, xlim=x_lim, ylim=y_lim, xaxs="i", yaxs="i", col=cl_col, pch=21);
        points(x=cl_pts$x, y=cl_pts$y, cex=cl_pt_size, xlim=x_lim, ylim=y_lim, xaxs="i", yaxs="i", col=cl_col, pch=21);
        
        if(SHOW_CLUSTER_BOUND_BOX) {
          xr <- range(cl_pts$x); yr <- range(cl_pts$y);
          rect(xr[1], yr[1], xr[2], yr[2], border=cl_col, lty="dashed");
        }
        
        # clip the points to the current view for calculating the centroid of the cluster
        #cl_pts <- cl_pts[]
        if(length(x_lim) > 0 && length(y_lim) > 0) { clipped_cl_pts <- clip_pts_from_view(cl_pts, c(x_lim, y_lim)); if(nrow(clipped_cl_pts)>0) {cl_pts <- clipped_cl_pts; } };
        
        cx <- mean(cl_pts$x); cy <- max(cl_pts$y) + 0.1;
        text(cx, cy, cl_name, col=cl_col);
      }
    }
  }
  
  # Highlight selected ctgs
  if(length(sel_ctgs) > 0) {
    print_debug_msg(paste("Number of selected points: ", length(sel_ctgs), sep=""))
    sel_pd_dat <- pd_dat[match(sel_ctgs, rownames(pd_dat)),];
    pt_size <- pd$pt_size[match(sel_ctgs, rownames(pd_dat))]
    print(points(x=sel_pd_dat$x, y=sel_pd_dat$y, cex=pt_size, xlim=x_lim, ylim=y_lim, col=sel_ctgs_col, pch=21));
  }    
  
  
  # Show markers
  if(marker_idx > 0) {
    marker_id <- names(d$markers)[marker_idx];
    
    print_debug_msg(paste("plot_cov().marker_id=", marker_id, sep=""));
    
    legend_marker_title <- marker_id;
    marker_pch <- 17; 
    print_debug_msg(paste("Preparing marker: ", marker_id, sep=""));
    
    # Retrieve a list of selected markers
    markers <- d$markers[[marker_id]];
    print_debug_msg(paste("Number of selected markers: ", length(markers), sep=""));
    
    # Subset the selected markers by their length
    markers <- markers[which(d$seq_len[match(names(markers), names(d$seq_len))] > MIN_CONTIG_LEN)];
    # Discard markers with "NA" value
    markers <- markers[!is.na(markers)];
    print_debug_msg(paste("Number of selected markers: ", length(markers), sep=""));
    print_debug_msg(paste("First few markers: ", paste(names(markers)[1:5], collapse=", "), sep=""));
    
    # Filtering top n marker
    # select the top N marker if a value is provided.
    if(length(top_n_marker) > 0) { top_n <- top_n_marker; } else { top_n <- DEFAULT_TOP_N_MARKER; }
    
    # Get the non-redundant values/names of top n markers
    nr_markers <- unlist(unique(markers));
    # Calculate the top N markers by their occurrence
    top_n_marker_ns <- sort(table(markers), decreasing=T)[1:min(top_n, length(nr_markers))];
    top_n_marker_list <- names(sort(table(markers), decreasing=T)[1:min(top_n, length(nr_markers))]);
    # Generate color palette for markers
    if(length(marker_col_palette) == 0) { 
      sel_marker_col_palette <- colors()[ceiling(runif(length(top_n_marker_list), 1, length(colors())))];  
    } else { 
      sel_marker_col_palette <- marker_col_palette(length(top_n_marker_list)); 
    }
    
    print_debug_msg(paste("Number of selected non-redundant markers: ", length(top_n_marker_list), sep=""));
    print_debug_msg(paste("Non-redundant markers: ", paste(top_n_marker_list, collapse=", "), sep=""));
    
    # Subset the markers with a value enlisted on the top_n_marker list
    sel_markers <- markers[which(markers %in% top_n_marker_list)]
    sel_marker_ctgs <- names(sel_markers);
    # Obtain the coordinates of the marker subset
    sel_marker_xy <- pd_dat[match(sel_marker_ctgs, rownames(pd_dat)),];
    # Assign color to the marker subset
    sel_marker_cols <- sel_marker_col_palette[match(sel_markers, top_n_marker_list)];
    
    print_debug_msg(paste("Number of contigs belonging to the selected markers: ", length(sel_marker_ctgs), sep=""));
    print_debug_msg(paste("Number of coordinates belonging to the selected markers: ", nrow(sel_marker_xy), sep=""));
    
    # Debug
    for(i in 1 : length(sel_markers)) {
      print_debug_msg(paste(names(sel_markers)[i], "=", sel_markers[i], "\n"))
    }
    
    marker_pt_size <- 0.8; 
    print_debug_msg(paste("Number of colors for the points belonging to the selected markers: ", length(sel_marker_col_palette), sep=""));
    
    #colorPalette <- colorRampPalette(c("darkred", "yellow", "steelblue"));
    #pt_col <- colorPalette(color_ramp_step)[as.numeric(cut(d$gc, breaks = color_ramp_step))]
    
    points(x=sel_marker_xy$x, y=sel_marker_xy$y, cex=marker_pt_size, xlim=x_lim, ylim=y_lim, xaxs="i", yaxs="i", col=sel_marker_cols, pch=marker_pch);
    
    if(SHOW_MARKER_LEGEND) {
      marker_lbls <- paste(top_n_marker_list, " (n=", top_n_marker_ns,")", sep="");
      
      #par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      legend("bottomright", marker_lbls, col=sel_marker_col_palette, pch=marker_pch, cex=0.6, box.col=adjustcolor("black", 0.1),bg=adjustcolor("white", 0.8), title=legend_marker_title);
    }
    
  }
  
  #   
  #   if(length(cl) > 0)
  #   {
  #     cl_ids <- names(cl);
  #     require(reshape2); cl <- melt(cl);
  #     cl_n <- length(cl_ids);
  #     sel_pd <- subset(pd, rownames(pd) %in% cl[, 1]);
  #     display_msg(paste("sel_pd: ", nrow(sel_pd), "\n", sep=""));
  #     
  #     sel_pd$cl_id <- cl[match(cl[,1 ], rownames(sel_pd)), 2];
  #     #sel_pd$pt_size <- d$cov$SEQ_LEN[match(cl[,1], rownames(d$cov))];
  #     with(sel_pd, points(x=x, y=y, cex=pt_size, xlim=x_lim, ylim=y_lim, col=sel_ctgs_col, pch=21));
  #     if(show_cl_id)
  #     {
  #       for(i in 1 : cl_n)
  #       {
  #         cl_id <- cl_ids[i];
  #         pts <- pd[match(rownames(sel_pd), rownames(pd)), ];
  #         cat(paste(cl_id, ": ", nrow(pts), "\n", sep=""));
  #         cx <- mean(pts$x); cy <- mean(pts$y);
  #         cat(paste(cl_id, ": ", cx, ", ", cy, "\n", sep=""))
  #         text(cx, cy, cl_id, col=sel_ctgs_col);
  #       }
  #     }
  #   }
  #   
  #   if(mode>=1) { zm(); } 
  #   
  #   if(mode==2) { 
  #     display_msg("Please select 6 points that form a region covering contigs of interest: ");
  #     bin_def <- NULL;
  #     bin_def <- ahull(locator(6, type="p", pch=20), alpha=100);
  #     print(plot(bin_def, add=T, col="black"));
  #     
  #     sel_pts <- {};
  #     for(i in 1:nrow(pd)) { if(inahull(bin_def, c(pd$x[i], pd$y[i]))) sel_pts <- rbind(sel_pts, rownames(pd)[i])};
  #     
  #     display_msg(paste("Number of points selected: ", length(sel_pts), sep=""));
  #   } else if(mode==3) { zm(); pick_cluster(pd); } else {}
  #   
  
  # Legend
  if(length(pdf_ofn) > 0) {
    dev.off();
  }
}

plot_cluster_diagram <- function(ds, x_lim=character(0), y_lim=character(0), log=LOG_ENABLED, pdf_ofn=character(0), pdf_w=PDF_WIDTH, pdf_h=PDF_HEIGHT)
{
  #dev.off();
  
  cl_pd <- data.frame(name=rep("", length(ds$clusters)), x=rep(0, length(ds$clusters)), y=rep(0, length(ds$clusters)), pt_size=rep(0, length(ds$clusters)), count=rep(0, length(ds$clusters)), pt_col=rep("black", length(ds$clusters)), stringsAsFactors = F);
  for(i in 1 : length(ds$clusters)) {
    cl_id <- names(ds$clusters)[i];
    cl_name <- ds$cluster_names[cl_id];
    
    print_debug_msg(paste(i, ": ", cl_name, sep=""));
    
    cl_ctgs <- unlist(ds$clusters[cl_id]);   
    cl_ctg_n <- length(cl_ctgs);
    
    cl_pts <- ds$pd$dat[match(cl_ctgs, rownames(ds$pd$dat)), ];
    cl_col <- "black";
    
    seq_lens <- ds$seq_len[match(cl_ctgs, rownames(ds$cov$dat))];
    mean_seq_len <- mean(seq_lens);
    
    # Centroid
    if(PRINT_DEBUG) {
      for(j in 1 : length(cl_ctgs)) { 
        print_debug_msg(paste(cl_ctgs[j], ", cl_pts$x=", cl_pts$x[j], ", cl_pts$y=", cl_pts$y[j], sep=""));
      }
    }
    
    cx <- mean(cl_pts$x[which(!is.na(cl_pts$x) & !is.na(cl_pts$y))]); cy <- mean(cl_pts$y[which(!is.na(cl_pts$x) & !is.na(cl_pts$y))]);
    
    cl_pd$name[i] <- cl_name;
    cl_pd$x[i] <- cx;
    cl_pd$y[i] <- cy;
    cl_pd$seq_lens[i] <- list(seq_lens);
    cl_pd$mean_seq_len[i] <- mean_seq_len;
    cl_pd$pt_size[i] <- sqrt(mean_seq_len)/50;
    cl_pd$count[i] <- cl_ctg_n;
    cl_pd$cl_col[i] <- adjustcolor("darkgrey", 0.6);
    
    print_debug_msg(paste("x=", cx, ", y=", cy, sep=""));
  }   
  
  #text(cx, cy, cl_name, col=cl_col);
  #points(x=cx, y=cy, cex=marker_pt_size, xaxs="i", yaxs="i", col=cl_col, pch=20);
  
  
  if(length(pdf_ofn) > 0) {
    pdf(pdf_ofn, width=pdf_w, height=pdf_h);
  }
  
  y_lab <- ds$pd$y_lab;
  x_lab <- ds$pd$x_lab;
  main_lab <- ds$pd$main_lab;
  print(plot(x=cl_pd$x, y=cl_pd$y, cex=cl_pd$pt_size, ylab=y_lab, xlab=x_lab, xaxt="n", yaxt="n", xlim=x_lim, ylim=y_lim, xaxs="i", yaxs="i", main=main_lab, col=cl_pd$cl_col, pch=19));
  
  paxy <- par("usr");
  abline(h=mean(paxy[3:4]), col="lightgrey", lty=2);
  abline(v =mean(paxy[1:2]), col="lightgrey", lty=2);
  
  for(i in 1 : nrow(cl_pd)) {
    #text(x=cl_pd$x[i], y=cl_pd$y[i], labels=paste(cl_pd$name[i], "\nn=", length(cl_pd$seq_lens[[i]]), "\nMean=", format(round(cl_pd$mean_seq_len[i]/1000, 1), nsmall = 2), "kb", sep=""), col="black", cex=0.7);
    text(x=cl_pd$x[i], y=cl_pd$y[i], labels=paste(cl_pd$name[i], " (", length(cl_pd$seq_lens[[i]]), ")", sep=""), col="black", cex=1);
  }
  
  tick_n <- 10;
  ticks <- seq(0, tick_n-1);
  if(log) { axis_labels <- (10^seq(0, tick_n, 0.5)[1:tick_n]); axis_labels[seq(2, tick_n, 2)] <- ""; } else { axis_labels <- ticks; }
  
  axis(1, at=ticks, labels=axis_labels);
  axis(2, at=ticks, labels=axis_labels);
  
  # Legend
  if(length(pdf_ofn) > 0) {
    dev.off();
  }
}




#
#
#
pick_points <- function(pd, msg=character(0))
{
  if(length(msg) == 0) { display_msg("Please select 6 points that form a region covering contigs of interest: "); } else { display_msg(msg); };
  
  sel_pts <- {};
  bin_def <- NULL;
  bin_def <- tryCatch({
    bf <- ahull(locator(6, type="p", pch=20), alpha=100);
    #return(bf);
  }, warning=function(w) {
    #display_msg(w);
    print("Warning: illegal operation.");
    NULL;
  }, error=function(e) {
    #display_msg(e);
    print("Error: Illegal operation.");
    NULL;
  }
  );
  
  #bin_def <- ahull(locator(6, type="p", pch=20), alpha=100);
  #plot(bin_def, add=T, col="black");
  
  if(!is.null(bin_def)) { 
    plot(bin_def, add=T, col="black");
    for(i in 1:nrow(pd)) { if(inahull(bin_def, c(pd$x[i], pd$y[i]))) sel_pts <- rbind(sel_pts, rownames(pd)[i])};
    print_debug_msg(paste("Number of points: ", length(sel_pts)));
  }
  
  return(sel_pts);
}


# Return the most left most, right most, lowest and highest points of selected points
get_ext_pts <- function(pd)
{
  # Not implement yet. For estimating the x,y,w and h of bound box
}

#
#
#
create_new_cluster <- function(ds, my_sel_ctgs)
{  
  sel_ctgs <- my_sel_ctgs;
  
  display_msg("", mode=-1);
  display_msg(SEPARATOR);
  display_msg("Create a new cluster", mode=1);
  display_msg(SEPARATOR);  
  display_msg("", mode=4);
  display_msg(paste("> Number of points picked: ", length(sel_ctgs), sep=""), mode=5);
  display_msg(paste("> Total size (bp): ", sum(ds$seq_len[match(sel_ctgs, rownames(ds$cov$dat))]), sep=""), mode=5);
  display_msg("", mode=4);
  display_msg("", mode=4);
  
  
  #   if(length(d$clusters) > 0)
  #   {
  #     cl_color_palette <- colors()[ceiling(runif(10, 1, length(colors()) -1))];
  #     for(i in 1 : length(d$clusters))
  #     {
  #       cl_id <- names(d$clusters)[i];
  #       cl_name <- d$cluster_names[cl_id];
  #       cl_ctgs <- unlist(d$clusters[cl_id]);
  #       #cat(paste("cl_ctgs: ", paste(cl_ctgs, sep=", ", collapse=", "), "\n", sep="")); 
  #       cl_pts <- d$pd$dat[match(cl_ctgs, rownames(d$pd$dat)), ];
  #       #cat(paste("cl_pts: ", nrow(cl_pts), "\n", sep=""));
  #       
  #       cl_pt_size <- d$pd$pt_size[match(cl_ctgs, rownames(d$pd$dat))];
  #       cl_col <- cl_color_palette[i];
  #       #cat(paste("cl_pt_size: ", nrow(cl_pt_size), "\n", sep=""));
  #       
  #       points(x=cl_pts$x, y=cl_pts$y, cex=cl_pt_size, xlim=x_lim, ylim=y_lim, xaxs="i", yaxs="i", col=cl_col, pch=21);
  #       
  #       cx <- mean(cl_pts$x); cy <- max(cl_pts$y) + 0.1;
  #       #print(paste(cl_name, ": ", cx, ", ", cy, sep=""))
  #       text(cx, cy, cl_name, col=cl_col);
  #     }
  #   }
  #  
  cluster_list <- {}
  cluster_list$clusters <- ds$clusters;
  cluster_list$cluster_names <- ds$cluster_names;
  removed_clusters <- {};
  removed_cluster_names <- {};
  
  # If clusters exist in the cluster list
  if(length(cluster_list$clusters) > 0) { 
    binned_ctgs <- unlist(cluster_list$clusters);
    ovlp_ctgs <- c(binned_ctgs, sel_ctgs)[duplicated(c(binned_ctgs, sel_ctgs))];
    #ovlp_ctgs <- sel_ctgs[which(binned_ctgs %in% sel_ctgs)];
    
    if(length(ovlp_ctgs) > 0) {
      opt <- -1;
      while(opt != "Y" && opt != "N") {
        display_msg("", mode=-1);
        display_msg(paste("Warning: ", paste(length(ovlp_ctgs), " points have been assigned to other clusters.", sep="")))
        display_msg("What do you want to do with the overlapping point(s)? ")
        display_msg("", mode=4);
        display_msg("> Remove the points from other clusters and add them to this new cluster [Y]")
        display_msg("> Exclude the points from the new cluster [N]")
        display_msg("", mode=4); 
        display_msg("", mode=4);
        display_msg("", mode=6);
        
        opt <- toupper(readline("Y/N? "));
      }
      
      # Remove selected contigs from other clusters
      if(opt == "Y") {
        # Start from new lists
        cluster_list$clusters <- {};
        cluster_list$cluster_names <- {};
        
        # Iterate all the existing clusters
        for(i in 1:length(ds$clusters)) {
          # Get the ID and name of the current cluster
          cl_id <- names(ds$cluster_names)[i];
          cl_name <- ds$cluster_names[cl_id];
          
          # Remove duplicates
          original_ctgs <- unlist(ds$clusters[[cl_id]]);
          updated_ctgs <- original_ctgs[! original_ctgs %in% sel_ctgs];
          
          print_debug_msg(paste("updated_ctgs=", length(updated_ctgs)));
          
          # Remove cluster with no contig assigned
          if(length(updated_ctgs)==0) { 
            removed_clusters[[cl_id]] <- unlist(ds$clusters[[i]]);
            removed_cluster_names <- c(removed_cluster_names, cl_name);
            print_debug_msg(paste(">> Deleting ", cl_name, "[", cl_id, "]", sep=""))
            next;
          } else {
            new_cl_n <- length(cluster_list$clusters);
            #new_cl_id <- paste("c", new_cl_n+1, sep="");
            new_cl_id <- get_unique_cl_id();
            
            print_debug_msg(paste(">> Updating ", cl_name, "[", cl_id, "] to ", cl_name, " [", cl_id, "]", sep=""))
            
            cluster_list$clusters[[new_cl_id]] <- unlist(updated_ctgs);
            cluster_list$cluster_names[[new_cl_id]] <- cl_name;
          }
        }
        
        print_debug_msg(paste("Size of new list of clusters: ", length(cluster_list$clusters), sep=""))
        #update_ds$clusters <- {}; update_ds$clusters <- new_clusters;
        #update_ds$cluster_names <- {}; update_ds$cluster_names <- new_cluster_names;
      }
      
      # Keep the selected contigs to their homes and remove them from the current selection
      if(opt == "N") {
        print_debug_msg(paste("sel_ctgs=", paste(sel_ctgs, collapse=", "), "\n",sep=""));
        print_debug_msg(paste("ovlp_ctgs=", paste(ovlp_ctgs, collapse=", "), "\n",sep=""));
        sel_ctgs <- sel_ctgs[which(! sel_ctgs %in% ovlp_ctgs)]; 
        print_debug_msg(paste("N: sel_ctgs=", length(sel_ctgs), "\n", sep=""))
        
        # If all contigs in sel_ctgs was excluded, we won't add any new cluster
        if(length(sel_ctgs) == 0) {
          display_msg("All contigs are excluded, no cluster is created.");
          return(cluster_list);
        }
      }
    } else {
      # No overlapping contigs
    }
  }
  
  print_debug_msg(paste("Pass: sel_ctgs=", length(sel_ctgs), "\n", sep=""));
  
  
  # Operation
  # opt_mode: 0 - new cluster, 1 - append to a cluster
  opt_mode <- 0;
  
  # Check the name of the new cluster
  cl_names <- unlist(cluster_list$cluster_names);
  new_cluster_name <- ""
  while(T) {        
    display_msg("", mode=-1);
    
    new_cluster_name <- readline("Please enter a name for this new cluster: "); 
    
    if(!check_name(new_cluster_name)) { display_msg(paste(new_cluster_name, " is not valid.", sep="")); next; }
    
    # new_cluster_name is identical to something in the cluster list
    #if(length(get_cl_id(ds, new_cluster_name) > 0))
    if(new_cluster_name %in% cl_names) {
      display_msg(paste("Cluster name, ", new_cluster_name, ", is identical to an existing cluster.", sep=""))
      
      opt <- -1;
      
      # Ask what to do
      while(opt != "Y" && opt != "N") { opt <- toupper(readline("Append to the existing cluster (Y) or enter another name (N)? (Y/N) ")); }
      
      # Append the selected contigs to the cluster specified by the name
      if(opt == "Y") { opt_mode <- 1; break; }
    } else {
      # The new name is acceptable
      break;
    }
  }
  
  # Perform operations
  if(opt_mode == 1) {
    # Append to an existing cluster
    #cl_id <- get_cl_id(update_ds, new_cluster_name);
    cl_id <- names(cluster_list$cluster_names)[which(cluster_list$cluster_names == new_cluster_name)];
    cluster_list$clusters[[cl_id]] <-  c(cluster_list$clusters[[cl_id]], sel_ctgs);
  } else {
    # Create a new id for the new cluster
    cl_n <- length(cluster_list$clusters);
    #cl_id <- paste("c", cl_n+1, sep="");
    cl_id <- get_unique_cl_id();
    
    
    cluster_list$clusters[[cl_id]] <- sel_ctgs;
    cluster_list$cluster_names[[cl_id]] <- new_cluster_name;
  }
  
  #eval.parent(substitute(ds<-update_ds));
  return(cluster_list);
}



get_unique_cl_id <- function(exclude_list=character(0))
{
  return(paste("cl", runif(1,1,90000), runif(1,1,200000), sep=""));
}

# Add contigs to a cluster specified by cl_id
append_to_cluster <- function(ds, cl_id, ctgs)
{
  #ds$clusters[[cl_id]] <- unique(c(ds$clusters[[cl_id]], ctgs));
  return(unique(c(ds$clusters[[cl_id]], ctgs)));
}

# Remove contigs from a cluster specified by cl_id
remove_from_cluster <- function(ds, cl_id, ctgs)
{
  return(unique(ds$clusters[[cl_id]][!ds$clusters[[cl_id]] %in% ctgs]));
}


#
#
# 
clip_pts_from_view <- function(pd_dat, paxy)
{
  x1<-paxy[1];x2<-paxy[2];y1<-paxy[3];y2<-paxy[4];
  print_debug_msg(paste("x1:", x1, " x2:", x2, " y1:", y1, " y2:", y2, sep=""));
  print_debug_msg(paste(colnames(pd_dat)));
  sel_pd_dat <- pd_dat[which(x1 <= pd_dat$x & x2 >= pd_dat$x),];
  print_debug_msg(paste("sel_pd_dat: ", nrow(sel_pd_dat), sep=""));
  sel_pd_dat <- sel_pd_dat[which(y1 <= sel_pd_dat$y & y2 >= sel_pd_dat$y),];
  
  return(sel_pd_dat);
}

# Return the cl_id for a cl_name
get_cl_id <- function(ds, cl_name)
{
  idx <- which(ds$cluster_names %in% cl_name);
  if(length(idx) > 0) { return(names(ds$cluster_names)[idx]); } else { return(character(0)); }
}



check_name <- function(name)
{
  if(nchar(name) == 0) { return(FALSE);}
  return(grepl("^([[:alpha:]]|[.][._[:alpha:]])[._[:alnum:]]*$", name));
}


display_msg <- function(msg, mode=0)
{
  if(mode == 0) { cat(paste(msg, "\n", sep="")); } 
  else if(mode == 1) { spacer_l <- floor((nchar(SEPARATOR) - nchar(msg)) / 2); SPACER_L <- paste(rep(" ", spacer_l), collapse=""); cat(paste(SPACER_L, msg, "\n", sep="")); }
  else if(mode == 2) { cat(paste(SEPARATOR, "\n", sep="")); }
  else if(mode == 6) { cat(paste(SEPARATOR2, "\n", sep="")); }
  else if(mode == 3) { cat(paste(paste(rep(" ", 3), collapse=""), msg, "\n", sep="")); }
  else if(mode == 4) { cat("\n"); }
  else if(mode == 5) { cat(paste("\t", msg, "\n", sep="")); }
  else if(mode == 7) { cat(paste("  version ", VERSION, "\n", sep="")); }
  #else if(mode == -1) { cat(paste(rep("\014\n", 200))); }
  else if(mode == -1) { cat(paste(rep("\n", 100))); }
  else {}
}

# Output debug msg
print_debug_msg <- function(msg) { if(PRINT_DEBUG) { cat(paste(msg, "\n", sep="")); } }

# Import data from input file 
# Format: tsv or csv
import_data <- function(ifn, format="tsv", header=T, comment.char="@", stringsAsFactors=F, convert_to_list=T)
{  
  print_debug_msg(paste("import_data().ifn=", ifn, sep=""));
  if(! file.exists(ifn)) { display_msg(paste("File, ", ifn, ", is not found", sep="")); return;}
  
  if(format == "tsv") { separator <- "\t"; };
  if(format == "csv") { separator <- ","; };
  
  d <- read.table(ifn, sep=separator, row.names=1, header=header, stringsAsFactors=stringsAsFactors, comment.char=comment.char);
  
  # Convert into list if the dataframe only contains one column
  if(ncol(d) == 1 && convert_to_list) { 
    nn <- rownames(d);d <- unlist(d[,1]); 
    #print_debug_msg(paste("Number of rows imported: ", length(d), sep=""));
    names(d) <- nn; 
    #print_debug_msg(paste("Number of rows imported: ", length(d), sep=""));
    display_msg(paste("Number of data imported: ", length(d), " from ", ifn, sep=""));
  } else {
    display_msg(paste("Number of data imported: ", nrow(d), ", ", ncol(d), " from ", ifn, sep=""));
  }
  return(d);
}


# This function exports binned contigs to list files that can be used in
# subsequent analyses
export_clusters_to_file <- function(ds, prefix=character(0), odir=getwd())
{
  if(length(prefix) == 0) { prefix <- ds$sample_id; }
  
  # Iterate all the clusters
  for(i in 1 : length(ds$clusters)) {
    cl_id <- names(ds$clusters)[i];
    cl_name <- ds$cluster_names[cl_id];
    ofn <- paste(odir, paste(prefix, cl_name, "list", sep="."), sep="/");
    display_msg(paste("Exporting ", cl_name, " (", cl_id, ") to ", ofn, sep=""));
    if(file.exists(ofn)) { file.rename(ofn, paste(ofn, ".bak", sep="")); }
    write.table(ds$clusters[[i]], ofn, row.names=F, quote=F, col.name=c("#Contig"));
  }
  #   ofn <- paste(odir, paste(bin_id, "list", sep="."), sep="/");
  #   if(file.exists(ofn)) { file.rename(ofn, paste(ofn, ".bak", sep="")); }
  #   display_msg(paste("Exporting contig ids of ", bin_id," to ", ofn, sep=""));
  #   write.table(rownames(binned_ctgs), ofn, row.names=F, quote=F, col.name=c("#Contig"));
}


# Quick save the current progress to *.cluster file
save_clusters_to_file <- function(ds, sel_ctgs=character(0), ofn=character(0), split=F)
{
  
  if(length(ofn) == 0) { ofn <- paste(ds$cov_src, ".clusters", sep=""); }
  
  display_msg("Exporting...")
  clusters <- {};
  
  # ctg_id, cluster_id, cluster_name, colors
  for(i in 1 : length(ds$clusters)) {
    cl_id <- names(ds$clusters)[i];
    
    cl_name <- ds$cluster_names[cl_id];
    cl_ctgs <- unlist(ds$clusters[cl_id]);
    display_msg(paste("cl_id:", cl_id, ", cl_name:", cl_name, "cl_ctgs: ", length(cl_ctgs), sep="\t"))
    
    new_df <- data.frame(CL_ID=rep(cl_id, length(cl_ctgs)), CL_NAME=rep(cl_name, length(cl_ctgs)), CL_CTG=cl_ctgs);
    if(i == 1) { clusters <- new_df; } else { clusters <- rbind(clusters, new_df, row.names=NULL); }
  } 
  
  # Saving current working selection
  if(length(sel_ctgs) > 0) { 
    cl_id <- TAG_WORKING_CTG;
    cl_name <- TAG_WORKING_CTG;
    new_df <- data.frame(CL_ID=rep(cl_id, length(sel_ctgs)), CL_NAME=rep(cl_name, length(sel_ctgs)), CL_CTG=sel_ctgs);
    clusters <- rbind(clusters, new_df, row.names=NULL);
  }
  
  
  display_msg(paste(length(ds$clusters), " clusters were exported to ", ofn, sep=""))
  
  if(file.exists(ofn)) { file.rename(ofn, paste(ofn, ".bak", sep="")); }
  write.table(clusters, ofn, quote=F, sep="\t", row.names=F, col.names =T);
}


export_mds_to_file <- function(pd, mds_ofn=character(0))
{
  if(length(mds_ofn) == 0) { mds_ofn <- "coords.mds"; }
  
  if(file.exists(ofn)) { file.rename(ofn, paste(ofn, ".bak", sep="")); }
}



# Import clusters from *.cluster file
import_cluster_file <- function(cluster_fn=character(0))
{
  #   if(length(cluster_fn) == 0)
  #   {
  #     cluster_fn <- ""
  #     while(!check_name(cluster_fn)) { 
  #       cluster_fn <- readline("Please enter a name for this new cluster: "); 
  #       if(toupper(new_cluster_name) == "C") { return(NULL); }
  #       if(file.exists(new_cluster))
  #     }
  #   }
  
  df_clusters <- {};
  
  clusters <- read.table(cluster_fn, header=T, stringsAsFactors = F, sep="\t");
  
  nr_cl_ids <- unique(clusters$CL_ID);
  cl_names <- {};
  cl_ctgs <- {};
  for(i in 1 : length(nr_cl_ids)) {
    cl_id <- nr_cl_ids[i];
    if(cl_id == TAG_WORKING_CTG) { next; }
    sel_ctgs <- subset(clusters, clusters$CL_ID == cl_id);
    cl_name <- sel_ctgs$CL_NAME[1];
    #cat(paste("cl_id=", cl_id, ", cl_name=", cl_name, ", sel_ctgs:", nrow(sel_ctgs), "\n", sep=""))
    cl_names[[cl_id]] <- cl_name;
    cl_ctgs[[cl_id]] <- unlist(sel_ctgs$CL_CTG);
  }
  
  # recover working selection
  df_clusters$working_clusters <- {};
  if(TAG_WORKING_CTG %in% cl_names) { df_clusters$working_clusters <- subset(clusters, clusters$CL_ID == TAG_WORKING_CTG); }
  
  df_clusters$clusters <- cl_ctgs;
  df_clusters$cluster_names <- cl_names;
  display_msg(paste(length(df_clusters$clusters), " clusters are imported.", sep=""));
  
  return(df_clusters);
}


# Main controller
main_ctr()

#
# Generate coverage plots
# 
plot_3d_cov <- function()
{
  
}



