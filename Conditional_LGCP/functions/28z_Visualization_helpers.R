# helper functions to work with the estimated dataset


# helper function

get_ns <- function(folder_name){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: within a folder, there are several datasets that end in a number then .RData
  #       get all possible n's 
  #
  # - example: 'CPGM_n_500.RData' --> 500
  #
  #
  #
  # inputs:
  #
  # - folder_name     (string)    'simu_results/banded_trig2/OG'
  # 
  #
  # outputs:
  #
  # - ns              (vector)
  #
  # ----------------------------------------------------------------------------  
  
  # 1) find the file in the folder and load it 
  files <- list.files(folder_name, full.names = TRUE)
  
  ns <- suppressWarnings({
    as.numeric(sub(".*_(.*)\\.RData$", "\\1", files))   # retrieve number before .RData and after recent underscore
  })
  
  ns <- ns[!is.na(ns)]
  
  return(ns)
}

# helper function
get_file_name_and_load <- function(folder_name, n){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: within a folder, there are several datasets that end in a number then .RData
  #       load that file
  #
  # inputs:
  #
  # - folder_name     (string)    'simu_results/banded_trig2/OG'
  # - n               (integer)   100
  # 
  #
  # outputs:
  #
  # - NULL, just loads that file
  #
  # ----------------------------------------------------------------------------
  
  # 1) find the file in the folder and load it 
  files <- list.files(folder_name, full.names = TRUE)
  
  ns <- suppressWarnings({
    as.numeric(sub(".*_(.*)\\.RData$", "\\1", files))   # retrieve number before .RData and after recent underscore
  })
  
  
  if(n %in% ns){
    idx <- which(ns == n)
  } else{
    stop('ERROR on 28z_get_file_name_and_load: n not available')
  }  
  
  env <- new.env()
  
  load(files[idx], envir = env) 
  
  return(as.list(env)[[1]])
  
}

# helper function
get_y_c_query <- function(folder_name, n){
  
  graph_results_i <- get_file_name_and_load(folder_name, n)
  
  return(names(graph_results_i$estimated_graphs_part_2))
  
}

get_method <- function(folder_name){
  
  
  # assume that method names lie at the beginning of the .RData file
  # OG_n_1000.RData
  # method <- 'OG'
  
  files <- list.files(folder_name, full.names = TRUE)
  
  # which files end in .RData?
  idx <- grep("\\.RData$", files)[1]
  
  method <- sub("_.*", "", sub(".*/", "", files[idx]))  # after last /, before next _
  
  return(method)
}
