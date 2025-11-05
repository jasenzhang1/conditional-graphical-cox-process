# helper functions to work with the estimated dataset


# helper function - return possible sample sizes from files

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

# helper function - load the file that ends with the number n
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

# helper function
old_to_new_graph_results_i <- function(graph_results_i){
  
  # old graph_results_i had 
  # - estimated_graphs_part_1
  # - estimated_graphs_part_2
  #
  #
  
  steps <- unique(unlist(lapply(graph_results_i$estimated_graphs_part_2, names)))
  reorganized <- setNames(lapply(steps, function(step) {
    sapply(graph_results_i$estimated_graphs_part_2, `[[`, step, simplify = FALSE)
  }), steps)
  
  # Return step_0, step_1 (shared) + reorganized per-subject steps

  return(c(graph_results_i$estimated_graphs_part_1, reorganized))  
  
}


load_all_results <- function(folder_name){
  
  # in a folder with multiple fitted datasets, load them all into a list
  
  files <- list.files(folder_name, full.names = TRUE)
  
  ns <- suppressWarnings({
    as.numeric(sub(".*_(.*)\\.RData$", "\\1", files)) # retrieve numbers
  })
  
  data_files <- files[!is.na(ns)]
  
  ns <- ns[!is.na(ns)]
  
  
  results_list <- lapply(data_files, function(f) {
    e <- new.env()          # create an isolated environment
    load(f, envir = e)      # load into that environment
    as.list(e)              # convert to a list (in case multiple objects)
  })
  
  names(results_list) <- ns    
  
  return(results_list)
  
}

# helper function - load a single RData file 
load_file <- function(file_name){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: with the base folder as `Conditional_LGCP`, locate an RData file and load it
  #
  # inputs:
  #
  # - file_name     (string)    'simu_results/banded_trig2/OG/temp.RData'
  # 
  #
  # outputs:
  #
  # - NULL, just loads that file
  #
  # ----------------------------------------------------------------------------
  
  # Create a new environment
  env <- new.env()
  
  # Load the file into the environment
  load(file_name, envir = env)
  
  # Return the first object in the environment as a list
  as.list(env)[[1]]
  
}
