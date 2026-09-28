# helper functions to work with the estimated dataset


# helper function - return possible sample sizes from files

get_ns <- function(folder_name){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: within a folder, there are several datasets that end in a number then .RData,
  #       get all possible n's
  #
  #       - example: 'CPGM_n_500.RData' --> 500
  #
  #
  # input:
  #
  # - folder_name      (string)               'simu_results/banded_trig2/OG'
  #
  #
  # output:
  #
  # - ns               (vector)               sample sizes
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

get_ns_with_rep_unsorted <- function(folder_name) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: extract 'n' from filenames like 'adj_type_n_500_rep_1.RData'
  #
  #       pattern logic:
  #         - look for "_n_"
  #         - capture digits (\\d+)
  #         - stop at "_rep_"
  #
  #
  # input:
  #
  # - folder_name      (string)               'simu_results/hub_block_v2/CPGM'
  #
  #
  # output:
  #
  # - ns               (vector)               unique sample sizes, in file order
  #
  # ----------------------------------------------------------------------------
  
  # 1) Get all file names in the directory
  files <- list.files(folder_name, full.names = FALSE)
  
  # 2) Use regmatches and regexec for more precise middle-string extraction
  # This regex looks for the digits immediately following "_n_"
  ns <- sapply(files, function(x) {
    match <- regmatches(x, regexec("_n_(\\d+)_rep_", x))
    return(match[[1]][2]) # Index 2 captures the group inside the parentheses
  })
  
  # 3) Clean up: Convert to numeric, remove NAs, and get unique values
  ns <- as.numeric(ns)
  ns <- unique(ns[!is.na(ns)])
  
  return(ns)
}


get_ns_with_rep <- function(folder_name) {
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: extract 'n' from filenames like 'adj_type_n_500_rep_1.RData'
  #
  #       pattern logic:
  #         - look for "_n_"
  #         - capture digits (\\d+)
  #         - stop at "_rep_"
  #
  #
  # input:
  #
  # - folder_name      (string)               'simu_results/hub_block_v2/CPGM'
  #
  #
  # output:
  #
  # - ns               (vector)               unique sample sizes, sorted
  #
  # ----------------------------------------------------------------------------
  
  # 1) Get all file names in the directory
  files <- list.files(folder_name, full.names = FALSE)
  
  # 2) Use regmatches and regexec for more precise middle-string extraction
  # This regex looks for the digits immediately following "_n_"
  ns <- sapply(files, function(x) {
    match <- regmatches(x, regexec("_n_(\\d+)_rep_", x))
    return(match[[1]][2]) # Index 2 captures the group inside the parentheses
  })
  
  # 3) Clean up: Convert to numeric, remove NAs, and get unique values
  ns <- as.numeric(ns)
  ns <- unique(ns[!is.na(ns)])
  
  return(sort(ns))
}

# helper function - load the file that ends with the number n
get_file_name_and_load <- function(folder_name, n){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: within a folder, there are several datasets that end in a number then .RData,
  #       load the one that ends with n
  #
  #
  # input:
  #
  # - folder_name      (string)               'simu_results/banded_trig2/OG'
  # - n                (integer)              e.g. 100
  #
  #
  # output:
  #
  # - object           (any)                  the first object stored in that file
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
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: get the queried y_c values of the dataset with sample size n
  #
  #
  # input:
  #
  # - folder_name      (string)               'simu_results/banded_trig2/OG'
  # - n                (integer)              e.g. 100
  #
  #
  # output:
  #
  # - y_c_query        (vector of strings)    names of estimated_graphs_part_2
  #
  # ----------------------------------------------------------------------------
  
  graph_results_i <- get_file_name_and_load(folder_name, n)
  
  return(names(graph_results_i$estimated_graphs_part_2))
  
}

get_method <- function(folder_name){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: get the estimation method name, assuming it lies at the beginning of the .RData file
  #
  #       - e.g. OG_n_1000.RData --> 'OG'
  #
  #
  # input:
  #
  # - folder_name      (string)               'simu_results/banded_trig2/OG'
  #
  #
  # output:
  #
  # - method           (string)               method name
  #
  # ----------------------------------------------------------------------------
  
  files <- list.files(folder_name, full.names = TRUE)
  
  # which files end in .RData?
  idx <- grep("\\.RData$", files)[1]
  
  method <- sub("_.*", "", sub(".*/", "", files[idx]))  # after last /, before next _
  
  return(method)
}

# helper function
old_to_new_graph_results_i <- function(graph_results_i){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: reorganize an old graph_results_i, which had estimated_graphs_part_1 and estimated_graphs_part_2,
  #       into the new step-first format
  #
  #
  # input:
  #
  # - graph_results_i  (list)                 old results with estimated_graphs_part_1 and estimated_graphs_part_2
  #
  #
  # output:
  #
  # - graph_results_i  (list)                 step_0, step_1 (shared) + reorganized steps --> y_c queries
  #
  # ----------------------------------------------------------------------------
  
  steps <- unique(unlist(lapply(graph_results_i$estimated_graphs_part_2, names)))
  reorganized <- setNames(lapply(steps, function(step) {
    sapply(graph_results_i$estimated_graphs_part_2, `[[`, step, simplify = FALSE)
  }), steps)
  
  # Return step_0, step_1 (shared) + reorganized per-subject steps

  return(c(graph_results_i$estimated_graphs_part_1, reorganized))  
  
}


load_all_results <- function(folder_name){
  
  # ----------------------------------------------------------------------------
  #
  # GOAL: in a folder with multiple fitted datasets, load them all into a list
  #
  #
  # input:
  #
  # - folder_name      (string)               'simu_results/banded_trig2/OG'
  #
  #
  # output:
  #
  # - results_list     (list)                 one list of loaded objects per file, named by n
  #
  # ----------------------------------------------------------------------------
  
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
  #
  # input:
  #
  # - file_name        (string)               'simu_results/banded_trig2/OG/temp.RData'
  #
  #
  # output:
  #
  # - object           (any)                  the first object stored in that file
  #
  # ----------------------------------------------------------------------------
  
  # Create a new environment
  env <- new.env()
  
  # Load the file into the environment
  load(file_name, envir = env)
  
  # Return the first object in the environment as a list
  as.list(env)[[1]]
  
}
