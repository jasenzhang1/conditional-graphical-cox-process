
# On 7/17, we decided to construct precision matrices based on:

# - E_{y_c, y_d} ground truths. A range of y_c values can belong to a single E_{y_c, y_d} graph
# - calculate h_ij to act as weights
# - previously, if h_ij is sufficiently low, this would correspond to a E_{y_c, y_d} graph that doesn't have a (i, j) edge


# borrowed functions

# - generate_covariance_matrix
# - 

find_yc_group <- function(y_c, y_c_borders){
  
  
  # ----------------------------------------------------------------------------
  # 
  # GOAL: we partition the euclidean space of y_c into regions.
  #
  # 
  # input:
  # 
  # - y_c         (q_c-dim vector) 
  # - y_c_borders (list of q_c vectors, each vector specifies bordering values)
  #
  #
  # output:
  #
  # - region (q_c vector), each value specifies region ID 
  #
  #
  # example:
  #
  # y_c = c(1, 5)
  # y_c_borders = list(c(0.5, 0.7), c(1, 5, 10))
  #
  # region = c(3, 2)
  #
  # 5 is in the 2nd region because of <= rule
  # ----------------------------------------------------------------------------
  
  region <- c()
  
  for(i in 1:length(y_c)){
    value <- y_c[i]
    borders <- y_c_borders[[i]]
  
    # Ensure borders are sorted
    borders <- sort(borders)
    
    # Define breaks including -Inf and Inf
    breaks <- c(-Inf, borders, Inf)
    
    # Find interval: which region `value` falls into
    region_i <- findInterval(value, vec = breaks, left.open = TRUE)
    
    region <- c(region, region_i)
  }
  
  return(region)
  
}

find_E_yc_yd_graph <- function(region, graphs){
  
  
  # graphs  (list) entries are border values separated by underscores
  #
  # - each item is a p x p adj_mat with 0's on diag
  #
  # region  (vector) vector of values that indicate region
  #
  # - i.e. c(1, 1)
  #
  #
  
  region_id <- paste(region, collapse = "_")
  
  return(graphs[[region_id]])
}

construct_precision_operator_E_yc_yd <- function(y_c, y_c_borders, graphs){

  
  # 
  # 
  # GOAL: we take in a list of graphs and construct a precision operator from the appropriate ground truth
  # 
  # y_c  (q_c dim vector)
  # graphs (list)   list of ground truth graphs E_{yc, yd}
  #
  # 
  
  # 1) for a single y_c vector, find its associated ground truth graph
  region_i <- find_yc_group(y_c, y_c_borders)
  graph_i <- find_E_yc_yd_graph(region_i, graphs)
  
  
  # 2) calculate the precision matrix from (graph_i = p x p adj_mat)
  
  
  
  
}



