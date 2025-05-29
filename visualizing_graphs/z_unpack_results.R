source('unpack_functions.R')


load('/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/spike_data/Brain_Region.RData')


n <- 50   
movements <- c(0, 1, 2)
VR <- 0
epoch_or_week <- 'week'
ew_nums <- c(18, 22, 26)
min_edges <- 1
results_root <- "/u/home/j/jasenzz/Graphical_Cox_Process/GMpp-main/results_alzheimers/"
# IDs <- c('346', '351', '366', '361', '362', '368')  # mouse ID
# ID2 <- c('Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3') # our name 
# 
# IDs <- c('346', '351', '366', '361')  # mouse ID
# ID2 <- c('Tau1', 'Tau2', 'Tau3', 'WT1') # our name 

for(movement in movements){
  for(ew_num in ew_nums){
    unpack(n, movement, VR, epoch_or_week, ew_num, min_edges, df_brain_region, results_root)
    
    print('done')
  }
}


