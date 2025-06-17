# introduce parallelization tools 

start_time <- Sys.time()
ffs = 1
set.seed(ffs)
source("GraphPP_FUN.R")
source("modified_funcs.R")
source("main_funcs.R")



#### simulate the data ####

# 1) Open 99_Create_Spiketrain_Dataset and load the data_df file

# 1a) Parameters that don't change

IDs <- c('346', '351', '366', '361', '362', '368')  # mouse ID
ID2 <- c('Tau1', 'Tau2', 'Tau3', 'WT1', 'WT2', 'WT3') # our name 

num_neurons <- c(169, 250, 240, 249, 235, 294)
names(num_neurons) <- ID2

# 1b) Parameters that change, double check!!! ---------------------------------

task_num <- '05'

movements <- c(0, 1, 2)
VR <- 0

ews <- c(17:29, 31, 33, 35, 38) # 1:4
ew_symb <- 'w'

tn_symb <- 't'
tn_vec <- c(5) # tns <- 10  # replicates or timescale

min_edges <- 1
ncores <- parallel::detectCores() - 1

print('Number of Cores')
print(ncores)

# -----------------------------------------------------------------------------

for(tns in tn_vec){
    main_algorithm(ews, ew_symb, movements, task_num, tn_symb, tns, min_edgges, ncores)
}


print('-----------------------------------------------------------')
end_time <- Sys.time()
total_time <- as.numeric(end_time - start_time, units = "mins") %>% round(2)
print('COMPLETE!!')
print(paste('total minutes taken: ', total_time, sep = ''))


