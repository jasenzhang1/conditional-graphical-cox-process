# mouse 346
# weeks 17, 18, 19

source('../../utils/0_Dependencies.R')
data_dir <- '/u/home/j/jasenzz/Alzheimers/Data/'

source('create_spiketrain_dataset_funcs.R')

weekly_dataset <- load_realigned_dataset(data_dir)[[1]]
neuron_df <- load_brain_region_dataset(data_dir)[[1]]

ID <- '346' #mouse_ID
age <- '17' #age
break_ties <- F # do we want to break ties?
n <- 600    #number of replicates

data_df2 <- get_spiketrain_dataset_multiple_weeks('346', c('17', '18', '19'), 600, weekly_dataset, neuron_df)
