#
#
# GOAL: analyze run time statistics from the .RData file
#
#

library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr) 

file_name <- 'run_time_results.RData'
folder_name <- '../results/runtime/'
load(paste0(folder_name, file_name)) # called strata_stats


# plot 1: active neurons vs runtime

t_cols <- paste0('t', 1:8)

strata_stats_long <- strata_stats %>% 
  pivot_longer(cols = all_of(t_cols),         
               names_to = "timepoint",
               values_to = "value") %>% 
  mutate(Movement = as.factor(Movement), VR = as.factor(VR)) %>% 
  mutate(MoveVR = interaction(Movement, VR, sep = "_")) %>% 
  mutate(Active_Neurons = factor(Active_Neurons, levels = c('5', '20', '50', '100'))) %>% 
  mutate(value = as.numeric(value)) %>% 
  mutate(x_group = interaction(Active_Neurons, MoveVR, sep = "_"))


g1 <- ggplot(strata_stats_long,
             aes(x = x_group, y = value, fill = timepoint)) +
  geom_bar(stat = "identity") +  # stacked by default
  labs(x = "Active Neurons + MoveVR", y = "Runtime", fill = "Timepoint") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

