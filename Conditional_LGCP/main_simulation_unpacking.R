# unpack
library(grid)
library(gridExtra)

setwd('simu_results')


load('n_100.RData')

grid.newpage()
#grid.draw(graph_results_i$estimated_graphs_part_1$g_71)
grid.draw(graph_results_100$estimated_graphs_part_1$g_21)