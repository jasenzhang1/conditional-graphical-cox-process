# unpack
library(grid)
library(gridExtra)


# 1) graphs that don't change over time 
# 2) graphs that change over time, but values below adjacency threshold are 0
# 3) ...

# 10) continuous covariates = (week, velocity, VR)


load('simu_results_banded_v2/n_1000.RData')


grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`5`$g_112)


grid.newpage()
pdf('temp2.pdf', width = 20, height = 14)
grid.draw(graph_results_i$estimated_graphs_part_2$`4`$g_84)
dev.off()