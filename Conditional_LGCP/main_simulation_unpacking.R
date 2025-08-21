# unpack
library(grid)
library(gridExtra)



load('simu_results_banded_v1/n_100.RData')

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`4`$g_112)


grid.newpage()
pdf('temp2.pdf', width = 20, height = 14)
grid.draw(graph_results_i$estimated_graphs_part_2$g_113)
dev.off()