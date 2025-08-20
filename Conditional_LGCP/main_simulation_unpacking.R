# unpack
library(grid)
library(gridExtra)



load('simu_results/n_10000.RData')

grid.newpage()
#grid.draw(graph_results_i$estimated_graphs_part_1$g_71)

pdf('temp2.pdf', width = 20, height = 14)
grid.draw(graph_results_i$estimated_graphs_part_2$g_113)
dev.off()