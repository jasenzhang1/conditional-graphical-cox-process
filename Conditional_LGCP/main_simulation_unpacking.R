# unpack
library(grid)
library(gridExtra)


# 1) graphs that don't change over time 
# 2) graphs that change over time, but values below adjacency threshold are 0
# 3) ...

# 10) continuous covariates = (week, velocity, VR)

# simu_results = constant




load('simu_results_banded_c1_3/n_300.RData')
grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_1$g_01) # ground truths

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_1$g_22) # rho_i

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2[[1]]$g_112) #HS norms

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2[[1]]$g_113) #ROC curves

grid.newpage()
visualize_V_cond_convergence('simu_results_banded_c0_4', 'V', 2, 3) %>% print()

# how does V_cond error change over time?



# values are -10, -7.5, -5.2, -1.95, 4.2
load('simu_results_banded_v2_3/n_100.RData')
grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_1$g_22)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`4.2`$g_112)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`4.2`$g_113)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-1.95`$g_112)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-1.95`$g_113)



grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-5.2`$g_112)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-5.2`$g_113)


grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-7.5`$g_112)

grid.newpage()
grid.draw(graph_results_i$estimated_graphs_part_2$`-7.5`$g_113)



grid.newpage()
pdf('temp2.pdf', width = 20, height = 14)
grid.draw(graph_results_i$estimated_graphs_part_2$`4`$g_84)
dev.off()