import problem as pr
import generate_graphs as gg
import importlib
import random
import osmnx as ox
import networkx as nx
import helper_functions as hf
import numpy as np
importlib.reload(gg)



seed = random.randint(0,1000)
seed = 204
random.seed(seed)
print("seed: ", seed)


# graph, demands = gg.generate_random_graph(n_max=20,t=2,k_num_max=2,cap_max=2,density_param=1)
# graph, demands = gg.generate_random_graph(n_max=10,t=2,k_num_max=2,cap_max=2,density_param=1,c_max=10)
# graph, demands = gg.generate_random_graph(n_max=1000,t=2,k_num_max=2,cap_max=2,density_param=1,c_max=10)



#################################################################################3

# graph = gg.place_to_nx("Castenaso")
graph = gg.place_to_nx("helsinki")


# demands = [(1,3,2),(13,6,3)]
# demands = [(1,0,2)]

demands = [(6,20,1),(6,100,3),(6,105,2)]

print(graph.number_of_nodes())
# hf.plot_multigraph(graph,with_labels=True, font_size=10)
# hf.plot_multigraph(graph,with_labels=True, font_size=10)
# graph_smooth, demands = gg.smooth_out_simple_nodes(graph,demands,min_c=200)
# print("demands_smooted: ", demands)
# print(graph_smooth.number_of_nodes())
# hf.plot_multigraph(graph_smooth,with_labels=True, font_size=10)

# gen = pr.Problem.Evolution.dfs_edges_random(graph,0,depth_limit=6)

# for e in gen:
#     print(e)

p = pr.Problem(graph,demands)
res = pr.Problem.Dai_solver.solve(p,MAX_ITER = 1,MAX_ITER_LR = 50)
hf.plot_solution_graph(graph,res.X,with_labels=True,font_size=10,figure_size=(20,20))
print(res)

# p = pr.Problem(graph,demands)
res = pr.Problem.Descent.solve(p,max_num_iter = 10)
print(res)
if res.X is not None:
    hf.plot_solution_graph(graph,res.X,with_labels=True,font_size=10,figure_size=(20,20))

