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
random.seed(seed)
print("seed: ", seed)


# graph, demands = gg.generate_random_graph(n_max=10,t=2,k_num_max=2,cap_max=2,density_param=1)
# p = pr.Problem(graph,demands)
# res = pr.Problem.Cplex_intlinprog.solve(p)
# print(res)


#################################################################################3

graph = gg.place_to_nx("Castenaso")


demands = [(1,3,2),(13,6,3)]


p = pr.Problem(graph,demands)
res = pr.Problem.Dai_solver.solve(p,MAX_ITER = 1,MAX_ITER_LR = 50)
hf.plot_solution_graph(graph,res.X,with_labels=True,font_size=5,figure_size=(20,20))
print(res)
