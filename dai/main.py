from problem import Problem
import generate_graphs as gg
import importlib
import random
import osmnx as ox
import networkx as nx
import helper_functions as hf
import numpy as np
importlib.reload(gg)
from solvers import Evolution, Descent, Optimal_cvxpy, Dai_solver, Dai2_solver, Keep_feasible


seed = random.randint(0,1000)
seed = 204
# seed = 205
random.seed(seed)
print("seed: ", seed)


# graph, demands = gg.generate_random_graph(n_max=20,t=2,k_num_max=2,cap_max=2,density_param=1)
# # graph, demands = gg.generate_random_graph(n_max=20,t=2,k_num_max=2,cap_max=2,density_param=1,c_max=10)
graph, demands = gg.generate_random_graph(n_max=10,t=2,k_num_max=10,cap_max=20,density_param=1,c_max=10)
# graph, demands = gg.generate_random_graph(n_max=1000,t=2,k_num_max=10,cap_max=2,density_param=1,c_max=10)

# hf.plot_multigraph(graph, with_labels=True,font_size=10,figure_size=(20,20))

#################################################################################3

# graph = gg.place_to_nx("Castenaso")
# graph = gg.place_to_nx("helsinki")


# demands = [(1,3,2),(13,6,3)]
# demands = [(1,0,2)]
# demands = [(13,15,2),(7,3,2)]

# demands = [(6,20,1),(6,100,3),(6,105,2)]
# demands = [(6,20,1),(60,100,3),(6,105,2)]

#############################################################################33

print(graph.number_of_nodes())
# hf.plot_multigraph(graph,with_labels=True, font_size=10)
# hf.plot_multigraph(graph,with_labels=True, font_size=10)
# graph_smooth, demands = gg.smooth_out_simple_nodes(graph,demands,min_c=200)
# print("demands_smooted: ", demands)
# print(graph_smooth.number_of_nodes())
# hf.plot_multigraph(graph_smooth,with_labels=True, font_size=10)

# gen = Evolution.dfs_edges_random(graph,0,depth_limit=6)

# for e in gen:
#     print(e)

import multiprocessing
import time



def all():

    p = Problem(graph,demands)


    # res = Optimal_cvxpy.solve(p)
    # if res.X is not None:
    #     hf.plot_solution_graph(graph,res.X,with_labels=True,font_size=10,figure_size=(20,20))
    # print(res)

    # res = Optimal_cvxpy.solve(p,SOLVER="CBC")
    # if res.X is not None:
    #     hf.plot_solution_graph(graph,res.X,with_labels=True,font_size=10,figure_size=(20,20))
    # print(res)

    # res = Optimal_cvxpy.solve(p,SOLVER="GUROBI")
    # if res.X is not None:
    #     hf.plot_solution_graph(graph,res.X,with_labels=True,font_size=10,figure_size=(20,20))
    # print(res)


    res = Dai_solver.solve(p,MAX_ITER = 1,MAX_ITER_LR = 50)
    Problem.verbose = False
    if res.X is not None:
        hf.plot_solution_graph(graph,res.X,with_labels=True,font_size=10,figure_size=(20,20))
    print(res)
    known_LB = int(res.message.split(".")[0])

    # Problem.verbose = True
    # res = Dai2_solver.solve(p,MAX_ITER = 10,MAX_ITER_LR = 50)
    # if res.X is not None:
    #     hf.plot_solution_graph(graph,res.X,with_labels=True,font_size=10,figure_size=(20,20))
    # print(res)

    Problem.verbose = True
    res = Descent.solve(p,max_num_iter = 100,known_LB=known_LB)
    print(res)
    if res.X is not None:
        hf.plot_solution_graph(graph,res.X,with_labels=True,font_size=10,figure_size=(20,20))
        
    
    # Problem.verbose = True
    # res = Keep_feasible.solve(p)
    # if res.X is not None:
    #     hf.plot_solution_graph(graph,res.X,with_labels=True,font_size=10,figure_size=(20,20))
    # print(res)

    print(p.results)
    
p = multiprocessing.Process(target=all, name="All")
p.start()

# Wait 10 seconds for foo
time.sleep(300)

print("terminate")
# Terminate foo
p.terminate()

# Cleanup
p.join()

