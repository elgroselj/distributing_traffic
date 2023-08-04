import pickle
from matplotlib import pyplot as plt
from problem import Problem
import generate_graphs as gg
import importlib
import random
import osmnx as ox
import networkx as nx
import helper_functions as hf
import numpy as np
importlib.reload(gg)
from solvers import Descent, Optimal_cvxpy, Dai_solver, Dai2_solver, Dai3_solver, Keep_feasible, Keep_feasible_shuffle
from stopit import threading_timeoutable as timeoutable

graphs_with_demands = []

# seed = random.randint(0,1000)
# seed = 204
# seed = 207
# random.seed(seed)
# print("seed: ", seed)

# graph, demands = gg.generate_random_graph(n_max=10,t=2,k_num_max=2,cap_max=3,density_param=1,c_max=4)
# # demands_infeasible = gg.random_demands(graph,10)
# graphs_with_demands.append((graph,[gg.random_demands(graph,1)]))

seed = random.randint(0,1000)
seed = 204
# seed = 205
random.seed(seed)
print("seed: ", seed)




graph, demands = gg.generate_random_graph(n_max=100,t=2,k_num_max=2,cap_max=3,density_param=1,c_max=4)
graphs_with_demands.append((graph,[demands]))

hf.plot_multigraph(graph, with_labels=True,font_size=10,figure_size=(10,10))



# graph, demands = gg.generate_random_graph(n_max=20,t=2,k_num_max=2,cap_max=2,density_param=1,c_max=10)
# # graphs_with_demands.append((graph,[demands,gg.random_demands(graph,1)]))
# graph, demands = gg.generate_random_graph(n_max=10,t=2,k_num_max=10,cap_max=20,density_param=1,c_max=10)
# graphs_with_demands.append((graph,[demands]))
# graph, demands = gg.generate_random_graph(n_max=1000,t=2,k_num_max=10,cap_max=2,density_param=1,c_max=10)
# # graphs_with_demands.append((graph,[demands]))
# graph, demands = gg.generate_random_graph(n_max=200,t=2,k_num_max=10,cap_max=2,density_param=1,c_max=10) # keep feasible takoj najde, ker je malo vrstnih redov
# graphs_with_demands.append((graph,[demands]))

# graph, demands = gg.generate_random_graph(n_max=200,t=10,k_num_max=10,cap_max=2,density_param=1,c_max=10)
# graphs_with_demands.append((graph,[demands]))

# hf.plot_multigraph(graph, with_labels=True,font_size=10,figure_size=(20,20))

#################################################################################3

graph = gg.place_to_nx("Castenaso")

list_of_demands = []
# demands = [(1,3,2),(13,6,3)]
# list_of_demands.append(demands)
# demands = [(1,0,2)]
# list_of_demands.append(demands)
# demands = [(13,15,2),(7,3,2)]
# list_of_demands.append(demands)

demands = gg.random_demands(graph,7)
list_of_demands.append(demands)

# demands = gg.random_demands(graph,20)
# list_of_demands.append(demands)

# demands = gg.random_demands(graph,50)
# list_of_demands.append(demands)

# demands = gg.random_demands(graph,100)
# list_of_demands.append(demands)


graphs_with_demands.append((graph,list_of_demands))




# graph = gg.place_to_nx("helsinki")

# list_of_demands = []
# demands = [(6,20,1),(6,100,3),(6,105,2)]
# list_of_demands.append(demands)
# # demands = [(6,20,1),(60,100,3),(6,105,2)] # infeasible
# # list_of_demands.append(demands)
# graphs_with_demands.append((graph,list_of_demands))








# # # graph = gg.place_to_nx("lj")
# graph = gg.place_to_nx("Ljubljana",save = True, mode = "point", point = (46.05407,14.52114))

# # # hf.plot_multigraph(graph, with_labels=False,font_size=5,figure_size=(30,30))

# list_of_demands = []

# # demands = [(0,1,1)]
# # list_of_demands.append(demands)

# brezovica = gg.point_to_closest_node((46.03333, 14.4),graph)
# medvode = gg.point_to_closest_node((46.1422, 14.41114),graph)

# demands = [(brezovica,medvode,5)]
# list_of_demands.append(demands)

# # # demands = [(10,21,3),(20,55,20)]
# # # list_of_demands.append(demands)


# graphs_with_demands.append((graph,list_of_demands))

# hf.plot_multigraph(graph, with_labels=True,font_size=5,figure_size=(20,20))

#############################################################################33

        
    
##########################################################################33

# hf.plot_multigraph(graph,with_labels=True, font_size=10)
# hf.plot_multigraph(graph,with_labels=True, font_size=10)
# graph_smooth, demands = gg.smooth_out_simple_nodes(graph,demands,min_c=200)
# print("demands_smooted: ", demands)
# print(graph_smooth.number_of_nodes())
# hf.plot_multigraph(graph_smooth,with_labels=True, font_size=10)

# gen = Evolution.dfs_edges_random(graph,0,depth_limit=6)

# for e in gen:
#     print(e)


import time
TIMEOUT = 100


Problem.verbose = True
p = Problem(graph,demands)
solvers = [
        #    lambda p,_: Optimal_cvxpy.solve(p,timeout=0.5*TIMEOUT),
           lambda p,_: Dai3_solver.solve(p,MAX_ITER=1,MAX_ITER_LR=50),
        
        #    lambda p, known_LB: Descent.solve(p, max_num_iter=1000, max_num_paths_generated=100,
                                    #   worsening_step_prob=0.1, max_depth=40, timeout=0.5*TIMEOUT,
                                    #   known_LB=known_LB),
        #    lambda p,known_LB: Keep_feasible.solve(p,timeout=0.5*TIMEOUT,known_LB=known_LB),
        #    lambda p,known_LB: Keep_feasible_shuffle.solve(p,max_cons_repetitions=100,timeout=0.5*TIMEOUT,known_LB=known_LB)
        ]

for graph, demands_list in graphs_with_demands:
    for di, demands in enumerate(demands_list):
        p = Problem(graph,demands)
        for solver in solvers:
            res = None
            known_LB = None
            
            l = [p.results[key] for key in list(p.results) if "Dai3" in key]
            if len(l) > 0:
                try:
                    known_LB = float(l[0].message.split(".")[0].strip("["))
                except:
                    Problem.print("no known LB")
                    pass
            
            @timeoutable()
            def f(solver, p):
                before_t = time.perf_counter()
                res = solver(p, known_LB)
                after_t = time.perf_counter()
                if res is not None:
                    res.time = after_t - before_t
                return res
            res = f(solver=solver, p=p, timeout=TIMEOUT)
            # hf.plot_solution_graph(graph,res.X,with_labels=False,font_size=10,figure_size=(20,20))
            
            if True and res is not None and res.X is not None and graph.number_of_nodes() < 20100:
                small = graph.number_of_nodes()<200
                hf.plot_solution_graph(graph,res.X,with_labels=small,font_size=10 if small else 0,figure_size=(20,20),with_arrows=small)
            
                
            
            
        demands_list[di] = [demands_list[di], p]

print("results: ")
for graph, demands_list in graphs_with_demands:
    print("-> ",end="")
    gg.about_graph(graph)
    for di, demands_data in enumerate(demands_list):
        demands, p = demands_data
        print("-> ",end="")
        gg.about_demands(demands)
        p.print_results()



