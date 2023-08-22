import time
from stopit import threading_timeoutable as timeoutable
from solvers import Biobjective_LP, Descent, Optimal_cvxpy, Dai3_solver, Keep_feasible, Keep_feasible_shuffle, Biobjective, DnC, LP_relaxation
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



    
    

graphs_with_demands = []

##########################################33



# graph, demands = gg.get_basic_example(UPPER_COST=1)
# graphs_with_demands.append((graph,[demands]))
# graph, demands = gg.get_basic_example(UPPER_COST=10)
# graphs_with_demands.append((graph,[demands]))
# graph, demands = gg.get_basic_example(UPPER_COST=100)
# graphs_with_demands.append((graph,[demands]))
# graph, demands = gg.get_basic_example(UPPER_COST=1000)
# graphs_with_demands.append((graph,[demands]))

# hf.plot_multigraph(graph, with_labels=True,font_size=12,figure_size=(20,20))

#######################################3

# seed = random.randint(0,1000)
# seed = 204
# seed = 207
# random.seed(seed)
# print("seed: ", seed)

# graph, demands = gg.generate_random_graph(n_max=10,t=2,k_num_max=2,cap_max=3,density_param=1,c_max=4)
# # demands_infeasible = gg.random_demands(graph,10)
# graphs_with_demands.append((graph,[gg.random_demands(graph,1)]))

seed = random.randint(0, 1000)
seed = 204
# seed = 205
random.seed(seed)
print("seed: ", seed)


# for graph, demands in gg.generate_random_graphs(5,n_max=100,t=3,k_num_max=2,cap_max=3,density_param=1,c_max=4):
#     graphs_with_demands.append((graph,[demands]))

# for graph, demands in gg.generate_random_graphs(4,n_max=10,t=3,k_num_max=1,cap_max=3,density_param=1,c_max=4):
#     graphs_with_demands.append((graph,[demands]))

# n_max=50, t = 15, k_num_max = 1, cap_max = 3 

# for n_max in [50,1000]:
#         for t in [5,20,100]:
#             if t < n_max:
#                 for k_num_max in [1,5]:
#                     for cap_max in [3,10]:
#                         for graph, demands in gg.generate_random_graphs(5,n_max,t,k_num_max,cap_max):
#                             graphs_with_demands.append((graph,[demands]))

# proti primer

# with open('LP_proti_primer.pkl', 'rb') as f:
#     graph, demands = pickle.load(f)
#     graphs_with_demands.append((graph,[demands]))

# NUM_SAMPLES = 10

# for graph, demands in gg.generate_random_graphs(2):
#     graphs_with_demands.append((graph,[demands]))
# ################################################################
# for graph, demands in gg.generate_random_graphs(NUM_SAMPLES):
#     # graphs_with_demands.append((graph,[demands]))
#     pass
    
# for graph, demands in gg.generate_random_graphs(NUM_SAMPLES,n_max=100):
#     # graphs_with_demands.append((graph,[demands]))
#     pass
    
# for graph, demands in gg.generate_random_graphs(NUM_SAMPLES,t=5, k_num_max=6):
#     # graphs_with_demands.append((graph,[demands]))
#     pass
    
# for graph, demands in gg.generate_random_graphs(NUM_SAMPLES,cap_max=10):
#     # graphs_with_demands.append((graph,[demands]))
#     pass

# for i, (graph, demands) in enumerate(gg.generate_random_graphs(NUM_SAMPLES,c_max=100)):
#     graphs_with_demands.append((graph,[demands]))
#     # if i == 1:
#     #     graphs_with_demands.append((graph,[demands]))
#     #     # with open('LP_proti_primer.pkl', 'wb') as f:
#     #     #     pickle.dump((graph,demands), f)
#     # else:
#     #     pass
#     pass

#########################################################################
# NUM_SAMPLES = 1
# for graph, demands in gg.generate_random_graphs(NUM_SAMPLES,t=5, k_num_max=6):
#     graphs_with_demands.append((graph,[demands]))
#     pass
#########################################################################

# # for n_max in [50,1000]:
# for n_max in [1000]:
#         # for t in [20,100]:
#         for t in [20]:
#             if t < n_max:
#                 for k_num_max in [5]:
#                     for cap_max in [3]:
#                         for graph, demands in gg.generate_random_graphs(1,n_max,t,k_num_max,cap_max):
#                             graphs_with_demands.append((graph,[demands]))

# hf.plot_multigraph(graph, with_labels=True,font_size=10,figure_size=(10,10))


# graph, demands = gg.generate_random_graph(n_max=20,t=2,k_num_max=2,cap_max=2,density_param=1,c_max=10)
# graphs_with_demands.append((graph,[demands,gg.random_demands(graph,1)]))
# graph, demands = gg.generate_random_graph(n_max=10,t=2,k_num_max=10,cap_max=20,density_param=1,c_max=10)
# graphs_with_demands.append((graph,[demands]))
# graph, demands = gg.generate_random_graph(n_max=1000,t=2,k_num_max=10,cap_max=2,density_param=1,c_max=10)
# graphs_with_demands.append((graph,[demands]))
# graph, demands = gg.generate_random_graph(n_max=200,t=2,k_num_max=10,cap_max=2,density_param=1,c_max=10) # keep feasible takoj najde, ker je malo vrstnih redov
# graphs_with_demands.append((graph,[demands]))

# graph, demands = gg.generate_random_graph(n_max=200,t=10,k_num_max=10,cap_max=2,density_param=1,c_max=10)
# graphs_with_demands.append((graph,[demands]))

# hf.plot_multigraph(graph, with_labels=True,font_size=10,figure_size=(20,20))

# 3

# graph = gg.place_to_nx("Castenaso",save=True)

# list_of_demands = []
# # demands = [(1,3,2),(13,6,3)]
# # list_of_demands.append(demands)
# # demands = [(1,0,1)]
# # list_of_demands.append(demands)
# # demands = [(1,0,2)]
# # list_of_demands.append(demands)
# # demands = [(13,15,2),(7,3,2)]
# # list_of_demands.append(demands)

# # demands = gg.random_demands(graph,7)
# # list_of_demands.append(demands)

# # demands = gg.random_demands(graph, 10)
# # list_of_demands.append(demands)

# demands = gg.random_demands(graph,50)
# list_of_demands.append(demands)

# # demands = gg.random_demands(graph,100)
# # list_of_demands.append(demands)


# graphs_with_demands.append((graph, list_of_demands))




#############HELSINKI######################
graph = gg.place_to_nx("Kamppi,Helsinki",save=True)

list_of_demands = []
# demands = [(6,20,1),(6,100,3),(6,105,2)]
# list_of_demands.append(demands)
# demands = [(6,20,1),(60,100,3),(6,105,2)]
# list_of_demands.append(demands)
# demands = [(6,20,15),(60,100,3),(6,105,2),(13,55,10)]
# list_of_demands.append(demands)
# demands = gg.random_demands(graph,10,repete_prob=0.3) # primer ko LP ne dela
# list_of_demands.append(demands)
# #demands = gg.random_demands(graph,20,repete_prob=0.3) # primer ko biobj ne dela
# #list_of_demands.append(demands)
demands = gg.random_demands(graph,5,repete_prob=0.8)
list_of_demands.append(demands)
graphs_with_demands.append((graph,list_of_demands))





# # # # graph = gg.place_to_nx("lj")
# graph = gg.place_to_nx("Ljubljana",save = True, mode = "point", point = (46.05407,14.52114))

# # # hf.plot_multigraph(graph, with_labels=False,font_size=5,figure_size=(30,30))

# list_of_demands = []

# # demands = [(0,1,1)]
# # list_of_demands.append(demands)

# # brezovica = gg.point_to_closest_node((46.03333, 14.4),graph)
# # medvode = gg.point_to_closest_node((46.1422, 14.41114),graph)

# # demands = [(brezovica,medvode,6)]
# # list_of_demands.append(demands)


# # demands = [(10,21,3),(20,55,20)]
# # list_of_demands.append(demands)


# demands = gg.random_demands(graph,20)
# list_of_demands.append(demands)

# # demands = gg.random_demands(graph,100)
# # list_of_demands.append(demands)


# graphs_with_demands.append((graph,list_of_demands))

# hf.plot_multigraph(graph, with_labels=True,font_size=5,figure_size=(20,20))



# graph = gg.place_to_nx("SanMarino")
# hf.plot_multigraph(graph, with_labels=True,font_size=7,figure_size=(20,20))

# 33


# 33

# hf.plot_multigraph(graph,with_labels=True, font_size=10)
# hf.plot_multigraph(graph,with_labels=True, font_size=10)
# graph_smooth, demands = gg.smooth_out_simple_nodes(graph,demands,min_c=200)
# print("demands_smooted: ", demands)
# print(graph_smooth.number_of_nodes())
# hf.plot_multigraph(graph_smooth,with_labels=True, font_size=10)

# gen = Evolution.dfs_edges_random(graph,0,depth_limit=6)

# for e in gen:
#     print(e)


TIMEOUT = 30
# TIMEOUT = 2*60

# seed = 1
# random.seed(seed)


Problem.verbose = True
p = Problem(graph, demands)
solvers = [
    ("Optimal_cvxpy", lambda p, _: Optimal_cvxpy.solve(
        p, timeout=0.5*TIMEOUT, SOLVER="CBC")),

    # ("LP_relaxation", lambda p, _: LP_relaxation.solve(p, SOLVER="GLOP")),
    # ("Dai3_solver", lambda p, _: Dai3_solver.solve(p, MAX_ITER=1, MAX_ITER_LR=50000, timeout=0.5*TIMEOUT)),

    # ("Keep_feasible_shuffle", lambda p, known_LB: Keep_feasible_shuffle.solve(
    #     p, max_cons_repetitions=3, timeout=0.5*TIMEOUT, known_LB=known_LB)),

    ("Descent", lambda p, known_LB: Descent.solve(p, max_num_iter=1000, max_num_paths_generated=100,
                                                  worsening_step_prob=0.1, max_depth=None, timeout=0.5*TIMEOUT,
                                                  known_LB=known_LB)),


    # ("Biobjective", lambda p, known_LB: Biobjective.solve(
    #     p, SOLVER="CBC", known_LB=known_LB)),
    # ("Biobjective_LP", lambda p, known_LB: Biobjective_LP.solve(
    #     p, SOLVER="GLOP", known_LB=known_LB)),
    # ("DnC", lambda p, known_LB: DnC.solve(p, known_LB=known_LB))
]

    # # ("dai2",lambda p,_: Dai2_solver.solve(p,MAX_ITER=50,MAX_ITER_LR=5)),
    # #    ("Keep_feasible",lambda p,known_LB: Keep_feasible.solve(p,timeout=0.5*TIMEOUT,known_LB=known_LB)),

ang_to_slo ={"Biobjective":"dvokrit.", "Optimal_cvxpy":"RiO", "LP_relaxation":"LP", "Dai3_solver":"LR","DnC":"'DiV'", "Descent":"sest.",
             "Keep_feasible_shuffle":"požrešna","Biobjective_LP":"dvokrit. LP"}

for graph, demands_list in graphs_with_demands:
    for di, demands in enumerate(demands_list):
        min_cost = None
        p = Problem(graph, demands)
        for solver_name, solver in solvers:
            res = None
            known_LB = None

            l = [p.results[key] for key in list(p.results) if "Dai3" in key or "LP" in key or "Opt" in key]
            if len(l) > 0:
                try:
                    LB = float(l[0].LB)
                    known_LB =  LB if known_LB is None else max(known_LB,LB) 
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
            if res is None:
                Problem.print("res is None")
            else:
                if res.cost is not None and res.status in [hf.Status.FEASIBLE, hf.Status.OPTIMAL]:
                    min_cost = res.cost if min_cost is None else min(res.cost,min_cost)

            if False and res is not None and res.X is not None and graph.number_of_nodes() < 201:
                small = graph.number_of_nodes() < 100
                s = np.array(p.cap - res.X.sum(axis=1).T).flatten()
                over_cap_edges_ids = list(np.where(s < 0)[0])
                over_cap_edges = [list(graph.edges())[ei] for ei in over_cap_edges_ids]
                print("over_cap_edges: ",[list(graph.edges())[ei] for ei in over_cap_edges_ids])
                hf.plot_solution_graph(graph, res.X, with_labels=small, font_size=7 if small else 0, figure_size=(
                    20, 20), with_arrows=small, over_cap_edges = over_cap_edges)
            if min_cost is None:
                l = [p.results[solver].cost for solver in p.results if ("Biobj" in solver or "DnC" in solver) and p.results[solver].cost is not None]
                min_cost = min(l) if len(l) > 0 else None
            
            if min_cost is not None:
                for solver in p.results:
                    if p.results[solver].cost is not None and p.results[solver].status in [hf.Status.FEASIBLE, hf.Status.OPTIMAL]:
                        p.results[solver].percent_of_min = p.results[solver].cost / min_cost

        demands_list[di] = [demands_list[di], p]

print("results: ")
for graph, demands_list in graphs_with_demands:
    print("-> ", end="")
    gg.about_graph(graph)
    for di, demands_data in enumerate(demands_list):
        demands, p = demands_data
        print("-> ", end="")
        gg.about_demands(demands)
        p.print_results()
        
# title_to_result = {graph.graph["graph_title"] : p.result for graph,(_ ,(_,p)) in graphs_with_demands}
titles_ps = []
for graph, demands_data in graphs_with_demands:
    for demands, p in demands_data:
        titles_ps.append((graph.graph["graph_title"], p))
# titles_results = [(graph.graph["graph_title"], demands_data) for demands_data in  [demands for graph,demands in graphs_with_demands]]
titles = {graph.graph["graph_title"] for graph,_ in graphs_with_demands}

if len(titles) < len(titles_ps) and any("_" in graph.graph["graph_title"] for graph,_ in graphs_with_demands):

    groups_info_with_results = []
    for title in titles:
        times = {}
        costs = {}
        LBs = {}
        over_caps = {}
        percents_of_min = {}
        
        status_count = {}
        
        nums_nodes = []
        densities = []
        for results, graph in [(p.results, p.graph) for title_,p in titles_ps if title_ == title]:
            for solver in results:
                result = results[solver]
                def update(li, val):
                    if val is None or np.isnan(val) or np.isinf(val):
                        return li[solver] if solver in li else []
                    
                    if solver not in li:
                        return [val]
                    else:
                        return li[solver] + [val]
                        
                # update = lambda li, val: [val] if solver not in li else li[solver] + [val]
                times[solver] = update(times,result.time)
                costs[solver] = update(costs,result.cost)
                LBs[solver] = update(LBs,result.LB)
                over_caps[solver] = update(over_caps,result.over_cap_count)
                percents_of_min[solver] = update(percents_of_min,result.percent_of_min)
                
                if solver not in status_count:
                    status_count[solver] = {hf.Status.BLANK:0,hf.Status.OPTIMAL:0,hf.Status.FEASIBLE:0,hf.Status.INFEASIBLE:0,hf.Status.OVER_CAP:0}
                status_count[solver][result.status] += 1
                
            nums_nodes.append(graph.number_of_nodes())
            densities.append(nx.density(graph))
                
        results = {}
        for solver in times:
            # results[solver] = Problem.Result_grouped(status_count = status_count[solver],
            #                                          mean_cost = np.mean([x for x in costs[solver] if x is not None and x != np.inf]),
            #                                          std_cost = np.std([x for x in costs[solver] if x is not None and x != np.inf]),
            #                                          mean_time = np.mean([x for x in times[solver] if x is not None]),
            #                                          std_time = np.std([x for x in times[solver] if x is not None]),
            #                                          mean_LB = np.mean([x for x in LBs[solver] if x is not None]))
            results[solver] = Problem.Result_grouped(status_count = status_count[solver],
                                                     mean_cost = round(np.mean(costs[solver]),3) if len(costs[solver]) > 0 else "-",
                                                     std_cost = round(np.std(costs[solver]),3) if len(costs[solver]) > 0 else "-",
                                                     mean_time = round(np.mean(times[solver]),3) if len(times[solver]) > 0 else "-",
                                                     std_time = round(np.std(times[solver]),3) if len(times[solver]) > 0 else "-",
                                                     mean_LB = round(np.mean(LBs[solver]),3) if len(LBs[solver]) > 0 else "-",
                                                     mean_over_cap = round(np.mean(over_caps[solver]),3) if len(over_caps[solver]) > 0 else "-",
                                                     mean_percent_of_min = round(np.mean(percents_of_min[solver]),3) if len(percents_of_min[solver]) > 0 else "-",
                                                     std_percent_of_min = round(np.std(percents_of_min[solver]),3) if len(percents_of_min[solver]) > 0 else "-")
        about_group = {"title":title, "mean_num_nodes":np.mean(nums_nodes), "mean_density":np.mean(densities)}
        groups_info_with_results.append((about_group, results))
            
    print("##############################")
    print(r"\begin{table}[]")
    print( r"\begin{tabular}{|c"+"|c"*(len(solvers)+1)+"|}")
    print(r"\hline")
    print("omrežja",end="")
    print("& polja", end="")
    print("".join(["& "+ ang_to_slo[solver_name] for solver_name,_ in solvers]))
    print(r"\\ \hline")

    for about_group, results in groups_info_with_results:
        # graph_cell = gg.graph_latex_cell(graph
        n_max, t, k_num_max, cap_max, _, c_max = about_group["title"].split("_")
        pairs = [("max($n$)", n_max), ("št. skupin",t), ("maks. vel. sk.", k_num_max),("maks. kapac.", cap_max),("maks. cena", c_max),
                 (r"$\overline{n}$", about_group["mean_num_nodes"]),(r"$\overline{\rho}$", round(about_group["mean_density"],3))]
        group_cell = hf.latex_cell_pairs(pairs)
        # ...
        print(group_cell)
        
        
        # print(r"&\begin{tabular}{l} $\overline{\text{cena}}$ \\ std(cena)      \\ $\overline{\text{rač. čas}}$ \\ std(rač. čas.) \\ -,opt,dop,nedop,preko  \\ $\overline{\text{sp. meja}}$ \\ $\overline{\text{preko}}$ \end{tabular}")
        print(r"&\begin{tabular}{l} $\overline{\text{cena}}$   \\   $\overline{\text{točnost}}$ \\ std(točnost)\\ $\overline{\text{rač. čas}}$ \\ std(rač. čas.) \\ statusi  \\ $\overline{\text{sp. meja}}$ \\ $\overline{\text{preko}}$ \end{tabular}")
        for solver_name, _ in solvers:
            l = [results[key] for key in results.keys() if solver_name in key]
            if len(l) == 0:
                print("&")
            else:
                res = l[0]
                result_cell = res.latex_cell()
                
                print("&",result_cell)
        print(r"\\ \hline")
    print("\end{tabular}")   
    print("\end{table}")
  
        
else:       
        
        
        
            
            
            
            
            
    print("##############################")
    print(r"\begin{table}[]")
    print( r"\begin{tabular}{|c"+"|c"*(len(solvers)+2)+"|}")
    print(r"\hline")
    print("omrežja",end="")
    print("& zahteve", end="")
    print("& polja", end="")
    print("".join(["& "+ ang_to_slo[solver_name] for solver_name,_ in solvers]))
    print(r"\\ \hline")

    for graph, demands_list in graphs_with_demands:
        graph_cell = gg.graph_latex_cell(graph)
        
        for di, demands_data in enumerate(demands_list):
            demands, p = demands_data
            demands_cell = gg.demands_latex_cell(demands)
            print(graph_cell)
            graph_cell = ""
            print("&",demands_cell)
            print(r"&\begin{tabular}{l}cena      \\točnost \\rač. čas  \\status   \\sp. meja  \\preko    \end{tabular}")
            for solver_name, _ in solvers:
                l = [p.results[key] for key in p.results.keys() if solver_name in key]
                if len(l) == 0:
                    print("&")
                else:
                    res = l[0]
                    result_cell = res.latex_cell()
                    
                    print("&",result_cell)
            print(r"\\ \hline")
    print("\end{tabular} ")
    print("\end{table}")
    
pass
