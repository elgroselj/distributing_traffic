import itertools
import multiprocessing
from problem import Problem

from collections import Counter
import sys
from matplotlib import pyplot as plt
import networkx as nx
import numpy as np
import random
import time
import copy

# from docplex.mp.model import Model
import cvxpy as cp

from scipy.optimize import minimize
from scipy import sparse
from scipy.optimize import linprog

import helper_functions as hf

import dai
import dai2
import dai3

class Descent:
        
    class Node:
        def init_class(p):
            Descent.Node.p = p
            Descent.Node.label = 0
            
            
        def __init__(self, paths, parent = None):
            self.paths = paths
            
            
            self.over_cap = None
            self.cost = None
            
            self.over_cap_count = None
            
            
            self.children = []
            self.parent = parent
            
            self.hash = None
            
            self.label = Descent.Node.label
            Descent.Node.label += 1
            
            self.level = 0 if self.parent is None else self.parent.level +1
            
        def calc_over_cap_count(mat,cap):
            edge_flow = mat.sum(axis=1).T
            diff = edge_flow-cap
            over_cap_count = diff[diff > 0].sum()
            return over_cap_count
        
        def evaluate(self):
            graph = Descent.Node.p.graph
            c = Descent.Node.p.c
            cap = Descent.Node.p.cap
            
        
            if self.cost is not None and self.over_cap is not None:
                return (self.over_cap,self.cost)
            
            r = sum([len(k) for k in self.paths]) # stevilo poti za določit (št posameznih avtov)
            t = len(self.paths) # stevilo skupin avtov
            m = len(graph.edges())
            self.Q = sparse.dok_matrix((m,r)) # Q je razprta različica X
            self.X = sparse.dok_matrix((m,t)) # Q je razprta različica X
            col = 0
            col_x = 0
            for k in self.paths:
                for path in k:
                    epath = hf.nodes_to_edges_path(path)
                    for edge in epath:
                        row = list(graph.edges()).index(edge)
                        self.Q[row,col] = 1
                        self.X[row,col_x] += 1
                        
                    col += 1
                col_x +=1
            self.Q = self.Q.tolil()
            self.X = self.X.tolil()
            
            
            edge_flow = self.Q.sum(axis=1).T
            
            _,over_cap_edges,_ = sparse.find(cap < edge_flow)
            mat = self.Q[over_cap_edges,:] > 0
            self.over_cap = {e : list(sparse.find(mat[ei,:]>0)[1]) for ei, e in enumerate(over_cap_edges)}
            
            self.cost = int((c @ edge_flow.T)[0,0])
            
            self.over_cap_count = Descent.Node.calc_over_cap_count(self.Q,cap)
            
            return self.over_cap, self.cost
        
                    
            
        def get_child(self,evaluated,max_num_paths_generated=1000, worsening_step_prob=0.001):
            graph = Descent.Node.p.graph
            demands = Descent.Node.p.demands
            # c = Descent.Node.p.c
            # cap = Descent.Node.p.cap
            
            alt_c = {e : graph.edges()[e]["c"] for e in graph.edges()}
                            
            
            
            over_cap_edges = random.sample(self.over_cap.keys(),len(self.over_cap)) # TODO
            
            for over_cap_edge in over_cap_edges:
                over_cap_paths_ids = self.over_cap[over_cap_edge]
                
                over_cap_paths_ids = random.sample(over_cap_paths_ids,len(over_cap_paths_ids)) # TODO
                
                
                e = list(graph.edges())[over_cap_edge]
                old_e_c = alt_c[e]
                alt_c[e] = np.inf
                nx.set_edge_attributes(graph,alt_c,"alt_c")
                
                for over_cap_path_id in over_cap_paths_ids:
                    id = 0
                    for ki, (_,_,num_k) in enumerate(demands):
                        id+=num_k
                        if id > over_cap_path_id:
                            pi = id - over_cap_path_id -1
                            break
                    # ki ... skupina od poti
                    # pi ... idx poti znotraj skupine
                    O, D, _ = demands[ki]
                    paths_generator = nx.shortest_simple_paths(graph,O,D,"alt_c")
                    for path_counter, path in enumerate(paths_generator):
                        if path_counter > max_num_paths_generated:
                            Problem.print(path_counter)
                            break
                        
                        paths = list(self.paths)
                        paths[ki] = list(paths[ki])
                        paths[ki][pi] = path
                        paths[ki] = tuple(paths[ki])
                        ch = Descent.Node(paths,parent=self)
                        ch.evaluate()
                        if ch in evaluated or (not random.random() < worsening_step_prob and ch.over_cap_count >= self.over_cap_count):
                            Descent.Node.label -= 1
                            continue
                        self.children.append(ch)
                        Problem.print("success")
                        yield ch
                        
                        
                Problem.print("**over_cap_paths exhausted")
                alt_c[e] = old_e_c
            Problem.print("***over_cap_edges exhausted")
            yield None
        
                                
                            
                            
            
        def __eq__(self, other):
            for ki in range(len(self.paths)):
                for path in self.paths[ki]:
                    if path not in other.paths[ki]:
                        return False
                for path in other.paths[ki]:
                    if path not in self.paths[ki]:
                        return False
            return True
            # return self.paths == other.paths
        
        def __hash__(self) -> int:
            if self.hash is not None: return self.hash
            self.hash = sum(self.paths[0][0])
            return self.hash
        
        def print_tree(self, level=0):
            ret = "\t"*level+repr(self)+"\n"
            if self.children is not None:
                for child in self.children:
                    ret += child.print_tree(level+1)
            return ret
        
        def __repr__(self) -> str:
            # return "({}):{} {} = {}".format(self.level, self.label, self.paths, self.over_cap)
            return "({}):{} = {} over_cap_count: {}".format(self.level, self.label, self.cost, self.over_cap_count)
        


    def lower_the_cost(node,feasible,allowed_over_cap_count=0):
        c = Descent.Node.p.c
        cap = Descent.Node.p.cap
        graph = Descent.Node.p.graph
        demands = Descent.Node.p.demands
        
        node_old = node

        while True:
            Problem.print("iter_of_lowering")
            glob_path_ind = -1
            mini = (node.cost, node.paths)
            for ki, k in enumerate(node.paths):
                for pi, path in enumerate(k):
                    mat = node.Q.copy()
                    glob_path_ind += 1
                    
                    O,D,_ = demands[ki]
                    paths_generator = nx.shortest_simple_paths(graph,O,D)
                    while True:
                        try:
                            path_lower = next(paths_generator)
                        except:
                            return node_old
                        
                        if path_lower == path:
                            break
                        
                        path_vec_lower = hf.edges_to_binary_vector(hf.nodes_to_edges_path(path_lower),list(graph.edges()))
                        mat[:,glob_path_ind] = path_vec_lower
                        
                        # edge_flow = np.sum(mat, axis=1).T
                        edge_flow = mat.sum(axis=1).T
                        diff = edge_flow-cap
                        over_cap_count = diff[diff > 0].sum()
                        if over_cap_count <= allowed_over_cap_count:
                            cost = (c @ edge_flow.T)[0,0]
                            if cost < mini[0]:
                                paths = list(node.paths)
                                paths[ki] = list(paths[ki])
                                paths[ki][pi] = path_lower
                                paths[ki] = tuple(paths[ki])
                                mini = (cost, paths)
                            break
            ch = Descent.Node(mini[1], parent=node)
            ch.evaluate()
            if ch.cost == node.cost:
                return node
            if ch in feasible:
                return ch
            node = ch     
                    
        
        


    def solve(p, max_num_iter = 10, timeout = None, known_LB = None, max_num_paths_generated=1000, worsening_step_prob=0.001,max_depth=10):
        before_t = time.perf_counter()
        
        Descent.Node.init_class(p)
        
        paths = [[nx.shortest_path(p.graph, O, D,"c") ] * num_k for O,D,num_k in p.demands]
        node0 = Descent.Node(paths)
        
        node0.evaluate()
        LB = max(known_LB, node0.cost) if known_LB is not None else node0.cost
        
        node = node0
        
        
        evaluated = set()
        feasible = set()
        
        st = 0
        while st < max_num_iter:
            after_t = time.perf_counter()
            if timeout is not None and after_t - before_t > timeout:
                break
            
            
            
            
            st += 1
            if node is None:
                Problem.print("node was None")
                break
            
            if node not in evaluated:
                Problem.print(("len_feasible:",len(feasible)))
                over_cap, _ = node.evaluate()
                evaluated.add(node)
                
                if len(over_cap) == 0: # dopustna rešitev
                    # izboljšaj kar se da TODO
                    # improved_node = node
                    improved_node = Descent.lower_the_cost(node, feasible)
                    if improved_node != node:
                        Problem.print("node improved")
                    if improved_node not in feasible:
                        feasible.add(improved_node)
                    if improved_node.cost <= LB:
                        print("optimal")
                        break
                    node = node.parent
                    continue
            # hf.plot_solution_graph(p.graph,node.X,with_labels=True,font_size=10,figure_size=(20,20))    ch = next(node.get_child(evaluated))
            if node.level > max_depth:
                Problem.print("max_depth")
                node = node.parent
                # node = node0
                continue
            
            Problem.print(st-1)
            ch = next(node.get_child(evaluated,max_num_paths_generated, worsening_step_prob))
            Problem.print(st-1)
            Problem.print((node,ch))
            if ch is None:
                Problem.print("no more ch")
                node = node.parent
                # node = node0
                continue
                
            node = ch
                
        Problem.print("len_evaluated:"+ str(len(evaluated)))               
        Problem.print("len_feasible:"+ str(len(feasible)))
        if Problem.verbose:
            tree = node0.print_tree()
            print(tree)
            
        message = ""
        over_cap_count = None
        if len(feasible) == 0:
            best_evaluated = min(evaluated, key= lambda node: (node.over_cap_count,node.cost))
            best_evaluated_improved = Descent.lower_the_cost(best_evaluated,feasible,allowed_over_cap_count=best_evaluated.over_cap_count)
            status = "no_feas"
            X = best_evaluated_improved.X
            cost = best_evaluated_improved.cost
            # message = "over_cap_count: {}".format(best_evaluated_improved.over_cap_count)
            over_cap_count = best_evaluated_improved.over_cap_count
        else:
            best_feasible_node = min(feasible, key= lambda node: node.cost)
            if best_feasible_node.cost == LB:
                status = "optimal"
            else:
                status = "suboptimal"
                # message = "at most {} percent of optimum ({})".format(best_feasible_node.cost/LB * 100,LB)
            X = best_feasible_node.X
            cost = best_feasible_node.cost
            over_cap_count = 0
        
        res = Problem.Result(status, X, cost, message, LB, over_cap_count)
        p.results[str(__class__)] = res
        return res
    


# class Dai_solver:
    
#     def solve(p,MAX_ITER,MAX_ITER_LR):
#         if not Problem.verbose:
#             save_stdout = sys.stdout
#             sys.stdout = open('trash', 'w')
#         ###############################
#         status, X, cost, message = dai.dai_solve(p.graph,p.demands,MAX_ITER,MAX_ITER_LR)
#         res = Problem.Result(status, X, cost, message)
#         p.results[str(__class__)] = res
        
#         #################################
#         if not Problem.verbose:
#             sys.stdout = save_stdout
        
#         return res
    
# class Dai2_solver:
#     def solve(p,MAX_ITER,MAX_ITER_LR):
#         if not Problem.verbose:
#             save_stdout = sys.stdout
#             sys.stdout = open('trash', 'w')
#         ###############################
#         status, X, cost, message = dai2.dai2_solve(p.graph,p.demands,MAX_ITER,MAX_ITER_LR)
#         res = Problem.Result(status, X, cost, message)
#         p.results[str(__class__)] = res
        
#         #################################
#         if not Problem.verbose:
#             sys.stdout = save_stdout
        
#         return res
    
class Dai3_solver:
    def solve(p,MAX_ITER,MAX_ITER_LR):
        if not Problem.verbose:
            save_stdout = sys.stdout
            sys.stdout = open('trash', 'w')
        ###############################
        status, X, cost, message, LB, over_cap_count = dai3.dai3_solve(p.graph,p.demands,MAX_ITER,MAX_ITER_LR)
        LB = LB if isinstance(LB, float) else int(LB[0,0])
        res = Problem.Result(status, X, cost, message,LB,over_cap_count)
        p.results[str(__class__)] = res
        
        #################################
        if not Problem.verbose:
            sys.stdout = save_stdout
        
        return res

class Keep_feasible:
      
    def permutation_generator(characters, counts, p):
        length = sum(counts)
        
        def generate_helper(pos, remaining):
            if pos == length:
                X_dict = {}
                code = [ki for ki,_ in remaining]
                Problem.print(code)
                for k in characters:
                    paths = [path for ki, path in remaining if ki==k]
                    for path in paths:
                        epath = hf.nodes_to_edges_path(path)
                        etoei = lambda e : list(p.graph.edges()).index(e)
                        path_dict = {(etoei(e),k):1 for e in epath}
                        X_dict = Counter(X_dict) + Counter(path_dict)
                t = len(p.demands)
                m = p.graph.number_of_edges()
                X = sparse.dok_matrix((m,t))
                for a,k in X_dict:
                    X[a,k] = X_dict[(a,k)]
                X = X.tolil()
                yield X
                
                
                return
            
            for idx, char in enumerate(characters):
                if counts[idx] > 0:
                    O,D,_ = p.demands[char]
                    try:
                        path = nx.dijkstra_path(p.graph,O,D,"alt_c")
                    except:
                        try:
                            code = [ki for ki,_ in remaining[:(pos)]]
                        except:
                            pass
                        Problem.print(("deadend, depth: ", pos, code))
                        return
                    
                    old_alt_c = []
                    epath = hf.nodes_to_edges_path(path)
                    for e in epath:
                        old_alt_c.append(p.graph.edges()[e]["alt_c"])
                        p.graph.edges()[e]["remaining_cap"] -= 1
                        if p.graph.edges()[e]["remaining_cap"] == 0:
                            p.graph.edges()[e]["alt_c"] = None
                    
                    remaining[pos] = (char,path)
                    counts[idx] -= 1
                    yield from generate_helper(pos + 1, remaining)
                    counts[idx] += 1
                    
                    for i,e in enumerate(epath):
                        p.graph.edges()[e]["alt_c"] = old_alt_c[i]
                        p.graph.edges()[e]["remaining_cap"] += 1
                        
        yield from generate_helper(0, [''] * length)
    
    def finish(p,feasible,LB,message=""):
        if len(feasible) == 0:
            status = "no_feas"
            X = None
            cost = None
            over_cap_count = None
        else:
            status = "optimal" if message == "optimal" else "feasible"
            cost, X = min(feasible, key= lambda x: x[0])
            cost = int(cost)
        
            over_cap_count = 0
    
        res = Problem.Result(status, X, cost,message,LB,over_cap_count)
        p.results[str(__class__)] = res
        return res
    
    def solve(p,timeout=None,known_LB=None):
        before_t = time.perf_counter()
        feasible = []
        
        t = len(p.demands)
        
        ks_org = []
        for ki, (_,_,num_k) in enumerate(p.demands):
            for i in range(num_k):
                ks_org.append(ki)
                

        remaining_cap = {e: p.graph.edges()[e]["cap"] for e in p.graph.edges()}
        
        counts = [num_k for _,_,num_k in p.demands]
        
        remaining_cap = {e: p.graph.edges()[e]["cap"] for e in p.graph.edges()}
        nx.set_edge_attributes(p.graph, remaining_cap, "remaining_cap")
        
        alt_c = {e: p.graph.edges()[e]["c"] for e in p.graph.edges()}
        nx.set_edge_attributes(p.graph, alt_c, "alt_c")
        
        for X in Keep_feasible.permutation_generator(range(t),counts,p):
            # print(X)
            
            after_t = time.perf_counter()
            Problem.print(after_t - before_t)
            if timeout is not None and after_t - before_t > timeout:
                Problem.print("timeouted")
                return Keep_feasible.finish(p,feasible,message="timeouted")
            
            if X is None:
                continue
            edge_flow = X.sum(axis=1).T
            cost = (p.c @ edge_flow.T)[0,0]
            feasible.append((cost,X))
            if known_LB is not None and cost <= known_LB:
                Problem.print("optimal")
                return Keep_feasible_shuffle.finish(p,feasible,known_LB,message="optimal")
            # Problem.print("*")
        Problem.print(feasible)
        return Keep_feasible.finish(p,feasible,known_LB,message="all_permutations")
    

class Keep_feasible_shuffle:
    
    def finish(p,feasible,LB,message=""):
        if len(feasible) == 0:
            status = message
            X = None
            cost = None
            over_cap_count = None
        else:
            status = "optimal" if message == "optimal" else "feasible"
            cost, X = min(feasible, key= lambda x: x[0])
            cost = int(cost)
        
            over_cap_count = 0
    
        res = Problem.Result(status, X, cost,"",LB,over_cap_count)
        p.results[str(__class__)] = res
        return res
    
    def solve(p,timeout=None,max_cons_repetitions = 100, known_LB = None):
        before_t = time.perf_counter()
        feasible = []
        
        t = len(p.demands)
        
        ks_org = []
        for ki, (_,_,num_k) in enumerate(p.demands):
            for i in range(num_k):
                ks_org.append(ki)
                

        remaining_cap = {e: p.graph.edges()[e]["cap"] for e in p.graph.edges()}
        
        # counts = [num_k for _,_,num_k in p.demands]
        
        # for ks in Keep_feasible.permutation_generator(range(t),counts,p):
        ks = ks_org
        seen = set()
        cons_repetitions = 0
        while True:
            random.shuffle(ks)
            if tuple(ks) in seen:
                cons_repetitions += 1
                if cons_repetitions >= max_cons_repetitions:
                    Problem.print("repetitions")
                    return Keep_feasible_shuffle.finish(p,feasible,known_LB,message="repetitions")
                Problem.print("seen")
                continue
            else:
                cons_repetitions = 0
                seen.add(tuple(ks))
            
            after_t = time.perf_counter()
            Problem.print(after_t - before_t)
            if timeout is not None and after_t - before_t > timeout:
                Problem.print("timeouted")
                return Keep_feasible_shuffle.finish(p,feasible,known_LB,message="timeouted")
            
            Problem.print(ks)
            
            X = hf.assemble(ks,p.graph,p.demands)
            if X is None:
                continue
            edge_flow = X.sum(axis=1).T
            cost = (p.c @ edge_flow.T)[0,0]
            feasible.append((cost,X))
            
            if known_LB is not None and cost <= known_LB:
                Problem.print("optimal")
                return Keep_feasible_shuffle.finish(p,feasible,known_LB,message="optimal")
            
class Optimal_cvxpy:
    def solve(p,timeout,SOLVER="GLPK_MI"):
        # if p.graph.number_of_nodes() > 30:
        #     return
        
        if not Problem.verbose:
            save_stdout = sys.stdout
            sys.stdout = open('trash', 'w')
        ###############################
        _, constraints, _, constraints_ex_additional, vp, obj_opt, _,_,_  = dai.init_from_graph(p.graph,p.demands)
        
        
        #################################
        if not Problem.verbose:
            sys.stdout = save_stdout
            
        
        # _, constraints, _, constraints_ex_additional, vp, obj_opt  = dai.init_from_graph(p.graph,p.demands)
        problem = cp.Problem(obj_opt, constraints + constraints_ex_additional)
        #####################################33
        if SOLVER == "CBC":
            problem.solve(verbose=True,solver=SOLVER, maximumSeconds=timeout)
        else:
            problem.solve(verbose=True,solver=SOLVER)
        ###############################333
        # def f(SOLVER):
        #     problem.solve(verbose=True,solver=SOLVER)
            
        # pr = multiprocessing.Process(target=f, name="F", args=(SOLVER,))
        # pr.start()

        # # Wait 10 seconds for foo
        # for i in range(int(timeout)):
        #     time.sleep(1)
        #     print("iter")
        #     if not pr.is_alive():
        #         print("nono")
        #         break
                
        # print("juhu")
        # # Terminate foo
        # pr.terminate()

        # # Cleanup
        # pr.join()
        ###########################################333
        status, X, cost = (problem.status, vp["X"].value, problem.value)
        try:
            cost = int(cost)
        except:
            pass
        res = Problem.Result(status, X, cost)
        p.results[str(__class__)] = res
        return res

class LP_relaxation:
    def solve(p,SOLVER="CBC"):
        
        if not Problem.verbose:
            save_stdout = sys.stdout
            sys.stdout = open('trash', 'w')
        ###############################
        _, _, _, _, vp, _, _, obj_LP, constraints_LP  = dai.init_from_graph(p.graph,p.demands)
        
        
        #################################
        if not Problem.verbose:
            sys.stdout = save_stdout
            
        
        # _, constraints, _, constraints_ex_additional, vp, obj_opt  = dai.init_from_graph(p.graph,p.demands)
        problem = cp.Problem(obj_LP, constraints_LP)
        #####################################33
        problem.solve(verbose=True,solver=SOLVER)       
        ###########################################333
        status_LP, X_real, val = (problem.status, vp["X_real"].value, problem.value)
        message = ""
        if status_LP == "infeasible":
            status = "infeasible"
            X = None
            cost = None
            LB = None
            over_cap_count = None
        else:
            X = np.round(X_real)
            LB = int(np.round(val))
            # message = "{} is LB".format(LB)
            s = vp["cap"] - X.sum(axis=1)
            cost = int(vp["c"].T @ X.sum(axis=1))
            if np.all(vp["B"] @ X == vp["H"]):
                if np.all(s >= 0):
                    if cost == LB:
                        status = "optimal"
                    else:
                        status = "feasible"
                    over_cap_count = 0
                
                else:
                    status = "over_cap"
                    over_cap_count = -np.sum(s[s<0])
                    # message += "over_cap: {}".format()
            else:
                status = "no_solution_found"
                cost = None
                over_cap_count = None
                
                
            
        res = Problem.Result(status, X, cost, message, LB, over_cap_count)
        p.results[str(__class__)] = res
        return res
    
class Biobjective:
    def solve(p, importance_factor = 0.5, SOLVER = "CBC",known_LB=None):
        if not Problem.verbose:
            save_stdout = sys.stdout
            sys.stdout = open('trash', 'w')
        ###############################
        _, constraints, _, _, vp, _, obj_biobj,_,_  = dai.init_from_graph(p.graph,p.demands)
        
        
        #################################
        if not Problem.verbose:
            sys.stdout = save_stdout
            
        problem = cp.Problem(obj_biobj, constraints)
        
        vp["gama"].value = importance_factor
        vp["zeta"].value = 1 - importance_factor
        problem.solve(verbose=True,solver=SOLVER)
        
        
        X = vp["X"].value
        cost = None
        message = ""
        over_cap_count = None
        
        if problem.status == "infeasible":
            status = problem.status
        else:
        
            cost = int(vp["c"].T @ X.sum(axis=1))
            
            over_flow = (X.sum(axis=1) - vp["cap"])
            over_cap_count = np.sum(np.max( [np.zeros(len(over_flow)), over_flow], axis=0))
            # over_cap_count = 100
            # message = "over_cap_count: {}".format(over_cap_count)
        
        
            if problem.status == "optimal" and over_cap_count == 0:
                if known_LB is not None and cost == known_LB:
                    status = "optimal"
                else:
                    status = "feasible"
            else:
                status = "over_cap"
            
        res = Problem.Result(status, X, cost, message,known_LB,over_cap_count)
        p.results[str(__class__)] = res
        return res
class DnC:
    def divide_demands(p,dim="x"):
        val = sum([ (p.graph.nodes()[O][dim] + p.graph.nodes()[D][dim])/2 for O,D,_ in p.demands]) / len(p.demands)
        demandsL = []
        demandsR = []
        demands_remaining = []
        for O,D,num_k in p.demands:
            Ox = p.graph.nodes()[O][dim]
            Dx = p.graph.nodes()[D][dim]
            if Ox < val and Dx < val:
                demandsL.append((O,D,num_k))
            elif Ox > val and Dx > val:
                demandsR.append((O,D,num_k))
            else:
                demands_remaining.append((O,D,num_k))
        return demandsL,demandsR,demands_remaining
                
            
        
        
    def solve(p,known_LB = None):
        try:
            demandsL,demandsR,demands_remaining = DnC.divide_demands(p,dim="x")
        except:
            i1 = int(len(p.demands)/4)
            i2 = i1*2
            demandsL= p.demands[:i1]
            demandsR= p.demands[i1:i2]
            demands_remaining = p.demands[i2:]
        if len(demandsL) == 0:
            X1 = None
        else:
            Problem.print((demandsL,demandsR,demands_remaining))
            subproblem1 = Problem(p.graph,demandsL)
            res = Biobjective.solve(subproblem1,importance_factor=10e-8)
            X1 = res.X
        if len(demandsR) == 0:
            X2 = None
        else:
            subproblem2 = Problem(p.graph,demandsR)
            res = Biobjective.solve(subproblem2,importance_factor=10e-8)
            X2 = res.X
            
        graph_ = nx.DiGraph(p.graph)
        if X1 is None and X2 is None:
            pass
        else:
            if X1 is None:
                X_sum = X2
            elif X2 is None:
                X_sum = X1
            else:
                X_sum = sparse.hstack([sparse.csr_matrix(X1),sparse.csr_matrix(X2)])
            new_cap = np.array(np.array(p.cap) - X_sum.sum(axis=1).T).flatten()
            print("new_cap.shape: ",new_cap.shape, new_cap[0])
            new_cap[new_cap < 0] = 0
            nx.set_edge_attributes(graph_, {e: float(new_cap[ei]) for ei,e in enumerate(graph_.edges())},"cap")
            
        final_problem = Problem(graph_,demands_remaining)
        res = Biobjective.solve(final_problem,importance_factor=10e-8)
        X3 = res.X
        
        l = []
        if X1 is not None: l.append(sparse.csr_matrix(X1))
        if X2 is not None: l.append(sparse.csr_matrix(X2))
        if X3 is not None: l.append(sparse.csr_matrix(X3))
        X = sparse.hstack(l)
        
        cost = np.array(p.c).T @ X.sum(axis=1)
        cost = cost if isinstance(cost, int) else int(cost[0,0])
        
        over_flow = np.array(X.sum(axis=1).T - np.array(p.cap)).flatten()
        over_cap_count = np.sum(np.max( [np.zeros(len(over_flow)), over_flow], axis=0))
        # message = "over_cap_count: {}".format(over_cap_count)
        message = ""
        
        if over_cap_count == 0:
            if known_LB is not None and known_LB == cost:
                status = "optimal"
            else:
                status = "feasible"
        else:
            status = "over_cap"
            
        res = Problem.Result(status, X, cost, message, known_LB, over_cap_count)
        p.results[str(__class__)] = res
        return res
                  
