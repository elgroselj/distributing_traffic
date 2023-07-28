
from matplotlib import pyplot as plt
import networkx as nx
import numpy as np
import random
import time
import copy

from docplex.mp.model import Model

from scipy.optimize import minimize
from scipy import sparse
from scipy.optimize import linprog

import helper_functions as hf

import dai

class Problem():
    def __init__(self, graph, demands):
        self.graph = graph
    
        self.demands = demands
        
        
        def unconnected_places(g,ZK):
            for z,k,_ in ZK:
                try:
                    nx.shortest_path(g, z, k)
                except Exception as e:
                    return(e)
            return None
    
        ups = unconnected_places(graph,demands)
        if ups is not None:
            raise Exception(unconnected_places)
        
        self.results = {}
        
    def __str__(self):
        return self.results
    
    
    def print(p):
        for k in p.results:
            print(k + ":")
            print(p.results[k]["message"], p.results[k]["fun"], p.results[k]["success"])
            
            
            
    class Result():
        def __init__(self, status, X, cost, message=""):
            self.status = status
            self.X = X # to je združene poti po komoditijih # sparse <3
            self.cost = cost
            self.message = message
        
        def __str__(self):
            return self.__repr__()
        
        def __repr__(self):
            return "status: {}, cost: {}, message: {}".format(self.status, self.cost, self.message)
  
        

    class Linprog_v1():
        def getB(num_edges, num_ZK, factors = None):
            B = None
            row_edge = []
            col_ZK = []
            data = []
            q = 0
            for i in range(num_edges):
                for j in range(num_ZK):
                    row_edge.append(i)
                    col_ZK.append(q)
                    if factors is None:
                        data.append(1)
                    else:
                        data.append(factors[i])
                    q += 1
            B = coo_array((data, (row_edge, col_ZK)), shape=(num_edges, num_ZK*num_edges)).tocsr()
            return B
        
        def MtoM2(M,num_ZK):
            num_nodes, num_edges = M.shape
            rows = []
            cols = []
            data = []
            
            shift = 0
            for k in range(num_ZK):
                
                for n in range(num_nodes):
                    for e in range(num_edges):
                        rows.append(k*num_nodes+n)
                        cols.append(e*num_ZK + shift)
                        data.append(M[n,e])
                shift += 1
            csr = coo_array((data, (rows, cols)), shape=(num_ZK*num_nodes, num_ZK*num_edges)).tocsr()
            return csr
        
        def nastavi_fbmm(g,ZK,t):
            f = np.repeat(t,len(ZK)) #min f*x
            
            B = Linprog_v1.getB(g.number_of_edges(),len(ZK)) # B*x <= c
            
            M = Problem.sparse_incidence_matrix(g.nodes(),g.edges()) # (M*X=M_ZK)
            
            a = [a for _,_,a in ZK]
            M_ZK = Problem.sparse_incidence_matrix(g.nodes(),[(z,k) for z, k, _ in ZK],factor=a)
            
            M2 = Linprog_v1.MtoM2(M.todense(), len(ZK)) # M2 * x = m_ZK
            #m_ZK = M_ZK.flatten('F')
            m_ZK = M_ZK.reshape((g.number_of_nodes()*len(ZK),1), order='F', copy=False)
            
            return (f,B,M2,m_ZK)
        
        def solve(p,integrality = 1):
            f,B,M2,m_ZK = Linprog_v1.nastavi_fbmm(p.g,p.ZK,p.t)
            m_ZK = m_ZK.todense()
            res = linprog(f, A_ub=B, b_ub=p.c, A_eq=M2, b_eq=m_ZK, integrality=integrality)
            
            X = None
            paths = None
            if res.x is not None:
                X = np.array(res.x).reshape(p.g.number_of_edges(),len(p.ZK))
                #paths = Problem.columns_to_paths(p.g,X) TODO
            p.results[str(__class__)] = {"X": X, "fun": res.fun, "message": res.message, "success": res.x is not None, "paths": paths}

    class Linprog_v2():
        def sestavi_QB(ZK, G, st_alternativ = 4):
            edges = G.edges()
            #edges_data = G.edges(data=True)
            Q = None
            B = None
            
            for i in range(len(ZK)):
                z,k,_ = ZK[i]
                X = nx.shortest_simple_paths(G, z, k, weight="t")
                try:
                    for counter, path in enumerate(X):
                        col = Problem.edges_to_binary_vector(Problem.nodes_to_edges_path(path), list(edges))
                        Q = col if Q is None else np.column_stack([Q,col])
                        if counter == st_alternativ-1:
                            break
                except Exception as e:
                    print(e)
                else:
                    b = np.repeat(np.eye(len(ZK),1,-i), counter+1, axis=1)
                    B = b if B is None else np.column_stack([B,b])
            
            # t = [d["t"] for _,_,d in edges_data]
            # c = [d["c"] for _,_,d in edges_data]
            
            # a = [ai for _,_,ai in ZK]
            # return(Q,B,t,c,a)
            return(Q,B)
        
        def solve(p,st_alternativ):
            Q,B = Linprog_v2.sestavi_QB(p.ZK, p.g, st_alternativ)
            f = np.dot(p.t, Q)
            
            res = linprog(f, A_ub=Q, b_ub=p.c, A_eq=B, b_eq=p.a, integrality=1)

            #print(res.fun, res.x, res.message)
            #print(res.fun, res.x, res.message)
            Q_used = Q[:,res.x > 0] # binarni vektorji uporabljenih poti
            paths = Problem.columns_to_paths(p.g,Q_used,edges_mode=False)
            #TODO naredi X iz Q_used (al pa obratno??) Q_used !== X !!!!
            p.results[str(__class__)+str(st_alternativ)] = {"X": Q_used, "fun": res.fun, "message": res.message, "success": res.x is not None, "paths": paths}

    class Linmin():
        def solve(p,x0_seed=None):
            t = np.array(p.t).T
            c = np.array(p.c).T
            
            # problematicen je izbor zacetnega priblizka, vecinoma se zatekne v lokalnem minimumu
            # kako vpeljat da gresta 2 avta po neki poti? TODO
            # bolj pameten zacetni priblizek
            num_ZK = len(p.ZK)
            num_edges = p.g.number_of_edges()
            
            # X0 = np.array([[1,0],[0,0],[1,1],[0,1],[0,0],[1,0],[0,1]])
            # x0 = X0.reshape(-1)
            x0_seed = x0_seed if x0_seed is not None else time.time()
            print("x0_seed: "+str(x0_seed))
            random.seed(x0_seed)
            x0 = np.random.rand(num_edges * num_ZK) #* max(a)
            # print(x0)
            
            # f,B,M2,m_ZK = Linprog_v1.nastavi_fbmm(p.g,p.ZK,p.t)
            # m_ZK = m_ZK.todense().reshape(-1)
            # #res = linprog(f, A_ub=B, b_ub=p.c, A_eq=M2, b_eq=m_ZK, integrality=integrality)
            # B = B.todense()
            # M2 = M2.todense()
            # def fun(x):
            #     r = np.dot(f,x)
            #     return r


            # constraints = []
            # constraints =  [{'type': 'ineq', 'fun': lambda x: x}]
            
            # def con1(x):
            #     r = c - np.dot(B,x)
            #     return r.reshape(-1)
            # constraints.append({'type': 'ineq', 'fun': con1})
            
            # def con2(x):
            #     r = np.dot(M2,x) - m_ZK
            #     return r
            # constraints.append({'type': 'eq', 'fun': con2})        
            
            
            M = Problem.sparse_incidence_matrix(p.g.nodes(),p.g.edges())
            M_ZK = Problem.sparse_incidence_matrix(p.g.nodes(),[(z,k) for z, k, _ in p.ZK],factor=p.a)
            
            def fun(x):
                X = x.reshape(num_edges,num_ZK)
                return np.sum(t.T@X)



            constraints = []
            constraints =  [{'type': 'ineq', 'fun': lambda x: x}]
            
            def con1(x):
                X = x.reshape(num_edges,num_ZK)
                unit = np.ones((num_ZK, 1))
                r = p.c - X@unit
                return np.repeat(r,x.size//r.size)
            constraints.append({'type': 'ineq', 'fun': con1})
            
            def con2(x):
                X = x.reshape(num_edges,num_ZK)
                r = M@X - M_ZK
                r = r.reshape(-1)
                return np.repeat(r,x.size//r.size)
            
            constraints.append({'type': 'eq', 'fun': con2})
            print(constraints)
            res = minimize(fun, x0, method='SLSQP', constraints=constraints, bounds=np.array((np.zeros(x0.size),np.ones(x0.size))).T)
            if not res.success:
                # Izpišemo optimalno rešitev
                # print('(NE)Optimalna rešitev:')
                # print(res.x, res.fun, res.message, res.success)

                #if res.x is not None:
                p.results[str(__class__)] = {"X": None, "fun": None, "message": res.message, "success": res.success, "paths": None}
                
            else:
                def smart_round(x):
                    x[x > 0.5] = 1
                    x[x < 0.5] = 0
                    return x
                
                round_x = smart_round(res.x)
                # print(round_x)
                # print(f(round_x))
                
                X = round_x.reshape(num_edges,num_ZK)
                paths = Problem.columns_to_paths(p.g,X)
                
                p.results[str(__class__)] = {"X": X, "fun": fun(round_x), "message": res.message, "success": res.success, "paths": paths}

    # class Lingen():
    #     def get_random_node(nodes):
    #             i = random.randint(0, len(nodes)-1)
    #             return (i, nodes[i])

    #     def get_random_path(g,z,k,path_seed=None,excluded_nodes=[]):
    #         path_seed = path_seed if path_seed is not None else time.time()
    #         #print("path_seed: "+str(path_seed))
    #         random.seed(path_seed)
            
    #         nodes_to_use = list(set(g.nodes()) - set(excluded_nodes))
    #         sub_g = g.subgraph(nodes_to_use)
    #         #while(True):
    #         for i in range(100):
    #             # print("*",end="")
    #             try:
    #                 _,node = Problem.Lingen.get_random_node(nodes_to_use)
    #                 path1 = nx.shortest_path(sub_g, z, node)
    #                 path2 = nx.shortest_path(sub_g, node, k)
    #             except:
    #                 continue
    #             path = path1+path2[1:]
    #             kvazi_path_g = nx.DiGraph(Problem.nodes_to_edges_path(path))
    #             path_correct = nx.shortest_path(kvazi_path_g, z, k)
    #             # print()
    #             return path_correct
    #         return None
            
    #     def mutate(g,pop,prob=0.1,path_seed=None):
    #         mutated_pop = []
            
    #         for i in range(int(np.ceil(len(pop)*prob))):
    #             ex = random.sample(pop,1)[0]
    #             col_id = random.randint(0,ex.shape[1]-1)
    #             path = Problem.nodes_to_edges_path(Problem.binary_vector_to_edges(ex[:,col_id],g.edges()), inverse=True)
    #             meja = random.randint(1,len(path)-1)
    #             i1, node1 = Problem.Lingen.get_random_node(path[:meja])
    #             i2, node2 = Problem.Lingen.get_random_node(path[meja:])
    #             i2 += meja
                
                
    #             new_path = Problem.Lingen.get_random_path(g,node1,node2,path_seed,excluded_nodes=path[:i1] + path[i2+1:])
    #             if new_path is None:
    #                 print("An alternative path wasn't found.")
    #                 continue
                
    #             mutated_path = path[:i1] + new_path + path[i2+1:]
    #             # print("path: ", path, path[:i1+1], path[i2:])
    #             # print("new_path: ", new_path)
    #             #print(len(mutated_path) == len(set(mutated_path)))
                    
                    
    #             mutated_ex = np.copy(ex)
    #             mutated_ex[:,col_id] = Problem.edges_to_binary_vector(Problem.nodes_to_edges_path(mutated_path), list(g.edges()))

    #             ids = map(id, pop + mutated_pop)
    #             if id(mutated_ex) not in ids:
    #                 mutated_pop.append(mutated_ex)
            
    #         return mutated_pop
        
    #     def crossover(pop,prob):
    #         crossovered_pop = []
    #         for i in range(int(np.ceil(len(pop)*prob))):
    #             ex1, ex2 = list(random.sample(pop,2))
    #             col_id = random.randint(0,ex1.shape[1]-1)
    #             v = ex1[:,col_id]
    #             ex1[:,col_id] = ex2[:,col_id]
    #             ex2[:,col_id] = v
                
    #             ids = map(id, pop + crossovered_pop)
    #             if id(ex1) not in ids:
    #                 crossovered_pop.append(ex1)
    #             if id(ex2) not in ids:
    #                 crossovered_pop.append(ex2)
                
    #         return crossovered_pop
        
        
        
    #     def best(pop, f, num):
    #         sorted_pop = sorted(pop, key = lambda ex: f(ex))
    #         return sorted_pop[:min(num,len(pop))]
                
            
    #     def solve(p,path_seed=None,num_iter=10,num_best=10,mutation_prob=0.5,crossover_prob=0.1):
            
    #         M = Problem.sparse_incidence_matrix(p.g.nodes(),p.g.edges())
    #         M_ZK = Problem.sparse_incidence_matrix(p.g.nodes(),[(z,k) for z, k, _ in p.ZK],factor=p.a)
    #         M_ZK_dupl = np.repeat(M_ZK.todense(), p.a, axis=1)
            
    #         f_t = lambda ex: np.sum(np.dot(ex.T,p.t))
    #         f_c = lambda ex: abs(np.sum(x for x in list(p.c - np.sum(ex, axis=1)) if x < 0))
                
    #         f_sums = lambda ex: abs(np.sum(M@ex - M_ZK_dupl))
    #         f = lambda ex: f_t(ex) + 1000 * f_c(ex) + 1000 * f_sums(ex)
            
    #         ex0 = None
    #         for z,k,a in p.ZK:
    #             #path = nx.shortest_path(p.g, z, k)
    #             path = Problem.Lingen.get_random_path(p.g, z, k, path_seed)
    #             col = Problem.edges_to_binary_vector(Problem.nodes_to_edges_path(path), list(p.g.edges()))
    #             for i in range(a):
    #                 ex0 = col if ex0 is None else np.column_stack([ex0,col])
            
    #         pop = [ex0]
            
    #         for i in range(num_iter):
    #             #print(len(pop))
                
    #             mutated_pop = Problem.Lingen.mutate(p.g,pop,prob=mutation_prob)
                
    #             crossovered_pop = []
    #             crossovered_pop = Problem.Lingen.crossover(pop+mutated_pop,prob=crossover_prob)
                
    #             best_pop = Problem.Lingen.best(pop+mutated_pop+crossovered_pop,f,num=num_best)
    #             pop = best_pop
                
    #         best_pop = Problem.Lingen.best(pop+mutated_pop+crossovered_pop,f,num=1)
    #         best_ex = best_pop[0]
    #         paths = Problem.columns_to_paths(p.g,best_ex)
    #         print("f_c_sums: ", f_c(best_ex), f_sums(best_ex))
    #         p.results[str(__class__)] = {"X": best_ex, "fun": f_t(best_ex), "message": "glej success", "success": f_c(best_ex) + f_sums(best_ex) == 0, "paths": paths}


    class Greedy():
        def update_path(paths, generators, path_ind):
            paths[path_ind] = next(generators[path_ind])


        def solve(p, max_num_iter = 10):
            generators = []
            for z,k,a in p.ZK:
                for i in range(a):
                    generator = nx.shortest_simple_paths(p.g, z, k, weight="t")
                    generators.append(generator)
            
            paths = [None] * len(generators) # init
            
            paths_to_update = list(range(0,len(generators))) # update all
            
            st = 0
            while len(paths_to_update) > 0 and st < max_num_iter:
                for path_ind in paths_to_update:
                    Greedy.update_path(paths, generators, path_ind)
                    
                Q = None
                for path in paths:
                    col = Problem.edges_to_binary_vector(Problem.nodes_to_edges_path(path), list(p.g.edges()))
                    Q = col if Q is None else np.column_stack([Q,col])
                    
                Q = Q.reshape((p.g.number_of_edges(), -1))
                edges_to_unload = np.where(np.sum(Q,axis=1) > p.c)[0]
                paths_in_jam = np.where(np.sum(Q[edges_to_unload,:], axis=0) > 0)[0]
                paths_to_update = [paths_in_jam[0]] if len(paths_in_jam) > 0 else [] # one per step
                
                st += 1
            
            f_t = lambda Q: np.sum(np.dot(Q.T,p.t))
            
            f_t(Q)
            p.results[str(__class__)] = {"X": Q, "fun": f_t(Q), "message": "glej success", "success": st < max_num_iter, "paths": paths}

    class Cplex_intlinprog():
        def cplex_solver_ex(c, cap, B, H):
            n,m = B.shape
            _,t = H.shape

            milp_model = Model(name = "myMILP")

            commodities = np.arange(t)

            edges = np.arange(m)#list(graph.edges())
            nodes = np.arange(n)#list(graph.nodes())

            X = milp_model.integer_var_matrix(edges, commodities, lb = 0, name="X")

            cost_expr = milp_model.sum(X[edge, com] * c[edge] for edge in edges for com in commodities)

            for node in nodes:
                for com in commodities:
                    milp_model.add_constraint(milp_model.sum(X[edge, com] * B[node, edge] for edge in edges)  == H[node, com])

            for edge in edges:
                milp_model.add_constraint(milp_model.sum(X[edge, com] for com in commodities)  <=  cap[edge])


            milp_model.set_objective("min", cost_expr)

            milp_model.print_information()
            
        
            sol = milp_model.solve()
            
            # print(sol.solve_details.status)
            # milp_model.print_solution()
            # X_dict = {}
            # for com in commodities:
            #     for edge in edges:
            #         if sol[X[edge,com]] != 0:
            #             X_dict[(edge,com)] = sol[X[edge,com]]
            # print(X_dict)
            
            X_res = sparse.dok_matrix((len(edges),len(commodities)))
            for com in commodities:
                for edge in edges:
                    if sol[X[edge,com]] != 0:
                        X_res[edge,com] = sol[X[edge,com]]
            return (sol.solve_details.status, X_res, sol.objective_value)
        
        def solve(p):
            c = np.round(np.array([data["c"] for _,_, data in p.graph.edges(data=True)])) # BŠS
            cap = np.floor(np.array([data["cap"] for _,_, data in p.graph.edges(data=True)])) # caps floor to int

            B = -1 * nx.incidence_matrix(p.graph,oriented=True)

            n,m = B.shape

            def demands_to_matrix_sparse(demands,n):
                H = sparse.dok_matrix((n,len(demands)))
                for k,(Ok,Dk,d) in enumerate(demands):
                    H[Ok,k] = d
                    H[Dk,k] = -d
                return H  
            
            H = demands_to_matrix_sparse(p.demands,n)
            
            status, X, cost = Problem.Cplex_intlinprog.cplex_solver_ex(c, cap, B, H)
            res = Problem.Result(status, X, cost)
            p.results[str(__class__)] = res
            return res
    
    class Dai_solver():
        
        def solve(p,MAX_ITER,MAX_ITER_LR):
            status, X, cost = dai.dai_solve(p.graph,p.demands,MAX_ITER,MAX_ITER_LR)
            res = Problem.Result(status, X, cost)
            p.results[str(__class__)] = res
            return res
            

    class Evolution():
        
        
        class Subject:
            id = 0
            
            
            
            def __init__(self,l,graph):
                self.id = Problem.Evolution.Subject.id
                Problem.Evolution.Subject.id += 1
                self.l = l
                self.graph = graph
            
                self.Q = None
                self.X = None
                
                self.over_cap = None
                self.cost = None
            
            def evaluate(self):
                if self.cost is not None and self.over_cap is not None:
                    return (self.cost, self.over_cap)
                
                r = sum([len(k) for k in self.l]) # stevilo poti za določit (št posameznih avtov)
                t = len(self.l) # stevilo skupin avtov
                m = len(self.graph.edges())
                self.Q = sparse.dok_matrix((m,r)) # Q je razprta različica X
                self.X = sparse.dok_matrix((m,t)) # Q je razprta različica X
                col = 0
                col_x = 0
                for k in self.l:
                    for path in k:
                        epath = hf.nodes_to_edges_path(path)
                        for edge in epath:
                            row = list(self.graph.edges()).index(edge)
                            self.Q[row,col] = 1
                            try:
                                self.X[row,col_x] += 1
                            except:
                                self.X[row,col_x] = 1
                            
                        col += 1
                    col_x +=1
                c = [self.graph.edges()[e]["c"] for e in self.graph.edges()]
                cap = [self.graph.edges()[e]["cap"] for e in self.graph.edges()]
                
                edge_flow = np.sum(self.Q, axis=1)
                self.over_cap = len(np.where(np.array(cap) < edge_flow))
                self.cost = int((np.array(c) @ edge_flow)[0])
                
                return self.over_cap, self.cost
                

            def __lt__(self, other):
                over_cap, cost = self.evaluate()
                over_cap_, cost_ = other.evaluate()
                # less je boljše
                if over_cap < over_cap_:
                    return True
                if over_cap == over_cap_:
                    if cost < cost_:
                        return True
                
                return False

            def __repr__(self) -> str:
                return "over_cap: {}, cost: {}, l: {}".format(self.over_cap, self.cost, self.l)
                
        
        # micro ops
        def mutate_path(graph_,path,depth_limit=3):
            if len(path) <= 2:
                return path
            old_path = list(path)
            node_i = random.sample(range(len(path)), 1)[0]
            node = path[node_i]
            
            rev = False
            if node_i > len(path)/2:
                graph = nx.reverse(graph_)
                path.reverse()
                node_i = len(path) - node_i -1
                rev = True
            else:
                graph = graph_
            
            skip_nodes = path[:node_i]
            sub_path = [node]
            for n0,n in hf.dfs_edges_random(graph,source=node,depth_limit=depth_limit,skip_nodes=skip_nodes):
                
                if n0 not in sub_path: continue
                sub_path = sub_path[:sub_path.index(n0)+1]
                
                
                # if n in path[:node_i]: # no loops
                #     continue
                
                
                     
                sub_path.append(n)
                
                # finish
                if n in path[node_i:]:
                    i1 = path.index(sub_path[0])
                    i2 = path.index(sub_path[-1])
                    path_ = path[:i1] + sub_path[:-1] + path[i2:]
                    if graph_ != graph:
                        path_.reverse()
                    return path_
            # if graph_ != graph:
            #     path.reverse()
            return old_path
        
        def crossover_paths_k(path1, path2):
            common = list(set(path1).intersection(set(path2)))
            
            if len(common) == 2: return (path1, path2)
            
            shift_node_i = random.sample(range(1,len(common)-1),1)[0]
            shift_node = common[shift_node_i]
            
            i1 = path1.index(shift_node)
            i2 = path2.index(shift_node)
            
            path1_ = path1[:i1] + path2[i2:]
            path2_ = path2[:i2] + path1[i1:]
            
            return path1_, path2_
        
        # macro ops
            
        def mutate(pop_,graph,prob=0.1):
            pop = list(pop_)
            mutated_pop = []
            for subj in pop:
                subj_ = None
                for ki, k in enumerate(subj.l):
                    for i, path in enumerate(k):
                        fl = random.random()
                        if fl < prob:
                            # new_path = Problem.Evolution.mutate_path(graph,path,depth_limit=graph.number_of_edges()/3)
                            new_path = Problem.Evolution.mutate_path(graph,list(path),depth_limit=graph.number_of_edges()/15)
                            if k[i] != new_path:
                                subj_ = Problem.Evolution.Subject(copy.deepcopy(subj.l),subj.graph) if subj_ is None else subj_
                                subj_.l[ki][i] = new_path
                if subj_ is not None:
                    mutated_pop.append(subj_)
            return mutated_pop
        
        def crossover(pop_,prob):
            pop = list(pop_)
            
            overcrossed_pop = []
            # prva sorta, znotraj dobrin
            for subj in pop:
                subj_ = None
                for ki, k in enumerate(subj.l):
                    if random.random() < prob and len(k) >= 2:
                        i1, i2 = random.sample(range(len(k)),2)
                        path1 = k[i1]
                        path2 = k[i2]
                        
                        if path1 != path2:
                            path1_, path2_ = Problem.Evolution.crossover_paths_k(list(path1),list(path2))
                            if path1_ != path1:
                                subj_ = Problem.Evolution.Subject(copy.deepcopy(subj.l),subj.graph) if subj_ is None else subj_
                                subj_.l[ki][i1] = path1_
                                subj_.l[ki][i2] = path2_
                if subj_ is not None:
                    overcrossed_pop.append(subj_)
            
            # druga sorta, izmenjava dobrin
            if random.random() < prob and len(pop) >= 2:
                subj1, subj2 = random.sample(pop,2)
                t = len(subj1.l)-1
                k = random.randint(0,t-1)
                if subj1.l[k] != subj2.l[k]:
                    v = subj1.l[k]
                    subj1.l[k] = subj2.l[k]
                    subj2.l[k] = v
                    overcrossed_pop += [subj1, subj2]
                    
                
            return overcrossed_pop
                  
            
        
        def selection(pop,max_pop_size):
            if len(pop) < max_pop_size:
                return sorted(pop)
            
            return sorted(pop)[:max_pop_size]
        
        def init_pop(graph, demands):
            l = []
            for O, D, num_k in demands:
                path = nx.dijkstra_path(graph,O,D)
                l.append(num_k*[path])
            subj = Problem.Evolution.Subject(l,graph)
            return [subj]

                
            
        def solve(p,num_iter=10,max_pop_size=1,mutation_prob=0.5,crossover_prob=0.5):
            
            pop = Problem.Evolution.init_pop(p.graph, p.demands)
            
            for subj in pop:
                    print(subj)
            
            for i in range(num_iter):
                #print(len(pop))
                
                mutated_pop = Problem.Evolution.mutate(pop,p.graph,prob=mutation_prob)
                print("num_mutants: ", len(mutated_pop))
                
                crossovered_pop = []
                crossovered_pop = Problem.Evolution.crossover(pop+mutated_pop,prob=crossover_prob)
                print("num_crossovered: ", len(crossovered_pop))
                
                
                best_pop = Problem.Evolution.selection(pop+mutated_pop+crossovered_pop,max_pop_size=max_pop_size)
                
                for subj in pop+mutated_pop+crossovered_pop:
                    print(subj)
                    
                pop = best_pop
                
                for subj in pop:
                    hf.plot_solution_graph(p.graph,subj.X,with_labels=True,font_size=5,figure_size=(20,20))
                    print(subj)
                
            # best_pop = Problem.Lingen.best(pop+mutated_pop+crossovered_pop,f,num=1)
            # best_ex = best_pop[0]
            # paths = Problem.columns_to_paths(p.g,best_ex)
            # print("f_c_sums: ", f_c(best_ex), f_sums(best_ex))
            # p.results[str(__class__)] = {"X": best_ex, "fun": f_t(best_ex), "message": "glej success", "success": f_c(best_ex) + f_sums(best_ex) == 0, "paths": paths}
    
    
    class Descent:
        class Node:
            def init_class(graph, demands, paths, generators):
                Problem.Descent.Node.generators = generators
                Problem.Descent.Node.paths = paths
                Problem.Descent.Node.graph = graph
                Problem.Descent.Node.demands = demands
                Problem.Descent.Node.label = 0
                
                
            def __init__(self,path_ids, parent = None):
                self.path_ids = [tuple(k) for k in path_ids]
                
                
                self.over_cap = None
                self.cost = None
                
                
                self.children = []
                self.parent = parent
                
                self.hash = None
                
                self.label = Problem.Descent.Node.label
                Problem.Descent.Node.label += 1
                
                self.level = 0 if self.parent is None else self.parent.level +1
            
            def evaluate(self):
                generators = Problem.Descent.Node.generators
                paths = Problem.Descent.Node.paths
                graph = Problem.Descent.Node.graph
                
                if self.cost is not None and self.over_cap is not None:
                    return (self.cost, self.over_cap)
                
                r = sum([len(k) for k in self.path_ids]) # stevilo poti za določit (št posameznih avtov)
                t = len(self.path_ids) # stevilo skupin avtov
                m = len(graph.edges())
                self.Q = sparse.dok_matrix((m,r)) # Q je razprta različica X
                self.X = sparse.dok_matrix((m,t)) # Q je razprta različica X
                col = 0
                col_x = 0
                for ki, k in enumerate(self.path_ids):
                    for path_id in k:
                        path = Problem.Descent.get_i_th_path(paths,generators,ki,path_id)
                        epath = hf.nodes_to_edges_path(path)
                        for edge in epath:
                            row = list(self.graph.edges()).index(edge)
                            self.Q[row,col] = 1
                            try:
                                self.X[row,col_x] += 1
                            except:
                                self.X[row,col_x] = 1
                            
                        col += 1
                    col_x +=1
                self.Q = self.Q.tolil()
                self.X = self.X.tolil()
                c = [graph.edges()[e]["c"] for e in graph.edges()]
                cap = [graph.edges()[e]["cap"] for e in graph.edges()]
                
                # edge_flow = np.array(np.sum(self.Q, axis=1)).flatten()
                edge_flow = self.Q.sum(axis=1).T
                # over_cap_edges = np.where(np.array(cap) < edge_flow)[0]
                _,over_cap_edges,_ = sparse.find(cap < edge_flow)
                mat = self.Q[over_cap_edges,:] > 0
                # self.over_cap = {e : np.where(mat[ei,:]>0) for ei, e in enumerate(over_cap_edges)}
                self.over_cap = {e : list(sparse.find(mat[ei,:]>0)[1]) for ei, e in enumerate(over_cap_edges)}
                # self.cost = int((np.array(c) @ edge_flow)[0])
                self.cost = (c @ edge_flow.T)[0,0]
                
                return self.over_cap, self.cost
                
                
            def get_child(self,evaluated):
                generators = Problem.Descent.Node.generators
                paths = Problem.Descent.Node.paths
                graph = Problem.Descent.Node.graph
                demands = Problem.Descent.Node.demands
                
                # c = np.array([graph.edges()[e]["c"] for e in graph.edges()])
                cap = np.array([graph.edges()[e]["cap"] for e in graph.edges()])
                
                edge_flow_before = self.Q.sum(axis=1).T
                diff_before = edge_flow_before-cap
                over_cap_count_before = diff_before[diff_before > 0].sum()
                
                over_cap_edges = random.sample(self.over_cap.keys(),len(self.over_cap)) # TODO
                for over_cap_edge in over_cap_edges:
                    over_cap_paths_ids = self.over_cap[over_cap_edge]
                    over_cap_paths_ids = random.sample(over_cap_paths_ids,len(over_cap_paths_ids)) # TODO
                    for over_cap_path_id in over_cap_paths_ids:
                        id = 0
                        for ki, (_,_,num_k) in enumerate(demands):
                            id+=num_k
                            if id > over_cap_path_id:
                                pi = id - over_cap_path_id -1
                                break
                        # ki ... skupina od poti
                        # pi ... idx poti znotraj skupine
                        mat = self.Q.copy()
                        # over_cap_count_before = sum([len(self.over_cap[e]) for e in self.over_cap])
                        paths_generator = Problem.Descent.generate_path(paths,generators,ki)
                        for path_id, path in enumerate(paths_generator):
                            path_vec = hf.edges_to_binary_vector(hf.nodes_to_edges_path(path),list(graph.edges()))
                            mat[:,over_cap_path_id] = path_vec
                            
                            # edge_flow = np.sum(mat, axis=1).T
                            edge_flow = mat.sum(axis=1).T
                            diff = edge_flow-cap
                            over_cap_count = diff[diff > 0].sum()
                            if over_cap_count < over_cap_count_before:
                                print(path)
                                path_ids = list(self.path_ids)
                                path_ids[ki] = list(path_ids[ki])
                                path_ids[ki][pi] = path_id
                                path_ids[ki] = sorted(path_ids[ki])
                                path_ids[ki] = tuple(path_ids[ki])
                                ch = Problem.Descent.Node(path_ids,parent=self)
                                if ch in evaluated:
                                    print("ch was in evaluated")
                                    Problem.Descent.Node.label -= 1
                                    continue
                                    
                                self.children.append(ch)
                                print("success")
                                yield ch
                            # elif over_cap_count == over_cap_count_before:
                            # korak v stran TODO
                            print("*alternative paths exhausted")
                    print("**over_cap_paths exhausted")
                print("***over_cap_edges exhausted")
                            
                                
                                
                
            def __eq__(self, other):
                return self.path_ids == other.path_ids  
            
            def __hash__(self) -> int:
                if self.hash is not None: return self.hash
                self.hash = 0
                for k in self.path_ids:
                    self.hash += tuple(k).__hash__()
                return self.hash
            
            def print_tree(self, level=0):
                ret = "\t"*level+repr(self)+"\n"
                if self.children is not None:
                    for child in self.children:
                        ret += child.print_tree(level+1)
                return ret
            
            def __repr__(self) -> str:
                return "({}):{} {} = {}".format(self.level, self.label, repr(self.path_ids),self.over_cap)
        
        def get_i_th_path(paths, generators, k, i):
            while len(paths[k]) <= i:
                next_path = next(generators[k])
                paths[k].append(next_path)
            return paths[k][i]
        
        def generate_path(paths,generators, k):
            i=0
            while True:
                path = Problem.Descent.get_i_th_path(paths, generators, k, i)
                yield path
                i += 1

        # def init_sol():
        #     X_dict = {}
        #         solution_found = True
        #         _,t = self.vp["H"].shape
        #         _,m = self.vp["B"].shape 
                
        #         remaining_cap = {e: graph.edges()[e]["cap"] for e in graph.edges()}
        #         for k in range(t):
        #             Ok, Dk, num_k = demands[k]
        #             # ubf
        #             ubf = get_lower_upper_bounds_on_flow(k,lower=False)
        #             remaining_cap_k = dict(remaining_cap)
        #             for e in ubf:
        #                 remaining_cap_k[e] = min(remaining_cap[e], ubf[e])
        #             nx.set_edge_attributes(graph, remaining_cap_k, "remaining_cap")
                
                    
        #             # lbf
        #             lbf = get_lower_upper_bounds_on_flow(k, lower=True)
                    
        #             c = {e: graph.edges()[e]["c"] for e in graph.edges()}
        #             for e in graph.edges():
        #                 if graph.edges()[e]["remaining_cap"] == 0:
        #                     c[e] = np.inf # nasičene / s kapaciteto 0 naj bodo neskončno drage
        #             alt_c = {e: (0 if e in lbf else c[e]) for e in c} # te k rabijo bit uporabljene naj bojo zastonj
                    
        #             nx.set_edge_attributes(graph, alt_c, "alt_c")
                    
        #             On = list(graph.nodes())[Ok]
        #             Dn = list(graph.nodes())[Dk]
                    
        #             for ki in range(num_k):
                        
        #                 path_of_nodes = nx.dijkstra_path(graph,On,Dn,"alt_c")
        #                 path = nodes_to_edges_path(path_of_nodes)
        #                 length = sum([graph.edges()[e]["alt_c"] for e in path])
        #                 print(length)
        #                 etoei = lambda e : list(graph.edges()).index(e)
        #                 path_dict = {(etoei(e),k):1 for e in path}
        #                 X_dict = Counter(X_dict) + Counter(path_dict)
                        
        #                 for e in path:
        #                     if e in lbf:
        #                         lbf[e] -= 1
        #                         if lbf[e] == 0:
        #                             del lbf[e]
        #                             graph.edges()[e]["alt_c"] = alt_c[e]
                        
        #                 for e in path:
        #                     graph.edges()[e]["remaining_cap"] -= 1
        #                     remaining_cap[e] -= 1
        #                     if graph.edges()[e]["remaining_cap"] == 0:
        #                         graph.edges()[e]["alt_c"] = np.inf
        #                     elif remaining_cap[e] < 0:
        #                         # not solution
        #                         solution_found = False
        #                         break
                
        #             if len(lbf) > 0:
        #                 # not solution    
        #                 print("lbf not empty")
        #                 solution_found = False
        #                 break
                             
        #         if solution_found:
        #             print(" heuristic success")
        #             sol["status"] += " heuristic "  
        #             X = sparse.dok_matrix((m,t))
        #             for a,k in X_dict:
        #                 X[a,k] = X_dict[(a,k)]
        #             X = X.todense()


        def solve(p, max_num_iter = 10):
            generators = []
            for O,D,num_k in p.demands:
                generator = nx.shortest_simple_paths(p.graph, O, D, weight="c")
                generators.append(generator)
            
            t = len(generators)
        
            paths = [[] for i in range(t)]
            
            Problem.Descent.Node.init_class(p.graph, p.demands, paths, generators)
            
            path_ids = [[0] * num_k for _,_,num_k in p.demands]
            node0 = Problem.Descent.Node(path_ids)
            node = node0
            
            evaluated = set()
            feasible = set()
            
            st = 0
            while st < max_num_iter:
                st += 1
                if node is None:
                    print("node was None")
                    break
                
                if node not in evaluated:
                    over_cap, _ = node.evaluate()
                    plt.show()
                    evaluated.add(node)
                    
                    if len(over_cap) == 0: # dopustna rešitev
                        # izboljšaj kar se da TODO
                        feasible.add(node)
                        node = node.parent
                # hf.plot_solution_graph(p.graph,node.X,with_labels=True,font_size=10,figure_size=(20,20))    
                try:
                    print(st-1)
                    ch = next(node.get_child(evaluated))
                    print(st-1)
                    print(node,ch)
                except:
                    print("no more ch")
                    node = node.parent
                    continue
                    
                # if ch in evaluated:
                #     print("already evaluated")
                #     continue
                node = ch
                 
            print("len_evaluated:", len(evaluated))               
            print("len_feasible:", len(feasible))
            tree = node0.print_tree()
            print(tree)
                
            message = ""
            if len(feasible) == 0:
                status = "no_feasible_found"
                X = None
                cost = None
            else:
                best_feasible_node = min(feasible, key= lambda node: node.cost)
                best_node = min(evaluated,  key= lambda node: node.cost)
                if best_feasible_node.cost <= best_node.cost:
                    status = "optimal"
                else:
                    status = "suboptimal"
                    message = "at most {} percent of optimum ({})".format(best_feasible_node.cost/best_node.cost * 100,best_node.cost)
                X = best_feasible_node.X
                cost = best_feasible_node.cost
            
            res = Problem.Result(status, X, cost, message)
            p.results[str(__class__)] = res
            return res
            

                
