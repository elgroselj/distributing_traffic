
import networkx as nx
import numpy as np
import random
import time

from docplex.mp.model import Model

from scipy.optimize import minimize
from scipy import sparse
from scipy.optimize import linprog

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

    class Lingen():
        def get_random_node(nodes):
                i = random.randint(0, len(nodes)-1)
                return (i, nodes[i])

        def get_random_path(g,z,k,path_seed=None,excluded_nodes=[]):
            path_seed = path_seed if path_seed is not None else time.time()
            #print("path_seed: "+str(path_seed))
            random.seed(path_seed)
            
            nodes_to_use = list(set(g.nodes()) - set(excluded_nodes))
            sub_g = g.subgraph(nodes_to_use)
            #while(True):
            for i in range(100):
                # print("*",end="")
                try:
                    _,node = Lingen.get_random_node(nodes_to_use)
                    path1 = nx.shortest_path(sub_g, z, node)
                    path2 = nx.shortest_path(sub_g, node, k)
                except:
                    continue
                path = path1+path2[1:]
                kvazi_path_g = nx.DiGraph(Problem.nodes_to_edges_path(path))
                path_correct = nx.shortest_path(kvazi_path_g, z, k)
                # print()
                return path_correct
            return None
            
        def mutate(g,pop,prob=0.1,path_seed=None):
            mutated_pop = []
            
            for i in range(int(np.ceil(len(pop)*prob))):
                ex = random.sample(pop,1)[0]
                col_id = random.randint(0,ex.shape[1]-1)
                path = Problem.nodes_to_edges_path(Problem.binary_vector_to_edges(ex[:,col_id],g.edges()), inverse=True)
                meja = random.randint(1,len(path)-1)
                i1, node1 = Lingen.get_random_node(path[:meja])
                i2, node2 = Lingen.get_random_node(path[meja:])
                i2 += meja
                
                
                new_path = Lingen.get_random_path(g,node1,node2,path_seed,excluded_nodes=path[:i1] + path[i2+1:])
                if new_path is None:
                    print("An alternative path wasn't found.")
                    continue
                
                mutated_path = path[:i1] + new_path + path[i2+1:]
                # print("path: ", path, path[:i1+1], path[i2:])
                # print("new_path: ", new_path)
                #print(len(mutated_path) == len(set(mutated_path)))
                    
                    
                mutated_ex = np.copy(ex)
                mutated_ex[:,col_id] = Problem.edges_to_binary_vector(Problem.nodes_to_edges_path(mutated_path), list(g.edges()))

                ids = map(id, pop + mutated_pop)
                if id(mutated_ex) not in ids:
                    mutated_pop.append(mutated_ex)
            
            return mutated_pop
        
        def crossover(pop,prob):
            crossovered_pop = []
            for i in range(int(np.ceil(len(pop)*prob))):
                ex1, ex2 = list(random.sample(pop,2))
                col_id = random.randint(0,ex1.shape[1]-1)
                v = ex1[:,col_id]
                ex1[:,col_id] = ex2[:,col_id]
                ex2[:,col_id] = v
                
                ids = map(id, pop + crossovered_pop)
                if id(ex1) not in ids:
                    crossovered_pop.append(ex1)
                if id(ex2) not in ids:
                    crossovered_pop.append(ex2)
                
            return crossovered_pop
        
        
        
        def best(pop, f, num):
            sorted_pop = sorted(pop, key = lambda ex: f(ex))
            return sorted_pop[:min(num,len(pop))]
                
            
        def solve(p,path_seed=None,num_iter=10,num_best=10,mutation_prob=0.5,crossover_prob=0.1):
            
            M = Problem.sparse_incidence_matrix(p.g.nodes(),p.g.edges())
            M_ZK = Problem.sparse_incidence_matrix(p.g.nodes(),[(z,k) for z, k, _ in p.ZK],factor=p.a)
            M_ZK_dupl = np.repeat(M_ZK.todense(), p.a, axis=1)
            
            f_t = lambda ex: np.sum(np.dot(ex.T,p.t))
            f_c = lambda ex: abs(np.sum(x for x in list(p.c - np.sum(ex, axis=1)) if x < 0))
                
            f_sums = lambda ex: abs(np.sum(M@ex - M_ZK_dupl))
            f = lambda ex: f_t(ex) + 1000 * f_c(ex) + 1000 * f_sums(ex)
            
            ex0 = None
            for z,k,a in p.ZK:
                #path = nx.shortest_path(p.g, z, k)
                path = Lingen.get_random_path(p.g, z, k, path_seed)
                col = Problem.edges_to_binary_vector(Problem.nodes_to_edges_path(path), list(p.g.edges()))
                for i in range(a):
                    ex0 = col if ex0 is None else np.column_stack([ex0,col])
            
            pop = [ex0]
            
            for i in range(num_iter):
                #print(len(pop))
                
                mutated_pop = Lingen.mutate(p.g,pop,prob=mutation_prob)
                
                crossovered_pop = []
                crossovered_pop = Lingen.crossover(pop+mutated_pop,prob=crossover_prob)
                
                best_pop = Lingen.best(pop+mutated_pop+crossovered_pop,f,num=num_best)
                pop = best_pop
                
            best_pop = Lingen.best(pop+mutated_pop+crossovered_pop,f,num=1)
            best_ex = best_pop[0]
            paths = Problem.columns_to_paths(p.g,best_ex)
            print("f_c_sums: ", f_c(best_ex), f_sums(best_ex))
            p.results[str(__class__)] = {"X": best_ex, "fun": f_t(best_ex), "message": "glej success", "success": f_c(best_ex) + f_sums(best_ex) == 0, "paths": paths}

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
            

