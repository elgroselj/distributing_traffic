from itertools import compress
import networkx as nx
import random
from scipy.optimize import linprog
import numpy as np
from scipy.sparse import coo_array
import osmnx as ox
import matplotlib.pyplot as plt
COLORS="brgymcbrgymc"
from scipy.optimize import minimize
import time

class Problem():
    def __init__(self, G=None, g=None, ZK=None, ZK_shape=None, ZK_seed = None, c=None, t=None, c_mode=None, t_mode=None):
        if G is None and g is None:
            raise Exception("Please provide G (ox) or g (nx).")
        elif g is None:
            g = nx.DiGraph(G)
        self.G = G
        self.g = g
        
        if ZK is None:
            if ZK_shape is None:
                raise Exception("Please provide ZK xor ZK_shape.")
            num_ZK,max_a = ZK_shape
            ZK = Problem.get_random_ZK(g,num_ZK,max_a,ZK_seed)
        self.ZK = ZK
        
        self.a = [ai for _,_,ai in ZK]
        
        unconnected_places = Problem.unconnected_places(g,ZK)
        if unconnected_places is not None:
            raise Exception(unconnected_places)
        
        if c is None and c_mode is None:
                raise Exception("Please provide c xor c_mode.")
        self.c = Problem.nastavi_c(g,c_mode,c)
        self.t = Problem.nastavi_t(g,t_mode,t)
        
        self.results = {}
    
    def unconnected_places(g,ZK):
        for z,k,_ in ZK:
            try:
                nx.shortest_path(g, z, k)
            except Exception as e:
                return(e)
        return None
    
    def get_random_ZK(g,num_ZK=2,max_a=2,ZK_seed=None):
        ZK_seed = ZK_seed if ZK_seed is not None else random.randint(1,1000)
        print("ZK_seed: " + str(ZK_seed))
        random.seed(ZK_seed)
        def get_random_node(g):
            nodes = list(g.nodes())
            i = random.randint(0, len(nodes)-1)
            return nodes[i]
        ZK = []
        for i in range(num_ZK):
            ZK.append((get_random_node(g),
                    get_random_node(g),
                    random.randint(1,max_a)))
        return ZK
    
    def fill_maxspeed(g):
        for ke in g.edges():
            e = g.edges()[ke]
            if "maxspeed" in e.keys():
                neki = e["maxspeed"]
                if not isinstance(neki,int):
                    e["maxspeed"] = int(neki) if isinstance(
                        neki, str) else int(list(neki)[0])
                else:
                    print("Maxspeed already fixed.")
                    return
            else:
                success = False
                for ke2 in g.edges():  # if nan set to neighbour
                    if (ke[0] in ke2 or ke[1] in ke2) and ("maxspeed" in g.edges()[ke2].keys()):
                        neki = g.edges()[ke2]["maxspeed"]
                        if not isinstance(neki,int):
                            neki = int(neki) if isinstance(neki, str) else int(list(neki)[0])
                        if neki != 999:
                            e["maxspeed"] = neki
                            success = True
                            break
                if not success:
                    e["maxspeed"] = 999

                # print(e["highway"], end="")
                # print(" is set to " + str(e["maxspeed"]))
    
    def nastavi_c(g, c_mode=None, c=None):
        if c is not None:
            cs = {e : c[i] for i, e in enumerate(g.edges())}
        elif c_mode == "maxspeed":
            Problem.fill_maxspeed(g)
            maxspeed_to_capacity = lambda maxspeed: int(maxspeed/10)
            cs = {e : maxspeed_to_capacity(g.edges()[e]["maxspeed"]) for e in g.edges()}
        elif isinstance(c_mode, int):
            cs = {e: c_mode for e in g.edges()}
        else:
            raise Exception("Invalid mode or c missing.")
        nx.set_edge_attributes(g, cs, name='c')
        
        return [d["c"] for _,_,d in g.edges(data=True)]


    def nastavi_t(g, t_mode=None, t=None):
        if t is not None:
            ts = {e : t[i] for i, e in enumerate(g.edges())}
        elif t_mode == "length":
            ts = {e: g.edges()[e]["length"] for e in g.edges()}
        elif t_mode == "time":
            Problem.fill_maxspeed(g)
            ts = {e: g.edges()[e]["length"]/g.edges()[e]["maxspeed"] for e in g.edges()}
        elif isinstance(t_mode, int):
            ts = {e: t_mode for e in g.edges()}
        else:
            raise Exception("Invalid mode or t missing.")
        nx.set_edge_attributes(g, ts, name='t')
        
        return [d["t"] for _,_,d in g.edges(data=True)]
    
    
    def sparse_incidence_matrix(nodes, edges, factor = None):
        nodes = list(nodes)
        edges = list(edges)
        row_nodes = []
        col_edges = []
        data = []
        for ei, e in enumerate(edges):
            if factor is not None:
                k = factor[ei]
            else:
                k = 1
            ni = nodes.index(e[0])
            row_nodes.append(ni)
            col_edges.append(ei)
            data.append(1*k)
            ni = nodes.index(e[1])
            row_nodes.append(ni)
            col_edges.append(ei)
            data.append(-1*k)

        csr = coo_array((data, (row_nodes, col_edges)), shape=(len(nodes), len(edges))).tocsr()
        return csr

      
    def nodes_to_edges_path(inp, inverse = False):
        # TODO to je dela za Q in le po cudezu za X, kaj ce mi ne pokaze vseh poti
        if not inverse:
            nl = inp
            el = []
            for i in range(1,len(nl)):
                el.append((nl[i-1], nl[i]))
            return el
        else:
            el = np.array(inp)
            
            prvi_mask = ~np.isin(el[:,0], el[:,1]) 
            node = el[:,0][prvi_mask]
            
            nl = [int(node)]
            for i in range(len(el)):
                for e in el:
                    if node == e[0]:
                        node = e[1]
                        nl.append(node)
                        break
            
            return nl

    def edges_to_binary_vector(el, edges):
        ev = np.zeros(len(edges))
        
        for e in el:
            ei = edges.index(e)
            if ei != -1:
                ev[ei] = 1
            
        return(ev)


    def binary_vector_to_edges(ev, edges):
        return list(compress(edges, ev))
        
    def columns_to_paths(g,X,edges_mode=False):
        # TODO a splitam dva avta al ne?
        paths = []
        for i in range(X.shape[1]):
            if edges_mode:
                paths.append(Problem.binary_vector_to_edges(X[:,i],g.edges()))
            else:
                paths.append(Problem.nodes_to_edges_path(Problem.binary_vector_to_edges(X[:,i],g.edges()), inverse=True))
        return paths
    
    def my_nx_draw(g,paths,with_labels = False,with_nodes = True):
        pos = nx.spring_layout(g)
        if with_nodes:
            nx.draw_networkx_nodes(g,pos, node_size = 100)
        if with_labels:
            nx.draw_networkx_labels(g, pos)
        nx.draw_networkx_edges(g, pos, edge_color='k', width=0.5)
        for ci, p in enumerate(paths):
            nx.draw_networkx_edges(g, pos, edgelist=p, edge_color=COLORS[ci], width=2) # highlight elist
            
    def draw(p,method=None):
        to_draw = p.results if method is None else {method: p.results[method]}
        for k in to_draw:
            print(k + ":")
            result = p.results[k]
            if result["paths"] is None:
                pass
                print(result["message"],result["success"])
            else:
                print(result["fun"],result["message"],result["success"])
                if p.G is not None:
                    fig, ax = ox.plot_graph_routes(p.G,result["paths"],route_colors=list(COLORS)[:len(result["paths"])])
                    # for ci, path in enumerate(paths):
                    #     fig, ax = ox.plot_graph_route(graph,path,route_color=COLORS[ci])
                else:
                    paths_edges = Problem.columns_to_paths(p.g,result["X"],edges_mode=True)
                    Problem.my_nx_draw(p.g,paths_edges,with_labels=True,with_nodes=False)
                    plt.show()
    
    def print(p):
        for k in p.results:
            print(k + ":")
            print(p.results[k]["message"], p.results[k]["fun"], p.results[k]["success"])

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
            paths = Problem.columns_to_paths(p.g,X)
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
        
            mutated_pop.append(mutated_ex)
        
        return mutated_pop
    
    def crossover(pop,prob):
        crossovered_pop = []
        for i in range(int(np.ceil(len(pop)*prob))):
            ex1, ex2 = list(random.sample(pop,2))
            col_id = random.randint(0,ex1.shape[1]-1)
            v = ex1[col_id]
            ex1[col_id] = ex2[col_id]
            ex2[col_id] = v
            crossovered_pop += [ex1,ex2]
            
        return crossovered_pop
    
    
    
    def best(pop, f, num):
        sorted_pop = sorted(pop, key = lambda ex: f(ex))
        return sorted_pop[:min(num,len(pop))]
            
        
    def solve(p,path_seed=None,num_iter=10):
        
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
            
            mutated_pop = Lingen.mutate(p.g,pop,prob=0.5)
            # crossovered_pop = Lingen.crossover(pop+mutated_pop,prob=0.1)
            crossovered_pop = []
            best_pop = Lingen.best(pop+mutated_pop+crossovered_pop,f,num=10)
            pop = best_pop
            
        best_pop = Lingen.best(pop+mutated_pop+crossovered_pop,f,num=10)
        best_ex = best_pop[0]
        paths = Problem.columns_to_paths(p.g,best_ex)
        print("f_c_sums: ", f_c(best_ex), f_sums(best_ex))
        p.results[str(__class__)] = {"X": best_ex, "fun": f_t(best_ex), "message": "glej success", "success": f_c(best_ex) + f_sums(best_ex) == 0, "paths": paths}
        
        