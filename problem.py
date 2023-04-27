from itertools import compress
import networkx as nx
import random
from scipy.optimize import linprog
import numpy as np
from scipy.sparse import coo_array
import osmnx as ox
import matplotlib.pyplot as plt
COLORS="brgymcbrgymc"

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

        coo = coo_array((data, (row_nodes, col_edges)), shape=(len(nodes), len(edges)))
        return coo

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
        B = coo_array((data, (row_edge, col_ZK)), shape=(num_edges, num_ZK*num_edges))
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
        coo = coo_array((data, (rows, cols)), shape=(num_ZK*num_nodes, num_ZK*num_edges))
        return coo
    
    def nastavi_fbcmm(g,ZK,c,t):
        f = np.repeat(t,len(ZK)) #min f*x
        
        B = Problem.getB(g.number_of_edges(),len(ZK)) # B*x <= c
        
        M = Problem.sparse_incidence_matrix(g.nodes(),g.edges()) # (M*X=M_ZK)
        
        a = [a for _,_,a in ZK]
        M_ZK = Problem.sparse_incidence_matrix(g.nodes(),[(z,k) for z, k, _ in ZK],factor=a)
        
        M2 = Problem.MtoM2(M.todense(), len(ZK)) # M2 * x = m_ZK
        m_ZK = M_ZK.todense().flatten('F')
        
        return (f,B,c,M2,m_ZK)
    
    def nodes_to_edges_path(inp, inverse = False):
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

    # def edge_to_incidence_vector(e, nodes):
    #     nv = np.zeros(len(nodes))
        
    #     ni = nodes.index(e[0])
    #     nv[ni] = 1
    #     ni = nodes.index(e[1])
    #     nv[ni] = -1
            
    #     return(nv)

    def binary_vector_to_edges(ev, edges):
        return list(compress(edges, ev))
        
    def columns_to_paths(g,X,edges_mode=False):
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
    
    def solve(self,mode,silent = True):
        if mode == "linprog_v1":
            
#                 c @ x
#                    such that:

            #     A_ub @ x <= b_ub
            #     A_eq @ x == b_eq
            #     lb <= x <= ub
            
            f,B,c,M2,m_ZK = Problem.nastavi_fbcmm(self.g,self.ZK,self.c,self.t)
            B = B.todense()
            res = linprog(f, A_ub=B, b_ub=c, A_eq=M2, b_eq=m_ZK, integrality=1)
            
            X = None
            if not silent:
                if res.x is None:
                    print(res.message)
                else:
                    print(res.x,res.fun,res.message)
                    X = np.array(res.x).reshape(self.g.number_of_edges(),len(self.ZK))
                    print(X)
                    # path0 = hf.nodes_to_edges_path(hf.binary_vector_to_edges(X[:,0],g.edges()), inverse=True)
                    # fig, ax = ox.plot_graph_route(graph,path0)
                    if self.G is not None:
                        paths = Problem.columns_to_paths(self.g,X)
                        fig, ax = ox.plot_graph_routes(self.G,paths,route_colors=list(COLORS)[:len(paths)])
                        # for ci, path in enumerate(paths):
                        #     fig, ax = ox.plot_graph_route(graph,path,route_color=COLORS[ci])
                    else:
                        paths_edges = Problem.columns_to_paths(self.g,X,edges_mode=True)
                        Problem.my_nx_draw(self.g,paths_edges,with_labels=True,with_nodes=False)
            self.results[mode] = {"X": X, "fun": res.fun, "message": res.message, "success": res.x is not None}
            
    