
import networkx as nx
import numpy as np
import helper_functions as hf


class Problem():
    verbose = False
    def __init__(self, graph, demands):
        self.graph = graph
    
        self.demands = demands
        
        self.c =  [self.graph.edges()[e]["c"] for e in self.graph.edges()]
        self.cap = [self.graph.edges()[e]["cap"] for e in self.graph.edges()]
        
        
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
        
        assert all(e >= 0 and isinstance(e, int) for e in self.c)
        # assert all(e > 0 for e in self.cap)
        assert all(e >= 0 for e in self.cap)
        
        self.results = {}
        
    def __str__(self):
        return self.results
    
    
    def print_results(self):
        for k in self.results:
            print(k + ":",end="")
            print(self.results[k])
            # print(self.results[k]["message"], self.results[k]["fun"], self.results[k]["success"])
            # if "time" in self.results[k]: print(self.results[k]["time"])
    
    def print(stri):
        if Problem.verbose:
            print(stri)
            
  
            
            
            
            
    class Result():
        def __init__(self, status, X, cost, message="",LB=None,over_cap_count=None):
            self.status = status
            self.X = X # to je združene poti po komoditijih # sparse <3
            self.cost = cost
            self.message = message
            self.time = None
            
            
            self.LB = LB
            self.over_cap_count = over_cap_count
            
        
        def __str__(self):
            return self.__repr__()
        
        def __repr__(self):
            stri = "cost: {}, time: {}, status: {}, message: {}".format(self.cost, self.time, self.status, self.message)
            if self.LB is not None: stri += ", LB: {}".format(self.LB)
            if self.over_cap_count is not None: stri += ", over_cap_count: {}".format(self.over_cap_count)
            return stri
        
        def latex_cell(self):
            # return hf.latex_cell([("cena", self.cost), ("rač. čas", round(self.time, 3)), ("status", self.status),
            #         ("sporočilo", self.message), ("sp. meja", self.LB), ("preko", self.over_cap_count)],with_field_names=False)
            return hf.latex_cell([self.cost, round(self.time, 3), self.status, self.message, self.LB, self.over_cap_count])
