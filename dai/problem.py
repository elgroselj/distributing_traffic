
import networkx as nx
import numpy as np
import helper_functions as hf


class Problem():
    verbose = False
    status_map = {hf.Status.BLANK:"-",hf.Status.OPTIMAL:"opt",hf.Status.FEASIBLE:"dop",hf.Status.INFEASIBLE:"nedop",hf.Status.OVER_CAP:"preko"}
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
            
            self.percent_of_min = None
            
        
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
            style_int = lambda x: int(x) if (isinstance(x,int) or isinstance(x,float)) and x is not None else "-"
            style_float = lambda x: round(x, 3) if (isinstance(x,int) or isinstance(x,float)) and x is not None else "-"
            return hf.latex_cell([style_int(self.cost), style_float(self.percent_of_min), round(self.time, 3), Problem.status_map[self.status], style_int(self.LB), style_int(self.over_cap_count)])
    
    class Result_grouped():
        def __init__(self, status_count, mean_cost, std_cost, mean_time, std_time, mean_LB, mean_over_cap, mean_percent_of_min, std_percent_of_min): # modus_status,LB=None,over_cap_count=None):
            self.status_count = status_count
            self.mean_cost = mean_cost
            self.std_cost = std_cost
            
            self.mean_time = mean_time
            self.std_time = std_time
            
            
            self.mean_LB = mean_LB
            self.mean_over_cap = mean_over_cap
            
            self.mean_percent_of_min = mean_percent_of_min
            self.std_percent_of_min = std_percent_of_min
            
        
        def __str__(self):
            return self.__repr__()
        
        # def __repr__(self):
        #     stri = "cost: {}, time: {}, status: {}, message: {}".format(self.cost, self.time, self.status, self.message)
        #     if self.LB is not None: stri += ", LB: {}".format(self.LB)
        #     if self.over_cap_count is not None: stri += ", over_cap_count: {}".format(self.over_cap_count)
        #     return stri
        
        def latex_cell(self):
            # return hf.latex_cell([("cena", self.cost), ("rač. čas", round(self.time, 3)), ("status", self.status),
            #         ("sporočilo", self.message), ("sp. meja", self.LB), ("preko", self.over_cap_count)],with_field_names=False)
            return hf.latex_cell([self.mean_cost, self.mean_percent_of_min, self.std_percent_of_min, self.mean_time, self.std_time,
                                  str(list(self.status_count.values())).strip(r"[|]"), self.mean_LB, self.mean_over_cap,
                                  ])
