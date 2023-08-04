
import networkx as nx


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
        
        assert all(e >= 0 for e in self.c)
        assert all(e > 0 for e in self.cap)
        
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
        def __init__(self, status, X, cost, message=""):
            self.status = status
            self.X = X # to je združene poti po komoditijih # sparse <3
            self.cost = cost
            self.message = message
            self.time = None
        
        def __str__(self):
            return self.__repr__()
        
        def __repr__(self):
            return "cost: {}, time: {}, status: {}, message: {}".format(self.cost, self.time, self.status, self.message)
