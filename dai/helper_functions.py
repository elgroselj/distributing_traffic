import cvxpy as cp
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
COLORS = "brgpy"

def plot_multigraph(graph, with_labels=True):
    G = nx.MultiDiGraph(graph)
    pos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=100, alpha=1)
    nx.draw_networkx_labels(G, pos)
    ax = plt.gca()
    for e in G.edges:
        edge_color = G.edges[e]["color"] if "color" in G.edges[e] else "k"
        ax.annotate("",
                    xy=pos[e[1]], xycoords='data',
                    xytext=pos[e[0]], textcoords='data',
                    arrowprops=dict(arrowstyle="->", color=edge_color,
                                    shrinkA=5, shrinkB=5,
                                    patchA=None, patchB=None,
                                    connectionstyle="arc3,rad=rrr".replace('rrr', str(0.1 * (e[2] + 1))
                                                                           ),
                                    ),
                    )
    if with_labels:
        digraph = nx.DiGraph(G)
        nx.draw_networkx_edge_labels(digraph, pos,
                                     edge_labels={e: "{},{}".format(digraph.edges[e]["c"], digraph.edges[e]["cap"])
                                                  for e in list(digraph.edges)})

    plt.axis('off')
    plt.show()

def plot_solution_graph(graph,X,with_labels=True):
    multi = nx.MultiDiGraph(graph)
    print("k\tCOLOR")
    for k in range(X.shape[1]):
        path = [e for i,e in enumerate(graph.edges) if X[i,k] != 0]
        multi.add_edges_from(path, color=COLORS[k])
        print(k,"\t",COLORS[k])
    # print(list(multi.edges()))
    plot_multigraph(multi,with_labels=with_labels)