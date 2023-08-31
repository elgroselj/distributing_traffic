from matplotlib import pyplot as plt
import generate_graphs as gg
import pandas as pd
import networkx as nx

def density_research():
    df = pd.DataFrame(columns=["n_max","t","k_num_max","density"])
    for n_max in [10,100,500,1000]:
        for t in [2,5,20,100]:
            if t < n_max:
                for k_num_max in [1,5]:
                    
                    graph,_ = gg.generate_random_graph(n_max=n_max,t=t,k_num_max=k_num_max,cap_max=2,density_param=1,c_max=10)
                    df.loc[len(df.index)] = [n_max,t,k_num_max, nx.density(graph)]
    print(df)
    print("mean: {}".format(df["density"].mean()))
    print("std: {}".format(df["density"].std()))
    
    for col in ["n_max","t","k_num_max"]:
        df_n_max = df[[col,"density"]].groupby(col).mean()
        df_n_max = df_n_max.reset_index()
        print(df_n_max)
        plt.plot(df_n_max[col],df_n_max["density"])
    plt.show()
density_research()