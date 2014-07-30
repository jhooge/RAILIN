'''
Created on Jul 11, 2012

@author: jhooge
'''

import pickle
import networkx as nx
import pylab as pl
from Pasta import Pasta
from networkx.algorithms import bipartite

if __name__ == '__main__':
    
    fn = "/is/ei/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/"\
         "src/Classification/Datasets/Ubiquitin/Results/"\
         "Ubq_Bax_Robustness_Tests_Max/NoisyData/"\
         "UbqBaxGaussPDFNoise=0.0.p"
         
    pasta_obj = pickle.load(open(fn,"r"))
    
    try:
        from networkx import graphviz_layout
    except ImportError:
        raise ImportError("This example needs Graphviz and either PyGraphviz or Pydot")
    
    B = nx.Graph()
    B.add_nodes_from([1,2,3,4], bipartite=0) # Add the node attribute "bipartite"
    B.add_nodes_from(['a','b','c'], bipartite=1)
    B.add_edges_from([(1,'a'), (1,'b'), (2,'b'), (2,'c'), (3,'c'), (4,'a')])
    bottom_nodes, top_nodes = bipartite.sets(B)
    
    pos = nx.graphviz_layout(B, prog="dot")
    fig = pl.figure(figsize=(10, 10))
    node_size = 300
    nx.draw(B, pos,
            with_labels=True,
            alpha=0.5,
            node_size=node_size)
    
    
    pl.show()