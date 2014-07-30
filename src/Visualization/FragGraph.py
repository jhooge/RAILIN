'''
Created on Jan 5, 2012

@author: jhooge
'''

from matplotlib.pyplot import draw, ion
import networkx as nx
import pickle
import pylab as pl

def plot_fragment_graph(G, labels=None, to_file=None, dpi=300):
    
#    try:
#        from networkx import graphviz_layout
#    except ImportError:
#        raise ImportError("This example needs Graphviz and either PyGraphviz or Pydot")
#    
#    if len(G.nodes()) == 0:
#        raise ValueError('No fragments to plot!')
#    
#    pos = nx.graphviz_layout(G, prog="dot")
    fig = pl.figure(figsize=(10, 10))
    node_size = 300
    nx.draw_circular(G)
#    nx.draw(G, pos,
#            with_labels=True,
#            alpha=0.5,
#            node_size=node_size)
    
    ## adjust the plot limits
#    offset = node_size / 2
#    xmax = offset + max(xx for xx, yy in pos.values())
#    ymax = offset + max(yy for xx, yy in pos.values())
#    xmin = min(xx for xx, yy in pos.values()) - offset
#    ymin = min(yy for xx, yy in pos.values()) - offset
#    pl.xlim(xmin, xmax)
#    pl.ylim(ymin, ymax)
    
    if to_file != None:
        filename = to_file.split('/')[-1]
        pl.savefig((filename + "_graph"), dpi=dpi)
    else:
        pl.show()
        return fig
    
if __name__ == '__main__':
    from numpy import eye
    from networkx import from_numpy_matrix
    A = eye(10)
    G = from_numpy_matrix(A)
    plot_fragment_graph(G)
    pl.show()
    
    
