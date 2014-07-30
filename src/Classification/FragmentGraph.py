'''
Created on Feb 17, 2011

@author: jhooge
'''

#def find_paths(start, graph, path=[]):
#    path = path + [start]
#    paths = []
#
#    for node in graph.neighbors(start):
#        if not node in path:
#            newpaths = find_paths(node, graph, path)
#            
#            for newpath in newpaths:
#                paths.append(newpath)
#        else:
#            paths.append(path)
#
#    if len(graph.neighbors(start)) == 0:
#        if len(path) > 1:
#            paths.append(path)
#    return paths

#def filter_fragment_list(fragments, min_length=2):
#
#    newlist = []
#    for fragment in fragments:
#        if len(fragment) >= min_length and fragment not in newlist:
#            newlist.append(fragment)
#    return newlist

#def residues_to_edges(residues):
#    all_edges = set([])
#    for r in residues:
#        curr_node = r
##        curr_node = residues.index(r)
##        curr_node = customIndex(r,residues)
#        successors = [p for p in r.succ]
##        successors = [residues.index(p) for p in r.succ]
##        successors = [customIndex(p,residues) for p in r.succ]
#        edges = [edge for edge in itertools.product([curr_node], successors)]
#        
#        all_edges = all_edges.union(set(edges))
#    all_edges = list(all_edges)
#    return all_edges

class Fragment(object):
    '''
    classdocs
    '''
#    def __init__(self, residues=[], seq_pos=[]):
#        self.residues = residues
#        self.errors = [[]]
    
    def __init__(self, residues=[]):
        self.residues = residues
        self.errors = [[]]
        
    def __str__(self):
        """Represents Fragment as string"""
        return '->'.join(map(str, self.residues))

    def __cmp__(self, other):
        """
        Compares two fragments
        Two fragments are considered equal, when the set of residues of 
        both are the same.
        """
        return cmp(self.residues, other.residues)

    def __len__(self):
        """
        Returns the length of a fragment.
        The length is defined as the number of residues it consists of.
        """
        return len(self.residues)
    
#    def remove_node(self, node):
#        """
#        Removes a node (residue) from fragment.
#        A node is considered as being removed, when its value is None
#        
#        @param node: A residue in fragment
#        @type node: object
#        """
#        while self.residues.count(node):
#            self.residues[self.residues.index(node)] = None
#        
#    def remove_nodes(self, nodes):
#        """
#        Removes a set of nodes from fragment
#        
#        @param nodes: A set of nodes to be removed from fragment
#        @type nodes: list
#        """
#        map(self.remove_node, nodes)
            
#    def split(self):
#        paths = [[]]
#
#        for node in self.residues:
#            if node is None:
#                paths.append([])
#            else:
#                paths[-1].append(node)
#
#        return filter_fragment_list([Fragment(path) for path in paths], 1)
#
#        
#    def ints_to_residues(self, residues):
#        '''
#        Converts residue list from residue indexes in residues
#        to Residue Objects
#        
#        MapToSeq runs on Residue objects
#        whereas
#        Fragment Splitting runs on integer residue lists
#        '''
#        residue_indexes = self.residues
#        self.residues = [residues[i] for i in residue_indexes]
#    
#    def residues_to_ints(self, residues):
#        '''
#        Converts residue list from Residue object to residue index
#        
#        MapToSeq runs on Residue objects
#        whereas
#        Fragment Splitting runs on integer residue lists
#        
#        '''
#        self.residues = [residues.index(residue) for residue in self.residues]
#    
#    def to_array(self, boolean=False):
#        
#        if boolean:
#            return np.array([r.profile for r in self.residues], bool)
#        else:
#            return np.array([r.profile for r in self.residues], float)
#    
#    def compare_to_seq(self, seq):
#        """
#        assumes len(frag1) >= len(frag2)
#        computes error score for all possible sequence positions
#        
#        m := sequence length
#        n := fragment length
#        p_i := aa_rec probability for amino acid at sequence position i
#        
#        error = sum(i=1,n) 1-p_i
#        """
#        
#        m = len(seq.residues)
#        n = len(self.residues)
#        assert  m >= n
#    
#        p = seq.to_array().astype(bool)
#        p = np.delete(p, 1, 1)
#        q = self.to_array()
#        
#        error = []
#        
#        for i in range(0, m - n + 1): # m length seq, n length fragment
#            error.append(np.abs(1. - q[p[i:i + n]]))
#        while len(error) < m:
#            error.append(np.ones(n) * 1.79769313e+308)
#            
#        self.errors = np.array(error)
#        
##    def fragment_mask(self, mask):
##        mask = deepcopy(mask)
##        m = len(self)
##        
##        i = m - 1
##        while i < len(mask):
##            if(mask[i] > 1.0e+200):
##                j = i - m + 1
##                mask[j:i] = 1.79769313e+308
##                i += m
##            else:
##                i += 1
##        return mask
#    
#    def get_min_error_pos(self, mask):
#        mask = self.fragment_mask(mask)
#        mask = mask.reshape(len(mask), 1)
#        masked_errors = np.multiply(mask, self.errors)
#        self.errors = masked_errors
#        i_min = np.argmin(masked_errors.sum(1))
#        
#        return i_min
#    
#    def update_mask(self, mask):
#        i = self.get_min_error_pos(mask)
#        mask[i:i + len(self)] = np.nan_to_num([np.inf])
#        
#        return mask

import networkx as nx
class FragmentGraph(nx.DiGraph):
    '''
    classdocs
    '''
    
#    def find_all_paths(self):
#        """
#        Returns all paths from nodes in node set of this graph
#        to each sink
#        
#        
#        """
#        paths = []
#
#        for node in self.nodes():
#            
#            paths += self.find_paths(node)
#
#        paths = set(map(tuple, paths))
#
#        return list(paths)
#    
#    def compute_fragment_pool(self):
#        paths = self.find_all_paths()
#        return FragmentPool(paths)
    
    '''NEW STUFF FOR RND_MAPPING'''
    
    def __get_sources(self):
        """
        Gets all sources, that is nodes with no incoming edges
        (indegree == 0)
        
        @return: list of nodes
        @rtype: list
        """
        
        sources = []
        for node in self.nodes():
            if self.in_degree(node) == 0:
                sources.append(node)
        return sources
    
    def __get_sinks(self):
        """
        Gets all sinks, that is nodes with no outgoing edges
        (outdegree == 0)
        
        @return: list of nodes
        @rtype: list
        """
        sinks = []
        for node in self.nodes():
            if self.out_degree(node) == 0:
                sinks.append(node)
        return sinks
    
    def __is_single_node(self, node):
        """
        Checks whether a node is a single node
        
        @param node: a node in networkx graph
        @type node: object
        
        @return: True if node has no incoming and no outgoing edges
        @rtype: bool
        """
        y = self.in_degree(node) == 0 and self.out_degree(node) == 0
        return y
    
    def __get_single_nodes(self):
        """
        Gets all nodes which have no incoming nor outgoing edges
        
        @return: single nodes
        @rtype: list of nodes
        """
        singles = [node for node in self.nodes() if self.__is_single_node(node)]
        return singles
    
    def find_paths(self, start, path=[]):
        """
        Find all paths from start node to every node that is a sink.
        
        @param start: node in graph
        @type start: object
        @param path: list of nodes on path
        @type path: list
        """
        path = path + [start]
        paths = []
    
        for node in self.neighbors(start):
            if not node in path:
                newpaths = self.find_paths(node, path)
                
                for newpath in newpaths:
                    paths.append(newpath)
            else:
                paths.append(path)
    
        if len(self.neighbors(start)) == 0:
            if len(path) > 1:
                paths.append(path)
        return paths
    
    def get_random_path(self, verbose=False):
        from numpy import random, array
        from Math import sample_from_histogram
        """
        Returns random path from random source to a sink in graph
        Samples from histogram that longest paths are chosen more 
        frequently.
        
        @param verbose: optional verbosity option
        @type verbose: bool
        
        @return: list of nodes on path
        @rtype: list
        """
        
        sources = self.__get_sources()
        if len(sources) == 0: ## if cycle then choose a random node
            for node in self.nodes():
                if self.out_degree(node) > 0:
                    sources.append(node)
                    
        if verbose:
            print 'Number Sources in Graph: ', len(sources), sources
            
        i = random.randint(len(sources))
        source = sources[i]
        if verbose:
            print 'Source: ', source
            
        paths = self.find_paths(source)
        if verbose:
            print 'Number paths: ', len(paths)
            
        p = array(map(float, map(len, paths)))
        p = p / sum(p)
        p = list(p)
        i = sample_from_histogram(p)
        path = paths[i]
        
        return path

    def remove_path_from(self, path):
        """
        Removes path from graph. If after the path has been removed there
        are single nodes in graph, they will be removed as well.
        
        @param path: list of nodes
        @type path: list
        """
        self.remove_nodes_from(path)
        self.remove_single_nodes_from()
    
    def remove_single_nodes_from(self):
        """Removes single nodes from graph."""
        for node in self.nodes():
            if self.__is_single_node(node):
                self.remove_node(node)
        
    def get_invalid_edges(self, invalid_nodes):
        from operator import xor
        """
        Takes a list of invalid nodes and removes all edges in graph
        that either end or (exclusive) start at an invalid node.
        A node is invalid if the posterior probability for a node at
        its assigned sequence position is below a certain threshold.
        
        @param invalid_nodes: list of nodes with posterior below threshold
        @type invalid_nodes: list
        
        @return: list of tuples representing an edge in graph
        @rtype: list
        """
        invalid_edges = []
        for edge in self.edges():
            e1, e2 = edge
            a = e1 in invalid_nodes
            b = e2 in invalid_nodes
            if xor(a, b):
                invalid_edges.append(edge)
        
        return invalid_edges
    
    def adj_matrix(self):
        from numpy import array
        """Returns this graph as adjacency matrix"""
        return array(nx.adj_matrix(self))
    
    def adj_list(self):
        """Returns this graph as adjacency list"""
        return self.adjacency_list()
    
#class FragmentPool(object):
#    
#    def __init__(self, paths):
#        '''
#        Constructor
#        
#        This class represents a list of fragment used for Sequence Mapping.
#        Every time a Fragment is mapped to the template sequence residues 
#        included in the fragment
#        '''
#        self.fragments = filter_fragment_list([Fragment(path) for path in paths]) # list of fragment objects
#        self.sort()
#    
#    def update(self, mapped_fragment, residues, invalidNodes=[]):
#        
#        '''
#        Sets mapped residues to None in all Fragments.
#        Splits Fragments in Subfragments.
#        
#        mapped_fragment = [4,5]
#        
#        [1,2,3,4,5] -> [1,2,3,None,None] -> [1,2,3]
#        [4,5,6,7,8] -> [None,None,6,7,8] -> [6,7,8]
#        [1,2,3,4,5,6,7,8] -> [1,2,3,None,None,6,7,8] -> [1,2,3], [6,7,8]
#        '''
#        if len(invalidNodes) > 0:
#            mapped_fragment.residues_to_ints(residues)
#            mapped_fragment.remove_nodes(invalidNodes)
#            split_frags = filter_fragment_list(mapped_fragment.split())
##            print 'Split to '
#            for f in split_frags:
#                f.ints_to_residues(residues)
##                print ''.join([r.name for r in f.residues])
#                self.add(f)
#            self.remove(mapped_fragment)
#            
#        else:
#            mapped = deepcopy(mapped_fragment)
#            
#            [frag.remove_nodes(mapped.residues) for frag in self.fragments]
#            newfrags = []
#            
#            for frag in self.fragments:
#                newfrags += frag.split()
#    
#            newfrags = filter_fragment_list(newfrags)
#            self.fragments = newfrags
#            self.sort()
#        
#    def remove(self, frag):
#        self.fragments.remove(frag)
#        self.sort()
#        
#    def add(self, frag):
#        self.fragments.append(frag)
#        self.sort()
#        
#    def sort(self):
#        '''
#        Sorts fragment pool by length (descending).
#        '''
#        self.fragments = sorted(self.fragments,
#                                key=lambda frag: len(frag),
#                                reverse=True)

if __name__ == '__main__':
    import pickle
    import pylab as pl
    
#    G = pickle.load(open('/kyb/agbs/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/multiple_test_files/31_NTD_fragment_graph_new.p'))
#    fig = plot_fragment_graph(G)
    
#    G1 = FragmentGraph()
#    G2 = FragmentGraph()
#    G3 = FragmentGraph()
#    G4 = FragmentGraph()
#    G5 = FragmentGraph()
#    G6 = FragmentGraph()
#    
#    G1.add_edges_from([(1, 2), (2, 3)])
#    G2.add_edges_from([(1, 2), (2, 3), (2, 4)])
#    G3.add_edges_from([(1, 2), (2, 1)])
#    G4.add_edges_from([(1, 2), (1, 3), (2, 4), (3, 4)])
#    G5.add_edges_from([(1, 1), (1, 2), (1, 2)])
#    G6.add_edges_from([(1, 2), (2, 3), (3, 4), (3, 5), (5, 6)])
#    
#    print type(G6.adj_matrix())
#    print G6.adj_matrix()
#    print type(G6.adjacency_list())
#    print G6.adjacency_list()
    
    
#    print G1.find_all_paths()
#    print G2.find_all_paths()
#    print G3.find_all_paths()
#    print G4.find_all_paths()
#    print G5.find_all_paths()
#    print G6.find_all_paths()
    
#    G6.remove_edges_from([(1, 2),(5,6)])
#    print G6.__get_single_nodes()
#    plot_fragment_graph(G6)
#    pl.show()
    
#    print G6.nodes()

#    while len(G.edges()) != 0:
#        print 'G(|V|,|E|) = G(%i,%i)' % (len(G.nodes()), len(G.edges()))
#        path = G.get_random_path()
##        plot_fragment_graph(G)
#        G.remove_path_from(path)
#        print 'Chosen Path: ', path
