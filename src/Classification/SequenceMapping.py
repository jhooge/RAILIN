'''
Created on Mar 6, 2012

@author: jhooge
'''
"""
This Module holds different algorithms for Sequence Mapping.
It is designed as Strategy Pattern (see: GoF)
"""

from abc import ABCMeta, abstractmethod
import pickle

class SequenceMapping(object):
    """
    The common interface of all sequence mapping strategies
    @requires: Typing and Linking has to be run on PASTA object, before
    Sequence Mapping can be run
    """
    __metaclass__ = ABCMeta
    
    def __repr__(self):
        return '{0}:{1}'.format(self.__class__.__name__, id(self))   
    
    @abstractmethod
    def generate_assignments(self, context):
        raise NotImplementedError('This is an abstract method !')
    
    def get_path(self, G, verbose=False):
        """
        Returns random path on FragmentGraph
        
        @param G: Fragment graph generated from linking matrix
        @type G: FragmentGraph
        @param verbose: verbosity option
        @type verbose: bool
        
        @return: list of residue indices on path
        @rtype: list
        """
        
        rnd_path = G.get_random_path(verbose)
        return rnd_path
    
    def P_from(self, pasta_obj):
        """
        Returns mxnxn posterior matrix from PASTA object, 
        where m is the number of residues and n is the number of
        amino acid labels
        
        @param pasta_obj: Pasta instance filled with typing and linking 
        information
        @type pasta_obj: Pasta
        @return: posterior matrix from pasta instance
        @rtype: numpy.ndarray
        """
        P = pasta_obj.P
        return P
    
    def S_from(self, pasta_obj):
        """
        Returns boolean mxn profile matrix for amino acid sequence,
        where m is the number of amino acids in template sequence and n
        is the number of amino acid class labels
        
        Rows in the matrix correspond to amino acid in sequence
        and columns correspond to an amino acid class label.
        S[i][j] = 1 if amino acid at sequence position i is equal to the
        class label j
        S[i][j] = 0, else
        
        @param pasta_obj: Instance of Pasta class
        @type pasta_obj: Pasta
        
        @return: Sequence profile matrix
        @rtype: numpy.ndarray
        """
        from numpy import zeros
        from copy import deepcopy
        
        seq = pasta_obj.seq
        aas = deepcopy(pasta_obj.amino_acids)
        
        # # remove CY* from amino acids
        aa = [a.three_let for a in pasta_obj.amino_acids]
        k = aa.index('CY*')
        del aas[k]
        
        S = zeros((len(seq), len(aas)))
        for i, s in enumerate(seq):
            for j, a in enumerate(aas):
                if s.one_let == a.one_let:
                    S[i][j] = 1
        
        return S
    
    def F_from(self, P, path):
        """
        Returns posterior probabilities in P for each residue on path as
        numpy.ndarray. In other words this method returns a slice of an 
        array of arbitrary dimensions by selecting indices in P provided
        by path parameter. 
        
        @attention: values in P are selected over axis 0, such that the 
        dimension of the returned array is the same than P.
        
        Example:
        path = [0,1]
        P.shape = (4,4,4)
        F_from(P, path).shape == (2,4,4)
        
        path = [0,1]
        P.shape = (4,4)
        F_from(P, path).shape == (2,4)
        
        @param P: Posterior Matrix
        @type P: numpy.ndarray
        @param path: index array with residue indices
        @type path: list
        
        @return: Posterior subarray only including posteriors for
        residue index on path
        @rtype: numpy.ndarray
        """
        
        from numpy import take
        F = take(P, path, axis=0)
        return F
    
    def frag_posteriors(self, F, S, i):
        """
        Computes posterior vector for residue 
        [i],[i-1]-joint posteriors at sequence position
        i to i+n, where i is sequence position where first fragment
        residue is mapped to and n is the fragment length.
        
        @todo: this method is only being used in SinglePosterior and might 
        as well be defined in this class
        @note: I'm unsure whether a multiplication here is correct 
        
        @param F: Fragment Matrix (len(F),len(amino_acids),len(amino_acids))
        @type F: array
        @param S: Sequence Matrix (len(S),len(amino_acids),len(amino_acids))
        @type S: array 
        @param i: Sequence position index
        @type i: int
        
        @return: 3d slice of fragment posteriors (len(F),len(amino_acids),len(amino_acids))
        @rtype: ndarray
        """
        from numpy import multiply
        
        n = F.shape[0]  # # Fragment length
        T = multiply(F, S[i:i + n])
        return T
   
class SinglePosterior(SequenceMapping):  # a concrete instrument
    
    def __init__(self, pasta_obj, no_assignments=100, verbose=False):
        from Linking import graph
        self.pasta = pasta_obj
        self.P = self.P_from(pasta_obj)  # # numpy array with posteriors
        self.L = pasta_obj.L  # # numpy array with linking
        self.S = self.S_from(pasta_obj)  # # list of amino acids in seq
        self.verbose = verbose
        self.fragment_graph = graph(self.L)
        self.no_assignments = no_assignments
    
    def P_from(self, pasta_obj):
        P = super(SinglePosterior, self).P_from(pasta_obj)
        P = P.sum(1)
        return P
    
    def generate_assignment(self):
        from Visualization.FragGraph import plot_fragment_graph
        from FragmentGraph import Fragment
        from Definitions import three2One
        from copy import deepcopy
        from numpy import ones
        
        P = self.P
        S = self.S
        m = ones(len(S))
        
        G = deepcopy(self.fragment_graph)
#        plot_fragment_graph(G)
        G.remove_single_nodes_from()
        assignment = []
        while len(G.edges()) != 0:
            if self.verbose:
                print 
                print 'G(|V|,|E|) = G(%i,%i)' % (len(G.nodes()), len(G.edges()))
                print 'Edges: ', G.edges()
            path = self.get_path(G, self.verbose)
            print "Got random path"
            F = self.F_from(P, path)
            i = self.max_pos(F, S, m)
#            print 'maximum: ', i
#            print 'mask at [%i:%i]: %s' % (i, i + len(F) - 1, reshape(m[i:i + len(F)], len(m[i:i + len(F)])))

            if not self.is_valid(F, S, m, i):
#                print "Invalid Mapping!"
                inval_nodes = self.invalid_nodes(F, S, m, i, path)
                if self.all_invalid(F, S, m, i, path):
                    G.remove_nodes_from(inval_nodes)
                    G.remove_single_nodes_from()
#                    print "All Residues Invalid!"
#                    print "Removing nodes: "
#                    print inval_nodes
                else:
                    inval_edges = G.get_invalid_edges(inval_nodes)
                    G.remove_edges_from(inval_edges)
                    G.remove_single_nodes_from()
#                    print "Removing edges:"
#                    print inval_edges
#                    
#                if len(G.nodes()) !=0:
#                    plot_fragment_graph(G, to_file='GeneratedDatasets/raw_data_%i'%c1)
            else:
                m = self.update_mask(F, m, i)
                
                for res_no in path:
                    assignment.append((res_no, i))
                    i += 1
#                assignment.append((Fragment(residues=path), i))
                G.remove_path_from(path)
                G.remove_single_nodes_from()
#                print "Valid Mapping!"
#                print "Removing Mapped Nodes:"
#                print path
#                print
                
            if self.verbose:
                posteriors = self.frag_posteriors(F, S, i).sum(1)
                print 'Posteriors:      ', posteriors
                print 'Valid?       ', self.is_valid(F, S, m, i, t=1. / 20)
                print 'Edges:       ', G.edges()
                
#        for frag, i in assignment:
#            frag.ints_to_residues(self.pasta.residues)
#            print i,''.join([three2One(r.name) for r in frag.residues])
#            frag.residues_to_ints(self.pasta.residues)
        return assignment
    
    def generate_assignments(self):
        assignments = []
        P = self.P
        S = self.S
        for i in range(0, self.no_assignments):
            conf = self.generate_assignment()
#            similarity = self.seq_similarity(conf, P, S)
            similarity = self.seq_similarity(conf, P)
            print 'Assignment %i %.4f' % (i, similarity)
            assignments.append((conf, similarity))
#        assignments = sorted(assignments, key=lambda conf: conf[1], reverse=True)
        return assignments
    
    def invalid_nodes(self, F, S, m, i, path, t=1. / 20):
        """
        @param F: Fragment Matrix
        @type F: numpy.ndarray
        @param S: Sequence Matrix
        @type S: numpy.ndarray
        @param i: Sequence position index
        @type i: int
        @param t: Minimal posterior for amino acid
        @type t: float
        
        @return: list of invalid nodes
        @rtype: list
        """
        from numpy import reshape
        
        posteriors = list(self.frag_posteriors(F, S, i).sum(1))
        inval = []
        local_mask = reshape(m[i:i + len(F)], len(m[i:i + len(F)]))
        j = 0
        for p in posteriors:
            if p <= t or local_mask[j] == 0:
                inval.append(path[j])
            j += 1
        return inval
    
    def all_invalid(self, F, S, m, i, path, t=1. / 20):
        return len(self.invalid_nodes(F, S, m, i, path, t)) == len(F)
    
    def is_valid(self, F, S, m, i, t=1. / 20):
        """
        A Mapping is valid if all posteriors are above a predefined threshold t.
        
        @param F: Fragment Matrix
        @param S: Sequence Matrix
        @param i: Sequence position index
        @param t: Minimal posterior for amino acid
        
        @return: Bool
        """
        from numpy import reshape
        
        
        T = self.frag_posteriors(F, S, i)
        local_mask = reshape(m[i:i + len(F)], len(m[i:i + len(F)])).sum()
        return (T.sum(1) > t).all() and local_mask == len(m[i:i + len(F)])
    
#    def frag_posteriors(self, F, S, i):
#        """
#        @param F: Fragment Matrix
#        @param S: Sequence Matrix
#        @param i: Sequence position index
#        
#        @return: 
#        """
#        from numpy import multiply 
#        
#    #    print 'F[%i:%i]'%(i,i+len(F)-1)
# #        print F
#    #    print 'S'
#    #    print S[i:i + len(F)]
#    #    print 'T = F*S'
#        T = multiply(F, S[i:i + len(F)])
# #        print T
#        return T
    
    def max_pos(self, F, S, m):
        """
        @param F: Fragment Matrix
        @param S: Sequence Matrix
        @param m: Mask column vector
        """
        from numpy import argmax, multiply, array
        
        k = []

        i = 0
        j = i + len(F)
        while j <= len(S):
            S = multiply(S.T, m).T
            T = self.frag_posteriors(F, S, i)
            k.append(T.sum())
            i += 1
            j += 1
        return argmax(array(k))

    def update_mask(self, F, M, i):
        """
        @param F: Fragment Matrix
        @param M: Mask Matrix
        @param i: Sequence position index
        """ 
        from numpy import zeros
        
        k = F.shape[0]  # # len fragment
        M[i:i + k] = zeros(k)
        return M
    
#    def seq_similarity(self, conf, P, S):
#        """
#        @param P: Posterior Matrix 
#        @param S: Sequence Matrix
#        """
#        from numpy import array
#        
#        errors = []
#        for frag, i in conf:
#            F = self.F_from(P, frag.residues)
#            T = self.frag_posteriors(F, S, i)
#            error = 1. - T.sum(1)
#            errors.append(error.sum())
#        
#        errors = array(errors)
#        m = len(S)
#        n = array([len(frag) for frag, i in conf]).sum()
#        gaps = m - n
#        similarity = 1. - (1. / m * (errors.sum() + gaps))
#        return similarity

    def seq_similarity(self, assignment, P):
        """
        @param assignment: Residue to Sequence assignment
        @param P: Posterior Matrix 
        """
        from numpy import array
        posteriors = []
        amino_acids = map(lambda x: x.one_let, self.pasta.amino_acids)
        del amino_acids[1]
        sequence = map(lambda x: x.one_let, self.pasta.seq)
        
        for res_no, seq_no in assignment:
            aa_no = amino_acids.index(sequence[seq_no])
            posteriors.append(P[res_no][aa_no])
        errors = (1 - array(posteriors)).sum()
        m = len(sequence)
        n = len(assignment)
        gaps = m - n
        similarity = 1. - (1. / m * (errors + gaps))
        return similarity
    
class JointPosteriorNew(SequenceMapping):  # a concrete instrument
    
    def __init__(self, pasta_obj, no_assignments=100, verbose=False):
        from Linking import graph
        self.pasta = pasta_obj
        self.P = self.P_from(pasta_obj)
        self.L = pasta_obj.L
        self.verbose = verbose
        self.fragment_graph = graph(self.L)
        self.no_assignments = no_assignments
    
    def S_from(self, pasta_obj):
        from Math import outer_2way
        from numpy import concatenate, ones, multiply, swapaxes
        import pylab as pl
        
        s = super(JointPosteriorNew, self).S_from(pasta_obj)
        s = concatenate([ones((1, s.shape[1])), s], 0)
        S = outer_2way(s[:-1].T, s[1:].T, multiply)
        S = swapaxes(S, 0, 1)
        return S
    
    def P_from(self, pasta_obj):
        from numpy import dot, reshape
        import pylab as pl
        
        P = pasta_obj.P  # # posterior matrix from aa recog
        n = len(pasta_obj.seq)  # # len sequence
        m = len(pasta_obj.residues)  # # len pasta list
        S = self.S_from(pasta_obj)
        P = dot(reshape(P, (m, -1)), reshape(S, (n, -1)).T)
        return P

    def generate_assignment(self):
        from FragmentGraph import Fragment
        from Definitions import three2One
        from copy import deepcopy
        from numpy import ones, reshape
        from Visualization.FragGraph import plot_fragment_graph
        
        P = self.P
        m = ones(P.shape[1])
        
        G = deepcopy(self.fragment_graph)
        G.remove_single_nodes_from()
        plot_fragment_graph(G)
        assignment = []
#        print '*' * 40
        while len(G.edges()) != 0:
            if self.verbose:
                print 
                print 'G(|V|,|E|) = G(%i,%i)' % (len(G.nodes()), len(G.edges()))
                print 'Edges: ', G.edges()
                
            path = self.get_path(G, self.verbose)
            F = self.F_from(P, path)
            i = self.max_pos(F, m)
#            print 'maximum: ', i
#            print 'mask at [%i:%i]: %s' % (i, i + len(F) - 1, reshape(m[i:i + len(F)], len(m[i:i + len(F)])))
            
            if not self.is_valid(F, m, i):
#                print 
#                print "Invalid Mapping!"
                frag_aas = ''.join(three2One(self.pasta.residues[r].name) for r in path)
                seq_aas = ''.join(self.pasta.seq[r].one_let for r in range(i, i + len(path)))
#                print 'Frag %s %s' % (frag_aas, path)
#                print 'Seq  %s' % (seq_aas)
#                print 'Posteriors: '
#                print self.frag_posteriors(F, m, i)
                inval_nodes = self.invalid_nodes(F, m, i, path)
#                print 'Invalid Nodes: %s' % inval_nodes
                if self.all_invalid(F, m, i, path):
                    G.remove_nodes_from(inval_nodes)
                    G.remove_single_nodes_from()
#                    print "All Residues Invalid!"
#                    print "Removing all invalid nodes from graph"
                else:
                    inval_edges = G.get_invalid_edges(inval_nodes)
                    G.remove_edges_from(inval_edges)
                    G.remove_single_nodes_from()
#                    print "Removing edges in graph"
#                    print inval_edges
                    
#                if len(G.nodes()) !=0:
#                    plot_fragment_graph(G, to_file='GeneratedDatasets/raw_data_%i'%c1)
            else:
#                frag_aas = ''.join([three2One(self.pasta.residues[i].name) for i in path])
#                print 'mask before: %s'%(''.join(map(str,m.astype(int))))
                m = self.update_mask(F, m, i)
#                print 'mask after:  %s %s %i-%i'%(''.join(map(str,m.astype(int))),frag_aas,i,i+len(frag_aas)-1)
#                print reshape(m, len(m))
                assignment.append((Fragment(residues=path), i))
                G.remove_path_from(path)
                G.remove_single_nodes_from()
#                frag_aas = ''.join(three2One(self.pasta.residues[r].name) for r in path)
#                seq_aas = ''.join(self.pasta.seq[r].one_let for r in range(i, i + len(path)))
#                print 
#                print "Valid Mapping!"
#                print 'Frag %s %s' % (frag_aas, path)
#                print 'Seq  %s' % (seq_aas)
#                print 'Posteriors: '
#                print self.frag_posteriors(F, ones(len(m)), i)
#                print "Removing Mapped Nodes:"
#                print path
#                print
            
            if self.verbose:
                posteriors = self.frag_posteriors(F, ones(P.shape[1]), i)
                print 'Posteriors:      ', posteriors
                print 'Valid?       ', self.is_valid(F, m, i, t=1. / 20)
                print 'Edges:       ', G.edges()

#        for frag, i in assignment:
#            if frag.residues[0].__class__.__name__ == 'PastaResidue':
#                frag.residues_to_ints(self.pasta.residues)
#            from Definitions import three2One
#            frag_aas = ''.join([three2One(self.pasta.residues[r].name) for r in frag.residues])
#            print i, frag_aas

        return assignment
    
    def generate_assignments(self):
        assignments = []
        P = self.P
#        S = self.S
        for i in range(0, self.no_assignments):
            assignment = self.generate_assignment()
            similarity = self.seq_similarity(assignment, P)
            simple_conf = []
            for frag, seq_no in assignment:
                for res_no in frag.residues:
                    simple_conf.append((res_no, seq_no))
                    seq_no += 1
            print 'Assignment %i %.4f' % (i, similarity)
            assignments.append((simple_conf, similarity))
        return assignments
    
    def frag_posteriors(self, F, m, i):
        from numpy import multiply, diag
        
        F = multiply(F, m)  # # mask fragment
        j = i + F.shape[0]
        p = diag(F.T[i:j].T)
        
        return p
    
    def max_pos(self, F, m):
        from numpy import argmax, prod
        
        k = F.shape[0]  # # len fragment
        n = F.shape[1]  # # len sequence
        
        p = []
        i = 0
        while i <= n - k:
            p.append(prod(self.frag_posteriors(F, m, i)))
            i += 1
        i = argmax(p)
        
        return i
    
    def is_valid(self, F, m, i, t=1. / 20):
        """
        A Mapping is valid if all posteriors are above a predefined threshold t.
        
        @param F: Fragment Matrix
        @param i: Sequence position index
        @param t: Minimal posterior for amino acid
        
        @return: Bool
        """
        
        p = self.frag_posteriors(F, m, i)
        return (p > t).all()
    
    def invalid_nodes(self, F, m, i, path, t=1. / 20):
        """
        @param F: Fragment Matrix
        @param m: mask vector
        @param i: Sequence position index
        @param t: Minimal posterior for amino acid
        """
        
        posteriors = list(self.frag_posteriors(F, m, i))
        inval = []
        j = 0
        for p in posteriors:
            if p <= t:
                inval.append(path[j])
            j += 1
        return inval
    
    def all_invalid(self, F, m, i, path, t=1. / 20):
        k = F.shape[0]  # # len fragment
        return len(self.invalid_nodes(F, m, i, path, t)) == k
    
    def update_mask(self, F, M, i):
        """
        Sets row i in mask array M to zero if sequence position i has been
        assigned a residue and returns the updated mask array. 
        
        @param F: Fragment Matrix
        @type F: numpy.ndarray
        @param M: Mask Matrix
        @type M: numpy.ndarray
        @param i: Start position of fragment in sequence
        @type i: int
        
        @return: Mask array
        @rtype: numpy.ndarray
        """
        
        from numpy import zeros
        
        k = F.shape[0]  # # len fragment
        M[i:i + k] = zeros(k)  # # set occupied positons to zero
        return M
    
    def seq_similarity(self, assignment, P):
        """
        @param assignment: Fragment to sequence assignment (Fragment, start pos)
        @param P: Posterior Matrix
        """
        from numpy import array, ones
        
        m = P.shape[1]  # # len sequence
        
        errors = []
        for frag, i in assignment:
            F = self.F_from(P, frag.residues)
            mask = ones(m)
            error = 1. - self.frag_posteriors(F, mask, i)
            errors.append(error.sum())
        
        errors = array(errors)
        n = array([len(frag) for frag, i in assignment]).sum()
        gaps = m - n
        similarity = 1. - (1. / m * (errors.sum() + gaps))
        return similarity
    
class JointPosterior(SequenceMapping):  # a second concrete instrument
    
    def __init__(self, pasta_obj, no_assignments=100, verbose=False):
        from Linking import graph
        """
        Constructor
        
        @param pasta_obj: Pasta instance filled with typing and linking
        information
        @type pasta_obj: Pasta
        @param no_assignments: Number of assignments to be generated
        @param no_assignments: int
        @param verbose: verbosity option
        @type verbose: bool
        """
        
        self.P = pasta_obj.P  # # numpy array with posteriors
        self.L = pasta_obj.L  # # numpy array with linking
        self.S = self.S_from(pasta_obj)  # # list of amino acids in seq
        self.verbose = verbose
        self.fragment_graph = graph(self.L)
        self.no_assignments = no_assignments
    
    def S_from(self, pasta_obj):
        from Math import outer_2way
        from numpy import concatenate, multiply, swapaxes, ones
        """
        Reformats sequence profile array and returns it
        
        @requires: Pasta object has to initialized with sequence
        
        @param pasta_obj: Pasta instance
        @type pasta_obj: Pasta
        
        @return S: reformatted sequence array
        @rtype S: numpy.ndarray
        """
        S = super(JointPosterior, self).S_from(pasta_obj)
        S = concatenate([ones((1, S.shape[1])), S], 0)
        S = outer_2way(S[:-1].T, S[1:].T, multiply)
        S = swapaxes(S, 0, 1)
        return S
    
    def generate_assignment(self):
        from FragmentGraph import Fragment
        from copy import deepcopy
        from numpy import ones_like
        """
        Generates a list of residue to sequence pos assignments and returns
        it as list of tuples. First element in tuple corresponds to residue
        index and second element corresponds to sequence pos.
        
        @return assignment: list of residue sequence pos tuples
        @rtype assignment: list 
        """
        
        G = deepcopy(self.fragment_graph)
        P = self.P
        S = self.S
        M = ones_like(S)
        G.remove_single_nodes_from()
        assignment = []
        
        while len(G.edges()) != 0:
            if self.verbose:
                print 
                print 'G(|V|,|E|) = G(%i,%i)' % (len(G.nodes()), len(G.edges()))
                print 'Edges: ', G.edges()
                
            path = self.get_path(G, self.verbose)
            F = self.F_from(P, path)
            i = self.max_pos(F, S, M)
#            print 'path: ', path
            if not self.is_valid(F, S, i):
#                print "Invalid Mapping!"
                inval_nodes = self.invalid_nodes(F, S, i, path)
                if self.all_invalid(F, S, i, path):
                    G.remove_nodes_from(inval_nodes)
                    G.remove_single_nodes_from()
#                    print "All Residues Invalid!"
#                    print "Removing nodes: "
#                    print inval_nodes
                else:
                    inval_edges = G.get_invalid_edges(inval_nodes)
                    G.remove_edges_from(inval_edges)
                    G.remove_single_nodes_from()
#                    print "Removing edges:"
#                    print inval_edges
            else:
                M = self.update_mask(F, M, i)
                assignment.append((Fragment(residues=path), i))
                G.remove_path_from(path)
                G.remove_single_nodes_from()
#                print "Valid Mapping!"
#                print "Removing Mapped Nodes:"
                
            if self.verbose:
                T = self.frag_posteriors(F, S, i)
                p = self.p_vec(T)
                print 'Posteriors:      ', p
                print 'Valid?       ', self.is_valid(F, S, i)
                print 'Edges:       ', G.edges()
                
        return assignment
    
    def generate_assignments(self):
        """
        Generates a number of assignments and returns it as a list of tuples
        where first element is the list of residue to sequence pos tuples and 
        the second element is a similarity score for the respective assignment.
        
        @return: list of assign
        @rtype: list
        """
        assignments = []
        for i in range(0, self.no_assignments):
            assignment = self.generate_assignment()
            similarity = self.seq_similarity(assignment, self.P, self.S)
            simple_conf = []
            for frag, seq_no in assignment:
                for res_no in frag.residues:
                    simple_conf.append((res_no, seq_no))
                    seq_no += 1
            print 'Assignment %i %.4f' % (i, similarity)
            assignments.append((simple_conf, similarity))
        return assignments
    
#    def frag_posteriors(self, F, S, i):
#        """
#        Computes posterior vector for residue 
#        [i],[i-1]-joint posteriors at sequence position
#        i to i+n, where i is sequence position where first fragment
#        residue is mapped to and n is the fragment length.
#        
#        @param F: Fragment Matrix (len(F),len(amino_acids),len(amino_acids))
#        @type F: 3d numpy array
#        @param S: Sequence Matrix (len(S),len(amino_acids),len(amino_acids))
#        @type S: 3d numpy array 
#        @param i: Sequence position index
#        @type i: int
#        
#        @return T: 3d slice of fragment posteriors (len(F),len(amino_acids),len(amino_acids))
#        @type T: 3d numpy array
#        """
#        from numpy import multiply
#        
#        n = F.shape[0] ## Fragment length
#        m = S.shape[0] ## Sequence length
#        j = i + n
#        T = multiply(F, S[i:j])
#    
#        return T

    def p_vec(self, T):
        """
        Computes a vector of posteriors from a 3d matrix.
        For proper input (see below) the index of the vector posterior
        coresponds to the index of the residue in the fragment T from
        which it has be computed from.
        
        @param T: Should be a Fragment matrix that has been masked by
        sequence matrix
        @type T: numpy.ndarray
        
        @return: vector of posteriors
        @rtype: numpy.ndarray
        """
        p = T.sum(1).sum(1)
        return p
    
    def max_pos(self, F, S, M):
        """
        Computes the maximum posterior sequence position i, using a 
        Fragment matrix, a mask matrix M and a sequence matrix S.
        
        Mijk = 0, if seq pos i is already occupied by a fragment 
        of a previous step
        Mijk = 1, else
             
        Sijk = 0, if amino acid pair [j,k] not at position [i-1,i] in sequence
        Sijk = 1, else
        
        M and S are multiplied element wise, resulting in a new sequence matrix 
        S' which is 0 at all positions that have been assigned to a residue.
        In a sliding window approach F then is multiplied with sequence positions
        [i, i+len(F)-1] for 0 <= i <= (len(S) - len(F) - 1), resulting in a position
        posterior matrix T.
        
        The maximum posterior position can then be computed by summing over k,j,i
        
        argmax{sum(k,j,i){T}}
        
        @param F: Fragment Matrix (len(F),len(amino_acids),len(amino_acids))
        @type F: 3d numpy array
        @param S: Sequence Matrix (len(S),len(amino_acids),len(amino_acids))
        @type S: 3d numpy array 
        @param M: Mask Matrix (len(S),len(amino_acids),len(amino_acids))
        @type M: 3d numpy array
        
        @return: position at which the sum over all residue posteriors 
                            is maximal.
        @rtype: int
        """
        from numpy import argmax, multiply, array
        
        k = []
        
        n = F.shape[0]  # # Fragment length
        m = S.shape[0]  # # Sequence length
        
        i = 0
        j = i + len(F)
        while j <= len(S):
            S = multiply(S, M)
            T = self.frag_posteriors(F, S, i)
            k.append(self.p_vec(T).sum())
            i += 1
            j += 1
        max_seq_pos = argmax(array(k))
        return max_seq_pos
    
    def update_mask(self, F, M, i):
        """
        Sets a slice of the mask matrix to zero on position i to len(F) 
        
        @param F: Fragment Matrix (len(F),len(amino_acids),len(amino_acids))
        @type F: 3d numpy array
        @param M: Mask Matrix (len(S),len(amino_acids),len(amino_acids))
        @type M: 3d numpy array 
        @param i: Sequence position index
        @type i: int
        
        @return: 3d slice of fragment posteriors (len(F),len(amino_acids),len(amino_acids))
        @rtype: 3d numpy array
        """ 
        from numpy import zeros_like
        
        k = F.shape[0]  # length of fragment
        M[i:i + k] = zeros_like(F)  # set matrix block to zero
        return M
    
    def is_valid(self, F, S, i, t=1. / 20):
        """
        A Mapping is valid if ALL posteriors are above a predefined threshold t.
        
        @param F: Fragment Matrix (len(F),len(amino_acids),len(amino_acids))
        @type F: 3d numpy array
        @param S: Sequence Matrix (len(S),len(amino_acids),len(amino_acids))
        @type S: 3d numpy array 
        @param i: Sequence position index
        @type i: int
        @param t: Minimal posterior for amino acid
        @type t: float
        
        @return: condition
        @rtype: bool
        """
        
        T = self.frag_posteriors(F, S, i)
        p = self.p_vec(T)
        condition = (p > t).all()
        return condition
    
    def invalid_nodes(self, F, S, i, path, t=1. / 20):
        """
        Returns a sublist of nodes in path that have a lower posterior than the 
        minimally required value of t.
        
        @param F: Fragment Matrix (len(F),len(amino_acids),len(amino_acids))
        @type F: 3d numpy array
        @param S: Sequence Matrix (len(S),len(amino_acids),len(amino_acids))
        @type S: 3d numpy array 
        @param i: Sequence position index
        @type i: int
        @param t: Minimal posterior for amino acid
        @type t: float
        @param path: Index list representing residues in posterior matrix
        @type path: list
        
        @return: list of invalid nodes
        @rtype: list
        """
        T = self.frag_posteriors(F, S, i)
        posteriors = self.p_vec(T)
        inval = []
        j = 0
        for p in posteriors:
            if p < t:
                inval.append(path[j])
            j += 1
            
        return inval
    
    def all_invalid(self, F, S, i, path, t=1. / 20):
        """
        Checks whether all residues have a posterior below threshold t
        for the sequence position they have been assigned to.
        
        @param F: Fragment Matrix (len(F),len(amino_acids),len(amino_acids))
        @type F: 3d numpy array
        @param S: Sequence Matrix (len(S),len(amino_acids),len(amino_acids))
        @type S: 3d numpy array 
        @param i: Sequence position index
        @type i: int
        @param t: Minimal posterior for amino acid
        @type t: float
        
        @return: condition
        @rtype: bool
        """
        condition = len(self.invalid_nodes(F, S, i, path, t)) == F.shape[0]
        return condition
    
    def seq_similarity(self, assignment, P, S):
        """
        Returns a sequence similarity score for residues in assignment
        
        @param assignment: list of residue sequence pos tuples
        @type assignment: list
        @param P: Posterior Matrix
        @type P: numpy.ndarray
        @param S: Sequence Matrix
        @type S: numpy.ndarray
        
        @return: similarity score for assigment
        @rtype: float
        """
        from numpy import array
        
        errors = []
        for frag, i in assignment:
            F = self.F_from(P, frag.residues)
            T = self.frag_posteriors(F, S, i)
            error = 1. - self.p_vec(T)
            errors.append(error.sum())
        
        n = len(S)
        m = array([len(frag) for frag, i in assignment]).sum()
        errors = array(errors)
        e = n - m + errors.sum()
        similarity = 1 - (e / n)
        
        return similarity
    
class ILP(SequenceMapping):
    def __init__(self, pasta_obj, no_assignments=10, verbose=False):
        self.pasta = pasta_obj
        self.C = self.P_from(pasta_obj)
        self.L = self.L_from(pasta_obj)
        self.no_assignments = no_assignments
        self.verbose = verbose
        
    def S_from(self, pasta_obj):
        from Math import outer_2way
        from numpy import concatenate, ones, multiply, swapaxes
        
        s = super(ILP, self).S_from(pasta_obj)
        s = concatenate([ones((1, s.shape[1])), s], 0)
        S = outer_2way(s[:-1].T, s[1:].T, multiply)
        S = swapaxes(S, 0, 1)
        return S
    
    def P_from(self, pasta_obj):
        from numpy import dot, reshape, zeros, log, nan_to_num, clip
        from copy import deepcopy
        import pylab as pl
        n = len(pasta_obj.seq)  # # len sequence
        m = len(pasta_obj.residues)  # # len pasta list
        
        P = pasta_obj.P  # # posterior matrix from aa recog
        S = self.S_from(pasta_obj)
        P = dot(reshape(P, (m, -1)), reshape(S, (n, -1)).T)
        
        # set probabilities to 0 for seq pos 0, when residue has [i-1] shifts
        for i, r in enumerate(pasta_obj.residues):
            if len(r.shifts_im1.values()) > 0:
                P[i][0] = P[i][0] * 0
#        pickle.dump(P,open("Datasets/Ubiquitin/Results/Pickled/posteriorMatrixILP.p","wb"))
        
        C = zeros((n, n))
        C[0:m] = -log(clip(P, 1e-308, 1.))
#        pickle.dump(C,open("Datasets/Ubiquitin/Results/Pickled/costMatrixILP.p","wb"))
        C = C.tolist()
        return C
    
    def L_from(self, pasta_obj):
        from numpy import ones, eye
        from copy import deepcopy
        import pylab as pl
        n = len(pasta_obj.seq)
        m = len(pasta_obj.residues)  # # len pasta list
        
#        pickle.dump(pasta_obj.L,open("Datasets/Ubiquitin/Results/Pickled/linkingMatrixILP.p","wb"))
        
        L = ones((n, n))
        L[0:m, 0:m] = pasta_obj.L  # # fill with linking constraints
        L = L - eye(len(L)) * L.diagonal()  # # set diagonal to zero
        L = L.astype('i')
        L = L.T  # transpose
        L = L.tolist()
        return L
    
    def max_prob(self):
        m = float(len(self.pasta.residues))
        n = float(len(self.pasta.seq))
        assert m <= n
        return m / n
    
    def rnd_prob(self, data_points=10e4):
        '''
        Returns estimated mean and standard deviation of score distribution
        for randomized amino acid recognition result
         
        n := sum of all fragment lengths
        m := length of sequence
        '''
        from numpy import mean, std, array
        from numpy.oldnumeric.random_array import random
        
        m = float(len(self.pasta.residues))
        n = float(len(self.pasta.seq))
        assert m <= n
        
        scores = []
        for i in range(int(data_points)):
        
            p = random(m)
            avg = 1 - ((n - m + p.sum()) / n)
            scores.append(avg)
            
        data = array(scores)
        mu = mean(data)  # # mean value
        sigma = std(data)  # # standard deviation
        
        return mu, sigma
    
    def generate_assignment(self, c, l):
        from numpy import reshape, array, exp
        from pymprog import beginModel, endModel, verbose, iprod, var, \
        par, st, minimize, solve
        
        n = len(c[0])  # # sequence length
        m = len(c)  # # number of residues
        N = range(n)
        M = range(m)
        
#        print 'Minimize:'
#        print '    sum_{i=1}^{n}sum_{i=1}^{n} c_{ij}x_{ij}'
#        print 'Subject to:'
#        print '(1) sum_{j=1}^{n} x_{ij} = 1, for all i in {1,2,...,n}'
#        print '    \"Every SeqPos can have max. one assigned Residue.\"'
#        print '(2) sum_{i=1}^{n} x_{ij} = 1, for all j in {1,2,...,n}'
#        print '    \"Every Residue can be assigned to max one SeqPos.\"'
#        print '(3) x_{ij+1} <= (L^TX)_{ij},'
#        print '    for all i in {1,2,...,n} and for all j in {1,2,...,n-1}'
#        print '    \"Assignment has to be consistent with linking, s.t.'
#        print '     IF   residue k is assigned to sequence position j,' 
#        print '     THEN position j+1 can only be assigned to a residue k\''
#        print '          that is a valid successor of residue k, that is'
#        print '          when L[k][k\'] = 1\"'
#        print 
        
        beginModel("assign")
        verbose(False)  # Turn on this for model output
        A = iprod(M, N)  # combine index
        # declare define_variables
        x = var(A, 'x', kind=bool)  # assignment decision vars
        # declare parameters: 
        tc = par(c, 'C')
        lt = par(l, 'L')
            
        minimize(sum(tc[i][j] * x[i, j] for i, j in A), 'totalcost')
        # subject to: each agent works on at most one task
        # one for each agent a name for this group of constraints, optional
        # subject to: each task must be assigned to somebody
        st([sum(x[k, j] for j in N) == 1 for k in M], 'residue')
        st([sum(x[i, k] for i in M) == 1 for k in N], 'sequence pos')
        st([x[i, j + 1] <= sum(lt[i][k] * x[k, j] for k in N) for i, j in A if j < (n - 1)], 'linking')
        
        solve('int')
        assign = [(i, j) for i in M for j in N if x[i, j].primal > .5]
        
        costs = []
        if len(assign) > 0:
#            print("Total Cost = %g" % vobj()) 
            for i, j in assign:
                costs.append(tc[i][j].value)
                if i < len(self.pasta.residues):
                    print "Residue %d assigned to SeqPos %d with Prob %g" % (i, j, exp(-tc[i][j].value))
                else:
                    print "(Dummy-)Residue %d assigned to SeqPos %d with Prob %g" % (i, j, exp(-tc[i][j].value))
        
        # # remove dummy assignments and dummy costs
        clean_assign = []
        for i, j in assign:
            if i <= len(self.pasta.residues) - 1:
                clean_assign.append((i, j))
        assign = clean_assign
        costs = array(costs)
        
        similarity = self.seq_similarity(costs)
        model = endModel()  # make sure the model is alive
        print '****************************'
        print "Total Probability = %.2f" % similarity
        print "Maximum Probability = %.2f" % self.max_prob()
        print "Random Probability: %.2f (+- %.2f)" % self.rnd_prob()
        
        # define_variables for solution
        X = reshape([x[i, j].primal for i in range(n) for j in range(n)], (n, n))
#        import pickle
#        pickle.dump(array(self.C), open('Datasets/Ubiquitin/CostMatrix.p', 'wb'))
#        pickle.dump(array(self.L).T, open('Datasets/Ubiquitin/LinkingMatrix.p', 'wb'))
#        pickle.dump(X, open('Datasets/Ubiquitin/Results/Pickled/AssignmentMatrixILP.p', 'wb'))
        
        return assign, similarity
    
    def generate_assignments(self):
        assignments = []
        assign, similarity = self.generate_assignment(self.C, self.L)
        assignments.append((assign, similarity))
        
        return assignments
    
    def seq_similarity(self, costs):
        from numpy import exp
        
        if costs == []:
            similarity = 0
        else:
            errors = 1 - exp(-costs)
            m = len(self.pasta.residues)
            n = len(self.pasta.seq)
            e = n - m + errors.sum()
            similarity = 1 - (e / n)
        return similarity
    
    def update_costs(self, cost_matrix, assignment, no_seq_pos, no_res):
        """
        Returns an updated costs for next ILP iteration
        Example:
        For the Assignment (i,j) C[i,:] and C[:,j] in costs will be set
        to the maximum value (which one) 
        
        @param cost_matrix: list of costs for residue to sequence assignment
        @type cost_matrix: list of lists
        @param assignment: tuple of residue to sequence assignments
        @type assignment: tuple
        @param no_seq_pos: Number of sequence positions
        @type no_seq_pos: int
        @param no_res: Number of residues
        @type no_res: int
        
        @return: Cost matrix with updated costs (list of lists)
        @rtype: list
        """
        from numpy import log, array, ones

        cost_matrix = array(cost_matrix)
        for i, j in assignment:
            cost_matrix[i, :] = ones(no_res) * -log(1e-308)
            cost_matrix[:, j] = ones(no_seq_pos) * -log(1e-308)
        
        return cost_matrix.tolist()
    
    def update_linking(self, linking_matrix, assignment, no_seq_pos, no_res):
        """
        Returns an updated linking matrix for next ILP iteration
        Example:
        For the Assignment (i,j) the row i and col j in costs will be set
        to the maximum value (which one) 
            
        @param linking_matrix: list of neighbor constraints for residues
        @type linking_matrix: list of lists
        @param assignment: tuple of residue to sequence assignments
        @type assignment: tuple
        @param no_seq_pos: Number of sequence positions
        @type no_seq_pos: int
        @param no_res: Number of residues
        @type no_res: int
        
        @return: Linking matrix with updated neighbor constraints (list of lists)
        @rtype: list
        """
        from numpy import array, zeros
        linking_matrix = array(linking_matrix)
        for i, j in assignment:
            linking_matrix[i, :] = zeros(no_res)
            linking_matrix[:, j] = zeros(no_seq_pos)
        
        return linking_matrix.tolist()

class CplexILP(SequenceMapping):
    
    def __init__(self, pasta_obj, no_assignments=10, verbose=False):
        self.pasta = pasta_obj
        self.C = self.P_from(pasta_obj)
        self.L = self.L_from(pasta_obj)
        self.no_assignments = no_assignments
        self.verbose = verbose
    
    def S_from(self, pasta_obj):
        from Math import outer_2way
        from numpy import concatenate, ones, multiply, swapaxes
        
        s = super(CplexILP, self).S_from(pasta_obj)
        s = concatenate([ones((1, s.shape[1])), s], 0)
        S = outer_2way(s[:-1].T, s[1:].T, multiply)
        S = swapaxes(S, 0, 1)
        return S
    
    def P_from(self, pasta_obj):
        from numpy import dot, reshape, zeros, log, clip
        n = len(pasta_obj.seq)  # # len sequence
        m = len(pasta_obj.residues)  # # len pasta list
        
        P = pasta_obj.P  # # posterior matrix from aa recog
        S = self.S_from(pasta_obj)
        P = dot(reshape(P, (m, -1)), reshape(S, (n, -1)).T)
        
        # set probabilities to 0 for seq pos 0, when residue has [i-1] shifts
        for i, r in enumerate(pasta_obj.residues):
            if len(r.shifts_im1.values()) > 0:
                P[i][0] = P[i][0] * 0
        
        C = zeros((n, n))
#        C[0:m] = P
        C[0:m] = -log(clip(P, 1e-308, 1.))
#        C = -log(clip(C, 1e-308, 1.))
        return C
    
    def L_from(self, pasta_obj):
        from numpy import ones, eye
        from copy import deepcopy
        import pylab as pl
        n = len(pasta_obj.seq)
        m = len(pasta_obj.residues)  # # len pasta list
        
        L = ones((n, n))
        L[0:m, 0:m] = pasta_obj.L  # # fill with linking constraints
        L = L - eye(len(L)) * L.diagonal()  # # set diagonal to zero
        L = L.astype('i')
#        L = L.T
        return L
    
    def define_variables(self, x, shape):
        """
        Returns a numpy.ndarray with symbolic identifiers.
        eg.: "x", "y",...
        
        @param x: Identifiers
        @type x: string
        @param shape: Number of rows and cols for symbolic matrix
        @type shape: tuple
        
        @return: Matrix with symbolic identifiers
        @rtype: numpy.ndarray
        """
        from numpy import ones, array
        X = ones(shape).tolist()
        for i in range(0, shape[0]):
            for j in range(0, shape[1]):
                X[i][j] = x + "[" + str(i + 1) + "," + str(j + 1) + "]"
        return array(X)
    
    def permutation_rows(self, X):
        """
        Takes a symbolic matrix and returns rows for
        permutation constraints.
        rows =
        [[[x11,x12,...],[1,1,...],
        [[x21,x22,...],[1,1,...]]],
        ...
        [[x11,x21,...],[1,1,...]],
        [[x12,x22,...],[1,1,...]]]
        
        @param X: symbolic matrix
        @type X: numpy.ndarray
        
        @return: define_variables with their assigned factors
        @rtype: list
        """
        from numpy import ones
        rows = []
        for row in X:
            rows.append([row.tolist(), ones(len(row)).tolist()])
        for row in X.T:
            rows.append([row.tolist(), ones(len(row)).tolist()])
        return rows
    
    def linking_rows1(self, X, L):
        """
        Takes a quadratic symbolic matrix and a 
        quadratic linking matrix and returns rows for
        linking constraints.
        
        @attention Linking matrix should NOT be L.T
        
        The Linking constraints can be written as
        
        x_ij-1 <= (LX)_ij, forall i in [0,n] and forall j in [1,n]
        x_ij-1 - (LX)_ij <= 0
        x_ij-1 - [sum(k to n) L_ik * X_kj] <= 0 with n := X.shape[0] = L.shape[0]
        x_ij-1 - L_i0*X_0j - ... - L_in*X_nj <= 0
        
        The left hand side then will be reformatted in the following notation:
        
        [[define_variables],[factors]]
        [[x_ij-1, x_0j, x_1j,...,x_nj],[1,-L_i0,-L_i1,...,-L_in]]
        
        @param X: Symbolic matrix
        @type X: numpy.ndarray
        @param L: Linking Matrix
        @type L: numpy.ndarray
        
        @return: Variables with their assigned factors
        @rtype: list
        """
        assert X.shape == L.shape
        
        rows = []
        for i in range(0, X.shape[0]):
            for j in range(1, X.shape[1]):
                a = [X[i, j - 1]]
                a.extend(X[:, j].tolist())
                b = [1]
                b.extend((L[i, :] * -1).tolist())
                rows.append([a, b])
        return rows
    
    def linking_rows2(self, X, L):
        """
        Takes a quadratic symbolic matrix and a 
        quadratic linking matrix and returns rows for
        linking constraints.
        
        @attention Linking matrix SHOULD BE the L.T
        
        The Linking constraints can be written as
        
        x_ij-1 <= (LX)_ij, forall i in [0,n] and forall j in [1,n]
        x_ij-1 - (LX)_ij <= 0
        x_ij-1 - [sum(k to n) L_ik * X_kj] <= 0 with n := X.shape[0] = L.shape[0]
        x_ij-1 - L_i0*X_0j - ... - L_in*X_nj <= 0
        
        The left hand side then will be reformatted in the following notation:
        
        [[define_variables],[factors]]
        [[x_ij-1, x_0j, x_1j,...,x_nj],[1,-L_i0,-L_i1,...,-L_in]]
        
        @param X: Symbolic matrix
        @type X: numpy.ndarray
        @param L: Linking Matrix
        @type L: numpy.ndarray
        
        @return: Variables with their assigned factors
        @rtype: list
        """
        assert X.shape == L.shape
        
        rows = []
        for i in range(0, X.shape[0]):
            for j in range(0, X.shape[1] - 1):
                a = [X[i, j + 1]]
                a.extend(X[:, j].tolist())
                b = [1]
                b.extend((L[i, :] * -1).tolist())
                rows.append([a, b])
        return rows
        
    def seq_similarity(self, costs):
        from numpy import exp
        
        if costs == []:
            similarity = 0
        else:
            errors = 1 - exp(-costs)
            m = len(self.pasta.residues)
            n = len(self.pasta.seq)
            e = n - m + errors.sum()
            similarity = 1 - (e / n)
        return similarity
    
    def generate_assignments(self):
        import cplex
        from numpy import array, ones_like, zeros_like, \
        hstack, ones, arange, zeros, eye, zeros_like 
        import pylab as pl
        C = self.C
        L = self.L
        self.pasta.C = C
        self.pasta.ILP_L = L
        assignments = []
        assignments2 = []
        assignment_matrices = []
        
        # # uncomment this line in to test whether linking constraint works
        # # use expected data from Data generator
        # # this means, that the diagonal has higher costs than any other
        # # assignment, but the linking constraints enforce the assignment
        # # matrix to be the identity matrix.
#        C = eye(C.shape[0])
        
        X = self.define_variables("x", C.shape)
        n = C.shape[0]
        C = C.flatten()
        
        problem = cplex.Cplex()  # # problem instance
#        problem.parameters.threads.set(4) ## limit number of threads
        objective = C  # # objective function
        # # Define variables (here "x00,...,xnn")
        x_up_bounds = ones_like(C)  # # upper bounds
        x_lo_bounds = zeros_like(C)  # # lower bounds
        x_ctypes = ['I'] * C.size  # # variable types (I = integer)
        
        # # permutation constraints
        r1 = self.permutation_rows(X)
        r2 = self.linking_rows1(X, L)
        r3 = self.linking_rows2(X, L.T)
        
        # # number of constraints
        no_p = len(r1)
        no_l = len(r2)
        
        # # Add variables to problem instance
        problem.variables.add(obj=objective,
                              lb=x_lo_bounds,
                              ub=x_up_bounds, types=x_ctypes,
                              names=hstack(X).tolist())
        
        # # Minimize:
        problem.objective.set_sense(problem.objective.sense.minimize)
        # Subject to:
        problem.linear_constraints.add(lin_expr=r1,
                                       senses=['E'] * no_p,
                                       rhs=ones(no_p),
                                       names=arange(0, no_p).astype("str"))
        problem.linear_constraints.add(lin_expr=r2,
                                       senses=['L'] * no_l,
                                       rhs=zeros(no_l),
                                       names=arange(0, no_l).astype("str"))
        problem.linear_constraints.add(lin_expr=r3,
                                       senses=['L'] * no_l,
                                       rhs=zeros(no_l),
                                       names=arange(0, no_l).astype("str"))
        problem.solve()
        pool = self.populate(problem)
        
        numsol = pool.get_num()  # # number of solutions
        for i in range(numsol):        
            objval_i = pool.get_objective_value(i)
            x_i = pool.get_values(i)
            A = array(x_i).reshape((n, n))
            assign = []
            costs = []
            for i, row in enumerate(A):
                for j, col in enumerate(row):
                    if i <= len(self.pasta.residues) - 1:
                        if col == 1:
                            assign.append((i, j))
                            costs.append(self.C[i][j])
            costs = array(costs)
            similarity = self.seq_similarity(costs)
            assignment_matrices.append(A)
            assignments.append((assign, similarity))
            assignments2.append((assign, objval_i))
        
        self.pasta.A = assignments
        self.pasta.B = assignments2
        self.pasta.Xs = assignment_matrices
#        assignments.append((assign, similarity))
        return assignments
    
    def populate(self, problem):
        """
        Generates solution pool and returns it
        """
        import cplex
        # set the solution pool relative gap parameter to obtain solutions
        # of objective value within 10% of the optimal 
        problem.parameters.mip.pool.relgap.set(0.1)
        
        try:
            problem.populate_solution_pool()
        except cplex.CplexSolverError:
            print "Exception raised during populate"
            return
    
        print
        # solution.get_status() returns an integer code
        print "Solution status = " , problem.solution.get_status(), ":",
        # the following line prints the corresponding string
        print problem.solution.status[problem.solution.get_status()]
    
        numcols = problem.variables.get_num()
    
        x = problem.solution.get_values()
    
        # Print information about other solutions
        print
        numsol = problem.solution.pool.get_num()
        print "The solution pool contains %d solutions." % numsol
    
        numsolreplaced = problem.solution.pool.get_num_replaced()
        print "%d solutions were removed due to the solution pool relative gap parameter." % numsolreplaced
        
        numsoltotal = numsol + numsolreplaced
        print "In total, %d solutions were generated." % numsoltotal
    
        meanobjval = problem.solution.pool.get_mean_objective_value()
        print "The average objective value of the solutions is %.10g." % meanobjval
    
        # write out the objective value of each solution and its
        # difference to the incumbent 
        names = problem.solution.pool.get_names()
        
        print
        print "Solution        Objective       Number of define_variables"
        print "                value           that differ compared to"
        print "                                the incumbent"
        
        epszero = 1e-10
        for i in range(numsol):        
    
            objval_i = problem.solution.pool.get_objective_value(i)
            
            x_i = problem.solution.pool.get_values(i)
            
            # compute the number of define_variables that differ in solution i
            # and in the incumbent
            numdiff = 0;
            for j in range (numcols):
                if abs (x_i[j] - x[j]) > epszero:
                    numdiff = numdiff + 1
            print "%-15s %-10g      %d / %d" % (names[i], objval_i, numdiff, numcols)
            
        return problem.solution.pool
    
#    def solution_to_assignment_matrix(self,x,n):
#        from numpy import array
#        A = array(x).reshape((n,n))
#        return A
#    
#    def solution_to_assignment_list(self, solution):
#        A = self.solution_to_assignment_matrix(solution,self.C.shape[0])
#        
#        from numpy import array
#        assign = []
#        costs = []
#        for i, row in enumerate(A):
#            for j, col in enumerate(row):
#                if i <= len(self.pasta.residues) - 1:
#                    if col == 1:
#                        assign.append((i, j))
#                        costs.append(self.C[i][j])
# #        print assign
#        costs = array(costs)
#        similarity = self.seq_similarity(costs)
#        return assign, similarity
    
    
#    def generate_assignments(self):
#        assignments = []
#        assign, similarity = self.generate_pool(self.C, self.L)
#        assignments.append((assign, similarity))
#        return assignments
    
class Context(object):
    """
    A class that utilizes a strategy to perform an operation
    """
    def __init__(self, pasta_obj, no_assignments, verbose=False, strategy=None):
        
        _strategies = {None:CplexILP(pasta_obj),
                       'Single':SinglePosterior(pasta_obj, no_assignments, verbose),
                       'JointNew':JointPosteriorNew(pasta_obj, no_assignments, verbose),
                       'Joint':JointPosterior(pasta_obj, no_assignments, verbose),
                       'ILP':ILP(pasta_obj),
                       'CplexILP':CplexILP(pasta_obj)}
        
        self.strategy = _strategies[strategy]
        
    @property
    def strategy(self):
        return self._strategy
    
    @strategy.setter
    def strategy(self, strategy):
        assert isinstance(strategy, SequenceMapping)
        self._strategy = strategy
        print 'Activated strategy: {0}'.format(strategy.__class__.__name__)
        
    def execute(self):
        """The operation being executed using the currently active strategy"""
        assignments = self.strategy.generate_assignments()
        return assignments
        
if __name__ == '__main__':
    from Pasta import Pasta
    from numpy import argmax, exp
    from Definitions import three2One
    from FileHandler import FileHandler
    import pylab as pl
    import pickle
    
    statsfile = 'shiftPresets/bmrb.shift'
    fh = FileHandler()
    
#    folder = "GeneratedDatasets/"
#    pastafile = folder + "exp_data.pasta"
#    seqfile = folder + "sequence.fasta"
#    result_folder = folder + "Results/"
    
#    folder = "tests/"
#    pastafile = folder + "residue_lists/all_pairs.pasta"
#    seqfile = folder + "sequences/all_pairs.fasta"
#    result_folder = folder + "Results/"
    
    folder = "Datasets/Ubiquitin/"
    pastafile = folder + "residue_lists/Ub_opt_unambiguous_unassigned.pasta"
    seqfile = folder + "sequences/Ub.fasta"
    result_folder = folder + "Results/"

    folder = "tests/"
    statsfile = folder + 'reference_lists/bmrb.shift'
    pastafile = folder + "residue_lists/all_singles.pasta"
    seqfile = folder + "sequences/all_singles.fasta"
    result_folder = folder + "Results/"

    pasta = Pasta(pastafile, statsfile, seqfile)
    pasta.typing()
    print "Typing done"
    pasta.linking("CplexILP", .6)
    print "Linking done"
    

#    strategies = ["Joint", "CplexILP", "ILP"]
    strategies = ["CplexILP"]
    for strat in strategies:
        context = Context(pasta, 100, strategy=strat)
        assignments = context.execute()
        
#        records = fh.to_fasta_records(assignments,
#                                      pasta.seq,
#                                      pasta.residues)
        
#        fh.write_fasta(result_folder + "Ub_bmrb_unassigned_costs=0_%s.fasta" % (strat), records, pasta.seq)
#        fh.conf_to_pasta(result_folder + "Ub_bmrb_unassigned_costs=0_%s.pasta" % (strat), assignments[0],
#                         pasta.residues, pasta.amino_acids)
        
    pickle.dump(pasta, open(result_folder + "FOO.p", "wb"))
    pasta = pickle.load(open(result_folder + "FOO.p", "r"))
    
    print pasta.A
    print pasta.B
    print pasta.C
    print pasta.ILP_L
    for X in pasta.Xs:
        pl.matshow(X)
    pl.show()
