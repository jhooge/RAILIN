'''
Created on Jan 21, 2010

@author: jhooge
'''

from Mapping import Mapping

def graph(A):
        from FragmentGraph import FragmentGraph
        import networkx as nx
        if A.sum() == 0:
            raise ValueError('No Residues Linked')
        else:
            G = nx.from_numpy_matrix(A, create_using=FragmentGraph())
        return G

class Linking(object):
    '''
    classdocs
    '''

    def __init__(self, residues):
        '''
        Constructor
        '''
        self.residues = residues
        
#    def linking_matrix(self,strategy,tolerance=.6):
#        """
#        Computes a nxn binary matrix of possible successors
#        rows := [i-1] Residues (predecessors)
#        cols := [i] Residues (successors
#        
#        Lij = 1, if [i] and [i-1] shift deviation is within tolerance
#        Lij = 0, else
#        
#        @bug: SequenceMapping with strategies other than ILP will 
#        only work without postProcessing
#        
#        @param tolerance: Instrument (or experimental) error
#        @type tolerance: float
#        @return: Matrix filled with ones when residue is a direct neighbor
#                   of another one
#        @rtype: numpy.ndarray
#        """
#        
#        from numpy import zeros
#        m = Mapping()
#        L = zeros((len(self.residues), len(self.residues)))
#        
#        i=0
#        for res1 in self.residues:
#            j=0
#            for res2 in self.residues:
#                mapping = m.compute_nucleus_mapping(res1, res2, tolerance)
#                if mapping != []:
#                    L[j][i] = 1.
#                j+=1
#            i+=1
#            
#        L = L.transpose()
#        if strategy == "CplexILP" or strategy == "ILP":
#                L = self.postProcessLinking(L)
#        
#        return L
    
    def linking_matrix(self, strategy, tolerance=.6, conservative=True):
        """
        Computes a nxn binary matrix of possible successors
        rows := [i-1] Residues (predecessors)
        cols := [i] Residues (successors)
        
        Lij = 1, if [i] and [i-1] shift deviation is within tolerance
        Lij = 0, else
        
        @bug: SequenceMapping with strategies other than ILP will 
        only work without postProcessing
        
        @param tolerance: Instrument (or experimental) error
        @type tolerance: float
        @return: Matrix filled with ones when residue is a direct neighbor
                 of another one
        @rtype: numpy.ndarray
        """
        
        from numpy import zeros, eye

        m = Mapping()
        L = zeros((len(self.residues), len(self.residues)), 'i')
        
#        if strategy != "ILP" or strategy != "CplexILP":
#            conservative = False
        
        for i, previous in enumerate(self.residues):
            for j, current in enumerate(self.residues):
                if m.linkable(previous, current, 
                              strategy, tolerance, 
                              conservative=conservative):
                    L[i, j] = 1
        return L * (1 - eye(len(L)))
    
    def postProcessLinking(self, L):
        from Mapping import Mapping
        from numpy import ones , eye
        """
        if len(shifts_i) == 0:
            residue could be predecessor for any other residue, but itself
            --> set row to 1
        if len(shifts_im1) == 0:
            residue could be successor for any other residue, but itself
            --> set col to 1
        set diagonal to zero
        """
        m = Mapping()
        n = len(self.residues)
        
        for i, res in enumerate(self.residues):
            
            shifts_im1, keys_im1 = m.get_carbons(res, True)
            shifts_i, keys_i = m.get_carbons(res, False)
            if len(shifts_im1) == 0:
                print 'Missing [i-1] shifts! Pasta No: ', res.pasta_no
                L[:, i] = 1
            if len(shifts_i) == 0:
                print 'Missing [i] shifts! Pasta No: ', res.pasta_no
                L[i, :] = 1
        
        return L * (1 - eye(len(L)))
        
    def diff(self, pred, succ, mapping):
        difference = {}
        
        for t in mapping:
            i = t[0]
            j = t[1]
            difference[t] = pred.shifts_i[j] - succ.shifts_im1[i]
            
        return difference
        
    def check_pos_constraints(self, linking_matrix, residues, tol):
        m = Mapping()
        L = linking_matrix
        for i in range(0, L.shape[0]):
            for j in range(0, L.shape[1]):
                if L[i][j] == 1:
                    mapping = m.compute_nucleus_mapping(residues[i],
                                                    residues[j],
                                                    tol)
                    d = self.diff(residues[i], residues[j], mapping)
                    print i, j
                    print "Any Mapping? ", len(mapping) > 0
                    for key in d.keys():
                        print key, d[key], "Within Tolerance? ", abs(d[key]) < tol
                        assert abs(d[key]) < tol
                    print residues[i]
                    print residues[j]
    
    def check_neg_constraints(self, linking_matrix, residues, tol):
        m = Mapping()
        L = linking_matrix
        for i in range(0, L.shape[0]):
            for j in range(0, L.shape[1]):
                if L[i][j] == 0:
                    mapping = m.compute_nucleus_mapping(residues[i],
                                                    residues[j],
                                                    tol)
                    print "Any Mapping? ", len(mapping) > 0
                    assert not len(mapping) > 0
                    d = self.diff(residues[i], residues[j], mapping)
                    print i, j
                    for key in d.keys():
                        print key, d[key], "Within Tolerance? ", abs(d[key]) < tol
                        assert abs(d[key]) < tol
                    print residues[i]
                    print residues[j]
    
if __name__ == '__main__':
    from FileHandler import FileHandler
    from Linking import Linking
    import pylab as pl
    
    fh = FileHandler()
    statsfile = 'shiftPresets/bmrb.shift'
    folder = "Datasets/Ubiquitin/"
    fn = folder + "residue_lists/Ub_opt_relabeled.pasta"
    seqfile = folder + "sequences/Ub.fasta"
    result_folder = folder + "Results/"
    
    tol = .6
    strat = "Joint"
    residues = fh.read_pasta(fn, statsfile)
    link = Linking(residues)
    
    L = link.linking_matrix(strat, tolerance=tol, conservative=True)
    L2 = link.linking_matrix(strat, tolerance=tol, conservative=False)
    pl.matshow(L)
    pl.matshow(L2)
    pl.show()
    
    
    
    
##    link.check_pos_constraints(L,residues,tol)
##    link.check_neg_constraints(L,residues,tol)
#    
#    A = fh.assignment_matrix_from_assigned(fn, seqfile, statsfile).astype("i")
##    L_expected = fh.linking_matrix_from_assigned(A).astype("i")
#    L_expected = fh.expected_linking(residues).astype("i")
#    L_estimated = link.linking_matrix(strat, tolerance=tol)
#    L_diff = L_expected - L_estimated
#    
#    link.check_pos_constraints(L_diff,residues,tol)
#    
#
##    for i in range(0,L.shape[0]):
##        for j in range(0,L.shape[1]):
##            if L[i][j] == 1:
##                mapping = m.compute_nucleus_mapping(residues[i],
##                                                residues[j],
##                                                tol)
##                d =  link.diff(residues[i],residues[j],mapping)
##                print i,j
##                for key in d.keys():
##                    print key,d[key], "Within Tolerance? ", abs(d[key])<tol
##                    assert abs(d[key])<tol
##                pred_no = residues[i].pasta_no
##                succ_no = residues[j].pasta_no
##                print residues[i]
##                print residues[j]
#
#    import pickle
#
#    pickle.dump(L_expected,open(folder+"Results/Pickled/L_expected.p","wb"))
#    pickle.dump(L_estimated,open(folder+"Results/Pickled/L_estimated.p","wb"))
#    pickle.dump(L_diff,open(folder+"Results/Pickled/L_diff.p","wb"))
#    pl.matshow(L_expected)
#    pl.title("Expected")
#    pl.matshow(L_estimated)
#    pl.title("Estimated")
#    pl.matshow(L_diff)
#    pl.title("L_expected - L_estimated")
#    pl.colorbar()
#    pl.show()
    
    
