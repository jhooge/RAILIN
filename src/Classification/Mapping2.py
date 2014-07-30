'''
Created on Jun 3, 2012

@author: jhooge
'''

from abc import ABCMeta, abstractmethod
from munkres import Munkres
from numpy import subtract, fabs, inf, array

class AbstractMapping(object):
    __metaclass__ = ABCMeta
    
    RULES = {} ## default rules for atom mapping, defined in subclass
    
    def init(self, rules=None):
        """
        Initializes Atom to Atom rules for two Residues
        
        @param rules: Dictionary of atom to atom mapping rules
        @type rules: dictionary
        
        @return: Dictionary of atom to atom mapping rules
        @rtype: dictionary
        """
        if rules is None:
            return self.RULES
        else:
            return rules
    
    @staticmethod
    def create(type, rules=None):
        if type == "on_amino_acid":
            return OnAminoAcid(rules=rules)
        if type == "on_residue":
            return OnResidue(rules=rules)
    
    @abstractmethod
    def find_mapping(self, residue1, residue2, previous=False):
        pass
    
    
class OnResidue(AbstractMapping):
    
    RULES = {'HN': ['H'],
             'N15':['N'],
             'CO' :['CO'],
             'CA' :['CA'],
             'CB' :['CB'],
             'CG' :['CG', 'CG1', 'CG2'],
             'CG1':['CG1', 'CG2'],
             'CG2':['CG1', 'CG2'],
             'CD' :['CD', 'CD1', 'CD2'],
             'CD1':['CD1', 'CD2'],
             'CD2':['CD1', 'CD2'],
             'CE' :['CE', 'CE1', 'CE2', 'CE3'],
             'CE1':['CE1', 'CE2', 'CE3'],
             'CE2':['CE1', 'CE2', 'CE3'],
             'CE3':['CE1', 'CE2', 'CE3'],
             'CH2':['CH2'],
             'CZ':['CZ', 'CZ2', 'CZ3'],
             'CZ2':['CZ2', 'CZ3'],
             'CZ3':['CZ2', 'CZ3'],
             'C1' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
             'C2' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
             'C3' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
             'C4' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
             'C5' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
             'C6' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
             'C7' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']}
    
    def __init__(self, rules=None):
        self.rules = super(OnResidue, self).init(rules)

    def isValid(self, residue1, residue2):
        shifts_i, keys_i = residue1.get_carbons(previous=False)
        shifts_im1, keys_im1 = residue2.get_carbons(previous=True)
        return residue1 != residue2 and len(shifts_im1) > 0 and len(shifts_i) > 0
    
    def compute_cost_matrix(self, residue1, residue2, tolerance):
        '''
        Computes a cost matrix C for Munkres for a predecessor/successor pair
        i := residue1
        j := residue2
        C_ij = 0, if deviation of i-shifts of residue1 and [i-1]-shifts of
                 residue2 is within tolerance
        C_ij = 1, else
        
        @todo: Double check the filling of cost matrix
        
        @param residue1: Possible predecessor of residue2
        @type residue1: Residue.PastaResidue
        @param residue2: Possible successor of residue1
        @type residue2: Residue.PastaResidue
        @return: C
        @rtype: numpy.ndarray
        '''
        
        shifts_i, keys_i = residue1.get_carbons(previous=False)
        shifts_im1, keys_im1 = residue2.get_carbons(previous=True)
        
        delta = fabs(subtract.outer(shifts_i, shifts_im1))
        C = 1 - (delta <= tolerance).astype('i')
        
        return C
    
    def find_mapping(self, residue1, residue2, tolerance):
        '''
        This method computes one possible mapping of [i] carbon
        nuclei of residue1 to [i-1] carbon nuclei of residue2.
        Valid mappings will be computed via Munkres algorithm.
        A mapping is valid if the sum of costs in the Munkres cost matrix
        is equal to zero.
        
        @todo: Better explanation for mapping if cost_matrix is rectangular
        
        @param residue1: Possible predecessor of residue2
        @type residue1: Residue.PastaResidue
        @param residue2: Possible successor of residue1
        @type residue2: Residue.PastaResidue
        @param tolerance: Instrument (or experimental) error
        @type tolerance: float
        @return mapping: list of tupels defnining a path through cost matrix
        @rtype: list
        '''
        costs = inf
        mapping = []
        
        shifts_i, keys_i = residue1.get_carbons(previous=False)
        shifts_im1, keys_im1 = residue2.get_carbons(previous=True)
        
        if self.isValid(residue1, residue2):
            m = Munkres()
            cost_matrix = self.compute_cost_matrix(residue1, residue2, tolerance)
            rows, cols = cost_matrix.shape
            if rows > cols:
                cost_matrix = cost_matrix.transpose()
            
            indices = m.compute(cost_matrix * 1)
            costs = sum([cost_matrix[i, j] for i, j in indices])
            
            if costs == 0.:
                if rows > cols:
                    mapping = [(keys_im1[j], keys_i[i]) for j, i in indices]
                else:
                    mapping = [(keys_im1[j], keys_i[i]) for i, j in indices]
                    
#            print 'Compute Nuceus Mapping for Residues'
#            print 'Residue %i [i] - shifts'%(residue1.pasta_no)
#            print 'Residue %i [i-1] - shifts'%(residue2.pasta_no)
#            print 'Cost Matrix: \n',cost_matrix
#            print 'Costs: ',costs
#            print 'Mapping: ', mapping
        
        return mapping
    
    def linkable(self, residue1, residue2, strategy, tolerance, conservative=True):

        from numpy import array, sum
        
        if not self.isValid(residue1, residue2):
            if strategy == "ILP" or strategy == "CplexILP":
                return True
            else: 
                return False
        else:
            m = Munkres()
            costs = self.compute_cost_matrix(residue1, residue2, tolerance)
            if conservative:
                return sum(costs == 0) > 0
            rows, cols = costs.shape
            if rows > cols:
                costs = costs.T
            
            i, j = array(m.compute(costs * 1)).T
            return costs[i, j].sum() == 0.
    
class OnAminoAcid(AbstractMapping):
    
    RULES = {'HN': ['H'],
             'N15':['N'],
             'CO' :['CO'],
             'CA' :['CA'],
             'CB' :['CB'],
             'CG' :['CG', 'CG1', 'CG2'],
             'CG1':['CG1', 'CG2'],
             'CG2':['CG1', 'CG2'],
             'CD' :['CD', 'CD1', 'CD2'],
             'CD1':['CD1', 'CD2'],
             'CD2':['CD1', 'CD2'],
             'CE' :['CE', 'CE1', 'CE2', 'CE3'],
             'CE1':['CE1', 'CE2', 'CE3'],
             'CE2':['CE1', 'CE2', 'CE3'],
             'CE3':['CE1', 'CE2', 'CE3'],
             'CH2':['CH2'],
             'CZ':['CZ', 'CZ2', 'CZ3'],
             'CZ2':['CZ2', 'CZ3'],
             'CZ3':['CZ2', 'CZ3'],
             'C1' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3'],
             'C2' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3'],
             'C3' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3'],
             'C4' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3'],
             'C5' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3'],
             'C6' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3'],
             'C7' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3']}
    
    def __init__(self, rules=None):
        self.rules = super(OnAminoAcid, self).init(rules)
    
#    def doMap(self, residue1, residue2, previous=False):
#        """
#        Generates a list of lists of tuples of possible assignments from 
#        atoms in residue1 to atoms of residue2
#        
#        @param residue1: Should be a Residue.PastaResidue object.
#        @type residue1: Residue.PastaResidue
#        @param residue2: Can be a Residue.PastaResidue object or a Residue.AminoAcid
#        object.
#        @type residue2: Residue
#        @param previous: If False, residue1.shifts_i() are mapped on residue2 shifts.
#        If True, residue1.shifts_im1() are mapped on residue2 shifts.
#        @type previous: bool
#        
#        @return: list of tuple lists representing nuclei mappings
#        @rtype: list
#        
#        @note: if previous == False:
#        residue1.shifts_i() are mapped on residue2 shifts
#        IF type(residue2) == Residue.AminoAcid:
#        shifts = residue2.shifts ## Residue.AminoAcid mean shifts
#        IF type(residue2) == Residue.PastaResidue:
#        shifts = residue2.shifts_i() ## Residue.PastaResidue [i]-shifts
#        ELSE:
#        residue1.shifts_im1() are mapped on residue2 shifts
#        IF type(residue2) == Residue.AminoAcid:
#        shifts = residue2.shifts ## Residue.AminoAcid mean shifts
#        IF type(residue2) == Residue.PastaResidue:
#        shifts = residue2.shifts_i() ## Residue.PastaResidue [i]-shifts
#        
#        Examples:
#        ---------
#        PastaResidue/PastaResidue Mapping (for Matching)
#        >>> m = Mapping(aa_mapping=False)
#        >>> r1 = PastaResidue(1, 0., 0, 'H', 0, {'C1':1.0,'C2':2.0})
#        >>> r2 = PastaResidue(2, 0., 0, 'H', 0, {'C1':3.0,'C2':4.0, 'C1i-1':1.0,'C2i-1':2.0})
#            
#        >>> print m.doMap(r2, r1, previous=True)
#        [[('C2i-1', 'C1'), ('C1i-1', 'C2')], [('C2i-1', 'C2'), ('C1i-1', 'C1')]]
#            
#        >>> print m.doMap(r1, r2, previous=True)
#        [[]]
#            
#        PastaResidue/AminoAcid Mapping (for Amino Acid Recognition)
#        >>> m = Mapping(aa_mapping=True)
#        >>> a = AminoAcid('X', 'XXX', 0, 'H', CA=[1.0,2.0], CB=[3.0,4.0])
#            
#        >>> print m.doMap(r2, a, previous=False) ## [i] shift recognition
#        [[('C2', 'CA'), ('C1', 'CB')], [('C2', 'CB'), ('C1', 'CA')]]
#            
#        >>> print m.doMap(r2, a, previous=True) ## [i-1] shift recognition
#        [[('C2i-1', 'CA'), ('C1i-1', 'CB')], [('C2i-1', 'CB'), ('C1i-1', 'CA')]]
#        """
#        import itertools
#        shifts, keys = residue1.get_carbons(previous)
#        res_shifts = dict(zip(keys, shifts))
#
#        _from = res_shifts.keys()
#        _to = [self.__strip_atom_label(atom) for atom in _from]
#        _to = [[atom for atom in atoms if atom in residue2.shifts] for atoms in _to]
#
#        ## get only injective mappings
#        mappings = [i for i in itertools.product(*_to)]
#        mappings = [m for m in mappings if len(set(m)) == len(m)]
#        mappings = [zip(_from, m) for m in mappings]
#        
#        return mappings
    
    def __strip_atom_label(self, atom_label):
        """
        Removes 'i-1' from PastaResidue atom label to get atom shift
        returns [] if label not in Mapping dictionary
        
        @param atom_key: String of atom label
        @type atom_key: string
        
        @return: Stripped 
        
        Examples:
        ---------
        >>> __strip_atom_label(self, 'CAi-1')
        ['CA']
        >>> __strip_atom_label(self, 'foobari-1')
        []
        """
        out = self.rules.get(atom_label.replace('i-1', ''), [])

        return out

    def compute_cost_matrix(self, shifts_i, shifts_j, previous=False):
        C = fabs(subtract.outer(shifts_i, shifts_j))
        return C
    
    def ambiguous(self, res_key, previous=False):
        if previous:
            return len(self.rules[res_key[:-3]]) > 1
        else:
            return len(self.rules[res_key]) > 1

    def atoms_valid(self, res_keys, aa_keys, previous=False):
        condition = []
        for k in res_keys:
            if not self.ambiguous(k, previous):
                if previous:
                    k = k.replace("i-1","")
                condition.append(k in aa_keys)
        return array(condition).all()

    def find_mapping(self, res, aa, previous=False):
        mapping = []
        
        res_keys = res.get_carbons(previous)[1]
        aa_keys = aa.get_carbons()[1]
        
#        print "res_keys ", res_keys
#        print "aa_keys ", aa_keys
        
        n = len(res_keys)
        m = len(aa_keys)
        if n <= m and self.atoms_valid(res_keys, aa_keys, previous): # assumes number of res shifts < number of aa shifts
            all_ambiguous = array(map(lambda x: self.ambiguous(x, previous),
                                      res_keys)).all()
            all_non_ambiguous = array(map(lambda x: not self.ambiguous(x, previous),
                                          res_keys)).all()
#            print "All Ambiguous ", all_ambiguous
#            print "All Non Ambiguous ", all_non_ambiguous
            
            namb_map = self.find_non_ambiguous_mapping(res, aa,
                                                       previous)
            amb_map = self.find_ambiguous_mapping(res, aa,
                                                  namb_map,
                                                  previous)
            
#            print res.name, aa.three_let, namb_map, amb_map, aa_keys
            
            if all_ambiguous:
                mapping = amb_map
            if all_non_ambiguous:
                mapping = namb_map
            else:
                amb_map.extend(namb_map)
                mapping = amb_map
        return mapping
    
    def find_ambiguous_mapping(self, res, aa,
                               non_ambiguous_mapping,
                               previous=False):
        from munkres import Munkres
        from numpy import fabs, sum, delete
        
        total_costs = None
        mapping = []
        ambiguous_keys = []
        ambiguous_shifts = []
        res_shifts, res_keys = res.get_carbons(previous)
        aa_shifts, aa_keys = aa.get_carbons()
        
        for i, j in non_ambiguous_mapping:
            if j in aa_keys:
                k = list(aa_keys).index(j)
                aa_shifts = delete(aa_shifts, k)
                aa_keys = delete(aa_keys, k)
        for i, key in enumerate(res_keys):
            if self.ambiguous(key, previous):
                ambiguous_keys.append(key)
                ambiguous_shifts.append(res_shifts[i])
                
        if len(aa_keys) > 0 and len(ambiguous_shifts) > 0:
            costs = fabs(subtract.outer(ambiguous_shifts, aa_shifts))
            munkres = Munkres()
            result = munkres.compute(costs * 1.)
            for i, j in result:
                mapping.append((ambiguous_keys[i], aa_keys[j]))
            
        return mapping
    
    def find_non_ambiguous_mapping(self, res, aa, previous=False):
#        non_ambiguous_keys = []
        mapping = []
        res_shifts, res_keys = res.get_carbons(previous)
        aa_shifts, aa_keys = aa.get_carbons()
        res_shifts = dict(zip(res_keys, res_shifts))
        aa_shifts = dict(zip(aa_keys, aa_shifts))
        
        for key in res_shifts.keys():
            if not self.ambiguous(key, previous):
#                non_ambiguous_keys.append(key)
                if previous:
                    mapping.append((key, self.rules[key[:-3]][0]))
                else:
                    mapping.append((key, self.rules[key][0]))
        return mapping
    

if __name__ == "__main__":
    import pylab as pl
    from numpy import max, mean, take, zeros
    from Residue import PastaResidue, AminoAcid
    from MaxLikelihood import Likelihood
    from FileHandler import FileHandler
    
    fh = FileHandler()
    fn = "/is/ei/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/"\
         "Classification/tests/reference_lists/bmrb.shift"
#    pastalist = 'tests/residue_lists/all_singles.pasta'
    pastalist = "Datasets/Ubiquitin/residue_lists/"\
                "incompleteData/Ub_bmrb_missing_shifts_0.50.pasta"
    
    statsfile = fn
    amino_acids = fh.read_preset(fn)
    residues = fh.read_pasta(pastalist, statsfile)
    toAA = AbstractMapping.create("on_amino_acid")
    toRes = AbstractMapping.create("on_residue")
    
    aas = amino_acids
    del aas[1]
    
    L = Likelihood()
    A = L.calc_likelihoods(residues, amino_acids, previous=False)
    
    pl.matshow(A)
    pl.colorbar()
    pl.show()
    
