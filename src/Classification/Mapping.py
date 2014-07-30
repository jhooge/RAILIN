'''
Created on Dec 9, 2010

@author: jhooge
'''
from munkres import Munkres
from numpy import subtract, fabs, inf, array

AA_ATOMS = \
            {'HN': ['H'],
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
                    'CE3', 'CH2', 'CZ2', 'CZ3']
             }
            
RES_ATOMS = \
            {'HN': ['H'],
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
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6','C7'],
             'C2' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6','C7'],
             'C3' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6','C7'],
             'C4' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6','C7'],
             'C5' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6','C7'],
             'C6' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6','C7'],
             'C7' :['CA', 'CB', 'CG', 'CD', 'CE', 'CZ',
                    'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2',
                    'CE3', 'CH2', 'CZ2', 'CZ3',
                    'C1', 'C2', 'C3', 'C4', 'C5', 'C6','C7']
             }

class Mapping(dict):
    """
    If Residue is mapped to AminoAcid then AA_ATOMS dict has to be used.
    If Residue is mapped to Residue then RES_ATOMS dict has to be used.
    """
    def __init__(self, aa_mapping=True):
        if aa_mapping:
            dict.__init__(self, AA_ATOMS)
        else:
            dict.__init__(self, RES_ATOMS)

    def __map_atom(self, pasta_name):
        '''
        Removes 'i-1' from PastaResidue atom label to get atom shift
        returns [] if label not in Mapping dictionary
        
        Parameters:
        -----------
        pasta_name : string
                String of nucleus label
        
        Returns:
        --------
        out : list
                If label in self.values(), list with stripped label
                If label not in self.values(), []
        
        Notes:
        ------
        
        Examples:
        ---------
        >>> __map_atom(self, 'CAi-1')
        ['CA']
        >>> __map_atom(self, 'foobari-1')
        []
        '''
        out = self.get(pasta_name.replace('i-1', ''), [])

        return out
    
    def generate_mappings(self, residue1, residue2, previous=False):
        '''
        Generates a list of lists of tuples of possible assignments.
        The mappings are all bijective and therefore valid.
        
        Parameters:
        -----------
        residue1 : Residue
                Should be a Residue.PastaResidue object.
        residue2 : Residue
                Can be a Residue.PastaResidue object or a Residue.AminoAcid
                object.
        previous : bool, optional
                If False, residue1.shifts_i() are mapped on residue2 shifts.
                If True, residue1.shifts_im1() are mapped on residue2 shifts.
        
        Returns:
        --------
        mappings : 2dlist
            list of tuple lists representing nuclei mappings
        
        Notes:
        ------
        if previous == False:
            residue1.shifts_i() are mapped on residue2 shifts
            if type(residue2) == Residue.AminoAcid:
                shifts = residue2.shifts ## Residue.AminoAcid mean shifts
            if type(residue2) == Residue.PastaResidue:
                shifts = residue2.shifts_i() ## Residue.PastaResidue [i]-shifts
        else:
            residue1.shifts_im1() are mapped on residue2 shifts
            if type(residue2) == Residue.AminoAcid:
                shifts = residue2.shifts ## Residue.AminoAcid mean shifts
            if type(residue2) == Residue.PastaResidue:
                shifts = residue2.shifts_i() ## Residue.PastaResidue [i]-shifts
        
        Examples:
        ---------
        PastaResidue/PastaResidue Mapping (for Matching)
        >>> m = Mapping(aa_mapping=False)
        >>> r1 = PastaResidue(1, 0., 0, 'H', 0, {'C1':1.0,'C2':2.0})
        >>> r2 = PastaResidue(2, 0., 0, 'H', 0, {'C1':3.0,'C2':4.0, 'C1i-1':1.0,'C2i-1':2.0})
            
        >>> print m.generate_mappings(r2, r1, previous=True)
        [[('C2i-1', 'C1'), ('C1i-1', 'C2')], [('C2i-1', 'C2'), ('C1i-1', 'C1')]]
            
        >>> print m.generate_mappings(r1, r2, previous=True)
        [[]]
            
        PastaResidue/AminoAcid Mapping (for Amino Acid Recognition)
        >>> m = Mapping(aa_mapping=True)
        >>> a = AminoAcid('X', 'XXX', 0, 'H', CA=[1.0,2.0], CB=[3.0,4.0])
            
        >>> print m.generate_mappings(r2, a, previous=False) ## [i] shift recognition
        [[('C2', 'CA'), ('C1', 'CB')], [('C2', 'CB'), ('C1', 'CA')]]
            
        >>> print m.generate_mappings(r2, a, previous=True) ## [i-1] shift recognition
        [[('C2i-1', 'CA'), ('C1i-1', 'CB')], [('C2i-1', 'CB'), ('C1i-1', 'CA')]]
        '''
        
        import itertools

        shifts, keys = self.get_carbons(residue1, previous)
        res_shifts = dict(zip(keys, shifts))

        _from = res_shifts.keys()
        _to = [self.__map_atom(atom) for atom in _from]
        _to = [[atom for atom in atoms if atom in residue2.shifts] for atoms in _to]

        ## get only injective mappings
        mappings = [i for i in itertools.product(*_to)]
        mappings = [m for m in mappings if len(set(m)) == len(m)]
        mappings = [zip(_from, m) for m in mappings]
        
        return mappings

    '''***************** RESIDUE NUCLEUS MAPPING *****************'''

    def get_carbons(self, residue, previous=False):
        shifts = []
        if not previous:
            keys = array(residue.shifts_i.keys())
            if len(keys) > 0:
                mask = array(map(lambda label: label.startswith('C'), keys))
                keys = keys[mask]
                shifts = array([residue.shifts_i[key] for key in keys])
        else:
            keys = array(residue.shifts_im1.keys())
            if len(keys) > 0:
                mask = array(map(lambda label: label.startswith('C'), keys))
                keys = keys[mask]
                shifts = array([residue.shifts_im1[key] for key in keys])
        return shifts, keys

    def isValid(self, residue1, residue2):
        shifts_i, keys_i = self.get_carbons(residue1, previous=False)
        shifts_im1, keys_im1 = self.get_carbons(residue2, previous=True)
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
        
        shifts_i, keys_i = self.get_carbons(residue1, previous=False)
        shifts_im1, keys_im1 = self.get_carbons(residue2, previous=True)
        
        delta = fabs(subtract.outer(shifts_i, shifts_im1))
        C = 1 - (delta <= tolerance).astype('i')
        
        return C
    
    def compute_nucleus_mapping(self, residue1, residue2, tolerance):
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
        @return mapping: list of tupels defining a path through cost matrix
        @rtype: list
        '''
        costs = inf
        mapping = []
        
        shifts_i, keys_i = self.get_carbons(residue1, previous=False)
        shifts_im1, keys_im1 = self.get_carbons(residue2, previous=True)
        
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
    
    def linkable(self, residue1, residue2, strategy, tolerance,conservative=True):

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
            return costs[i,j].sum() == 0.
    
if __name__ == "__main__":
    pass
#    from FileHandler import FileHandler
#    from Definitions import one2Three
#    from numpy import argmin
#    fh = FileHandler()
#    m = Mapping(aa_mapping=True)
#    
#    folder = "Datasets/Ubiquitin/"
#    pastafile = folder + "residue_lists/Ub_opt_unambiguous.pasta"
#    statsfile = "shiftPresets/bmrb.shift"
#    seqfile = folder + "sequences/Ub.fasta"
#    
#    residues = fh.read_pasta(pastafile, statsfile)
#    amino_acids = fh.read_preset(statsfile)
#    del amino_acids[1]
#    sequence = fh.read_seq(seqfile, statsfile)
#    aas = [aa.one_let for aa in amino_acids]
#    
#    for res in residues:
#        if res.name_im1 != "NAA":
#            k = aas.index(res.name_im1[0])
#            aa = amino_acids[k]
#            mappings = m.generate_mappings(res, aa,previous=True)
#            print res.pasta_no, aa.one_let
#            diffs = []
#            for i, mapping in enumerate(mappings):
#                
#                d = [abs(res.shifts_im1[a] - aa.shifts[b][0])
#                     for a, b in mapping]
#                diffs.append(d)
#            diffs = array(diffs)
#            if diffs != []:
#                s = argmin(diffs.sum(1))
#                min_map = mappings[s]
#                min_delta = diffs[s]
#                    
#                for a, b in min_map:
#                    res.shifts_im1[b+"i-1"] = res.shifts_im1.pop(a)
#    
#    for res in residues:
#        res.name = "NAA"
#        res.name_im1 = "NAA"
#    fh.write_pasta(folder + "residue_lists/Ub_opt_unambiguous_unassigned.pasta", residues)
