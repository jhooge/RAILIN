from Mapping2 import AbstractMapping

class Likelihood(object):

    def __init__(self):
        pass
        
    def __call__(self, residue1, residue2, previous=False):
        """
        """
        m = AbstractMapping.create("on_amino_acid")
        mapping = m.find_mapping(residue1, residue2, previous)
        prob = self.calc_likelihood(residue1, residue2, mapping)
        if prob == []:
            prob = [0.]
        # TODO
#        print "Prob: ", prob
#        print "Max Prob: ", max(prob)
#        print "Mean Prob: ", mean(prob)
        return prob
#        return nan_to_num(summarize(prob))

    def calc_likelihood(self, residue1, residue2, mapping):
        """
        Calculates likelihood for particular mapping
        Depending on the type of residue2, calc_likelihood is called either on a 
        Residue.PastaResidue or on Residue.AminoAcid
        Likelihood is only calculated for shifts of a given mapping
        
        IF type(residue2) is Residue.AminoAcid:
        THEN return log-likelihood that residue1 one is the corresponding amino acid of residue2
        ELSE:
        return log-likelihood that residue2 is the successor of residue1
        
        @param residue1: A Residue.PastaResidue which shifts are to be evaluated
        @type residue1: Residue.PastaResidue
        @param residue2: A Residue.PasaResidue or a Residue.AminoAcid against residue1 is being evaluated
        @type residue2: Residue
        
        @return: Likelihood for residue1
        @rtype: float
        """
        like = 0
        if len(mapping) > 0:
            shifts = residue1.get_shifts(mapping)
            like = residue2.calc_likelihood(shifts)
        return like
    
    def calc_likelihoods(self,residue_list, reference_list, previous=False):
        '''
        Takes a residue list and reference list and computes the log-likelihood for corresponding
        amino acids or the log-likelihood for an [i],[i-1] residue pair.
        log-likelihoods for all possible res,aa- or res,res-mappings are 
        then summarized by their mean value.
        Returns ndarray with log-likelihoods
        
        @param residue_list: a list of Residue.PastaResidue objects
        @type residue_list: list
        @param reference_list: a list of Residue.AminoAcid objects
        @type reference_list: list
        @param previous: option to set which shifts are used to compute likelihood
        @type previous: bool
        @param summarize: function to summarize likelihoods from differrent atom mappings
        @type summarize: numpy.mean
        
        @return: numpy.ndarray with likelihoods for different amino_acids for each residue
        @rtype: array
        '''
        
        from numpy import array, ones_like
        
        x = [[self(r, a, previous) for a in reference_list] for r in residue_list]
        x = array(x)
        
        # set likelihood to 1 if residue has no [i-1] shifts
        if previous:
            i = 0
            for r in residue_list:
                shifts_im1, keys_im1 = r.get_carbons(previous=True)
                if len(shifts_im1) == 0:
                    x[i] = ones_like(x[i])
                i += 1
        
        # set likelihood to 1 if residue has no [i] shifts
        else:
            i = 0
            for r in residue_list:
                shifts_i, keys_i = r.get_carbons(previous=False)
                if len(shifts_i) == 0:
                    x[i] = ones_like(x[i])
                i += 1
            
        return x
    
if __name__ == '__main__':
    
    from FileHandler import FileHandler
    from numpy import argmax, mean, max, zeros, ones
    import pylab as pl
    fh = FileHandler()
    L = Likelihood()
#    pastafile = 'multiple_test_files/Ubiquitin/Ub_new.pasta'
    pastafile = 'tests/residue_lists/all_singles.pasta'
    statsfile = 'tests/reference_lists/bmrb.shift'
    residues = fh.read_pasta(pastafile, statsfile)
    amino_acids = fh.read_preset(statsfile)
    aas = amino_acids
    print ' '.join([aa.three_let for aa in amino_acids])
#    r = residues[1]
#    print r.name
    
    A = L.calc_likelihoods(residues, amino_acids, previous=False)
    
    pl.matshow(A)
    pl.colorbar()
    pl.show()
