from numpy import mean, max, nan_to_num
import numpy as np

class Likelihood(object):

    def __init__(self):
        '''
        Initialize mapping dictionaries.
            - aa_mapping: Map atom labels from residue list to amino acid reference list
            - res_mapping: Map atom labels from residue list to residue list
        '''
        from Mapping import Mapping
        self.aa_mapping = Mapping(aa_mapping=True)
        self.res_mapping = Mapping(aa_mapping=False)
        
    def __call__(self, residue1, residue2, previous=False, summarize=mean):
        """
        Evaluates total probability of all atom mappings from residue1 to residue2
        
        Tests on residue2 type: 
        If it is an amino acid then use AA_MAPPING.
        If it is a residue then use RES_MAPPING
        
        @param residue1: A Residue.PastaResidue which shifts are to be evaluated
        @type residue1: Residue.PastaResidue
        @param residue2: A Residue.PasaResidue or a Residue.AminoAcid against residue1 is being evaluated
        @type residue2: Residue
        @param previous: option to set which shifts are used to compute likelihood
        @type previous: bool
        @param summarize: function to summarize likelihoods from differrent atom mappings
        @type summarize: numpy.mean
        
        @return: Likelihood for residue1
        @rtype: float
        """
        if type(residue2.shifts.values()[0]) == list:
            mapping = self.aa_mapping
        if type(residue2.shifts.values()[0]) == float:
            mapping = self.res_mapping
        
        mappings = mapping.generate_mappings(residue1, residue2, previous)
        prob = [self.calc_likelihood(residue1, residue2, m) for m in mappings]
        if prob == []:
            prob = [0.]
        # TODO
#        print "Prob: ", prob
#        print "Max Prob: ", max(prob)
#        print "Mean Prob: ", mean(prob)
        return summarize(prob)
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
        
        @todo: Check whether it really does this !!!
        
        @param residue1: A Residue.PastaResidue which shifts are to be evaluated
        @type residue1: Residue.PastaResidue
        @param residue2: A Residue.PasaResidue or a Residue.AminoAcid against residue1 is being evaluated
        @type residue2: Residue
        
        @return: Likelihood for residue1
        @rtype: float
        """
        
        shifts = residue1.get_shifts(mapping)
        return residue2.calc_likelihood(shifts)
    
    def calc_likelihoods(self,residue_list, reference_list, previous=False, summarize=mean):
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
        from Mapping import Mapping
        m = Mapping()
        
        x = [[self(r, a, previous, summarize) for a in reference_list] for r in residue_list]
        x = array(x)
        
        # set likelihood to 1 if residue has no [i-1] shifts
        if previous:
            i = 0
            for r in residue_list:
                shifts_im1, keys_im1 = m.get_carbons(r, previous=True)
                if len(shifts_im1) == 0:
                    x[i] = ones_like(x[i])
                i += 1
        
        # set likelihood to 1 if residue has no [i] shifts
        else:
            i = 0
            for r in residue_list:
                shifts_i, keys_i = m.get_carbons(r, previous=False)
                if len(shifts_i) == 0:
                    x[i] = ones_like(x[i])
                i += 1
            
        return x
    
if __name__ == '__main__':
    
    from FileHandler import FileHandler
    from numpy import argmax, mean, max
    import pylab as pl
    fh = FileHandler()
    L = Likelihood()
#    pastafile = 'multiple_test_files/Ubiquitin/Ub_new.pasta'
    pastafile = 'tests/residue_lists/all_singles.pasta'
    statsfile = 'shiftPresets/bmrb.shift'
    residues = fh.read_pasta(pastafile, statsfile)
    amino_acids = fh.read_preset(statsfile)
    print ' '.join([aa.three_let for aa in amino_acids])
#    r = residues[1]
#    print r.name
    
#    p11 = L.calc_likelihoods(residues, amino_acids, previous=False, summarize=mean)
#    p12 = L.calc_likelihoods(residues, amino_acids, previous=True, summarize=mean)
    p21 = L.calc_likelihoods(residues, amino_acids, previous=False, summarize=max)
    p22 = L.calc_likelihoods(residues, amino_acids, previous=True, summarize=max)
    
    print p21
    print p22
    
#    pl.matshow(p11)
#    pl.colorbar()
#    pl.matshow(p21)
#    pl.colorbar()
    pl.matshow(p21)
    pl.colorbar()
    pl.matshow(p22)
    pl.colorbar()
    pl.show()
#    print argmax(p_i[0])
#    print argmax(p_i.sum())
#    print argmax(p_im1), p_im1
