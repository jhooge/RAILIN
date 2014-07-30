'''
Created on Feb 28, 2012

@author: jhooge
'''
from Residue import PastaResidue
from FileHandler import FileHandler
from numpy import random, arange
from numpy.random import uniform, shuffle
from copy import deepcopy
import os

class DataGenerator(object):
    '''
    This class can be used to generate random datasets.
    It includes methods to add uniform noise to chemical shifts
    instrument specific instrument error.
    
    Example:
    map(type,generate_dataset(2))
    >>> [<type 'list'>, <type 'list'>, <type 'list'>]
    '''

    def __init__(self, statsfile):
        '''
        Constructor
        '''
        fh = FileHandler()
        ## list of Residue.AminoAcid objects
        self.amino_acids = fh.read_preset(statsfile)
        del self.amino_acids[1] ## Remove Cy*
        
    def generate_dataset(self, size, noise=0., instrument_error=0., fn=None):
        """
        This method generates input files for use in Pasta
        A dataset in general consists of ASCII formatted
        sequence template file, an unassiged pseudo residue
        list and a fully assigned pseudo residue list for validation
        purposes.
        
        @param size: number of amino acids in generated sequence 
        and pasta residues 
        @type size: int
        @param noise: uniform noise on [i] and [i-1] chemical shifts in sigma-folds
        @type noise: float
        @param instrument_error: instrument (or experiment) specific 
        error in ppm
        @type instrument_error: float
        @param fn: filename where dataset will be saved
        @type fn: string
        
        @return: complete dataset: (sequence,raw_data,expected_data)
        @rtype: tuple(list,list,list)
        """
        
        seq = self.generate_rnd_seq(size)
        raw = self.generate_raw_data(seq, noise, instrument_error)
        exp = self.generate_expected_data(seq, noise, instrument_error)
        
        if fn is not None:
            if not os.path.exists(fn):
                os.makedirs(fn)
            os.chdir(fn)
        
            fh = FileHandler()
            fh.write_seq('sequence.fasta',"Template", seq)
            fh.write_pasta('raw_data.pasta', raw)
            fh.write_pasta('exp_data.pasta', exp)
        
        dataset = (seq, raw, exp)
        return dataset
        
    def generate_rnd_seq(self, length):
        '''
        Generates a random list of Residue.AminoAcid objects.
        
        @param length: desired sequence length
        @type length: int
        
        @return: list of random Amino Acid objects
        @rtype: list
        '''
        
        seq = []
        i = 0
        while i < length:
            j = random.randint(len(self.amino_acids))
            seq.append(self.amino_acids[j])
            i += 1
            
        return seq
    
    def generate_expected_data(self, seq, noise=0., instrument_error=0.):
        '''
        Generates a list of Residue.PastaResidue objects
        from Residue.AminoAcid object list. The shift values are equal
        to the reference values and every PastaResidue comes with an
        assigned three letter name.
        
        @param seq: list of Residue.AminoAcid objects
        @type seq: list
        @param noise: noise level in factors of standard dev 
        @type noise: float
        @param instrument_error: NMR specific measurement error
        @type instrument_error: float
        
        @return: list of Residue objects
        @rtype: list
        '''
        
        data = []
        i = 0
        while i < len(seq):
            res = self.__convert(seq[i])
            data.append(res)
            i += 1
            
        ## add uniform noise [i] and [i-1 shifts]
        data = self.add_uniform_noise_to(data, seq, noise)
        
        i = len(data) - 1
        j = i - 1
        while j >= 0: ## set [i-1] shifts
            shifts = deepcopy(data[j].shifts_i)
            keys = shifts.keys()
            values = shifts.values()
            
            k = 0
            while k < len(keys):
                keys[k] += 'i-1'
                k += 1
                
            shifts = dict(zip(keys, values))
            data[i].shifts_im1 = shifts
            i -= 1
            j = i - 1
        
        data = self.add_instrument_error(data, instrument_error)
        self.__assign_attributes(data, seq)
        
        return data
    
    def generate_raw_data(self, seq, noise, instrument_error):
        """
        Generates a list of Residue.PastaResidue objects
        from Residue.AminoAcid object list adds uniform noise and 
        instrumental error to the shift values and shuffles the dataset.
        
        
        @param seq: list of Residue.AminoAcid objects
        @type seq: list
        @param noise: noise level in factors of standard dev 
        @type noise: float
        @param instrument_error: NMR specific measurement error
        @type instrument_error: float
        
        @return: list of Residue.PastaResidue objects 
        with noisy shift data
        @rtype: list
        """
        
        data = self.generate_expected_data(seq, noise, instrument_error)
        self.unassign_attributes(data, seq)
        data = self.__shuffle_residues(data)
        return data
    
    def add_uniform_noise_to(self, data, seq, noise):
        """
        Adds random uniform noise to data. 
        Amount of noise depends of expected amino acid standard deviation.
        
        @attention: assumes len(data) == len(seq)
    
        @param data: list of Residue.PastaResidue objects
        @type data: list
        @param seq: list of Residue.AminoAcid objects
        @type seq: list
        @param noise: noise level in factors of standard dev 
        @type noise: float
        
        @return: list of Residue.PastaResidue objects with noisy shift data
        @rtype: list
        """
        
        i = 0
        while i < len(data): ## add noise to i
            for key in data[i].shifts_i.keys():
                mu = seq[i].shifts[key][0]
                std = seq[i].shifts[key][1]
                if noise != 0.:
                    a = mu - (std * noise)
                    b = mu + (std * noise)
                    data[i].shifts_i[key] = uniform(a, b)
                    
    #                print '%.2f <= %.2f <= %.2f'%(a,data[i][key],b)
                    assert data[i].shifts_i[key] >= a and data[i].shifts_i[key] <= b

            i += 1
            
        return data
    
    def add_instrument_error(self, data, error):
        '''
        Adds random uniform noise to data. This is independent
        of expected amino acid and therefore should only effect
        linking performance.
        The error is defined as the deviation of 
        measurements of corresponding [i] and [i-1] shifts in ppm.
        
        @attention: assumes len(data) == len(seq)
        
        @param data: list of Residue.PastaResidue objects
        @type data: list
        @param error: NMR specific measurement error in ppm
        @type error: float
        
        @return: list of Residue.PastaResidue objects with added instrumental error.
        @rtype: list
        '''
        from copy import deepcopy
        data = deepcopy(data)
        for res in data:
            shifts_i = res.shifts_i
            shifts_im1 = res.shifts_im1
            for key_i, key_im1 in zip(shifts_i.keys(), shifts_im1.keys()):
                a = shifts_i[key_i] - error
                b = shifts_i[key_i] + error
                c = shifts_im1[key_im1] - error
                d = shifts_im1[key_im1] + error
                shifts_i[key_i] = uniform(a, b)
                shifts_im1[key_im1] = uniform(c, d)
                
                assert shifts_i[key_i] >= a and shifts_i[key_i] <= b 
                assert shifts_im1[key_im1] >= c and shifts_im1[key_im1] <= d
                
        return data
    
    def __assign_attributes(self, data, seq):
        '''
        For each measurement in data
        sets name and sequence position 
        attributes to expected values.
        
        @attention: assumes len(data) == len(seq)
        
        @param data: list of Residue.PastaResidue objects
        @type data: list
        @param seq: list of Residue.AminoAcid objects
        @type seq: list
        '''
        
        i = 0
        while i < len(seq):
            data[i].name = seq[i].three_let
            data[i].pasta_no = i + 1
            data[i].original_no = i + 1
            data[i].seq_pos = i
            i += 1
        i = len(seq) - 1
        j = i - 1
        while j >= 0:
            data[i].name_im1 = seq[j].three_let
            i -= 1
            j = i - 1
    
    def unassign_attributes(self, data, seq):
        '''
        For each measurement in data
        sets name and sequence position 
        attributes to initial values.
        
        @param data: list of Residue.PastaResidue objects
        @type data: list
        @param seq: list of Residue.AminoAcid objects
        @type seq: list
        '''
        
        i = 0
        while i < len(data):
            data[i].name = 'NAA'
            data[i].seq_pos = -1
            i += 1
        i = len(seq) - 1
        j = i - 1
        while j >= 0:
            print i
            data[i].name_im1 = 'NAA'
            i -= 1
            j = i - 1
    
    def __shuffle_residues(self, data):
        '''
        Shuffles a list of residues and returns the result
        
        @param data: list of residue objects
        @type data: list
        @return: Shuffled residue list
        @rtype: list
        '''
        
        ind = arange(0, len(data))
        shuffle(ind)
        shuffled = [data[i] for i in ind]
        return shuffled
    
    def __is_aa(self, obj):
        '''
        Tests whether an object is a Residue.AminoAcid
        
        @param obj: any object
        @type obj: object
        
        @return: True if obj is an amino acid, else False
        @rtype: bool
        '''
        y = type(obj).__name__ == 'AminoAcid'
        return y
    
    def __convert(self, aa):
        '''
        Converts an Residue.AminoAcid to a Residue.PastaResidue object.
        [i] Shifts of PastaResidue correspond to amino acid reference shifts
        
        @raise ValueError: if aa is not a Residue.AminoAcid object
        
        @param aa: AminoAcid object with reference shifts
        @type aa: Residue.AminoAcid
        
        @return: Residue object with reference shifts
        @rtype: Residue.PastaResidue
        '''
        
        if not self.__is_aa(aa): 
            raise ValueError('Value is not an Amino Acid object!')
        
        keys = aa.shifts.keys()
        values = [mu for mu, std in aa.shifts.values()]
        res = PastaResidue(shifts_i=dict(zip(keys, values)))
        return res

def add_noise_to_residues(data, noise):
    i = 0
    while i < len(data): ## add noise to i
        
        for key in data[i].shifts_i.keys():
            mu = data[i].shifts_i[key]
            std = seq[i].shifts[key][1]
            if noise != 0.:
                a = mu - (std * noise)
                b = mu + (std * noise)
                data[i].shifts_i[key] = uniform(a, b)
                
#                print '%.2f <= %.2f <= %.2f'%(a,data[i][key],b)
                assert data[i].shifts_i[key] >= a and data[i].shifts_i[key] <= b
        for key in data[i].shifts_im1.keys():
            mu = data[i].shifts_im1[key]
            print seq[i-1].three_let, seq[i].shifts
            std = seq[i-1].shifts[key.replace("i-1","")][1]
            if noise != 0.:
                a = mu - (std * noise)
                b = mu + (std * noise)
                data[i].shifts_im1[key] = uniform(a, b)
                
#                print '%.2f <= %.2f <= %.2f'%(a,data[i][key],b)
                assert data[i].shifts_im1[key] >= a and data[i].shifts_im1[key] <= b

        i += 1
        
    return data

if __name__ == "__main__":
    from numpy import linspace
    import pickle
    
    pastafile = "Datasets/Ubiquitin/residue_lists/Ub_bmrb.pasta"
    statsfile = 'tests/reference_lists/bmrb.shift'
    seqfile = "Datasets/Ubiquitin/sequences/Ub_bmrb.fasta"
    fh = FileHandler()
    dg = DataGenerator(statsfile)
    
    residues = fh.read_pasta(pastafile, statsfile)
    seq =  fh.read_seq(seqfile, statsfile)
    
    for noise in linspace(0,4,41): 
        noisy_residues = add_noise_to_residues(residues, noise)
        dg.unassign_attributes(noisy_residues, seq)
        fn = 'Datasets/Ubiquitin/residue_lists/NoisyData/UbqBaxNoise=%.1f.pasta'%noise
        fh.write_pasta(fn, noisy_residues)
    
    
#    seq, raw, exp = dg.generate_dataset(10,
#                                        noise=2.,
#                                        instrument_error=0.,
#                                        fn='GeneratedDatasets')

#    seq = fh.read_seq(seqfile,statsfile)
#    filenames = []
#    for n in linspace(0,4,50):
#        data = dg.generate_expected_data(seq,n)
#        dg.unassign_attributes(data, seq)
#        fn = 'tests/residue_lists/all_pairs_noise=%.4f.pasta'%n
#        fh.write_pasta(fn, data)
#        print fn, "WRITTEN"
