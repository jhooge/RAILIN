'''
Created on Mar 2, 2012

@author: jhooge
'''

from abc import abstractmethod, ABCMeta
from copy import deepcopy
from Math import gaussian_pdf as pdf
#from Math import laplace_pdf as pdf
from numpy import sum, exp, zeros, argmax, array
from numpy.ma.core import prod, mean
import Definitions as DEF

__metaclass__ = type

class Residue(object):
    '''
    classdocs
    '''
    
    __metaclass__ = ABCMeta
    
    def __init__(self,
                 seq_pos= -1,
                 sec_struc='RC',
                 profile=None,
                 shifts={}):
        self.seq_pos = seq_pos
        self.sec_struc = sec_struc
        self.profile = profile
        self.shifts = shifts
        
    def set_attribute(self, key, value):
        attributes = self.__dict__
        if key not in attributes:
            raise AttributeError('Attribute does not exist in %s!' 
                                 % (self.__class__.__name__))
        attributes[key] = value
    
    def get_attribute(self, key):
        attributes = self.__dict__
        if key not in attributes:
            raise AttributeError('Attribute does not exist in %s!' 
                                 % (self.__class__.__name__))
        return attributes[key]
    
    @abstractmethod
    def calc_likelihood(self, shifts):
        pass
    
    @abstractmethod
    def set_profile(self, param):
        pass

class PastaResidue(Residue):
    
    def __init__(self,
                 pasta_no= -1,
                 original_no= -1,
                 energy=0,
                 assigned=False,
                 shifts_i={},
                 shifts_im1={},
                 name='NAA',
                 name_im1='NAA',
                 comment_i='',
                 comment_im1=''):
        super(PastaResidue, self).__init__()
        self.pasta_no = pasta_no
        self.original_no = original_no
        self.energy = energy
        self.assigned = assigned
        self.shifts_i = shifts_i
        self.shifts_im1 = shifts_im1
        self.name = name
        self.name_im1 = name_im1
        self.comment_i = comment_i
        self.comment_im1 = comment_im1
        
    def get_shifts(self, mapping):
        shifts = dict(self.shifts_i.items() + 
                      self.shifts_im1.items())
        
        a = array([item in shifts.items() for item in self.shifts_i.items()]).all()
        b = array([item in shifts.items() for item in self.shifts_im1.items()]).all()
        assert a and b
        
        return dict([(b, shifts[a]) for a, b in mapping])
    
    def get_carbons(self,previous=False):
        carbon_shifts = []
        if previous:
            keys = array(self.shifts_im1.keys())
            if len(keys) > 0:
                mask = array(map(lambda label: label.startswith('C'), keys))
                keys = keys[mask]
                carbon_shifts = array([self.shifts_im1[key] for key in keys])
        else:
            keys = array(self.shifts_i.keys())
            if len(keys) > 0:
                mask = array(map(lambda label: label.startswith('C'), keys))
                keys = keys[mask]
                carbon_shifts = array([self.shifts_i[key] for key in keys])
        return carbon_shifts, keys
    
    def calc_likelihood(self, shifts):
        if len(shifts): ## self has no [i-1] shifts
            sigma = .5 * .3
            return exp(sum([pdf(x, self[n], sigma, log=True) for n, x in shifts.items()]))
        else: ## self has [i-1] shifts
            return .0

    def set_profile(self, param):
        '''
        @param param: list of posteriors
        @type param: list
        '''
        self.profile = param
    
    def assign_name(self,amino_acids):
        '''
        @todo: check whether this method is used anywhere
        '''
#        tmp = DEF.THREE_LET[:]
        amino_acids = deepcopy(amino_acids)
        tmp = map(lambda a: a.three_let, amino_acids)
        del tmp[1] # remove CY*
        if self.profile is not None:
            self.name = tmp[argmax(self.profile)]
        else:
            raise ValueError('Residue profile has not been set in residue with pasta no %i !' % self.pasta_no)
    
    def assign_name_im1(self, amino_acids):
        '''
        @todo: check whether this method is used anywhere
        '''
        tmp = DEF.THREE_LET[:]
        amino_acids = deepcopy(amino_acids)
        tmp = map(lambda a: a.three_let, amino_acids)
        del tmp[1] # remove CY*
        if self.profile is not None:
            self.name_im1 = tmp[argmax(self.profile)]
        else:
            raise ValueError('Residue profile has not been set in residue with pasta no %i !' % self.pasta_no)
        
    def __str__(self):

        s = '#PASTA no.: %d' % self.pasta_no
        s += '\tPseudo energy: %d' % self.energy
        s += '\tSequence: %d' % self.seq_pos
        s += '\tSecondary: %s\n' % self.sec_struc
        
        cs = self.shifts_i
        for key in sorted(cs.keys()):
            s += '%s\t%s\t%.3f\t(%d)\n' % (self.name, key, 
                                           cs[key], self.original_no)

        cs = self.shifts_im1
        for key in sorted(cs.keys()):
            s += '%s\t%s\t%.3f\t(%d)\n' % (self.name_im1, key, 
                                           cs[key], self.original_no)

        s += self.comment_i
        s += self.comment_im1
        
        return s
        
class AminoAcid(Residue):
    
    def __init__(self, one_let='X', three_let='XXX'):
        super(AminoAcid, self).__init__()
        self.one_let = one_let
        self.three_let = three_let
        
    def calc_likelihood(self, shifts):
        return exp(sum([pdf(x, self.shifts[n][0], self.shifts[n][1], log=True) 
                        for n, x in shifts.items()]))

    def set_profile(self, param):
        '''
        @param param: a list of amino acids
        @type param: list
        '''

        p = zeros(len(param))
        i = [i for i, a in enumerate(param) if a.three_let == self.three_let]
        if not len(i) == 1: raise
        i = i[0]
        p[i] = 1.

        self.profile = p
    
    def get_carbons(self):
        
        carbon_shifts = []
        keys = array(self.shifts.keys())
        if len(keys) > 0:
            mask = array(map(lambda label: label.startswith('C'), keys))
            keys = keys[mask]
            carbon_shifts = array([self.shifts[key] for key in keys])
            carbon_shifts = carbon_shifts.T[0]
            
        return carbon_shifts, keys
    
if __name__ == '__main__':
    pasta_res = PastaResidue(23,
                             shifts_i={'CA':0},
                             shifts_im1={'CAi-1':0,'CBi-1':0})
    amino_acid = AminoAcid()
    amino_acid.set_attribute('shifts', {'CA':[0, 1], "H":[100,100]})
    print amino_acid.get_carbons()
    
