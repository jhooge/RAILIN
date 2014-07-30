'''
Created on Aug 27, 2010

@author: jhooge
'''
import Definitions as DEF
import numpy as np
from FileHandler import FileHandler

from abc import ABCMeta, abstractmethod
class AbstractPrior(object):
    """
    This abstract class represents Priors for the Bayesian Model 
    that computes the posterior probability for the typing step.
    """
    
    __metaclass__ = ABCMeta
    
    def init(self, seq, amino_acids):
        """
        Formats sequence and class labels to one letter codes
        and returns it as a list of parameters, namely a list of
        amino acid one letter codes in sequence and amino acid one letter
        codes for sequence alphabet (class labels)
        
        @attention: Class labels from BMRB file include Cy* (B), this
        label will be removed as it can not occur in a template sequence.
        In a later step when the posterior will be computed, it will be
        averaged over the two oxidization states of Cystein.
        
        @param seq: List of Residue.AminoAcid objects in sequence
        @type seq: list
        @param amino_acids: List of Residue.AminoAcid objects in Reference file
        @type amino_acids: list
        
        @return: List of model parameters consisting of sequence 
        @rtype: list
        """
        seq_labels = [aa.one_let for aa in seq]
        aa_labels = [aa.one_let for aa in amino_acids]
        
        # remove Cy* (B) as this code never occurs in a sequence
        k = aa_labels.index('B')
        del aa_labels[k]
        
        return [seq_labels, aa_labels]
    
    @abstractmethod
    def evaluate(self):
        """Abstract method to evaluate the prior"""
        pass
    
    @staticmethod
    def create(seq, amino_acids, type=None):
        """
        Factory method to create priors
        Valid options:
         - indy: for Independent AbstractPrior
         - joint for Joint AbstractPrior (default)
        
        @param seq: List of Residue.AminoAcids in sequence
        @type seq: list
        @param amino_acids: List of Residue.AminoAcids in sequence alphabet
        @type amino_acids: list
        @param type: AbstractPrior type to create
        @type type: string
        
        @return: Prior of given type
        @rtype: Prior
        """
        if type == None:
            type = "joint"
        if type == "indy":
            return IndependentPrior(seq, amino_acids)
        if type == "joint":
            return JointPrior(seq, amino_acids)
        if type=="uniform":
            return UniformPrior(seq, amino_acids)
        
class IndependentPrior(AbstractPrior):
    """
    Concrete Prior class which represents independent pair frequencies 
    of amino acids in template sequence.
    Evaluation of this prior results in the product of the single frequencies
    of a pair of amino acids in the template sequence
    """
    
    def __init__(self, seq, amino_acids):
        """
        Constructor
        
        @param seq: List of sequence one letter labels
        @type seq: list
        @param class_labels: List of amino acid class labels in one letter codes
        @type class_labels: list
        """
        params = super(IndependentPrior, self).init(seq, amino_acids)
        self.seq_labels = params[0] ## sequence
        self.aa_labels = params[1] ## class labels

    def evaluate(self):
        """
        Computes all pair frequencies of aas in given sequence.
        The pair frequency is the product of the single priors
        of the amino acids in a pair.
        
        @return: Normalized frequency of (independent) amino acid pairs in sequence
        @rtype: numpy.ndarray
        """
        aas = [self.seq_labels[i] for i in range(len(self.seq_labels))] # all aas occurring in sequence

        f = np.zeros(len(self.aa_labels))
        for i in range(f.shape[0]):
            f[i] = aas.count(self.aa_labels[i])
        f = f / f.sum()

        # first residue has unknown predecessor
        f = np.multiply.outer(f, f)
        return f

class JointPrior(AbstractPrior):
    """
    Concrete Prior class which represents joint frequencies 
    of amino acids in template sequence.
    Evaluation of this prior results in joint frequencies
    of a pair of amino acids in the template sequence
    """
    def __init__(self, seq, amino_acids):
        """
        Constructor
        
        @param seq: List of sequence one letter labels
        @type seq: list
        @param class_labels: List of amino acid class labels in one letter codes
        @type class_labels: list
        """
        params = super(JointPrior, self).init(seq, amino_acids)
        self.seq_labels = params[0] ## sequence
        self.aa_labels = params[1] ## class labels

    def evaluate(self):
        """
        Computes all pair frequencies of aas in given sequence.
        returns f(i-1,i)/f.sum()

        @return: Normalized frequency of amino acid pairs in sequence
        @rtype: numpy.ndarray
        """
        dipeptides = zip(self.seq_labels[:-1], self.seq_labels[1:])
        f = np.zeros((len(self.aa_labels), len(self.aa_labels)))

        for i in range(f.shape[0]):
            for j in range(f.shape[1]):
                f[i, j] += dipeptides.count((self.aa_labels[i], self.aa_labels[j]))

        # first residue has unknown predecessor
        j = self.aa_labels.index(self.seq_labels[0])
        f[:, j] += 1 / float(len(self.aa_labels))

        return f / f.sum()

class UniformPrior(AbstractPrior):
    """
    Concrete Prior class which represents a uniform prior
    """
    def __init__(self, seq, amino_acids):
        """
        Constructor
        
        @param seq: List of sequence one letter labels
        @type seq: list
        @param class_labels: List of amino acid class labels in one letter codes
        @type class_labels: list
        """
        params = super(UniformPrior, self).init(seq, amino_acids)
        self.seq_labels = params[0] ## sequence
        self.aa_labels = params[1] ## class labels

    def evaluate(self):
        """
        Computes uniform prior

        @return: Normalized frequency of amino acid pairs in sequence
        @rtype: numpy.ndarray
        """
        f = np.ones((len(self.aa_labels), len(self.aa_labels)))
        return f / f.sum()

if __name__ == "__main__":
    fh = FileHandler()
    
    from Pasta import Pasta
    import pylab as pl
    
    folder = "tests/"
    statsfile = folder+'reference_lists/bmrb.shift'
#    pastafile = folder+"residue_lists/all_singles.pasta"
#    seqfile = folder+"sequences/all_singles.fasta"
#    seq = fh.read_seq(seqfile, statsfile)
#    amino_acids = fh.read_preset(statsfile)

    folder = "Datasets/Ubiquitin/"
    pastafile = folder+"Ubiquitin_newVT.pasta"
    seqfile = folder+"Ub.fasta"
    result_folder = folder+"Results/"
    
    pasta = Pasta(pastafile, statsfile, seqfile, verbose=True)
    P1 = pasta.typing(type="indy")
    P2 = pasta.typing()
    P1_i = P1.sum(1)
    P1_im1 = P1.sum(2)
    P2_i = P2.sum(1)
    P2_im1 = P2.sum(2)
    
    pl.matshow(P1_i)
    pl.matshow(P1_im1)
    pl.matshow(P2_i)
    pl.matshow(P2_im1)
    pl.show()
#    single = AbstractPrior.create(seq, amino_acids, "single")
#    indy = AbstractPrior.create(seq, amino_acids, "indy")
#    pair = AbstractPrior.create(seq, amino_acids, "joint")
#    priors = [single,indy, pair]
    
