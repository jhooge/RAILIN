from numpy import mean
from Priors import AbstractPrior

class Posterior(object):
    """
    Callable class that returns a numpy.ndarray with dimensions
    mxnxn populated  with posterior probabilities from amino acid
    recognition.
    
    Summing up this object along axis 1 results in posterior scores 
    for amino acid types given the [i] shifts of the residues.
    Summing up this object along axis 2 results in posterior scores 
    for amino acid types given the [i-1] shifts of the residues.
    In both cases the resulting numpy.ndarray has mxn dimensions.
    Where m is the number or residues and n the number of amino acids
    labels.
    """
    def __init__(self, residues, amino_acids, seq, type=None):
        """
        @param residues: PastaResidue objects generated from pasta file
        @type residues: list
        @param amino_acids: Amino_acid objects generated from bmrb file
        @type amino_acids: list
        @param seq: Amino_acid objects generated from sequence
        @type seq: list
        """
        
        self.residues = residues
        self.amino_acids = amino_acids
        self.seq = seq
        self.prior = AbstractPrior.create(seq, amino_acids, type)
        
    def __call__(self, cache_likelihood=True, summarize=mean):
        
        from MaxLikelihood import Likelihood
        from numpy import take
        
        L = Likelihood()
        
#        A = L.calc_likelihoods(self.residues, self.amino_acids, False,summarize)
#        B = L.calc_likelihoods(self.residues, self.amino_acids, True, summarize)
        A = L.calc_likelihoods(self.residues, self.amino_acids, previous=False)
        B = L.calc_likelihoods(self.residues, self.amino_acids, previous=True)
        C = self.prior.evaluate()
        
        ## average over cysteines
        aa = [a.three_let for a in self.amino_acids]

        k = aa.index('CY*')
        l = aa.index('CYS')
    
        if False:
            weight = 0.5
            A[:, k] = (A[:, k] + A[:, l]) * weight
            B[:, k] = (B[:, k] + B[:, l]) * weight
        else:
            a,b = A[:, k],A[:, l]
            m = a > b
            A[:, k] = m * a + (1-m) * b
            a,b = B[:, k],B[:, l]
            m = a > b
            B[:, k] = m * a + (1-m) * b
                
        i = [i for i, a in enumerate(aa) if a != 'CYS']
    
        A = take(A, i, 1)
        B = take(B, i, 1)
        
        if cache_likelihood:
            self.likelihood = (A, B)
        return self.calc_posterior(A, B, C)
    
    def calc_posterior(self, L_curr, L_prev, prior):
        """
        Combine likelihoods for i and i-1 recognition with pair
        prior.
        """
        from Math import outer_2way, normalize_rows
        from numpy import multiply, reshape
    
        p = outer_2way(L_prev.T, L_curr.T, multiply)
        p = p.swapaxes(0, 1)
        p *= prior
            
        p = reshape(p, (p.shape[0], -1))
        normalize_rows(p)
        
        return reshape(p, (p.shape[0], prior.shape[0], prior.shape[1]))
    
if __name__ == "__main__":
    import pickle
    from Definitions import three2One
    from Pasta import Pasta
    import pylab as pl
    
    folder = "tests/"
    statsfile = folder+'reference_lists/bmrb.shift'
    pastafile = folder + "residue_lists/all_singles.pasta"
    seqfile = folder + "sequences/all_singles.fasta"
    pasta = Pasta(pastafile, statsfile, seqfile)
    P = pasta.typing(type="joint")
    
    pl.matshow(P.sum(1))
    pl.show()
