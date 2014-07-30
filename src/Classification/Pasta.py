from SequenceMapping import Context
from Posterior import Posterior
from numpy import array, mean, arange
from FileHandler import FileHandler
import datetime
import time

class Pasta(object):

    def __init__(self, pastafile, reffile, seqfile, verbose=False):
        """
        Constructor
        
        @param pastafile: filename for residue list
        @type pastafile: string
        @param reffile: filename for reference list (eg. BMRB ASCII file)
        @type reffile: string
        @param seqfile: filename for FASTA sequence file
        @type seqfile: string
        """
        
        fh = FileHandler()
        ## list of Residue.PastaResidue objects
        self.residues = fh.read_pasta(pastafile, reffile)
        ## list of Residue.AminoAcid objects 
        self.amino_acids = fh.read_preset(reffile)
        ## list of Residue.AminoAcid objects
        self.seq = fh.read_seq(seqfile, reffile)
        self.P = None ## numpy.ndarray for typing posterior probabilities
        self.L = None ## numpy.ndarray for linking constraints
        self.S = None ## numpy.ndarray for aa type in sequence
        self.A = None ## list of assignments and respective similarity score
        ## ILP STUFF
        self.B = None ## list of assignments and respective costs
        self.C = None ## cost matrix ILP
        self.Xs = None ## assignment matrices from solution pool,
        self.ILP_L = None ## Linking Matrix of ILP
        self.typing_elapsed = 0
        self.linking_elapsed = 0
        self.mapping_elapsed = 0
        self.full_running_time = 0
        
    def run(self, no_assignments=100, tolerance=.6, verbose=False, map_strat=None):
        """
        A run includes three steps:
            - Typing: Amino acid recognition
            - Linking: Estimation of potential direct neighbors for each residue
            - Sequence Mapping: Assignment of residues to sequence positions
        After a full run a number of assignments with a similarity score is
        returned, which is represented as a list of tuples. Each tuple consists 
        of a list of Residue to sequence position tuples and a similarity score.
        
        Example:
        assignments = [([(0,1), (1,4)], 0.99),([(0,2),(1,3)], 0.89)]
        This list consists of two possible assignments.
        The first resulted in an similarity score of 99%, in which 
        Residue 0 has been assigned to sequence position 1 and Residue 1
        has been assigned to sequence position 4.
            
        @see: SequenceMapping.Context for valid map_strat options
            
        @param no_assignments: Number of assignments to be generated
        @type no_assignments: int
        @param verbose: verbosity option
        @type verbose: bool
        @param map_strat: Mapping strategy
        @type map_strat: string
        
        @return: list of possible assignments with 
        respective similarity scores
        @rtype: list
        """
        
        import pylab as pl
        
        print '**** Typing Running ****'
        self.typing(type="joint")
        print '**** Typing Complete! ****'
        print '**** Time Elapsed: %s****\n' % str(datetime.timedelta(seconds=self.typing_elapsed))
        print '**** Linking Running ****'
        self.linking(strategy=map_strat, tolerance=tolerance)
        print '**** Linking Complete! ****'
        print '**** Time Elapsed: %s ****\n' % str(datetime.timedelta(seconds=self.linking_elapsed))
        print '**** Sequence Mapping Running ****'
        assignments = self.sequence_mapping(no_assignments=no_assignments, verbose=verbose, strategy=map_strat)
        print '**** Mapping Completed ****'
        print '**** Time Elapsed: %s ****\n' % str(datetime.timedelta(seconds=self.mapping_elapsed))
        proc_times = array([self.typing_elapsed,
                            self.linking_elapsed,
                            self.mapping_elapsed])
        print '****************** PROCESSING TIME ******************'
        print 'Processing Times for %i Residues' % len(self.residues)
        print '-' * 20
        print 'Typing: %s' % str(datetime.timedelta(seconds=self.typing_elapsed))
        print 'Linking: %s' % str(datetime.timedelta(seconds=self.linking_elapsed))
        print 'Mapping: %s' % str(datetime.timedelta(seconds=self.mapping_elapsed))
        print '-' * 20
        self.full_running_time = proc_times.sum()
        print 'Full Run: %s' % str(datetime.timedelta(seconds=self.full_running_time))
        print '**** DONE ****'
        
        return assignments
        
    def typing(self, cache_posterior=True, type=None):
        """
        Runs an amino acid recognition on Residue.PastaResidues, and returns
        the result in the form of posterior probabilities in a 
        mxnxn numpy.ndarray. Where m is the number of residues and n is the
        number of amino acids in the sequence.
        
        @attention: As a side effect, sets posterior attribute of this PASTA 
        instance
        
        @param cache_posterior: option to store posterior object as PASTA state
        @type cache_posterior: bool
        
        @return: mxnxn array of posterior probabilities
        @rtype: numpy.ndarray
        """
        t = time.clock()
        posterior_matrix = Posterior(self.residues, self.amino_acids,
                                     self.seq, type=type)
        self.prob = posterior_matrix(summarize=mean)
        if cache_posterior: self.posterior = posterior_matrix
        
        self.__set_res_profiles()
        self.P = posterior_matrix()
        self.__assign_res_names()
        posteriors = self.P
        self.typing_elapsed = (time.clock() - t)
        return posteriors
        
    def linking(self, strategy, tolerance=.6):
        """
        Runs linking on chemical shifts and returns a numpy.ndarray
        which is 1 if residue at current row is the predecessor of
        residue at current column.
        The tolerance value defines how large the deviation of the shifts
        between the chemical shifts of two residues can be, such that they
        are considered direct neighbors in the protein sequence.
        This tolerance value strongly depends on the underlying
        NMR experiment, that generates the chemical shift data.

        @attention: Ask your local NMR specialist about the optimal value!
        @attention: As a side effect, sets linking attribute of this PASTA instance
        @see: Linking.linking_matrix for information how it is filled
        
        @param tolerance: maximum chemical shift deviation for residue linking
        @type tolerance: float
        
        @return: matrix with boolean linking constraints
        @rtype: numpy.ndarray
        """
        from Linking import Linking
#        from Visualization.MatrixPlot import mat_fig
        
        if strategy == "ILP" or strategy == "CplexILP" or "None":
            conservative = True
        if strategy == "Joint":
            conservative = False
        
        t = time.clock()
        link = Linking(self.residues)
        L = link.linking_matrix(strategy, tolerance=tolerance, conservative=conservative)
        self.L = L
        linking_matrix = self.L
        self.linking_elapsed = (time.clock() - t)
        return linking_matrix
        
    def sequence_mapping(self, no_assignments, verbose=False, strategy=None):
        """
        Runs sequence mapping for residues and template sequence and returns
        a number of different assignment depending on certain mapping strategy.
        
        @see: SequenceMapping.Context for valid strategy parameters
        @see: SequenceMapping.generate_assignments for return value of this 
        method
        
        @param no_assignments: Number of assignments to be generated by mapping
        strategy
        @type no_assignments: int
        @param verbose: verbosity option
        @type verbose: bool
        @param strategy: Mapping strategy
        @type strategy: string
        
        @return: list of residue to sequence assignments
        @rtype: list
        """
        from operator import itemgetter
        
        t = time.clock()
        context = Context(self, no_assignments, verbose, strategy=strategy)
        assignments = context.execute()
        assignments = sorted(assignments, key=itemgetter(1), reverse=True)
        self.A = assignments
        self.mapping_elapsed = (time.clock() - t)
        return assignments

    def __assign_res_names(self):
        from numpy import unravel_index, max
        
        P = self.P
        if P is None:
            raise ValueError("Typing has not run yet!")
            
        aas = self.amino_acids
        aa_names = map(lambda a: a.three_let, aas)
        k = aa_names.index("CY*")
        del aa_names[k]
        
        for i, p in enumerate(P):
            res = self.residues[i]
            k, l = unravel_index(p.argmax(), p.shape)
            if sum(p[k, :] == max(p)) == 1:
                res.set_attribute("name", aa_names[l])
            if sum(p[:, l] == max(p)) == 1:
                res.set_attribute("name_im1", aa_names[k])

#    def __assign_res_names(self):
#        """Assigns residue names after typing has been run on PASTA instance."""
#        for r in self.residues:
#            r.assign_name(self.amino_acids)
        
    def __set_res_profiles(self):
        """
        Sets profile attribute in Residue with posterior probabilities 
        for each amino acid index residue list of this PASTA instance
        
        @raise ValueError: If Typing has not been run yet
        """
        if self.posterior is not None:
            posterior = self.posterior().sum(1)
            i = 0
            for r in self.residues:
                r.set_profile(posterior[i])
                i += 1
        else:
            raise ValueError("Typing has not run yet!")
    
    def violations(self, assignment):
        """
        @param A: assignment
        @type: list
        @rtype: int
        
        a := sum(dot(x, A) > dot(L.T, x))
        b := sum(dot(x, A.T) > dot(L, x))
        #violations = a+b
        
        """
        
        from numpy import dot, eye, zeros, ones
        import pylab as pl
        
        n = len(self.residues)
        X = zeros((len(self.seq), len(self.seq)))
        l = self.L
        L = ones((len(self.seq), len(self.seq)))
        L[:n, :n] = l
        L *= (1 - eye(len(L)))
        A = eye(len(L), k= -1)
        for i, j in assignment[0]:
            X[i][j] = 1.
        
        v1 = sum(dot(X, A) > dot(L.T, X))
        v2 = sum(dot(X, A.T) > dot(L, X))
        
        return sum(v1) + sum(v2)
    
    def assignment_score(self, assignment):
        return assignment[1]

if __name__ == '__main__':
    import pickle
    import pylab as pl
    from numpy import linspace
    fh = FileHandler()
    statsfile = 'shiftPresets/bmrb.shift'

######################## BMRB DATASETS ########################
#    from Visualization.MatrixPlot import mat_fig
#    from Visualization.FragGraph import plot_fragment_graph
#    from networkx import from_numpy_matrix
#    from FragmentGraph import FragmentGraph

    ## Ubiquitin (6457) 76 RESIDUES
    name = "Ubiquitin"

#    ### Sgs1 (4445) 91 RESIDUES
#    name = "Sgs1"

    ### MOZ (18142) 110 RESIDUES
#    name = "MOZ"
    
    ### CCL2 (17890) 153 RESIDUES
#    name = "CCL2"

#    ### SUPEROXIDE DISMUTASE (4341) 192 RESIDUES
#    name = "SuperoxideDismutase"

#    ### HUMAN PRION PROTEIN (4402) 210 RESIDUES
#    name = "HumanPrionProtein"

    ## Thiamine (15063) 224 RESIDUES
#    name = "Thiamine"

#    ### ENZYME 1 (4106) 259 RESIDUES
#    name = "EIN"

    seqfile = "Datasets/Ubiquitin/sequences/Ub_bmrb.fasta"
    pastafile = "Datasets/Ubiquitin/residue_lists/Ub_bmrb_unassigned.pasta"
    pasta = Pasta(pastafile, statsfile, seqfile)
    pasta.run()

    
    
    
    
#    result_folder = "/agbs/cluster/jhooge/pasta/%s/" % name
#    result_folder = "Datasets/%s/Results/"%name
#    outfile = result_folder + "%s_Gauss.p" %name
#    
#    print "File ", pastafile
#    print "Out ", outfile
#    pasta = Pasta(pastafile, statsfile, seqfile)
#    print "RUNNING: ", pastafile
#    pasta.run(500, verbose=False, map_strat="CplexILP")

    
################################################################

# MISSING SHIFTS
#    percentages = arange(0, .40, .05)
#    seqfile = "Datasets/Ubiquitin/sequences/Ub_bmrb.fasta"
#    folder = "Datasets/Ubiquitin/residue_lists/incompleteData/"
    
#    result_folder = "/agbs/cluster/jhooge/pasta/UbqAnalysisMaxLikelihood/"
#    pastafiles = [folder+"Ub_bmrb_missing_shifts_2_%.2f.pasta"%p for p in percentages]
##    dumpfiles = [result_folder + "Ub_bmrb_missing_shifts_Gauss_%.2fTLM.p"%p for p in percentages]
#    dumpfiles = [result_folder + "Ub_bmrb_missing_shifts_2_Gauss_%.2fTLM.p"%p for p in percentages]
#    
#    for pastafile in pastafiles:
#        i = pastafiles.index(pastafile)
#        print "RUNNING: ",pastafile
#        pasta = Pasta(pastafile, statsfile, seqfile)
#        pasta.run(500, verbose=False,map_strat="CplexILP")
#        print "DUMPED: ",dumpfiles[i]
#        pickle.dump(pasta, open(dumpfiles[i], "wb"))
    
## MISSING RESIDUES
#    percentages = arange(0, .55, .05)
#    seqfile = "Datasets/Ubiquitin/sequences/Ub_bmrb.fasta"
#    folder = "Datasets/Ubiquitin/residue_lists/incompleteData/"
#    result_folder = "/agbs/cluster/jhooge/pasta/UbqAnalysisMaxLikelihood/"
#    pastafiles = [folder + "Ub_bmrb_missing_residues_%.2f.pasta"%p for p in percentages]
##    dumpfiles = [result_folder + "Ub_bmrb_missing_residues_Gauss_%.2fTLM.p"%p for p in percentages]
#    dumpfiles = [result_folder + "Ub_bmrb_missing_residues_Laplace_%.2fTLM.p"%p for p in percentages]
#    
#    for pastafile in pastafiles:
#        i = pastafiles.index(pastafile)
#        print "RUNNING: ",pastafile
#        pasta = Pasta(pastafile, statsfile, seqfile)
#        pasta.run(500, verbose=False,map_strat="CplexILP")
#        print "DUMPED: ",dumpfiles[i]
#        pickle.dump(pasta, open(dumpfiles[i], "wb"))

# AMBIGUOUS DATA
#    percentages = arange(0, .55, .05)
#    seqfile = "Datasets/Ubiquitin/sequences/Ub_bmrb.fasta"
#    folder = "Datasets/Ubiquitin/residue_lists/ambiguousData/"
#    result_folder = "/agbs/cluster/jhooge/pasta/UbqAnalysisMaxLikelihood/"
#    pastafiles = [folder + "Ub_bmrb_ambiguous_residues_%.2f.pasta"%p for p in percentages]
##    dumpfiles = [result_folder + "Ub_bmrb_ambiguous_residues_Gauss_%.2fTLM.p"%p for p in percentages]
#    dumpfiles = [result_folder + "Ub_bmrb_ambiguous_residues_Laplace_%.2fTLM.p"%p for p in percentages]
#    
#    for pastafile in pastafiles:
#        i = pastafiles.index(pastafile)
#        print "RUNNING: ",pastafile
#        pasta = Pasta(pastafile, statsfile, seqfile)
#        pasta.run(500, verbose=False,map_strat="CplexILP")
#        print "DUMPED: ",dumpfiles[i]
#        pickle.dump(pasta, open(dumpfiles[i], "wb"))

## TOLERANCE TEST BAX
#    folder = "Datasets/Ubiquitin/"
#    pastafile = folder + "residue_lists/Ub_bmrb_unassigned.pasta"
#    seqfile = folder + "sequences/Ub_bmrb.fasta"
#    result_folder = "/agbs/cluster/jhooge/pasta/UbqAnalysisTolerance/"
#    
#    tolerances = linspace(0,2,41)
#    dumpfiles = [result_folder + "UbqBaxGauss_tol=%.2f.p"%t for t in tolerances]
##    dumpfiles = [result_folder + "UbqBaxLaplace_tol=%.2f.p"%t for t in tolerances]
#    strat = "CplexILP"
#    for i,tol in enumerate(tolerances):
#        print dumpfiles[i]
#        print "RUNNING WITH TOL=",tol
#        pasta = Pasta(pastafile, statsfile, seqfile)
#        pasta.typing()
#        pasta.linking(strat, tolerance=tol)
#        pasta.sequence_mapping(500, verbose=False, strategy=strat)
#        pickle.dump(pasta, open(dumpfiles[i], "wb"))

## TOLERANCE TEST TRUFFAULT
#    folder = "Datasets/Ubiquitin/"
#    pastafile = folder + "residue_lists/Ub_opt_unambiguous_unassigned.pasta"
#    seqfile = folder + "sequences/Ub_bmrb.fasta"
#    result_folder = "/agbs/cluster/jhooge/pasta/UbqAnalysisTolerance/"
#    
#    tolerances = linspace(2,0,41)
##    dumpfiles = [result_folder + "UbqMPIGauss_tol=%.2f.p"%t for t in tolerances]
#    dumpfiles = [result_folder + "UbqMPILaplace_tol=%.2f.p"%t for t in tolerances]
#    strat = "CplexILP"
#    for i,tol in enumerate(tolerances):
#        print dumpfiles[i]
#        print "RUNNING WITH TOL=",tol
#        pasta = Pasta(pastafile, statsfile, seqfile)
#        pasta.typing()
#        pasta.linking(strat, tolerance=tol)
#        pasta.sequence_mapping(500, verbose=False, strategy=strat)
#        pickle.dump(pasta, open(dumpfiles[i], "wb"))

## TYPING TEST 
#    noise = linspace(0,3,31)
#    seqfile = "Datasets/Ubiquitin/sequences/Ub_bmrb.fasta"
#    folder = "Datasets/Ubiquitin/residue_lists/NoisyData/"
#    result_folder = "/agbs/cluster/jhooge/pasta/UbqAnalysisNoise/"
#    pastafiles = [folder + "UbqBaxNoise=%.1f.pasta"%n for n in noise]
#    
##    dumpfiles = [result_folder + "UbqBaxGaussFullNoise=%.1f.p"%n for n in noise]
##    dumpfiles = [result_folder + "UbqBaxGaussCACBNoise=%.1f.p"%n for n in noise]
##    dumpfiles = [result_folder + "UbqBaxLaplaceFullNoise=%.1f.p"%n for n in noise]
#    dumpfiles = [result_folder + "UbqBaxLaplaceCACBNoise=%.1f.p"%n for n in noise]
#    
#    strat = "CplexILP"
#    for i,n in enumerate(noise):
#        print "RUNNING %s "%pastafiles[i]
#        pasta = Pasta(pastafiles[i], statsfile, seqfile)
#        atoms = ["CA","CB","CO"]
#        atoms_im1 = ["CAi-1","CBi-1","COi-1"]
#        for r in pasta.residues:
#            for key in r.shifts_i.keys():
#                if key not in atoms:
#                    del r.shifts_i[key]
#            for key in r.shifts_im1.keys():
#                if key not in atoms_im1:
#                    del r.shifts_im1[key]
#        pasta.typing()
#        pickle.dump(pasta, open(dumpfiles[i], "wb"))

#    result_folder = "/agbs/cluster/jhooge/pasta/"
#    i = 10
#    pastafile = dumpfiles[i]
#    print "RUNNING ",pastafile
#    pasta = pickle.load(open(pastafile,"r"))
#    pasta.sequence_mapping(no_assignments=1, verbose=False, strategy="CplexILP")
#    records = fh.to_fasta_records(pasta.A, pasta.seq, pasta.residues)
#    fh.write_fasta(result_folder+"Ubq_Bax_%i_clustermapped.fasta"%i, records, pasta.seq)
#    fh.conf_to_pasta(result_folder+"Ubq_Bax_%i_clustermapped.pasta"%i, pasta.A[0], 
#                     pasta.residues, pasta.amino_acids)
#    pickle.dump(pasta, open(result_folder+"Ubq_Bax_%i_clustermapped.p"%i, "wb"))
#    print 
#    print "File: ", pastafile
#    print "Score: ",pasta.assignment_score(pasta.A[0])
#    print "No of violations: ",pasta.violations(pasta.A[0])
    

