'''
Created on Jun 10, 2012

@author: jhooge
'''
if __name__ == '__main__':
    from FileHandler import FileHandler
    from Mapping import Mapping
    from Mapping2 import AbstractMapping
    from Definitions import three2One
    from MaxLikelihood import Likelihood
    from numpy import array,mean, max
    import pylab as pl
    
    fh = FileHandler()
    pairs = 'tests/residue_lists/all_pairs.pasta'
    singles = 'tests/residue_lists/all_singles.pasta'
    seqfile_pairs = 'tests/sequences/all_pairs.fasta'
    seqfile_singles = 'tests/sequences/all_singles.fasta'
    statsfile = 'tests/reference_lists/bmrb.shift'
    
    single_residues = fh.read_pasta(singles, statsfile)
    amino_acids = fh.read_seq(seqfile_singles, statsfile)
    atoms = {'CO' :0, 'CA' :1,
             'CB' :2, 'CG' :3,
             'CG1':4, 'CG2':5,
             'CD' :6, 'CD1':7,
             'CD2':8, 'CE' :9,
             'CE1':10, 'CE2':11,
             'CE3':12, 'CZ':13, 'CZ2':14,
             'CZ3':15, "CH2":16}
    aa_names = ''.join([aa.one_let for aa in amino_acids])
#    k = aa_names.index("F")
#    print "Residue: ",res.name_im1, res.get_carbons(previous=True)
#    print "Amino Acid: ",aa.three_let, aa.get_carbons()
    
    mapper1 = Mapping(aa_mapping=True)
    L = Likelihood()
#    mapping1 = mapper1.generate_mappings(res, aa, previous=False)
#    mapper2 = AbstractMapping.create("on_amino_acid")
#    mapping2 = mapper2.doMunkresMap(res, aa, previous=True)
    
    expected_mappings = {}
    for res in single_residues:
        for aa in amino_acids:
            key = ''.join([three2One(res.name), aa.one_let])
            expected_mappings[key] = mapper1.generate_mappings(res, aa, previous=False)
            mean_L = L.calc_likelihoods(single_residues, amino_acids,summarize=mean)
            max_L = L.calc_likelihoods(single_residues, amino_acids,summarize=max)
            mean_L = L.calc_likelihoods(single_residues, amino_acids)
            max_L = L.calc_likelihoods(single_residues, amino_acids)
            
    
#    for k in expected_mappings.keys():
#        print k, expected_mappings[k]

#    mean_L = L.calc_likelihoods(single_residues, amino_acids,summarize=mean)
#    for k, aa in enumerate(amino_acids):
#        res = single_residues[k]
#        aa = amino_acids[k]
#        mappings = expected_mappings[aa.one_let*2]
##        print res.name
#        test_mean_L = []
#        for m in mappings:
#            test_mean_L.append(L.calc_likelihood(res, aa, m))
#        print mean(array(test_mean_L)) == mean_L[k][k]
    
#    pl.matshow(mean_L)
#    pl.show()
