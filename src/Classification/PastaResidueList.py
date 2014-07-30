'''
Created on Jun 28, 2012

@author: jhooge
'''

from FileHandler import FileHandler
from copy import deepcopy
from numpy import random, linspace, arange
from random import choice

def remove_rnd_residues(residues, percentage):
    n = 76
    m = float(len(residues))
    while 1. - m / n < percentage:
        print n,m
        j = random.randint(len(residues))
#        print "Removing Res ",residues[j].pasta_no
        del residues[j]
        m = float(len(residues))
    return residues

def add_noise(data, noise):
    """
    data
    noise in ppm
    """
    from numpy.random import uniform
    i = 0
    while i < len(data): ## add noise to i
        
        for key in data[i].shifts_i.keys():
            shift = data[i].shifts_i[key]
            if noise != 0.:
                a = shift - noise
                b = shift + noise
                data[i].shifts_i[key] = uniform(a, b)
                
#                print '%.2f <= %.2f <= %.2f'%(a,data[i][key],b)
                assert data[i].shifts_i[key] >= a and data[i].shifts_i[key] <= b
        for key in data[i].shifts_im1.keys():
            shift = data[i].shifts_im1[key]
            if noise != 0.:
                a = shift - noise
                b = shift + noise
                data[i].shifts_im1[key] = uniform(a, b)
                
#                print '%.2f <= %.2f <= %.2f'%(a,data[i][key],b)
                assert data[i].shifts_im1[key] >= a and data[i].shifts_im1[key] <= b

        i += 1
        
    return data
    
def get_no_of_carbon_shifts(residues):
    n = 0
    for r in residues:
        shifts_i, keys_i = r.get_carbons(previous=False)
        shifts_im1, keys_im1 = r.get_carbons(previous=True)
        
        n += len(shifts_i) + len(shifts_im1)
    return n

def get_no_of_ambiguous_keys(residues):
    n = 0
    ambiguous_keys = ["C1", "C2", "C3",
                      "C4", "C5", "C6",
                      "C7"]
    for r in residues:
        shifts_i, keys_i = r.get_carbons(previous=False)
        shifts_im1, keys_im1 = r.get_carbons(previous=True)
        
        for key in keys_i:
            if key in ambiguous_keys:
                n += 1
        for key in keys_im1:
            key = key.replace("i-1", "")
            if key in ambiguous_keys:
                n += 1
        
    return n

def get_non_empty_indices(residues, previous):
    indices = []
    
    for i,r in enumerate(residues):
        shifts, atoms = r.get_carbons(previous)
        if len(atoms) > 1:
            indices.append(i)
    return indices

def all_empty(residues):
    ind_i = get_non_empty_indices(residues, False)
    ind_im1 = get_non_empty_indices(residues, True)
    empty = len(ind_i) == 0 and len(ind_im1) == 0
    
    return empty

def remove_rnd_shifts(residues, percentage):
    from random import choice
    n = float(get_no_of_carbon_shifts(residues))
    m = float(get_no_of_carbon_shifts(residues))
    while not all_empty(residues) and 1.-m/n < percentage:
        previous = choice([True,False])
        ind = get_non_empty_indices(residues, previous)
        if len(ind) > 0:
            r = get_rnd_residue(residues, ind)
            shifts, keys = r.get_carbons(previous=previous)
            j = random.randint(0, len(keys))
            if previous:
                print "in ",r.pasta_no
                print keys[j],"i-1 deleted"
                del r.shifts_im1[keys[j]]
            else: 
                print "in ",r.pasta_no
                print keys[j]," deleted"
                del r.shifts_i[keys[j]]
        m = float(get_no_of_carbon_shifts(residues))
    return residues
    
def get_indices(residues, previous):
    from numpy import array
    indices = []
    ambiguous = ["C1", "C2", "C3",
                "C4", "C5", "C6",
                 "C7", 
                 "C1i-1", "C2i-1", "C3i-1",
                 "C4i-1", "C5i-1", "C6i-1",
                 "C7i-1"]
    if not previous:
        for i,r in enumerate(residues):
            shifts, atoms = r.get_carbons(previous)
            if array(map(lambda a: a not in ambiguous, atoms)).any():
                indices.append(i)
    if previous:
        for i,r in enumerate(residues):
            shifts, atoms = r.get_carbons(previous)
            if array(map(lambda a: a not in ambiguous, atoms)).any():
                indices.append(i)
    return indices

def get_rnd_residue(residues, ind):
    from random import choice
    return residues[choice(ind)]

def swap_rnd_atom(residue, previous):
    from random import choice
    ambiguous = ["C1", "C2", "C3",
                "C4", "C5", "C6",
                 "C7", 
                 "C1i-1", "C2i-1", "C3i-1",
                 "C4i-1", "C5i-1", "C6i-1",
                 "C7i-1"]
    k = 1
    shifts, atoms = residue.get_carbons(previous)
    non_ambiguous_atoms = list(set(atoms).difference(ambiguous))
    if not previous:
        old_atom = choice(non_ambiguous_atoms)
        new_atom = "C" + str(k)
        if new_atom not in atoms:
            new_atom = "C" + str(k)
#            print "Swapping in ", residue.pasta_no
#            print new_atom, old_atom
            residue.shifts_i[new_atom] = residue.shifts_i.pop(old_atom)
        else:
            while new_atom in atoms and k <= len(atoms):
                k += 1
                new_atom = "C" + str(k)
#            print "Swapping in ", residue.pasta_no
#            print new_atom, old_atom
            residue.shifts_i[new_atom] = residue.shifts_i.pop(old_atom)
    else:
        old_atom = choice(non_ambiguous_atoms)
        new_atom = "C" + str(k)+"i-1"
        if new_atom not in atoms:
            new_atom = "C" + str(k)+"i-1"
#            print "Swapping in ", residue.pasta_no
#            print new_atom, old_atom
            residue.shifts_im1[new_atom] = residue.shifts_im1.pop(old_atom)
        else:
            while new_atom in atoms and k <= len(atoms):
                k += 1
                new_atom = "C" + str(k)+"i-1" 
#            print "Swapping in ", residue.pasta_no
#            print new_atom, old_atom
            residue.shifts_im1[new_atom] = residue.shifts_im1.pop(old_atom)

def make_ambiguous(residues, percentage):
    from random import choice
    n = float(get_no_of_carbon_shifts(residues))
    m = float(get_no_of_ambiguous_keys(residues))

    while m/n < percentage:
        print m/n
        previous = choice([True,False])
        i = get_indices(residues, previous)
        res= get_rnd_residue(residues,i)
        swap_rnd_atom(res,previous)
        i = get_indices(residues, previous)
        m = float(get_no_of_ambiguous_keys(residues))
    return residues

if __name__ == '__main__':
    fh = FileHandler()
    
    folder = "Datasets/Ubiquitin/"
    pastafn = "residue_lists/Ub_bmrb_unassigned.pasta"
    pastafn2 = "residue_lists/Ub_opt_unambiguous_unassigned.pasta"
    presetfn = 'shiftPresets/bmrb.shift'
    seqfn = folder + "sequences/Ub_bmrb.fasta"
    result_folder = folder + "Results/"
    residues = fh.read_pasta(folder + pastafn, presetfn)
    residues2 = fh.read_pasta(folder + pastafn2, presetfn)

## ADD NOISE
#    noise = linspace(0,3,31)
#    for n in noise:
#        print residues[0].shifts_i.values()
#        new_res = add_noise(deepcopy(residues), n)
#        outfolder = folder + "UbqBaxNoise=%.1f.pasta" %n
#        fh.write_pasta(outfolder, new_res)
#        print outfolder, "written"

## INTRODUCE AMBIGUOUS SHITS
#    old_res = deepcopy(residues)
#    for p in arange(.0, .55, .05):
#        new_res = make_ambiguous(old_res,p)
#        old_res = deepcopy(new_res)
#        outfolder = folder + pastafn[:-17] + "_ambiguous_residues_%.2f.pasta" % p
#        fh.write_pasta(outfolder, new_res)
#        print outfolder, "written"

## TEST CORRECTNESS OF GENERATED LISTS
#    folder = "Datasets/Ubiquitin/residue_lists/ambiguousData/"
#    for p in arange(.0, .55, .05):
#        pastafn = "Ub_bmrb_ambiguous_residues_%.2f.pasta"% p
#        print folder+pastafn
#        residues = fh.read_pasta(folder+pastafn, presetfn)
#        for res in residues:
#            res.pasta_no
#        m = get_no_of_ambiguous_keys(residues)
#        n = get_no_of_carbon_shifts(residues)
#        print m,n,float(m)/n

# REMOVE RESIDUES
#    old_res = deepcopy(residues)
#    for p in arange(.0,.55,.05):
#        new_res = remove_rnd_residues(old_res,p)
#        old_res = deepcopy(new_res)
#        outfolder = folder+pastafn[:-17]+"_missing_residues_%.2f.pasta"%p
#        fh.write_pasta(outfolder, new_res)
#        print outfolder+" written"

## REMOVE SHIFTS
#    print arange(.025,.20,.05)
#    old_res = deepcopy(residues)
#    for p in arange(.025,.20,.05):
#        new_res = remove_rnd_shifts(old_res, p)
#        old_res = deepcopy(new_res) 
#        outfolder = folder + pastafn[:-17] + "_missing_shifts_2_%.2f.pasta" % p
#        fh.write_pasta(outfolder, new_res)
#        print outfolder + " written"
#        assert len(new_res) == len(residues)
