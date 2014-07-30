'''
Created on Nov 11, 2010

@author: jhooge
'''

ONE_LET = [ 'A', 'B', 'C', 'D', 'E', 'F', 'G',
            'H', 'I', 'K', 'L', 'M', 'N', 'P',
            'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
            'X' ]

THREE_LET = [ "ALA", "CY*", "CYS", "ASP", "GLU", "PHE", "GLY",
              "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO",
              "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
              "NAA" ]

ATOMS = ['CA', 'CB', 'CG', 'CD', 'CE', 'CZ']

def type2Index(list):
    '''
    returns a list of index tuple at which array position
    the type is stored
    '''
    indexes = []
        
    for k in list:
        if len(k) < 3:
            i = 0
        else:
            i = int(k[2])
        j = ATOMS.index(k[0:2])
        indexes.append((i, j))
            
    return indexes

def index2Type(indices):

    '''
    Takes a list of index tuples and returns a corresponding list
    of nucleus identifiers.
    '''
    types = []
    for i, j in indices:
        if i == 0: 
            types.append(ATOMS[j])
        else:
            types.append(ATOMS[j] + str(i))
        
    return types

def three_let2Index(list):
    indices = map(THREE_LET.index, list)
    return indices

def one_let2Index(list):
    indices = map(ONE_LET.index, list)
    return indices

def index2Three_let(indices):
    return [THREE_LET[i] for i in indices]

def index2One_let(indices):
    return [ONE_LET[i] for i in indices]

def three2One(three_let):
    three_let = three_let.upper()
    convert = dict(zip(THREE_LET, ONE_LET))
    return convert[three_let]

def one2Three(one_let):
    one_let = one_let.upper()
    convert = dict(zip(ONE_LET, THREE_LET))
    return convert[one_let]
