'''
Created on Feb 15, 2012

@author: jhooge
'''
import numpy as np
import pickle
import pylab as plt
import Definitions as DEF
from numpy.core.numeric import zeros_like

def plot_conf_scores(conf_space, to_file=None, dpi=300):
    
    conf_scores = 1. - np.array([c[1] for c in conf_space])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    n, bins, patches = ax.hist(conf_scores, 10, normed=0, alpha=.3)
    
    ax.set_xlabel("Sequence Similarity")
    ax.set_ylabel("Number of Configurations")
    '''TODO: Legend will be clipped if anchored'''
    ax.set_yticks(np.arange(0, max(n) + 1))
    
    if to_file != None:
        filename = to_file.split('/')[-1]
        fig.set_size_inches(18.5,10.5)
        fig.savefig(filename + "_confspace", dpi=dpi)
    else:
        return fig

# a stacked bar plot with errorbars

def conf_matrix(conf_space, n):
    from Classification.Pasta import read_seq, read_preset
    '''
    m = length of sequence
    '''
    M = []
    for conf, error in conf_space:
        line = list('*'*n)
        for frag, i in conf:
            line[i:i + len(frag)] = [DEF.three2One(r.name) for r in frag.residues]
        M.append(line)
    return M

def transpose_2Dlist(l):
    l_T = zip(*l)
    l_T = map(list, l_T)
    return l_T
    
def count_labels(l, aa_labels):
    counts = np.zeros(len(aa_labels))
    i = 0
    for label in aa_labels:
        counts[i] = l.count(label)
        i += 1
    return counts

def frequency_matrix(conf_space, n, aa_labels):
    '''
    cols := sequence position
    rows := aa_labels
    '''
    F = []
    
    M = conf_matrix(conf_space, n)
    for row in M:
        print ''.join(row)
#    print 'dim(M)=(%i,%i)' % (len(M), len(M[0]))
    M_T = transpose_2Dlist(M)
#    print 'dim(M_T)=(%i,%i)' % (len(M_T), len(M_T[0]))
    for row in M_T:
#        print ''.join(row)
        F.append(count_labels(row, aa_labels).astype(int))
#    print 'dim(F)=(%i,%i)' % (len(F), len(F[0]))
    F_T = transpose_2Dlist(F)
#    print 'dim(F_T)=(%i,%i)' % (len(F_T), len(F_T[0]))
    
    ## normalize
    F_T = np.array(F_T)
    
    return F_T
    

def stacked_bar_plot(conf_space, seqfile, statsfile, to_file=None, dpi=300):
    from Pasta import read_seq, read_preset
    
    '''
    n := sequence_length
    '''
    amino_acids = read_preset(statsfile)
    seq = read_seq(seqfile)
    seq_labels = [aa.one_let for aa in seq]
    aa_labels = [aa.one_let for aa in amino_acids]
    del aa_labels[1] ## remove Cy*
    seq_len = len(seq_labels)
    
    print '\nConfiguration Space'
    F = frequency_matrix(conf_space, seq_len, aa_labels)
    
    ind = np.arange(seq_len)    # the x locations for the groups
    width = 1.       # the width of the bars: can also be len(x) sequence
    
    colors = ['#0000FF',
          '#00FF00',
          '#FFFF00',
          '#FFA500',
          '#FF0000',
          '#A0522D',
          '#EE82EE',
          '#A020F0',
          '#191970',
          '#6495ED',
          '#7FFFD4',
          '#556B2F',
          '#FFD700',
          '#B22222',
          '#CDB79E',
          '#C1CDC1',
          '#8B7D7B',
          '#7A8B8B',
          '#008B45',
          '#8B864E']
    
    fig = plt.figure()

    n_per_subplot = 50
    n = len(F[0]) % n_per_subplot ## rest
    m = (len(F[0]) - n) / n_per_subplot + 1 ## number of subplots
    
    i = 0
    k = 1
    while k <= m and i < len(seq):
        j = i + n_per_subplot
        if j >= len(seq):
            j = i + (len(seq) - i)
        ax = fig.add_subplot(m, 1, k)
        l = 0
        for row, color in zip(F, colors):
            ax.bar(ind[i:j], row[i:j], width, color=color, label=aa_labels[l])
            l += 1
        
        ax.set_ylabel("%i-%i" % (i, j))
        ax.set_xticks(ind[i:j] + width / 2.)
        ax.set_xticklabels(seq_labels[i:j])
        ax.set_xlim(ind[i], ind[j - 1] + (width))
        
        i = j + 1
        k += 1

        
    # loop through all patch objects and collect ones at same x
    all_patches = []
    for ax in fig.get_axes():
        all_patches.extend(ax.patches)
    patch_at_x = {}
    for patch in all_patches:
        if patch.get_x() not in patch_at_x: patch_at_x[patch.get_x()] = []
        patch_at_x[patch.get_x()].append(patch)
    
    # custom sort function, in reverse order of height
    def yHeightSort(i, j):
        if j.get_height() > i.get_height(): return 1
        else: return - 1
    
    # loop through sort assign z-order based on sort
    for x_pos, patches in patch_at_x.iteritems():
        if len(patches) == 1: continue
        patches.sort(cmp=yHeightSort)
        [patch.set_zorder(patches.index(patch)) for patch in patches]
        
    fig.get_axes()[0].legend(bbox_to_anchor=(1.01, 0, 1, 1), loc=2, borderaxespad=0.)
    fig.get_axes()[0].set_title(to_file + '\nConfiguration Space Alignment')
    fig.get_axes()[-1].set_xlabel('Sequence')
    
    if to_file != None:
        filename = to_file.split('/')[-1]
        fig.set_size_inches(18.5,10.5)
        fig.savefig(filename + '_confalign', dpi=dpi)
    else:
        return fig

if __name__ == '__main__':
    
    conf_space = pickle.load(open('/kyb/agbs/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/multiple_test_files/31C/31_CTD_conf_space.p'))
    seqfile = '/kyb/agbs/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/multiple_test_files/31C/31_CTD.fasta'
    statsfile = '/kyb/agbs/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/shiftPresets/bmrb_ascii.shift'
    
    n = len(read_seq(seqfile))
#    print n
    amino_acids = read_preset(statsfile)
    aa_labels = [aa.one_let for aa in amino_acids]
    
#    F = frequency_matrix(conf_space, n, aa_labels)
#    for row in F:
#        print row
        
    stacked_bar_plot(conf_space, seqfile, statsfile, to_file='FOOBAR')
        
#    print list(list('*'*10)*5)
#    print ['A','A','B'].count(['A'])

