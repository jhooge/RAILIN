'''
Created on Dec 21, 2011

@author: jhooge
'''
from Math import gaussian_pdf
from matplotlib.figure import Figure
from numpy import mean, std, array, arange, sqrt, exp, linspace
from numpy.oldnumeric.random_array import random
from scipy import special
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pickle

def plot_mapping(conf, seqfile, statsfile, n=50, to_file=None,dpi=300):
    
    from Classification.Pasta import read_seq, read_preset
    from numpy import array, zeros, zeros_like, argmax, mean
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    
    amino_acids = read_preset(statsfile)
    seq = read_seq(seqfile)
    seq_labels = [aa.one_let for aa in seq]
    aa_labels = [aa.one_let for aa in amino_acids]
    fragment_labels = ['*'] * len(seq)
    del aa_labels[1] ## remove Cy*

    data = zeros(len(seq)) ## scores of aa in sequence
    data2 = zeros(len(seq)) ## maximal scores in residue
    conf, similarity = conf
    
    for frag, i in conf:
        
        j = i + len(frag.residues)
#        print i,''.join([r.name for r in frag.residues]), frag.errors[i]
        data[i:j] = 1 - frag.errors[i]
        data2_maxima_at = [argmax(res.profile) for res in frag.residues]
        data2[i:j] = [max(res.profile) for res in frag.residues]
        fragment_labels[i:j] = [aa_labels[k] for k in data2_maxima_at]
    
    data_incorrect = zeros_like(data)
    data_correct = zeros_like(data)
    data2_incorrect = zeros_like(data2)
    data2_correct = zeros_like(data2)
    
    i = 0
    for a, b in zip(seq_labels, fragment_labels):
        if a == b:
            data_correct[i] = data[i]
            data2_correct[i] = data2[i]
        else:
            data_incorrect[i] = data[i]
            data2_incorrect[i] = data2[i]
        i += 1
    
    data_mean = mean(data)
    data2_mean = mean(data2)
    
    n_per_subplot = n
    n = len(data) % n_per_subplot ## rest
    m = (len(data) - n) / n_per_subplot + 1 ## number of subplots
    
    frag_posteriors = plt.figure()
    seq_posteriors = plt.figure()
    i = 0
    k = 1
    print '\nTop Configuration'
    while k <= m and i < len(seq):
        j = i + n_per_subplot
        if j >= len(seq):
            j = i + (len(seq) - i)
            
        width = .5
        xlocations = array(range(len(seq)))
        
        print '%i\t%s\t%i' % (i, ''.join(fragment_labels[i:j + 1]), j)
        
        ax = seq_posteriors.add_subplot(m, 1, k)
        ax.bar(xlocations[i:j] - (width / 2.), data_correct[i:j],
               color='black', linewidth=0, width=width, aa=True,
               label='Match')
        ax.bar(xlocations[i:j] - (width / 2.), data_incorrect[i:j],
                color='r', linewidth=0, width=width, aa=True,
                label='No Match')
        ax.plot([xlocations[i], xlocations[j - 1]], [data_mean, data_mean], color='r', label='Mean')
        ax.set_xticks(xlocations[i:j])
        ax.set_xticklabels(seq_labels[i:j])
        ax.set_xlim(xlocations[i] - (width / 2.), xlocations[j - 1] + (width / 2.))
        ax.set_ylim(0, 1)
        ax.set_ylabel("%i-%i" % (i, j))
        
        ax2 = frag_posteriors.add_subplot(m, 1, k)
        ax2.bar(xlocations[i:j] - (width / 2.), data2_correct[i:j],
               color='black', linewidth=0, width=width, aa=True,
               label='Match')
        ax2.bar(xlocations[i:j] - (width / 2.), data2_incorrect[i:j],
                color='r', linewidth=0, width=width, aa=True,
                label='No Match')
        ax2.plot([xlocations[i], xlocations[j - 1]], [data2_mean, data2_mean], color='r', label='Mean')
        ax2.set_xticks(xlocations[i:j])
        ax2.set_xticklabels(fragment_labels[i:j])
        ax2.set_xlim(xlocations[i] - (width / 2.), xlocations[j - 1] + (width / 2.))
        ax2.set_ylim(0, 1)
        ax2.set_ylabel("%i-%i" % (i, j))
        
        i = j + 1
        k += 1
    
    print 'Sequence Similarity: %.4f'%similarity
    
    top_ax = frag_posteriors.get_axes()[0]
    bottom_ax = frag_posteriors.get_axes()[-1]
    top_ax.legend(bbox_to_anchor=(1.01, 0, 1, 1), loc=2, borderaxespad=0.)

    bottom_ax.set_xlabel("Frag Residues")
#    top_ax.set_title(to_file + ' Posteriors')
    
    top_ax2 = seq_posteriors.get_axes()[0]
    bottom_ax2 = seq_posteriors.get_axes()[-1]
    top_ax2.legend(bbox_to_anchor=(1.01, 0, 1, 1), loc=2, borderaxespad=0.)

    bottom_ax2.set_xlabel("Sequence Residues")
#    filename = to_file.split('/')[-1]
#    top_ax2.set_title(filename + ' Posteriors')

    if to_file != None:
        filename = to_file.split('/')[-1]
        frag_posteriors.set_size_inches(18.5,10.5)
        seq_posteriors.set_size_inches(18.5,10.5)
        frag_posteriors.savefig(filename + '_fraglabeled', dpi=dpi)
        seq_posteriors.savefig(filename + '_seqlabeled', dpi=dpi)
    else:
#        plt.show()
        return seq_posteriors

def plot_conf_score(conf, seqfile, to_file=None,dpi=300):
    from Pasta import read_seq
    
    c, similarity = conf ## configuration and similarity tuple
    n = array([len(frag) for frag, pos in c]).sum()
    m = len(read_seq(seqfile))
    mu, sigma = avg_score_distribution(n, m, data_points=10e4, bins=100, visualize=False)
    bound_min = mu
    bound_max = max_score(n, m)
    c_score = similarity
    X = linspace(bound_min - .1, 1, 1000, endpoint=True)
    Z = array([z_score(x, mu, sigma) for x in X]) ## Z Scores
    P = array(map(p_value, Z)) ## P Values
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(X, P, label='Err. Func.')
    ax.plot(c_score, p_value(z_score(c_score, mu, sigma)), 'o', color='r')
    ax.plot(1.96 * sigma + mu, p_value(1.96), 'o', color='black')
    ax.plot(bound_max, p_value(z_score(bound_max, mu, sigma)), 'o', color='black')
    
    ax.vlines(c_score, 0., p_value(z_score(c_score, mu, sigma)), colors='r', linestyles='--', label='Conf. Score')
#    ax.vlines(1.96 * sigma + mu, 0., p_value(1.96), colors='black', linestyles='-.', label='Crit. Val.')
#    ax.vlines(bound_max, 0., p_value(z_score(bound_max, mu, sigma)), colors='black', linestyles='--', label='Upper Bound')
    
    ax.text(1.96 * sigma + mu, p_value(1.96) + .02, '%.4f' % (1.96 * sigma + mu))
    ax.text(c_score, p_value(z_score(c_score, mu, sigma)) - .04, '%.4f' % c_score)
    ax.text(bound_max, p_value(z_score(bound_max, mu, sigma)) + .02, '%.4f' % bound_max)
    
    area = linspace(1.96 * sigma + mu, bound_max, 100)
    ax.fill_between(area, 0., p_value(z_score(area, mu, sigma)), color='lightgrey')
    
    ax.set_xlabel("Sequence Similarity")
    ax.set_ylabel("P-Value")
    ax.set_xticks(arange(0., 1.1, .1))
    ax.set_yticks(arange(0., 1.1, .1))
    ax.set_ylim(0., 1.1)
    '''TODO: Legend will be clipped if anchored'''
#    ax.legend(bbox_to_anchor=(1.01, 0, 1, 1), loc=2, borderaxespad=0.)
    ax.legend(loc='lower right', borderaxespad=0.)
    
    if to_file != None:
        filename = to_file.split('/')[-1]
        fig.set_size_inches(18.5,10.5)
        fig.savefig(filename + "_p_val", dpi=dpi)
    else:
        return fig
#        plt.show()

def avg_score_distribution(n, m, data_points=10e4, bins=100, visualize=False):
    '''
    Returns estimated mean and standard deviation of score distribution
    for randomized amino acid recognition result
     
    n := sum of all fragment lengths
    m := length of sequence
    '''
    
    assert n <= m
    
    scores = []
    for i in range(int(data_points)):
    
        p = random(n)
        avg = 1 - ((m - n + p.sum()) / m)
        scores.append(avg)
        
    data = array(scores)
    mu = mean(data) ## mean value
    sigma = std(data) ## standard deviation
    
    if visualize:
        n, bins, patches = plt.hist(data, bins, normed=1, alpha=.3)
        y = mlab.normpdf(bins, mu, sigma)
        plt.plot(bins, y, 'r-', linewidth=1)
        plt.vlines(mu, 0, mlab.normpdf([mu], mu, sigma), colors='r')
        plt.show()
    
    return mu, sigma
    
def max_score(n, m):
    '''
    Returns maximal score if no error in mapped fragments
    n := sum of all fragment lengths
    m := length of sequence
    '''
    assert n <= m
    
    n = float(n)
    m = float(m)
    return 1 - ((m - n) / m)

def z_score(x, mu, sigma):
    '''
    It represents the number of standard deviations from the mean. 
    
    z = (x - mu) / sigma
    
    So, for a distribution that has a mean of 100 and a standard 
    deviation of 10, the value 90 is represented by a z-value of:
    
    z = (90 - 100) / 10 = -1

    '''
    return (x - mu) / sigma

def p_value(z):
    '''
    p value is large when z is large
    if a value is very far from the mean value
    p value becomes large
    '''
    return special.erf(z / sqrt(2))

if __name__ == '__main__':
    from Pasta import read_assignments, read_seq, read_pasta
    from FragmentGraph import Fragment
    
    pastafile = '/kyb/agbs/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/PastaTest/pastatest_full.txt'
    assi_file = '/kyb/agbs/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/PastaTest/CH-SFRRM.longlist'
    seqfile = '/kyb/agbs/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/PastaTest/pastatest.fasta'
    statsfile = '/kyb/agbs/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/shiftPresets/bmrb_ascii.shift'
    best_conf = pickle.load(open('/kyb/agbs/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/PastaTestResults/best_conf.p'))
    path = '/kyb/agbs/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/PastaTestResults/CH-SFRRM'
    
    fig1 = plot_mapping(best_conf, seqfile, statsfile, to_file=path)
    
#    fig2 = plot_conf_score(best_conf, seqfile, to_file=path)
#    plt.show()
    
    
