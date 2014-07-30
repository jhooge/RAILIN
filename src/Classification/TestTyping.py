
'''
Created on Apr 22, 2012

@author: jhooge
'''

class TestTyping(object):
    '''
    classdocs
    '''


    def __init__(self, pasta_obj):
        '''
        Constructor
        '''
        self.pasta = pasta_obj
        self.P = self.P_from(pasta_obj)
        
def S_from(pasta_obj):
    from mulch.math import outer_2way
    from numpy import array, concatenate, multiply, swapaxes, ones
    from copy import deepcopy
    import pylab as pl
    
    seq = pasta_obj.seq
    aas = deepcopy(pasta_obj.amino_acids)
    del aas[1]
    
    S = []
    for aa in seq:
        aa.set_profile(aas)
        S.append(aa.profile)
    S = array(S)
    S = concatenate([ones((1, S.shape[1])), S], 0)
    S = outer_2way(S[:-1].T, S[1:].T, multiply)
    S = swapaxes(S, 0, 1)
    return S
    
def P_from(pasta_obj):
    from numpy import dot, reshape
    from copy import deepcopy
    n = len(pasta_obj.seq) ## len sequence
    m = len(pasta_obj.residues) ## len pasta list
    
    P = pasta_obj.P ## posterior matrix from aa recog
    S = S_from(pasta_obj)
    P = dot(reshape(P, (m, -1)), reshape(S, (n, -1)).T)
    ## set probabilities to 0 for seq pos 0, when residue has [i-1] shifts
    for i, r in enumerate(pasta_obj.residues):
        if len(r.shifts_im1.values()) > 0:
            P[i][0] = P[i][0] * 0
            
    return P
        
def savefig(self, fn, P, xlabel, ylabel, title, colorbar_label, cbticks=None, dpi=300):
    from Visualization.MatrixPlot import mat_fig
    fig = mat_fig(P, xlabel, ylabel,
                  title, colorbar_label=colorbar_label,
                  cbticks=cbticks)
    fig.savefig(fn, dpi=dpi)

#def result(pasta_obj):
#    P = P_from(pasta_obj)
#    sequence = pasta_obj.seq
#    residues = pasta_obj.residues
#    i = 0
#    correct = 0
#    incorrect = 0
##    print "pasta_no|expected pair|seq(argmax)|argmax|max(posterior)|expected posterior"
#    while i < len(residues):
#        
#        res = residues[i]
#        expected = ''.join([three2One(res.name_im1),
#                            three2One(res.name)])
##        expected = ''.join([res.name_im1,
##                            res.name])
#        
#        if argmax(P[i]) == 0:
#            classified = ''.join(['X',
#                                  sequence[argmax(P[i])].one_let])
#        else:
#            classified = ''.join([sequence[argmax(P[i]) - 1].one_let,
#                                  sequence[argmax(P[i])].one_let])
#        
#        if expected != classified:
##            print res.pasta_no, expected, classified, argmax(P[i]), max(P[i]), P[i][i]
#            incorrect += 1
#        else:
#            correct += 1 
#        i += 1
#    
##    print "Correct:   ", correct
##    print "Incorrect: ", incorrect
#    
#    return float(correct) / (correct + incorrect)

def entropy_plot(P, labels=None):
    from numpy import log, arange, take, array
    import pylab as pl
    from scipy.stats.distributions import entropy
    
    sorted_labels = sorted(labels)
    e = []
    for label in sorted_labels:
        k = labels.index(label)
        e.append(entropy(P[k]) / log(2))
    
    k = sorted_labels.index("XA")
    sorted_labels.remove("XA")
    del e[k]
    
    e = array(e)
    
    for k, label in enumerate(sorted_labels):
        print label, e[k]
    
    xlocations = arange(0, len(sorted_labels), 20)
    labels = take(array(sorted_labels), xlocations)
    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.plot(arange(0, len(e)), e, "blue")
    ax.set_ylabel("Bit")
    ax.set_xticks(xlocations)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_title("Entropie mit fehlendem CA")
    ax.axhline(mean(e), color="red",
               label=r"$-\frac{1}{n \cdot \log(2)} H(p_1,\ldots,p_m)$"\
               "\n $= %.3f Bit$" % mean(e))
    ax.legend()
    pl.rcParams.update({'font.size': 22})
    fig.savefig("doc/Diplomarbeit/bilder/Figures/Entropy_missing_CA.eps", dpi=600)

def accuracy(title, x, y, y2):
    print x
    print y
    print y2
    y = array(y) * 100
    y2 = array(y2)*100
    z = polyfit(x, y, deg=2)
    p = poly1d(z)
    z2 = polyfit(x, y2, deg=2)
    p2 = poly1d(z2)
    print p
    print p2
    pl.scatter(x, y, s=75, color="black", label="Daten (vollstaendig)")
    pl.plot(x, p(x),
            linestyle='--', color="black",
            label=r"$f(x) = 0.6004 x^2 + 1.276 x + 4.404$")
    pl.scatter(x, y2, s=75, color="gray",label=r"Daten ($C_\alpha,C_\beta$)")
    pl.plot(x, p2(x),
            linestyle='-.',color="gray",
            label=r"$f(x) = 2.401 x^2 - 2.454 x + 32.27$")
    pl.title(title, fontsize=20)
    pl.ylabel("Inkorrekt klassifiziert [%]", fontsize=20)
    pl.xlabel("Streuungsinterval [ppm]", fontsize=20)
    pl.ylim(min(y) - 5, 100)
    pl.xlim(0, 3)
    pl.legend(loc=2)
    pl.rcParams.update({'font.size': 22})
    pl.show()

def test_typing(typing_files):
    acc = []
    run_times = []
    ratios = []
    for fn in typed_files:
        print fn
        p = pickle.load(open(fn, "r"))
        correct = 0
        incorrect = 0
        for i, r in enumerate(p.residues):
            dipep = ''.join(map(three2One, [r.name_im1, r.name]))
            if dipep == labels[i]:
                correct += 1
            else:
                print i, r.pasta_no, labels[i], dipep
                incorrect += 1
            ratios.append((correct, incorrect))
        acc.append(float(incorrect) / len(labels))
        print "Accuracy = %.4f" % acc[-1], ratios[-1]
#    pickle.dump(acc, open("UbqBaxLaplaceCACB.p", "wb"))

if __name__ == "__main__":
    import pylab as pl
    from Pasta import Pasta
    from FileHandler import FileHandler
    from numpy import argmax, arange, log, mean, array, linspace
    from Definitions import three2One, one2Three
    from scipy.stats.distributions import entropy
    from numpy import polyfit, poly1d
    from copy import deepcopy
    import pickle
    
    fh = FileHandler()
    pastafile = 'Datasets/Ubiquitin/residue_lists/Ub_bmrb.pasta'
    seqfile = 'Datasets/Ubiquitin/sequences/Ub_bmrb.fasta'
    statsfile = 'tests/reference_lists/bmrb.shift'
    sequence = fh.read_seq(seqfile, statsfile)
    orig_res = fh.read_pasta(pastafile, statsfile)
    labels = map(lambda x: ''.join([three2One(x.name_im1), 
                                    three2One(x.name)]), 
                                    orig_res)
    labels.reverse()
    noise = linspace(0, 3, 31)
    typed_files = ["Datasets/Ubiquitin/Results/Ubq_Bax_Robustness_Tests_Max/"\
                   "NoisyData/UbqBaxGaussFullNoise=%.1f.p" % n for n in noise]
#    typed_files = ["Datasets/Ubiquitin/Results/Ubq_Bax_Robustness_Tests_Max/"\
#                   "NoisyData/UbqBaxLaplaceCACBNoise=%.1f.p" % n for n in noise]
    
    test_typing(typed_files)
    
#    acc = pickle.load(open("typing_accuracy_noisy_pairs.p", "r"))
#    acc2 = pickle.load(open("typing_accuracy_noisy_pairsCACB.p", "r"))

#    acc = pickle.load(open("UbqBaxGaussFull.p", "r"))
#    acc2 = pickle.load(open("UbqBaxGaussCACB.p", "r"))
#    title = "Typisierungsgenauigkeit (Gauss-PDF)\nUbiquitin (Cornilescu et al.)"
    
    acc = pickle.load(open("UbqBaxLaplaceFull.p", "r"))
    acc2 = pickle.load(open("UbqBaxLaplaceCACB.p", "r"))
    title = "Typisierungsgenauigkeit (Laplace-PDF)\nUbiquitin (Cornilescu et al.)"
    
    accuracy(title, noise, acc, acc2)
#    accuracy(title, noise, acc2)


