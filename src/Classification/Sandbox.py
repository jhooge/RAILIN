from Linking import Linking
from Pasta import Pasta
from numpy import linspace, zeros
from DataGenerator import DataGenerator
import time

import pylab as pl
import pickle
from numpy import meshgrid

def no_constraints_matrix(pasta):
    tol = linspace(0., 1., 100)
    no_neighbors = []
    for t in tol:
        linking_mat = pasta.linking(t)
        n = linking_mat.size - linking_mat.sum().astype('i')
        print "tol=%.4f | constraints=%i"%(round(t,4), n)
        no_neighbors.append(n)
    return no_neighbors

def neighbors(tol,s_len):
    no_neighbors = zeros((len(tol),len(s_len)))
    for i,t in enumerate(tol):
        for j,l in enumerate(s_len):
            sequence = dg.generate_rnd_seq(l)
            residues = dg.generate_expected_data(sequence, noise=2., 
                                                 instrument_error=0.3)
            linking_mat = Linking(residues).linking_matrix(t)
            n = linking_mat.size - linking_mat.sum().astype('i')
            print "S_len=%i | tol=%.4f | n=%i"%(len(residues),t,n)
            no_neighbors[i][j] = n
    print no_neighbors
    return no_neighbors

if __name__ == "__main__":
    
    statsfile = 'shiftPresets/bmrb.shift'
#    folder = "Datasets/Ubiquitin/"
#    pastafile = folder + "Ubiquitin_newVT.pasta"
#    seqfile = folder + "Ub.fasta"
#    result_folder = folder + "Results/"
    
    dg = DataGenerator(statsfile)
    
#    pasta = Pasta(pastafile, statsfile, seqfile)
    
#    no_neighbors = neighbors(tol, s_len)
    
#    data = []
#    seq_lengths = linspace(1,30,30).astype("i")
#    tol = linspace(0., 1., 20)
#    sequences = map(dg.generate_rnd_seq, seq_lengths)
#    for s in sequences:
#        data.append(dg.generate_expected_data(s, noise=2., 
#                                             instrument_error=0.3))
#        print "Data generated [%i] "%len(s)
#        
#    no_neighbors = zeros((len(tol),len(data)))
#    for i,t in enumerate(tol):
#        for j,residues in enumerate(data):
#            n = Linking(residues).linking_matrix(t).sum()
#            print "S_len=%i | tol=%.4f | n=%i"%(len(residues),t,n)
#            no_neighbors[i][j] = n
            
#    pl.matshow(no_neighbors)
#    pickle.dump(no_neighbors, open("neighbors.p","wb"))
    from numpy import polyfit, poly1d
    tol = linspace(0., 1., 20)
    s_len = linspace(1,250,250)
    no_neighbors = pickle.load(open("neighbors2d.p","r"))
    
    colors = ["r","g","b"]
    for i,k in enumerate([0,9,19]):
        x = s_len
        y = no_neighbors[k]
        z = polyfit(x, y, deg=2)
        p = poly1d(z)
        print "Tol=",tol[k]
        print p
        print "len = 50 | No Constraints=",no_neighbors[k][50]
        print "len = 100 | No Constraints=",no_neighbors[k][100]
        print "len = 250 | No Constraints=",no_neighbors[k][249]
        pl.plot(s_len[0:20],no_neighbors[k][0:20],marker='.',linestyle='',color=colors[i],label='tol=%.2f'%tol[k])
        pl.plot(s_len[0:20],p(s_len[0:20]),linestyle='-',color=colors[i],label='fit(%.2f)'%tol[k])
        pl.ylabel("Anzahl aktiver Constraints")
        pl.xlabel("Sequenzlaenge")
        pl.legend()
    pl.show()
    
#    X = tol
#    Y = s_len
#    Z = pickle.load(open("neighbors2d.p","wb"))
#    pl.plot_surface()
    
#    pl.plot(seq_lengths,neighbors)
#    pl.show()
