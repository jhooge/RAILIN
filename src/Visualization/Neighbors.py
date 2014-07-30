'''
Created on May 16, 2012

@author: jhooge
'''

from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from numpy import linspace,meshgrid
import pickle
import pylab as pl

if __name__ == '__main__':
    
    ## plot up to 250 residues
#    seq_lengths = linspace(1,20,20).astype("i")
#    tol = linspace(0., 1., 100)
#    neighbors = pickle.load(open("pickled/n_seq_len_Ubiquitin.p","r"))
#    times = pickle.load(open("pickled/n_seq_len_time.p","r"))
    
#    fig = pl.figure()
#    ax = fig.add_subplot(111)
#    ax.plot(tol,neighbors)
#    pl.title('Linking Constraints Ubiquitin\n73 Reste',fontsize=20)
#    pl.ylabel('Anzahl aktiver Constraints',fontsize=20)
#    pl.xlabel('Toleranz [ppm]',fontsize=20)
#    fig.savefig("Figures/nConstraintsUbiquitin", dpi=600)
#    pl.show()


    X = linspace(0., 1., 20)
    Y = linspace(1,250,250)

#    X = linspace(0., 1., 20)
#    Y = linspace(1,10,10)
    Z = pickle.load(open("pickled/neighbors.p","r"))
    X,Y = meshgrid(X,Y)
    fig = pl.figure()
    ax = axes3d.Axes3D(fig)
    ax.plot_surface(X,Y,Z.T,rstride=1,cstride=1,
                    cmap=cm.jet,
                    linewidth = 1)
    ax.set_title('Linking Constraints',fontsize=20)
    ax.set_zlabel("Anzahl aktiver Constraints",fontsize=20)
    ax.set_ylabel('Sequenzlaenge',fontsize=20)
    ax.set_xlabel('Toleranz [ppm]',fontsize=20)
    fig.savefig("Figures/nConstraints2d", dpi=300)
    pl.show()
    