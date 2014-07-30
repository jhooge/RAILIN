'''
Created on Jun 19, 2011

@author: jhooge
'''
from Classification.Pasta import read_preset
import pylab as pl

statsfile = 'bmrb_ascii.shift'

amino_acids = read_preset(statsfile)
labels = [aa.three_let for aa in amino_acids] 
ca = [aa['CA'] for aa in amino_acids if aa.one_let != 'G']
cb = [aa['CB'] for aa in amino_acids if aa.one_let != 'G']

N = 100
x = pl.linspace(40., 80, N)
y = pl.linspace(0., 80., N)
X, Y = pl.meshgrid(x, y)
Z = []

z = pl.bivariate_normal(x, y, ca[0][1], cb[0][1], ca[0][0], cb[0][0])

#for i in range(0, len(ca)):
#    mu_a = ca[i][0]
#    sig_a = ca[i][1]
#    mu_b = cb[i][0]
#    sig_b = cb[i][1]
#    Z.append(pl.bivariate_normal(X, Y, sig_a, sig_b, mu_a, mu_b))
#Z = array(Z)
#
#print Z

axContour = pl.subplot(1,1,1)
pl.title(r'Verteilungen der $C_\alpha$ und $C_\beta$ Verschiebungen')
pl.xlabel(r'$C_\alpha$ (in ppm)')
pl.ylabel(r'$C_\beta$ (in ppm)')

#for z in Z:
#    pcolor(X, Y, z)
#    contour(X, Y, z, zdir='x')
#    pl.contour(X, Y, z, zdir='z')
#    pl.scatter(X, Y)

#help(pl.contour)
pl.show()
