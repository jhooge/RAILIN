'''
Created on Apr 14, 2012

@author: jhooge
'''
from matplotlib.axes import Subplot

def nucleus_mean(nuc, amino_acids):
    '''
    computes the mean over all amino acid means for a specific nucleus
    '''
    from numpy import array, mean
    means = array([aa.shifts[nuc][0] for aa in amino_acids if nuc in aa.shifts.keys()])
    return mean(means)
            
def data_array(amino_acids, nuclei):
    from numpy import zeros
    
    means = zeros((len(amino_acids),len(nuclei)))
    stds = zeros((len(amino_acids),len(nuclei)))
    
    for i,aa in enumerate(amino_acids):
        for j,nuc in enumerate(nuclei):
            if nuc in aa.shifts.keys():
                means[i,j] = aa.shifts[nuc][0] 
                stds[i,j] = aa.shifts[nuc][1]
            else:
                means[i,j] = None
                stds[i,j] = None
                
    return means,stds

def nan_mean(ndarray):
    from numpy.ma import masked_array
    from numpy import isnan, mean
    mndarray = masked_array(ndarray,isnan(ndarray))
    return mean(mndarray)

def plot_chem_shift_dist(title,residues, amino_acids):
    three_lets = [aa.three_let for aa in amino_acids]
    nuclei = ['CO', 'CA', 'CB', 'CG', 'CG1', 'CG2',
              'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2',
              'CE3', 'CZ', 'CZ2', 'CZ3', 'CH2']
    colors = ['b','g','r','c','c','c',
              'm','m','m','y','y','y',
              'y','gray','gray','gray','k']
    
    means, stds = data_array(amino_acids, nuclei)
    fig = pl.figure()
    ax = Subplot(fig,111)
    ax = pl.gca()
    fig.add_subplot(ax)
    
    x = means
#    print nan_to_num(x)
    y = ones(x.shape)
    for i in range(0,x.shape[0]):
        y[i]*=i
        
    for j in range(0,x.shape[1]):
        ax.errorbar(x.T[j], y.T[j], xerr=nan_to_num(stds.T[j]), 
                    linestyle='None',marker='o', color=colors[j], 
                    label=(nuclei[j]), linewidth=4, capsize=20)
    
    for r in residues:
        x = ones_like(means)*nan
        name = r.name
        shifts, atoms = r.get_carbons(previous=False)
        for shift,atom in zip(shifts,atoms):
            if atom in nuclei:
                j = nuclei.index(atom)
#                if name != "NAA":
#                    name = one2Three(r.name[0])
                i = three_lets.index(name)
                x[i][j] = shift
        for j in range(0,x.shape[1]):
            ax.plot(x.T[j],y.T[j], marker="x",linestyle="None", color=colors[j],markersize=15)
    
    pl.yticks(arange(0,20), three_lets)
    pl.ylim(-1,20)
    pl.title('Verteilung chemischer Verschiebungen\n %s'%title)
    pl.ylabel('Aminosaeure')
    pl.xlabel('ppm')
    pl.rcParams.update({'font.size': 22})
    # Shink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    pl.show()

if __name__ == '__main__':
    from FileHandler import FileHandler
    import pylab as pl
    from numpy import ones, array, mean, nan_to_num, arange, ones_like,nan
    from Definitions import one2Three
    
    statsfile = '/is/ei/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/tests/reference_lists/bmrb.shift'
    folder = "/is/ei/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/Datasets/Ubiquitin/residue_lists/"
    pastafn = folder+"Ub_bmrb.pasta"
#    pastafn = folder+"Ub_opt_unambiguous.pasta"
    fh = FileHandler()
    amino_acids = fh.read_preset(statsfile)
    del amino_acids[1]
    name = "Ubiquitin"
#    pastafn = "/is/ei/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/Datasets/"\
#             "%s/residue_lists/%s.pasta"%(name,name)
    residues = fh.read_pasta(pastafn, statsfile)
    plot_chem_shift_dist(name, residues, amino_acids)