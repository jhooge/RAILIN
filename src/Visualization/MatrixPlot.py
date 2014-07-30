'''
Created on May 23, 2012

@author: jhooge
'''

def mat_fig(A,title, xlabel, ylabel,
            xticklabels,
            yticklabels,
            colorbar_label=None,
            cbticks=None):
    """
    Returns a matshow figure for a linking matrix
    
    @param A: Some matrix to generate figure from
    @type A: numpy.ndarray
    @param xlabel: Label for x axis
    @type xlabel: string
    @param ylabel: Label for y axis
    @type ylabel: string
    @param title: Title for figure
    @type title: string
    @param cbticks: Ticks for colorbar
    @type cbticks: list
    
    @return: Figure
    @rtype: pylab.figure()
    """
    import pylab as pl
    import numpy as np
    
    fig = pl.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(A, interpolation='nearest', cmap=pl.cm.Greys)
    cb = fig.colorbar(cax, ticks=cbticks)
    if colorbar_label is not None:
        cb.set_label(colorbar_label,fontsize=20) 
    ax.set_xticklabels(xticklabels)
    ax.set_xticks(np.arange(0,len(xticklabels)))
    ax.set_yticklabels(yticklabels)
    ax.set_yticks(np.arange(-1,len(yticklabels)+1))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    pl.rcParams.update({'font.size': 16})
    return fig

if __name__ == '__main__':
    from numpy import eye
    import pylab as pl
    from pylab import savefig, show
    fig = mat_fig(eye(21),"bla","blub", "foo")
#    fig.savefig("TestFig.png")
    pl.show()
    pass