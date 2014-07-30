'''
Created on Feb 15, 2012

@author: jhooge
'''
import pylab as plt
import numpy as np

def plot_info_content(array, title, xlabels, ylabel, dpi=300):
    
    xlocations = np.arange(len(xlabels))
    width = .5
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.bar(xlocations - (width / 2.), array,
               color='black', linewidth=0, width=width, aa=True,)
    
    ax.set_ylabel("Some label")
    ax.set_xticks(xlocations)
    ax.set_xticklabels(xlabels)
    ax.set_title(title)
    '''TODO: Legend will be clipped if anchored'''
#    ax.legend(bbox_to_anchor=(1.01, 0, 1, 1), loc=2, borderaxespad=0.)
    ax.legend(loc=2, borderaxespad=0.)
    fig.autofmt_xdate()
    
    return ax
    

if __name__ == '__main__':
    from numpy import array, prod, multiply
    import pylab as pl
#    proc_times = [1,2,3,4]
#    plot_times(proc_times)

a = array([[1, 2, 1],
          [2, 2, 2],
          [3, 3, 3]]).astype("f")
          
b = 1. / a.sum(1)
for i, ele in enumerate(b):
    a[i] = ele * a[i]
print a
fig = plot_info_content(a[0],"Some title",["a","b","c"],"Some label")
pl.show()

