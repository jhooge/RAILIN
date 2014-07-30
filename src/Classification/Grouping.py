'''
Created on Mar 29, 2012

@author: jhooge
'''

def Load(filename, gzip = 0, lock = None):
    import cPickle, os

    ## file locked?
    if lock is not None:
        import os

        path, file = os.path.split(filename)
        lockfile = '%s/%s.lock' % (__lockPath, file)
        
        if os.path.exists(lockfile):
            raise '%s is locked (%s)' % (filename, lockfile)

        ## lock file
        open(lockfile, 'w').close()

    filename = os.path.expanduser(filename)
        
    if gzip:
        import gzip
        try:
            f = gzip.GzipFile(filename)
        except:
            return
        
    f = open(filename)

    try:
        this = cPickle.load(f)
    except:
        import Numeric

        f.close()
        try:
            f = gzip.GzipFile(filename)
        except:
            f = open(filename)
        
        u = Numeric.Unpickler(f)
        this = u.load()

    f.close()

    if lock is not None:
        ## release lock
        import os
        try:
            os.remove(lockfile)
        except:
            raise 'missing lockfile %s' % lockfile

    return this

class Grouping(object):
    '''
    classdocs
    '''

    def __init__(selfparams):
        '''
        Constructor
        '''
        pass 
    
    def cluster_points(self, x, tol=(0.4, 0.04)):

        from numpy import mean, fabs, array, ones
    
        clusters = []
        tol = array(tol)
        labels = -ones(len(x)).astype('i')
        counter = -1
        
        for i in range(len(x)):
    
            for counter, cluster in enumerate(clusters):
    
                y = mean(cluster,0)
                if (fabs(x[i]-y) < tol).all():
                    cluster.append(x[i])
                    labels[i] = counter
                    break
            else:
                clusters.append([x[i]])
                labels[i] = counter + 1
    
        return clusters, labels

    def assign_amides(self,x, y, tol=(0.4,0.04), n=3):
    
        from numpy import subtract, fabs, argmin, array, compress, nonzero
    
        ## assign based on closest Mahalanobis distance
    
        tol = array(tol)
        
        d1 = fabs(subtract.outer(x[:,0],y[:,0])) / tol[0]
        d2 = fabs(subtract.outer(x[:,1],y[:,1])) / tol[1]
    
        labels = argmin(d1**2 + d2**2,0)
    
        ## throw out wrong assignments (i.e. assignments that are
        ## more than n*tol away from cluster center)
    
        newidx = 0 * labels - 1 
    
        for i in range(len(y)):
    
            d = fabs(x[labels[i]]-y[i])
    
            ## correct assignment
    
            if (d < n*tol).all():
                newidx[i] = labels[i]
    
        labels = newidx
    
        ## collect wrongly assigned data points and cluster them without a reference
        ## (i.e. unsupervised clustering)
    
        invalid = labels==-1
    
        clusters_invalid, idx_invalid = self.cluster_points(compress(invalid,y,0), tol)
    
        labels[nonzero(invalid)[0]] = idx_invalid + len(x)
    
        ## calculate new cluster centers 
    
        centers = []
        for i in range(labels.max()+1):
            z = compress(labels==i,y,0)
            if len(z):
                centers.append(z.mean(0))
            else:
                centers.append(None)
        
        return labels, centers
    
    def create_spin_systems(self,labels, data):
        from numpy import compress
    
        systems = []
        for i in set(labels.tolist()):
            systems.append(compress(labels==i,data,0))
    
        return systems
    
    def define_spin_system(self,peaks, tol=0.8):
    
        from numpy import mean, subtract, fabs, argmin, array, nonzero, argmax, sum, std, \
             argsort, concatenate, compress
    
        ## clusters CA shifts
    
        ## TODO: CA dimension hard coded
    
        clusters, labels = self.cluster_points(peaks[:,:1],tol=[tol])
        ## two distinct CA resonances (regular case)
        
        leftover_peaks = None
        if len(clusters) == 2:
    
            ## detect CA-CA peak
    
            ca_shifts = array(map(mean,clusters))
            ca_peaks = peaks[argmin(fabs(subtract.outer(ca_shifts,peaks[:,1])),1)]
    
            ## more intense peak is CA(i)
    
            if fabs(ca_peaks[0][-1]) > fabs(ca_peaks[1][-1]):
    
                shifts_i = peaks[nonzero(labels==0)[0]]
                shifts_im= peaks[nonzero(labels==1)[0]]
    
            else:
                
                shifts_i = peaks[nonzero(labels==1)[0]]
                shifts_im= peaks[nonzero(labels==0)[0]]
    
        ## only one CA resonance (related to prolines?)
    
        elif len(clusters) == 1:
    
            ## TODO: discuss with Vincent whether this is correct
    
            shifts_im= peaks[nonzero(labels==0)[0]]
            shifts_i = None
    
        ## more than two CA resonances (i.e. maybe the spin-system is contaminated with
        ## a resonance from another spin-system)
        ## keep maximally populated clusters
                
        else:
            print 'more than two CA-carbons'
            
            index  = argmax(map(len,clusters))
            shifts = compress(labels==index,peaks,0)
            mean_nh= take(shifts,(2,3),1).mean(0)
            std_nh = std(take(shifts,(2,3),1),0)
            
            ## select peaks with amide shifts that are closest to reference
    
            distances = []
            for index in range(len(clusters)):
                shifts = compress(labels==index,peaks,0)
                distances.append(sum((mean_nh-take(shifts,(2,3),1).mean(0))**2/std_nh**2)**0.5)
    
            index1, index2 = argsort(distances)[:2]
            
            ca_shifts = array([mean(clusters[i],0)[0] for i in [index1, index2]])
            ca_peaks = peaks[argmin(fabs(subtract.outer(ca_shifts,peaks[:,1])),1)]
    
            ## more intense peak is CA(i)
    
            if fabs(ca_peaks[0][-1]) > fabs(ca_peaks[1][-1]):
    
                shifts_i = peaks[nonzero(labels==index1)[0]]
                shifts_im= peaks[nonzero(labels==index2)[0]]
    
            else:
                
                shifts_i = peaks[nonzero(labels==index2)[0]]
                shifts_im= peaks[nonzero(labels==index1)[0]]
    
            indices = range(len(clusters))
            indices.remove(index1)
            indices.remove(index2)
    
            leftover_peaks = concatenate([compress(labels==index,peaks,0) for index in indices],0)
        
        return shifts_im, shifts_i, leftover_peaks
    
    def make_carbon_dict(self, peaks):

        from numpy import argsort, fabs, take, array
    
        if peaks is None: return {}
        ca_shift = peaks[:,0].mean()
        shifts = array(list(set(peaks[:,1].tolist())))
        shifts = take(shifts,argsort(fabs(ca_shift-shifts)))
        labels = ['C%d'%i for i in range(len(peaks))]
        labels[0] = 'CA'

        return dict(zip(labels,shifts.tolist()))
    
    def make_pasta_residue(self, shifts_im, shifts_i, number):

        from Classification.Residue import PastaResidue
    
        shifts_im = self.make_carbon_dict(shifts_im)
        shifts_im = dict([(k+'i-1',v) for k,v in shifts_im.items()])
        shifts_i  = self.make_carbon_dict(shifts_i)
    
        residue = PastaResidue()
        residue.set_attribute('pasta_no', number)
        residue.set_attribute('original_no', number)
        residue.set_attribute('energy', -1)
        residue.set_attribute('shifts_i', shifts_i)
        residue.set_attribute('shifts_im1', shifts_im)

        return residue

    def add_centroids_to_shift_dicts(self,centroids,labels,residues):
        labels = list(set(labels))
        for res,id in zip(residues,labels):
            if centroids is not None:
                res.shifts_i['N15'] = centroids[id][0]
                res.shifts_i['HN'] = centroids[id][1]

    def get_negative_intensities(self, systems):
        idx = []
        for i, system in enumerate(systems):
            if (system < 0).any():
                idx.append(i)
        return idx
    
if __name__ == '__main__':
    from FileHandler import FileHandler
#    from Visualization.Grouping import plot_2d_clustering
    from Pasta import Pasta
    from numpy import take
    group = Grouping()
    fh = FileHandler()
    
#    data4d = Load('/tmp/ubq4d.pkl')
#    data2d = Load('/tmp/ubq2d.pkl')
    identifiers, data2d = fh.pasta_from_2D_Nhsqc('/is/ei/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/multiple_test_files/Ubiquitin/Nhsqc-Ub.list')
    identifiers, data4d = fh.pasta_from_4D_Nhsqc('/is/ei/jhooge/EclipseWorkspaces/PASTA/PyPASTA/GMM/src/Classification/multiple_test_files/Ubiquitin/CCCANH-4D-ref.list')
    


#    x = data2d
#    y = take(data4d,(2,3),1)
#
#    tol = (0.3, 0.03)
#    n = 1
#
#    labels, centers = group.assign_amides(x, y, tol, n=n)
#    
#    systems = group.create_spin_systems(labels, data4d)
#    residues = map(group.define_spin_system, systems)
#    pasta_residues = [group.make_pasta_residue(residue[0],residue[1],i+1)
#                      for i, residue in enumerate(residues)]
#    pasta_residues = [residue for residue in pasta_residues
#                      if len(residue.shifts_i) or len(residue.shifts_im1)]
#    group.add_centroids_to_shift_dicts(centers,labels, pasta_residues)
#    fh.create_residue_list('multiple_test_files/Ubiquitin/Ub.pasta',pasta_residues)
#    
#    plot_2d_clustering(labels,centers,y)
