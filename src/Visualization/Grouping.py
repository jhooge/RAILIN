'''
Created on Apr 3, 2012

@author: jhooge
'''

def cluster_colors(cluster_dict, noseed=1):
    """
    n:= number of labels
    """
    import random
    
    colors = {}
    for key in cluster_dict.keys():
        noseed = noseed + 1
        random.seed(noseed)
        colors[key] = [float(random.randint(0, 255)) / 255,
                       float(random.randint(0, 255)) / 255,
                       float(random.randint(0, 255)) / 255]
    return colors

def clusters(data,centroid_ids):
    """
    @param data: mxn array with m number of observations and n number of features
    @type data: ndarray
    @param centroid_ids: indexes from kmeans2 clustering
    @param centroid_ids: list
    
    @return dic: dictionary where keys are class labels and valuesa are data values
    @rtype dic: dict
    """
    
    from numpy import where,array
    dic ={}
    labels = list(set(centroid_ids))
    for label in labels:
#        print label, where(centroid_ids == label)[0]
        dic[label] = array([data[i] for i in where(centroid_ids == label)[0]])
    
    return dic

def plot_2d_clustering(labels, centroids, data):
    """
    Visualizes 
    
    @param labels: label[i] is the index of the centroid the i'th observation is closest to.
    @type labels: ndarray
    @param centroids: A 'k' by 2 array of centroids found at the last iteration of k-means.
    @type centroids: ndarray
    @param data: A 'm' by 2 array of m observations in 2 dimensions.
    @type data: ndarray
    """
    
    from numpy import take, array
    import pylab as pl
    x = centroids
    y = data
#    idx, centroids = assign_amides(centroids, data2d)
#    print len(list(set(idx)))
#    print idx
    cluster_dict = clusters(y,labels)
    color_dict = cluster_colors(cluster_dict, noseed=1)
    
#    for clustid,clustered_data_points in cluster_dict.items():
#        print 'Centroid ',clustid, res[clustid]
#        print 'Data in cluster'
#        print clustered_data_points
        
    
    ## plot centroids
    clustered_data_points =[]
    colors = []
    for id in cluster_dict.keys():
#        print len(clustered_data_points), id
        clustered_data_points.append(centroids[id])
        colors.append(color_dict[id])
    clustered_data_points=array(clustered_data_points)
    pl.scatter(clustered_data_points[:, 0], clustered_data_points[:, 1], marker='o',alpha=.6, c=colors, s=50)
    clustered_data_points =[]
    
    ## plot clustered_data_points points colored by 
    for id in cluster_dict.keys():
        colors = []
        clustered_data_points = cluster_dict[id]
        i=0
        while i < clustered_data_points.shape[0]:
            colors.append(color_dict[id])
            i+=1
        pl.scatter(clustered_data_points[:,0], clustered_data_points[:,1], s=100, marker='+',alpha=.6, edgecolor=colors)
        
    pl.savefig('SpinSysCustering.png',dpi=500)

if __name__ == '__main__':
    pass