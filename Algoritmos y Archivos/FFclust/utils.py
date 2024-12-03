import collections
import numpy as np
import clustering

def save_clusters(dataset, clusters, filename):
    """
    Generate output file for clusters, the output file consists of
    the index to which cluster each element in the dataset(streamlines)
    belongs to, an index of -1 indicates that the element isnt assigned anywhere
    """
    output_array = -1*np.ones(len(dataset), np.int32)
    for cluster in clusters:
        output_array[cluster.indices] = cluster.id
    np.savetxt(X=output_array, fname=filename, delimiter='\n', fmt='%d')

def save_clusters_fibers(clusters, filename):
    text=""
    file = open(filename,"w")
    for cluster in clusters:
        for fiber_index in cluster.indices:
            text += str(fiber_index)+" "
        file.write(text+"\n")
        text =""
    file.close()

def save_clusters_centroids(clusters, filename):
    text=""
    file = open(filename,"w")
    for cluster in clusters:
        text += str(cluster.centroid)+" "
        file.write(text+"\n")
        text =""
    file.close()

def save_join_clusters_fibers(results, filename):
    text=""
    file = open(filename,"w")
    for clusters_names in results:
        for name in clusters_names:
            text = text + name+ "-"
            for fiber_index in clusters_names[name].indices:
                text = text + str(fiber_index)+" "
            file.write(text+"\n")
            text =""
    file.close()      

def load_clusters(dataset, filename, not_assigned=False):
    """
    From file where each line is the id of the cluster to which each streamline belongs
    generate a dipy ClusterMapCentroid to use
    """

    cls_dict = collections.defaultdict(list)
    with open(filename, 'r') as f:
        for i,id in enumerate(f):
            id = int(id)
            cls_dict[id].append(i)
            
    if -1 in cls_dict:
        no_clusters = cls_dict[-1]
        del cls_dict[-1]


    clusters = clustering.clusters_to_clustermap(list(cls_dict.values()), dataset, cls_dict.keys())
    
    if not_assigned:
        return clusters, no_clusters

    else:
        return clusters
