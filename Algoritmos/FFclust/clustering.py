import itertools
import functools

import time
import threading
import string


import numpy as np
import sklearn.cluster
import dipy.segment.clustering

from joblib import Parallel, delayed

import bundleTools
import processing
import metric
import os
import subprocess
import segmentation as seg
import logging
import utils
import copy


def cluster_kmeans(x, k,random_state=0):
    kmeans = sklearn.cluster.MiniBatchKMeans(n_clusters=k,
                                             random_state=random_state)
    return kmeans.fit_predict(x), kmeans

def parallel_points_clustering(X, ks, n_jobs=-1, backend='threading'):
    """
    Compute MiniBatchKMeans for each point according to cluster size
    Do in parallel, conserve clusterers object.
    """
    results = Parallel(n_jobs=n_jobs, backend=backend)(delayed(cluster_kmeans)(x, k) for x, k in zip(X, ks))
    labels, clusterers = list(map(list, zip(*results)))
    labels = np.array(labels).T
    return labels, clusterers


def map_clustering(labels):
    """
    Create groups based on the labels that composed each element

    @arg labels: List of lists, each element in the dataset
    has m features, then
    if each feature has an associated label each element is the
    collection of such labels, in this case the labels come from
    a previous clustering step performed in each feature.
    """

    #Hash tuple instead of a string? should be better
#    string_labels_it = map(lambda i: ''.join(map(str, i)), labels)
    labels_it = map(tuple, labels)
    clusters = []
    d = {}
    
    #For every element check its label
    #if its new then create a new cluster
    #else assign to the cluster that has all elements with the
    #same label
    for i, label in enumerate(labels_it):
        if label in d:
            clusters[d[label]].append(i)
        else:
            clusters.append([i])
            d[label] = len(d)
    return clusters

#Maybe should be moved to processing?
def split_fibers(fibers,points):
    return np.array(np.split(fibers, len(points), axis=1))[:,:,0]

#Wrapper for dipy given that we have clusters and not centroids
def clusters_to_clustermap(clusters, fibers, ids=None):
    """
    Wrapper for dipy clustering interface
    @arg clusters : list of lists representing the clusters, each
    inner list contains the index of the elements in the dataset that
    conform this cluster.
    
    @arg fibers : dataset
    """
    
    cluster_map = dipy.segment.clustering.ClusterMapCentroid(fibers)
    if ids:
        for i, cluster in zip(ids, clusters):
            #compute a centroid as a mean of its elements

            centr = fibers[cluster].mean(axis=0) 
            c = dipy.segment.clustering.ClusterCentroid(centroid=centr,
                                                        id=i, indices=cluster,
                                                        refdata=fibers)
            cluster_map.add_cluster(c)
        return cluster_map
    for i, cluster in enumerate(clusters):
        #compute a centroid as a mean of its elements
        centr = fibers[cluster].mean(axis=0) 
        c = dipy.segment.clustering.ClusterCentroid(centroid=centr,
                                                    id=i, indices=cluster,
                                                    refdata=fibers)
        cluster_map.add_cluster(c)

    return cluster_map

def get_groups(labels, ngroups):
    groups = [[] for i in range(ngroups)]
    for i, l in enumerate(labels):
        groups[l].append(i)
    return groups
def new_clusters_in_group(group, refdata, centroids,
                          actual_clusters, max_dist, better_join=False):
    """
    group: indices of clusters that conform a group
    max_dist: maximum distance allowed for joining two clusters inside a group
    centroids: all centroids where to index and perform computation to
               determine wether or not join two clusters based on its centroids
               done after filtering out small clusters
    actual_clusters: similar to centroids, but has all the elements in a cluster
    """
    new_clusters_in_group = []
    processed = set()
    for i in range(len(group)):
        #New cluster with objects from first cluster c1
        c1 = group[i]
        if c1 in processed:
            continue
        #Maybe id=ind doesnt matter, can be changed accordingly in the
        #main thread to be ascending without the use of the ind variable.
        #id=ind
        ob = dipy.segment.clustering.ClusterCentroid(centroid=actual_clusters[c1].centroid,
                                                     indices=actual_clusters[c1].indices,
                                                     refdata=refdata)

        processed.add(c1)
        centroids_ids = []#no considerar c1
        for j in range(i+1, len(group)):
            c2 = group[j]
            if c2 in processed:
                continue
            dist1 = metric.direct_dist(centroids[c1], centroids[c2])
            dist2 = metric.direct_dist(centroids[c1], centroids[c2][::-1])
            #add all elements from c2 into ob
            #max distance

            #Add option to check for distance with other elements in the
            #joined group and join in that case
            if (dist1 <= dist2) and (dist1 <= max_dist):
                if better_join:
                    brk = False
                    for cid in centroids_ids:
                        if metric.max_dist(centroids[c2], centroids[cid]) > max_dist:
                            brk = True
                            break
                    if brk: continue
                processed.add(c2)
                centroids_ids.append(c2)
                for index in actual_clusters[c2].indices:
                    ob.assign(id_datum=index, features=ob.refdata[index])
                    ob.update()

            elif (dist2 < dist1) and (dist2 <= max_dist):
                if better_join:
                    brk = False
                    for cid in centroids_ids:
                        if metric.max_dist(centroids[c2], centroids[cid]) > max_dist:
                            brk = True
                            break
                    if brk: continue

                processed.add(c2)
                centroids_ids.append(c2)
                for index in actual_clusters[c2].indices:
                    ob.assign(id_datum=index, features=ob.refdata[index][::-1])
                    ob.update()
                
        #Should return clusters to be added at the end.
        new_clusters_in_group.append(ob)

    return new_clusters_in_group

#Write .bundles .bundlesdata of small (< 1 fiber) and long centroids (> 1 fiber) of clusters
def small_clusters_reassignment(clusters,min_size_filter, max_size_filter,input_dir, output_dir, reassignment_dir, threshold,refdata):

    small_clusters = clusters.get_small_clusters(max_size_filter) #Cambiar filtro a 5
    large_clusters = clusters.get_large_clusters(min_size_filter) #Cambiar filtro a 1
    small_centroids = np.asarray([x.centroid for x in small_clusters])
    large_centroids = np.asarray([x.centroid for x in large_clusters])

    #Create folders
    small_bundles_dir = output_dir+"/small_clusters"
    large_bundles_dir = output_dir+"/large_clusters"
    os.makedirs(small_bundles_dir)
    os.makedirs(large_bundles_dir)

    #Write bundles before reassignment
    for i,cluster in enumerate(small_clusters):
        small_bundles_file = small_bundles_dir+"/map_bundles_"+str(i)+".bundles"
        c = np.asarray(cluster[:])
        #bundleTools.write_bundle(small_bundles_file,c)

    for i,cluster in enumerate(large_clusters):
        large_bundles_file = large_bundles_dir+"/map_bundles_"+str(i)+".bundles"
        c = np.asarray(cluster[:])
        #bundleTools.write_bundle(large_bundles_file,c)        

    print('Número total de clusters: ',str(len(clusters)))
    print('Número total de clusters con ',str(max_size_filter),' fibras o menos: ',str(len(small_clusters)))
    print('Número total de clusters con más de ',str(min_size_filter),' fibras: ',str(len(large_clusters)))

    #Call segmentation method
    nPoints = 21
    print('...... Comenzando la segmentación: ...... ')
    reassignment = seg.segmentation(nPoints,threshold, large_centroids,small_centroids,len(small_centroids), len(large_centroids))

    count = 0
    num_fibers_reass = 0
    num_discarded = 0
    #Reassign small clusters to large clusters 
    for small_index,large_index in enumerate(reassignment):
        fibers = small_clusters[small_index].indices
        reassignment_file = reassignment_dir+"/reass_bundles_"+str(small_index)+"s"+str(large_index)+"l"+".bundles"
        if int(large_index)!=-1:
            #Hacer bucle para reasignar todas las fibras
            for fiber in fibers:
                centroid = large_clusters[large_index].centroid
                f = small_clusters[small_index].refdata[fiber]
                large_clusters[large_index].assign(id_datum=fiber, features= small_clusters[small_index].refdata[fiber])
                #large_clusters[large_index].update()
                num_fibers_reass += 1
            metric.recalc_centroid(large_clusters[large_index].indices, refdata)
            count+=1
            #bundleTools.write_bundle(reassignment_file,large_clusters[large_index])     
        else:
            if len(fibers)>2:
                recover_cluster = small_clusters[small_index]
                large_clusters.append(recover_cluster)
                recover_file = reassignment_dir+"/recover_bundles_"+str(small_index)+"to"+str(len(large_clusters)-1)+".bundles"
                #bundleTools.write_bundle(recover_file, recover_cluster) 
            else:
                num_discarded +=1
            #else sacar el cluster de small_clusters y añadirselo a large_clusters, SI len(fibers>1)

                c = np.asarray(cluster[:])

    print('Número de clusters reasignados: '+str(count)+' con un total de ',str(num_fibers_reass))
    print('Número de clusters descartados: '+str(num_discarded)+' con un total de ',str(num_discarded))

    return large_clusters




def join_by_group(clusters, actual_clusters, centroids, groups, max_dist=8.0, better_join=False):
    """
    Use a point-base clusterer to group centroids of clusters according to a certain point
    after doing this construct clusters that are the union of previous ones,
    based on wether their
    max distance(considering inverse too) is less than or equal to a
    threshold.
    Pass groups as arguments?
    """
    # actual_clusters = clusters.get_large_clusters(min_size_filter)
    # centroids = np.asarray([x.centroid for x in clusters.get_large_clusters(min_size_filter)])
    # centroids_points = centroids[:,point_index]
    # labels = clusterer.predict(centroids_points)
    # groups = get_groups(labels, ngroups=ngroups)
#    groups = [[] for _ in range(ngroups)]
 #   for i, c in enumerate(labels):
  #      groups[c].append(i)

    final_clusters = dipy.segment.clustering.ClusterMapCentroid(clusters.refdata)

    part_func = functools.partial(new_clusters_in_group, refdata=clusters.refdata, centroids=centroids,
                                  actual_clusters=actual_clusters, max_dist=max_dist, better_join=better_join)
    results = Parallel(n_jobs=-1, backend='threading')(delayed(part_func)(group) for group in groups)
    
    ind = 0
    for result in results:
        for x in result:
            x.id = ind
            ind += 1
            final_clusters.add_cluster(x)
            
    return final_clusters
    
#Classes for use with dipy ClusteringMap objects
class MapClustering(dipy.segment.clustering.Clustering):
    #Must return ClusterMap object
    #So we wrap our function around this
    def cluster(self, fibers, labels):
        clusters = map_clustering(labels)
        return clusters_to_clustermap(clusters, fibers)
def new_merged_clusters(to_merge, clusters_map, refdata, filtered_size):
    mapping = {}
    for i, cluster in enumerate(clusters_map.get_large_clusters(filtered_size)):
        mapping[i] = cluster.id
    #use indices into mapping to get the clusters_map correct indices

    cluster_map = dipy.segment.clustering.ClusterMapCentroid(refdata)
    for i, cluster in enumerate(to_merge):
        ind = [mapping[i] for i in cluster.indices]
        #clusters[ind]#clusters to be joined
        gen = (clusters_map[ind][j].indices for j in range(len(cluster)))
        indices = list(itertools.chain.from_iterable(gen))

        #indices = list(itertools.chain.from_iterable(clusters_map[ind]))
        #Improve centroid!
        c = dipy.segment.clustering.ClusterCentroid(centroid=cluster.centroid,
                                                    id=i, indices=indices, refdata=refdata)

        cluster_map.add_cluster(c)
    return cluster_map


def clusters_join_remove(clusters, n_clusters, max_dist=10, min_size=30):
    import collections
    """
    Take as argument the clusters, the number of kmeans clusters, the
    maximum distance alowed inside the kmeans clusters
    and the number of fibers required to separate between large and
    small clusters
    """

    new_clusters = dipy.segment.clustering.ClusterMapCentroid(clusters.refdata)
    
    midpoint = clusters.centroids[0].shape[0]//2
    #get middle point of every centroid
    midpoint_centroids = np.asarray(clusters.centroids)[:,midpoint]

    clusterer = sklearn.cluster.MiniBatchKMeans(n_clusters)
    labels = clusterer.fit_predict(midpoint_centroids)
    groups = get_groups(labels, n_clusters)

    new_clusters = clusters
    reagrupables = np.zeros(n_clusters)
    no_reagrupables = np.zeros(n_clusters)
    count = np.zeros(n_clusters)
    clusters_array = np.asarray(clusters)

    for i, group in enumerate(groups):
        not_joinable_ingroup = []
        diff = midpoint_centroids[group] - clusterer.cluster_centers_[i]
        midpoint_distances = np.sqrt((diff**2).sum(axis=1))
        mask = (midpoint_distances <= max_dist)
        reagrupables[i] = mask.sum()
        no_reagrupables[i] = (~mask).sum()
        joinable = np.asarray(group)[mask] #Clusters del grupo que pueden ser unidos
        #Add clusters that remain the same to list
        for idx in np.asarray(group)[~mask]:
            not_joinable_ingroup.append(idx)
        #numero de clusters con menos de 30 fibras que son agrupables
        count[i] = len(list(filter(lambda x: len(x) < min_size, (clusters_array)[joinable])))
                       
        #separar entre pequeños y grandes dentro del grupo
        mask_small = np.asarray(list(map(len, clusters_array[joinable]))) < min_size

        small_clusters_group = clusters_array[joinable][mask_small]
        large_clusters_group = clusters_array[joinable][~mask_small]
        
        small_centroids = np.asarray([x.centroid for x in small_clusters_group])
        large_centroids = np.asarray([x.centroid for x in large_clusters_group])
        merge_clusters = collections.defaultdict(list)
        #No large clusters to join in this group skip to next group
        #Or join among the small group
        if large_centroids.shape[0] == 0:
            for idx in joinable:
                not_joinable_ingroup.append(idx)
            continue

        #Compute which of the larges cluster is closer to the small ones
        for i, c in enumerate(small_centroids):
            #index returned is the cluster among the large_clusters_group which we join
#            print(c)
            dist, index = metric.max_dist_one_to_res(c, large_centroids)
            if dist <= max_dist:
                merge_clusters[index].append(i)#then we create clusters from dict
            else:
                not_joinable_ingroup.append(i)#Or just add data

        ind = 0
        #WRONG, MODIFYING ORIGINAL DATA
#        for idx in not_joinable_ingroup:
#            c = dipy.segment.clustering.ClusterCentroid(centroid=clusters_array[idx].centroid,
#                                                        id=ind, indices=clusters_array[idx].indices, refdata=clusters_array[idx].refdata)
#            ind += 1
#            new_clusters.add_cluster(c)
            
        for key, added_clusters in merge_clusters.items():
            c = dipy.segment.clustering.ClusterCentroid(centroid=large_clusters_group[key].centroid,
                                                    id=ind, indices=large_clusters_group[key].indices, refdata=large_clusters_group[key].refdata)
            for small_cluster in added_clusters:
                for index in small_clusters_group[small_cluster].indices:
                    c.assign(id_datum=index, features=c.refdata[index])
                    c.update()
            new_clusters.add_cluster(c)

        #Por cada cluster de los pequeños buscar en los grandes
        #Unir al cluster grande cuya maxima distancia usando centroides es la menor
        
        #seleccionar id clusters reagrupables
        #sum(filter(lambda x: x < 30, map(len, (clusters_array)[np.asarray(group)[mask]])))

    return new_clusters
#    return reagrupables, no_reagrupables, count
       # idx = np.asarray(group)[np.where(midpoint_distances > max_dist)]
       #new_clusters = np.delete(new_clusters, idx)
def joinable_clusters(clique, visited):
    """
    Visits the nodes of a clique and returns nodes that are not visited.
    In the clustering, the nodes represent the centroids of the clusters.
    """
    clusters_ids = []
    for node in clique:
        #some clusters may be joined already
        if not visited[node]:
            clusters_ids.append(node)
            visited[node] = True
    return clusters_ids

def clique_join(clusters, refdata, joined_bundles_dir, ident_clusters, object_dir, threshold):
    """
    Returns a list of new clusters(dipy.clustering.ClusterCentroid container)
    based on computing the distance between centroids and using this matrix as input
    to finding cliques, distances greater than the specified threshold are ignored
    """
    # cl = []

    # for c in clusters:
    #     ident = int(ident_clusters[str(c.indices)])
    #     if (ident == 6953 or ident == 7113):
    #         cl.append(c)

    import networkx
    centroids = np.asarray([x.centroid for x in clusters])
 
     #use numpy matrix dist
    max_dists = metric.matrix_dist(centroids, get_max=True)
    # if (len(cl)>0): 
    #     centroids2 = np.asarray([x.centroid for x in cl])
    #     max_dists2 = metric.matrix_dist(centroids2, get_max=True)
    #     max_dists2[max_dists2 > threshold] = 0
    #     print(max_dists2)
    #     network2 = networkx.from_numpy_matrix(max_dists2)
    #     cliques2 = sorted(networkx.find_cliques(network2), key=len, reverse=True)
    #     print(cliques2)

    max_dists[max_dists > threshold] = 0
    network = networkx.from_numpy_array(max_dists)
    cliques = sorted(networkx.find_cliques(network), key=len, reverse=True)


    visited = np.zeros(len(centroids), np.bool_)
    new_clusters = []
    cluster_names = []
    map_cluster_names = {}
    for clique in cliques:
        #each clique is potentially a cluster
        clusters_ids = joinable_clusters(clique, visited)

        if not clusters_ids:
            continue

        #need to copy indices to preserve number of streamlines, no idea why
        ob = dipy.segment.clustering.ClusterCentroid(centroid=clusters[clusters_ids[0]].centroid,
                                                     indices=clusters[clusters_ids[0]].indices[:],
                                                     refdata=clusters[clusters_ids[0]].refdata)


        #node value indexes into centroids
        if (len(clusters_ids)) > 1:
            for idx in clusters_ids[1:]:
                for index in clusters[idx].indices:
                    ob.assign(id_datum=index, features=refdata[index])
            ob.centroid = metric.recalc_centroid(ob.indices,refdata)

        new_clusters.append(ob)

        string_ids = ""
        for idx in clusters_ids:
            index_cluster = ident_clusters[str(clusters[idx].indices)]
            string_ids = string_ids + str(index_cluster)+"c"        

        cluster_names.append(string_ids)

        #c = np.asarray(ob[:])
        #bundleTools.write_bundle(join_bundles_file,c)                

    for i,cluster in enumerate(new_clusters):
        map_cluster_names[cluster_names[i]] = cluster           


    return map_cluster_names

def parallel_group_join_clique(clusters, groups, refdata,joined_bundles_dir,final_centroids_dir,ident_clusters,object_dir,thr_join):
    """
    Same as group_join_clique, using multiprocessing.dummy.Pool for threading
    Return clusters after distributing them in groups and joining with both distances and
    cliques inside those groups.
    
    groups: Each group if a list of indices of clusters which are in that group.
    """
    import multiprocessing.dummy
    
    final_clusters = dipy.segment.clustering.ClusterMapCentroid(refdata)
    ind = 0
    pool_input = []
    for group in groups:
        cls = [clusters[i] for i in group]
        pool_input.append(cls)
        
    p = multiprocessing.dummy.Pool()
    func = functools.partial(clique_join, refdata=refdata, joined_bundles_dir = joined_bundles_dir, ident_clusters = ident_clusters, object_dir = object_dir,threshold = thr_join)
    results = p.map(func, pool_input)

    final_output_filename = object_dir + '/final_clusters.txt' 
    utils.save_join_clusters_fibers(results=results, filename=final_output_filename) 

     
     #ESCRIBIR LOS CLUSTERS CON TRAZABILIDAD"
    # for cluster_names in results:
    #     for name in cluster_names:
    #         cluster = cluster_names[name]
    #         join_bundles_file = joined_bundles_dir+"/final_bundles_"+name+".bundles"
    #         c = np.asarray(cluster[:])
    #         bundleTools.write_bundle(join_bundles_file,c) 

     #ESCRIBIR LOS CLUSTERS SIN TRAZABILIDAD"
    final_centroids = []
    i = 0
    for cluster_names in results:
        for name in cluster_names:
            cluster = cluster_names[name]
            join_bundles_file = joined_bundles_dir+"/"+str(i)+".bundles"
            i+=1
            c = np.asarray(cluster[:])
            final_centroids.append(cluster.centroid)
            bundleTools.write_bundle(join_bundles_file,c) 

    bundleTools.write_bundle(final_centroids_dir + '/centroids.bundles',np.asarray(final_centroids))

    for new_clusters in results:
        for k in new_clusters:
            new_clusters[k].id = ind
            ind += 1
            final_clusters.add_cluster(new_clusters[k])

    return final_clusters

def group_join_clique(clusters, groups, refdata):
    """
    Return clusters after distributing them in groups and joining with both distances and
    cliques inside those groups.
    
    groups: Each group if a list of indices of clusters which are in that group.
    """
    final_clusters = dipy.segment.clustering.ClusterMapCentroid(refdata)
    ind = 0
    for group in groups:
        cls = [clusters[i] for i in group]
        new_clusters = clique_join(cls, refdata)
        for c in new_clusters:
            c.id = ind
            ind += 1
            final_clusters.add_cluster(c)

    return final_clusters
    
def dist_filter(clusters, thr_size, thr_dist, fibers):
	final_clusters = clusters.get_large_clusters(thr_size+1)
	count_discarded = 0
	small_clusters = clusters.get_small_clusters(thr_size)
	for c in small_clusters:
		bad_cluster = False
		matrix = metric.matrix_dist(np.asarray(c))
		for row in matrix:
			for element in row:       
				if (element>=thr_dist):
					bad_cluster = True
					count_discarded = count_discarded + 1
					break
			if bad_cluster == True:
				break
		if bad_cluster == False:
			final_clusters.append(c)
	print("Número de clusters descartados por el filtro de distancia máxima: ",count_discarded)
	return final_clusters
  