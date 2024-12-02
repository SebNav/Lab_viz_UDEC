import numpy as np
import math
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
from dipy.tracking.streamline import length as streamline_length
from random import randint
import matplotlib.patches as mpatches
import operator


def euclidean_distance(x1,y1,z1,x2,y2,z2):
    return (((x1-x2)*(x1-x2))+((y1-y2)*(y1-y2))+((z1-z2)*(z1-z2)))

def is_flipped(fiber1,fiber2):
    p0_f1 = fiber1[0]
    p0_f2 = fiber2[0]
    p20_f2 = fiber2[20]

    ed_direct = euclidean_distance(p0_f1[0],p0_f1[1],p0_f1[2],p0_f2[0],p0_f2[1],p0_f2[2])
    ed_flip = euclidean_distance(p0_f1[0],p0_f1[1],p0_f1[2],p20_f2[0],p20_f2[1],p20_f2[2])
    if (ed_direct>ed_flip):
        return True
    else:
        return False

def recalc_centroid(indices, refdata):
    fibers = {}
    for i in indices:
        f = refdata[i]
        f_lenght = streamline_length(f)
        fibers[i] = f_lenght
    sorted_items = sorted(fibers.items(), key=operator.itemgetter(1))
    sorted_fibers = []
    for item in sorted_items:
        sorted_fibers.append(item[0])
    index_60length = math.ceil(len(sorted_fibers)*0.6) - 1
    index_80length = math.ceil(len(sorted_fibers)*0.8) - 1
    longest_fibers = []
    while index_60length <= index_80length:
        longest_fibers.append(sorted_fibers[index_60length])
        index_60length+=1
    number_fibers_10percent = math.ceil(len(longest_fibers)*0.1)
    fibers_10percent = []
    for i in range(number_fibers_10percent):
        rand_index = randint(1, number_fibers_10percent)
        fibers_10percent.append(longest_fibers[rand_index])

    fibers_data = []
    for index in fibers_10percent:
        fiber = refdata[index]
        fibers_data.append(fiber)
    matrix = matrix_dist(streamlines = np.asarray(fibers_data))

    min_sum = 0
    selected_index = 0
    for i in range(len(matrix)):
        sum = 0
        for j in range(len((matrix)[0])):
            sum += matrix[i][j]
        if (sum < min_sum):
            min_sum = sum
            selected_index = i

    centroid = fibers_data[i]
    return centroid


def recalc_centroid_inv(fiber,fiber2):
    point_inv = 20
    new_fiber = fiber
    for i in range(21):
        x = fiber[i][0]
        x_inv = fiber2[point_inv][0]
        y = fiber[i][1]
        y_inv = fiber2[point_inv][1]
        z = fiber[i][2]
        z_inv = fiber2[point_inv][2]
        x_mean = (x + x_inv) / 2
        y_mean = (y + y_inv) / 2
        z_mean = (z + z_inv) / 2
        new_fiber[i][0] = x_mean
        new_fiber[i][1] = y_mean
        new_fiber[i][2] = z_mean
        point_inv -= point_inv
    return new_fiber

def recalc_centroid_direct(fiber,fiber2):
    new_fiber = fiber
    for i in range(21):
        x1 = fiber[i][0]
        x2 = fiber2[i][0]
        y1 = fiber[i][1]
        y2 = fiber2[i][1]
        z1 = fiber[i][2]
        z2 = fiber2[i][2]
        x_mean = (x1 + x2) / 2
        y_mean = (y1 + y2) / 2
        z_mean = (z1 + z2) / 2
        new_fiber[i][0] = x_mean
        new_fiber[i][1] = y_mean
        new_fiber[i][2] = z_mean
    return new_fiber


def direct_dist(v1, v2):
    """
    Returns the maximum direct distance between two streamlines
    """
    return math.sqrt(((v1-v2)**2).sum(axis=1).max())
def max_dist(v1, v2):
    """
    Returns the maximum distancia between two streamlines independent of orientation
    by obtaining the direct distance two times, one of them with a flipped streamline, finally choosing the minimum  between the two.
    """
    return min(direct_dist(v1, v2), direct_dist(v1, v2[::-1]))

def intra_cluster_distances(cluster, dataset, filename, metric=max_dist, save=False):
    """
    Calculate distance between all elements in a cluster using
    a metric (by default max distance), then save the resulting matrix in 
    """
    
    #callable should be tested for metric.
        
    c = cluster.indices
    n = len(c)
    dist_matrix = np.zeros((n, n))
        #Upper triangular matrix only
    for i in range(n):
        for j in range(i,n):
            dist_matrix[i][j] = metric(dataset[c[i]], dataset[c[j]])
    #Since metric(x,y) = metric(y,x) we can make it so that the lower
    #triangular is the same as the upper
    y = np.triu(dist_matrix) + np.triu(dist_matrix,k=1).T
    # for i in range(n):
    #     for j in range(n):
    #         dist_matrix[i][j] = metric(dataset[c[i]], dataset[c[j]])

    if save:
        np.savetxt(X=y, fname=filename, delimiter=' ', fmt='%f', newline='\n')
    return y
def clusters_metrics(clusters, dataset, output_dir='./', metric=max_dist,n_jobs=1,save=False):

    if output_dir[-1] != '/':
        output_dir.append('/')
    results = Parallel(n_jobs=-1, backend='threading')(delayed(intra_cluster_distances)
                                             (cluster, dataset, output_dir+str(cluster.id)+'-intra-cluster-distance.txt', max_dist,save) for
                                             cluster in clusters)
    return results
    #should return the matrices per cluster for further analysis
    
def max_dist_one_to_res(fiber, cluster):
    """
    Compute squared maximum distance from a fiber to all the fibers in a cluster
    this maximum distance considers both direct and inverse calculation

    Returns maximum_distances, index where it was found
    """
    f2 = np.asarray([fiber, fiber[::-1]])
    #    return ((f2 - cluster)**2).sum(axis=3).max(axis=2).min(axis=0).max()
    #retorna índice del menor
    distances = ((f2[:,np.newaxis] - cluster)**2).sum(axis=3).max(axis=2).min(axis=0)
    max_dist_ind = np.argmax(distances)
    return np.sqrt(distances[max_dist_ind]), max_dist_ind


def all_pairs_max_dist(cluster):
    """
    Returns upper triangular squared max distance matrix between every pair
    avoids extra computation by considering simmetry in matrix.

    Maximum in the whole cluster can be found by applying max in the output and then
    taking the square root

    """
    import itertools
    if len(cluster) < 2:
        return np.zeros(1)
    #Compute difference between fiber a and every other pair, just upper triangular.
    d1 = np.asarray(list(itertools.chain.from_iterable(((cluster - cluster[i])[i+1:]**2
                                                        for i in range(cluster.shape[0]-1)))), dtype=np.float32)
    d2 = np.asarray(list(itertools.chain.from_iterable(((cluster[:,::-1] - cluster[i])[i+1:]**2
                                                        for i in range(cluster.shape[0]-1)))), dtype=np.float32)
    dist_dir = (d1).sum(axis=2).max(axis=1)
    dist_inv = (d2).sum(axis=2).max(axis=1)
    return np.asarray([dist_dir, dist_inv]).min(axis=0)

def dask_matrix_dist(streamlines, get_max=True, get_mean=False):
    """
    Return distance matrix between streamlines
    By default returns only maximum distances, can also return
    the mean distance matrix.

    Uses dask library to avoid memory errors.
    Input streamlines should be numpy array with shape (num_streamlines, num_points, point_dimension)
    """
    import dask.array as da
    num_streamlines, num_points, dim = streamlines.shape
    x = da.from_array(streamlines, chunks=(500, num_points, dim))
    distances = ((da.stack((x,x[:,::-1]))[:,None]-x[:,None])**2).sum(axis=4)
    max_distances = da.sqrt(distances.max(axis=3).min(axis=0))
    mean_distances = da.sqrt(distances.mean(axis=3).min(axis=0))

    if get_max and get_mean:
        max_dists, mean_dists = da.compute(max_distances, mean_distances)
        return max_dists, mean_dists
    elif get_max:
        return max_distances.compute()
    elif get_mean:
        return mean_distances.compute()
    else:
        print('error, should return atleast one matrix')
        
def matrix_dist(streamlines, get_max=True, get_mean=False):
    """
    Compute distance between all streamlines, can return both max/min distance matrices
    """
    x = streamlines
    #output of this computation is a (num_streamlines, num_streamlines, num_points) array
    #first stack both direct and reversed in an array -> (2, num_streamlines, num_points, dim)
    #using [:,None] create new axis to broadcast computation and perform difference between
    #every pair of streamlines(both direct and reverse), then perform x^2 + y^2 + z^2
    distances = ((np.stack((x,x[:,::-1]))[:,None]-x[:,None])**2).sum(axis=4)

    if get_max and get_mean:
        max_distances = np.sqrt(distances.max(axis=3).min(axis=0))
        mean_distances = np.sqrt(distances.mean(axis=3).min(axis=0))
        return max_distances, mean_distances
    elif get_max:
        return np.sqrt(distances.max(axis=3).min(axis=0))
    elif get_mean:
        return  np.sqrt(distances.mean(axis=3).min(axis=0))
    else:
        print('error, should return atleast one matrix')

