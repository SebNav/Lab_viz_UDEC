import itertools
import string
import numpy as np


def slope(curve):
    """
    Slope between points of a 2D curve, if a curve has
    n points returns an array with n - 1 points.
    """
    a = curve[:-1]
    b = np.roll(curve, -1, axis=0)[:-1]

    dxdy = b.T - a.T
    return dxdy[1] / dxdy[0]


def string_prefix(data):
    return np.array([string.ascii_lowercase[i] + str(d)
                     for i, d in enumerate(data)])
def string_slopes(slopes):
    return np.array([string.ascii_lowercase[i] + str(slope)
                     for i, slope in enumerate(np.nditer(slopes))])

def get_slopes(curves, sign=True):
    """
    Get multidimensional array of slopes
    given an array of curves.
    """
    if sign:
        return np.array([np.sign(slope(x)) for x in curves])
    else:
        return np.array([slope(x) for x in curves])


#turn into generator?
def str_prefix_slopes(slopes):
    """
    Returns slopes with a prefix based on it's position, example:
    [3,2,1] -> ['a3','b2','c1']
    """
    return np.fromiter((string.ascii_lowercase[i] + str(value)
                        for i,value in enumerate(slopes)), dtype=(str, len(slopes)))


#transform curve to slopes with a string prefix
def curve_to_str_slopes(curve, sign=True):
    if sign:
        return str_prefix_slopes(np.sign(slope(curve)))
    else:
        return str_prefix_slopes(slope(curve))

#general function to transform an array of curves into
#data for minhashing
#here curves is just an iterable
def encode_curves(curves, encoding_function):
    #yield np.array([encoding_function(curve) for curve in curves])
   for curve in curves:
       yield encoding_function(curve)

def flatten_lists(clusters_streamlines):
    return list(itertools.chain.from_iterable(clusters_streamlines))

def streamline_len(s):
    #aproximate arc length
    return sum(np.linalg.norm(s[:-1] - s[1:], axis=1))

def get_long_clusters(clusters, min_len):
    """
    Return clusters where the centroids has at least min_len length
    """
    return [cluster for cluster in clusters
            if streamline_len(cluster.centroid) >= min_len]
def get_short_clusters(clusters, max_len):
    """
    Return clusters where the centroids has at most max_len length
    """
    return [cluster for cluster in clusters
           if streamline_len(cluster.centroid <= max_len)]
def get_clusters_len(clusters, min_len, max_len):
    """
    Clusters where length(centroids) is between min and max len.
    """
    lens = [streamline_len(cluster.centroid) for cluster in clusters]
    return [cluster for i, cluster in enumerate(clusters)
           if lens[i] >= min_len and lens[i] < max_len]
