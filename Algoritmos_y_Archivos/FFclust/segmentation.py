import shutil
import os
import bundleTools
import ctypes

_seg = ctypes.CDLL('segmentation_clust_v1.1/segmentation.so')

_seg.segmentation.argtypes = (ctypes.c_uint,ctypes.POINTER(ctypes.c_char),ctypes.POINTER(ctypes.c_char),ctypes.c_float, ctypes.POINTER(ctypes.c_char),ctypes.c_uint,ctypes.c_uint)



#The result value is a list whose indices are the indices of small clusters, each index contains the index of it's large_cluster reassignment
def segmentation(nPoints,threshold, large_clusters_fibers,small_clusters_fibers,nfibers_subject,nfibers_atlas):
    global _seg

    nBundles = len(large_clusters_fibers)
    n_small_clusters_fibers = len(small_clusters_fibers)

    #Create folders
    bundles_dir = "segmentation_clust_v1.1/bundles"
    if os.path.exists(bundles_dir):
        shutil.rmtree(bundles_dir)
        os.makedirs(bundles_dir)
    else:
        os.makedirs(bundles_dir)

     #Write file with the largest cluster's 
    large_centroids_file = bundles_dir+"/large_clusters.bundles"
    bundleTools.write_bundle(large_centroids_file, large_clusters_fibers)

    #Write file with the smallest cluster's 
    small_centroids_file = bundles_dir+"/small_clusters.bundles"
    bundleTools.write_bundle(small_centroids_file, small_clusters_fibers)

    ouputWorkDirectory = bundles_dir+'/result'

    s_subject =small_centroids_file.encode()
    s_atlas= large_centroids_file.encode()
    s_output= ouputWorkDirectory.encode()
    _seg.segmentation.restype = ctypes.POINTER(ctypes.c_int)
    result = _seg.segmentation(ctypes.c_uint(nPoints), s_subject, s_atlas,
                                                                          ctypes.c_float(threshold),s_output,ctypes.c_uint(nfibers_subject),ctypes.c_uint(nfibers_atlas))
    result_list = [result[i] for i in range(n_small_clusters_fibers)]
    _seg.freeme(result)
    shutil.rmtree(bundles_dir)
    return(result_list)