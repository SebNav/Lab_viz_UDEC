import argparse
import logging
import time
import os

import numpy as np
import dipy.segment.clustering
import joblib
import matplotlib.image as mpimg
import bundleTools
import clustering 
import metric
import shutil
from utils import save_clusters, save_clusters_fibers, save_clusters_centroids


parser = argparse.ArgumentParser(description='Perform clustering on a dataset of streamlines')

parser.add_argument('--points', nargs='*', type=int,
                    help='Points to be used in map clustering')

parser.add_argument('--ks', nargs='*', type=int,
                    help='Number of clusters to be used for each point in K-Means for map')

parser.add_argument('--thr-segmentation', type=float, default=6,
                    help='Minimum threshold for segmentation')

parser.add_argument('--thr-join', type=float, default=6,
                    help='Minimum threshold for join')

parser.add_argument('--output-directory',
                    help='Directory where to place all output, including images, and pickled files')

parser.add_argument('--infile', help='Input streamlines file')
args = parser.parse_args()

# =============================================================================
# 
# args.points = [0, 3, 10, 17, 20]
# args.ks = [300, 200, 200, 200, 300]
# args.infile = "../../019/019/output/019_21p_filtered_trm.bundles"
# args.output_directory = "../../019/result/FFclust/th_6/"
# args.thr_segmentation = 6
# args.thr_join = 6
# 
# =============================================================================
if not args.points:
    points = np.arange(21)
else:
    points = args.points
    
if len(args.ks) == len(points):
    ks = args.ks
else:
    ks = [args.ks[0] for _ in range(len(points))]

    
fibers = np.asarray(bundleTools.read_bundle(args.infile))

#Root work directory creation, output should be placed here
if args.output_directory:
    work_dir = args.output_directory
else:
    work_dir = 'output_points{}_ks{}_'.format(str(points), str(ks)) + args.infile

if os.path.exists(work_dir):   
  shutil.rmtree(work_dir)

os.mkdir(work_dir)   

final_bundles_dir = work_dir + '/FinalBundles'
os.makedirs(final_bundles_dir, exist_ok=True)

final_centroids_dir = work_dir + '/FinalCentroids'
os.makedirs(final_centroids_dir, exist_ok=True)

map_bundles_dir = work_dir + '/MapBundles'
os.makedirs(map_bundles_dir, exist_ok=True)
individual_map_bundles_dir =  map_bundles_dir + '/IndividualBundles'
os.makedirs(individual_map_bundles_dir, exist_ok=True)

reassigned_bundles_dir = work_dir + '/ReassignedBundles'
os.makedirs(reassigned_bundles_dir, exist_ok=True)


object_dir = work_dir + '/output'
os.makedirs(object_dir, exist_ok=True)

map_output_filename = object_dir + '/clusters_map.txt'

#New logging file each time, changing filemode to a should append instead
logging.basicConfig(filename=work_dir+'/info.log', filemode='w', level=logging.INFO)


if __name__ == '__main__':
    colormap = mpimg.imread('colors256.jpg')[0]
    inicio_tiempo = time.time()
    t1 = time.time()
    X = clustering.split_fibers(fibers[:,points,:],points)
    labels, clusterers = clustering.parallel_points_clustering(X=X, ks=ks)
    logging.info('Tiempo Kmeans: {}'.format(time.time() - t1))

    t1 = time.time()
    m = clustering.MapClustering()
    map_clusters = m.cluster(fibers, labels)
    logging.info('Tiempo Map: {}'.format(time.time() - t1))
    

    final_output_filename = object_dir + '/map_clusters.txt' 
    save_clusters(dataset=fibers, clusters=map_clusters, filename=map_output_filename)

    with open(work_dir+'/stats.txt', 'w') as f:
        f.write('Number of clusters in map_clusters: {}\n'.format(len(map_clusters)))
        f.write('Number of fibers in map_clusters: {}\n'.format(sum(map_clusters.clusters_sizes())))

    for p,k, clusterer in zip(points, ks, clusterers):
        joblib.dump(clusterer, object_dir + '/clusterer-Point{}-k{}.pkl'.format(p,k))

    """Write .bundles .bundlesdata of small (< 1 fiber) and long centroids (> 1 fiber) of clusters
     and reassign small clasters to large clusters.
     Input:
      ouput: path of output file directory to write the clusters
      min_size_filter = minimum size to get the largest clusters
      max_size_filter = maximum size to get the clusters with only 1 fiber
      input_dir = path of the input dir that contains the results of segmetation/reassignation
      threshold = minimum distance for the segmentation's method
     """
    #Compenzando la reasignación
    t1 = time.time()
    actual_clusters = clustering.small_clusters_reassignment(clusters=map_clusters,
                                                          min_size_filter=6,
                                                          max_size_filter=5,
                                                          input_dir = 'segmentation/bundles/result/parallelFastCPU',
                                                          output_dir=individual_map_bundles_dir,
                                                          reassignment_dir=reassigned_bundles_dir,
                                                          threshold = args.thr_segmentation,
                                                          refdata = fibers)
	
    #actual_clusters = map_clusters.get_large_clusters(6)
 
    ident_clusters={}
    for i,c in enumerate(actual_clusters):
        ident_clusters[str(c.indices)] = i

    logging.info('Tiempo Segmentacion {}'.format(time.time() - t1))

    t1 = time.time()
    point_index = 10
    ngroups = ks[len(ks)//2]
    
    #actual_clusters = map_clusters.get_large_clusters(args.filter_min_size)
    centroids = np.asarray([x.centroid for x in actual_clusters])
    centroids_points = centroids[:,point_index]
    clusterer=clusterers[len(points)//2]
    labels = clusterer.predict(centroids_points)
    
    groups = clustering.get_groups(labels, ngroups=ngroups)
    joined_clusters = clustering.parallel_group_join_clique(actual_clusters, groups, fibers,final_bundles_dir,final_centroids_dir,ident_clusters,object_dir,args.thr_join)
    

    final_clusters = joined_clusters
    tiempo_final = time.time() - inicio_tiempo
    with open(work_dir+'/tiempo.txt', 'w') as f:
        f.write(str(tiempo_final))
   # logging.info('Tiempo total {}'.format(time.time() - inicio_tiempo))
    #calculo de tiempo
    final_centroids_filename = object_dir + '/final_centroids.txt' 
    save_clusters_centroids(clusters=final_clusters, filename=final_centroids_filename)

    # final_centroids = np.asarray([x.centroid for x in final_clusters])
    # bundleTools.write_bundle(final_centroids_dir + '/centroids.bundles',final_centroids)
