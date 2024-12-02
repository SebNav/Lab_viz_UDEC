"""
Created on Wed Oct 23 16:00:22 2024

@author: Sebastián Navarrete

Funciones creadas inicialmente por Claudio Roman Godoy; Paralelizadas por Cristobal Mendoza;
Optimizadas y reescritas por Sebastián Navarrete
"""


import numpy as np
import read_write_bundle as bt
from dipy.tracking.streamline import transform_streamlines
from nibabel.streamlines.tractogram import Tractogram
import nibabel.streamlines.tck as TF
import os

T_inv=np.linalg.inv(np.load('tck2bundles.npy'))


def file_bundles2tck(bundle_path,bundle_output_path):

    """
    Comando para transformar streamlines en formato bundles a tck
    Parametros
    ----------
    bundle_path : str
        Path al archivo .bundles.
        ejemplo = "Pepito/tractografia/Tractografía.bundles"
        ejemplo: .
    bundle_output_path :str
        Path donde se guardara el archivo .tck.
        ejemplo = "Pepito/tractografia/Tractografía_2.tck"
        
    Returns
    None.
    
    """
    
    streamlines        = bt.read_bundle(bundle_path)
    streamlines        = np.array(streamlines)
    streamlines_tck    = transform_streamlines(streamlines, T_inv)
    
    centroids_tractogram_file = Tractogram(streamlines = streamlines_tck)
    centroids_tractogram_file.affine_to_rasmm = np.eye(4)
    centroids_tck = TF.TckFile(centroids_tractogram_file,header = {'timestamp':0})
    centroids_tck.save(bundle_output_path)



def folder_bundle2tck(bundle_path,bundle_output_path,crear_carpetas=False):

    """
    Parametros
    ----------
    bundle_path : str
        Path a la carpeta con los archivos .bundles.
        ejemplo = "Pepito/fasciculos/"

    bundle_output_path :str
        Path a la carpeta donde se guardaran los distintos archivos en .tck.
        ejemplo = "Pepito/fasciculos/bundlestck/"
    
    crear_carpetas : bool
        Si crear_carpetas es False no se creara la carpeta donde se guardaran los fascículos.
        Por el contrario si es True la carpeta especificada se creara si no existe
        
    Returns
    None.
    
    """

    if not os.path.exists(bundle_output_path):
        
        if not crear_carpetas:
            assert  crear_carpetas ,"Error la carpeta donde se quieren guardar los fascículos no existe, modifique crear_carpetas a True para que el comando cree la carpeta"
        else:
            os.makedirs(bundle_output_path)
            

    archivos_bundles = os.listdir(bundle_path)
    archivos = [name.rstrip('.bundles').rstrip('.bundlesdata') for name in archivos_bundles]
    archivos = sorted(list(set(archivos)))
    
    for fasc in archivos:
        streamlines        = bt.read_bundle(bundle_path+fasc+".bundles")
        streamlines        = np.array(streamlines)
        streamlines_tck    = transform_streamlines(streamlines, T_inv)
        
        centroids_tractogram_file = Tractogram(streamlines = streamlines_tck)
        centroids_tractogram_file.affine_to_rasmm = np.eye(4)
        centroids_tck = TF.TckFile(centroids_tractogram_file,header = {'timestamp':0})
        centroids_tck.save(bundle_output_path+fasc+".tck")


def file_tck2bundles(bundle_path,bundle_output_path):

    """
    Comando para transformar streamlines en formato tck a bundles
    Parametros
    ----------
    bundle_path : str
        Path al archivo .bundles.
        ejemplo = "Pepito/tractografia/Tractografía.tck"
        ejemplo: .
    bundle_output_path :str
        Path donde se guardara el archivo .bundles.
        ejemplo = "Pepito/tractografia/Tractografía_2.bundles"
        
    Returns
    None.

    """
    
    streamlines        = bt.read_bundle(bundle_path)
    streamlines        = np.array(streamlines)
    streamlines_tck    = transform_streamlines(streamlines, T_inv)
    
    centroids_tractogram_file = Tractogram(streamlines = streamlines_tck)
    centroids_tractogram_file.affine_to_rasmm = np.eye(4)
    centroids_tck = TF.TckFile(centroids_tractogram_file,header = {'timestamp':0})
    centroids_tck.save(bundle_output_path)




folder_bundle2tck("test_data/fibras_bundles/", "test_data/fibras_tck_prueba/",crear_carpetas=True)

