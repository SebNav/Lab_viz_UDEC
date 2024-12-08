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
from joblib import Parallel, delayed
from dipy.io.streamline import load_tractogram,save_tractogram

"""
Faltan añadir documentacion sobre los header y algunos problemas relacionados con el cambio de tck y trk al visuaizarlos en dsi studio

"""


T_inv=np.linalg.inv(np.load('tck2bundles.npy'))
T=np.load('tck2bundles.npy')


#funciones para aplicar el cambio de Tck a bundles

def apply_aff_bundle_parallel(bunIn,bunOut,t,nthreads):

	newPoints = Parallel(n_jobs=nthreads)(delayed(apply_aff_fiber)(f,t) for f in bunIn)
	return newPoints

def apply_aff_fiber(f,t):
	#print(idx)
	newfib=[]
	for p in f:
		pt=apply_aff_point(p,t)
		newfib.append((pt))
	return np.asarray(newfib,dtype=np.float32)


def apply_aff_point(inPoint,t):
	#Tfrm = N.array([[t[1,0], t[1,1], t[1,2], t[0,0]],[t[2,0], t[2,1], t[2,2], t[0,1]],[ t[3,0], t[3,1], t[3,2], t[0,2]],[0, 0, 0, 1]])
	Tfrm = t
	#Tfrm = N.array([[0.6, 0, 0, 0],[0, -0.6, 0, 5],[ 0, 0, -0.6, 0],[0, 0, 0, 1]])
	tmp = Tfrm * np.transpose(np.matrix(np.append(inPoint,1)))
	outpoint = np.squeeze(np.asarray(tmp))[0:3]
	return outpoint



def file_fchange(bundle_path,bundle_output_path,format_change,header="MNI152_T1_1mm.nii.gz",nthreads=-1):

    """
        Comando para transformar formatos de tractos para archivos
        
        Parámetros  
        ----------  
        bundle_path : str  
            Ruta al archivo .tck.  
            Ejemplo: "Pepito/tractografia/Tractografía.tck"  
        
        bundle_output_path : str  
            Ruta donde se guardará el archivo .bundles.  
            Ejemplo: "Pepito/tractografia/Tractografía_2.bundles"  
            
        format_change : str
            string que especifica el cambio de formato que se realizara
        
        header : str (opcional)  
            Los archivos .tck y .trk requieren un header para ser guardados.  
            Ruta a una imagen en formato .nii o .nii.gz de la cual se utilizará el header para guardar el archivo .tck.  
        
            *Para la visualización y uso general de tractos en formato .tck, el header no es relevante y se puede usar cualquier imagen (por ejemplo, la imagen MNI especificada).  
            Sin embargo, en casos más específicos es necesario ser más preciso, por lo que se recomienda utilizar una imagen T1 o de difusión del mismo sujeto en el mismo espacio que la tractografía.*  
        
        nthreads : int (opcional)  
            Número de núcleos del procesador a utilizar. Por defecto, se usan todos los núcleos disponibles.  
        
        Returns  
        -------  
        None.  
    """
    
    
    
    if   format_change == "bundles2tck":
        
        streamlines        = bt.read_bundle(bundle_path)
        streamlines        = np.array(streamlines)
        streamlines_tck    = transform_streamlines(streamlines, T_inv)
        
        centroids_tractogram_file = Tractogram(streamlines = streamlines_tck)
        centroids_tractogram_file.affine_to_rasmm = np.eye(4)
        centroids_tck = TF.TckFile(centroids_tractogram_file,header = {'timestamp':0})
        centroids_tck.save(bundle_output_path)
        
    elif format_change == "tck2bundles":
    
        #bbox_valid_check, si es True entrega un error si la posicion del tracto esta fuera del header especificado.
        #Si se quiere asegurar que esten "alineados" especifcar el header de una imagen del mismo sujeto y bbox_valid_check =True
        fibs_tck = load_tractogram(bundle_path,header,bbox_valid_check=False).streamlines 
        fibs_aff = apply_aff_bundle_parallel(fibs_tck,'',T,nthreads)
        bt.write_bundle(bundle_output_path,fibs_aff)
        
    
    elif format_change == "trk2tck":
        
        fibs_trk = load_tractogram(bundle_path,'same',bbox_valid_check=False)
        save_tractogram(fibs_trk, bundle_output_path,bbox_valid_check=False)
    
    elif format_change == "tck2trk":
        
        fibs_tck = load_tractogram(bundle_path,header,bbox_valid_check=False)
        save_tractogram(fibs_tck, bundle_output_path,bbox_valid_check=False)
    

    else:
        
        assert  format_change in ["tck2trk","bundles2tck","tck2bundles","trk2tck"],"Cambio de formato especificado no existe"
        


def folder_fchange(bundle_path,bundle_output_path,format_change,header="MNI152_T1_1mm.nii.gz",nthreads=-1,crear_carpeta=False):

    """
       Comando para transformar formatos de tractos para carpetas
        
        Parámetros  
        ----------  
        bundle_path : str  
            Ruta al archivo .tck.  
            Ejemplo: "Pepito/tractografia/Tractografía.tck"  
        
        bundle_output_path : str  
            Ruta donde se guardará el archivo .bundles.  
            Ejemplo: "Pepito/tractografia/Tractografía_2.bundles"  
            
        format_change : str
            string que especifica el cambio de formato que se realizara
        
        header : str (opcional)  
            Los archivos .tck y .trk requieren un header para ser guardados.  
            Ruta a una imagen en formato .nii o .nii.gz de la cual se utilizará el header para guardar el archivo .tck.  
        
            *Para la visualización y uso general de tractos en formato .tck, el header no es relevante y se puede usar cualquier imagen (por ejemplo, la imagen MNI especificada).  
            Sin embargo, en casos más específicos es necesario ser más preciso, por lo que se recomienda utilizar una imagen T1 o de difusión del mismo sujeto en el mismo espacio que la tractografía.*  
        
        nthreads : int (opcional)  
            Número de núcleos del procesador a utilizar. Por defecto, se usan todos los núcleos disponibles.  
        
        crear_carpeta : bool  
            Si es False, no se creará la carpeta donde se guardarán los fascículos.  
            Si es True, la carpeta especificada se creará si no existe.  

        Returns  
        -------  
        None.  
    """
    
        
    if not os.path.exists(bundle_output_path):
        
        if not crear_carpeta:
            assert  crear_carpeta ,"Error la carpeta donde se quieren guardar los fascículos no existe, modifique crear_carpetas a True para que el comando cree la carpeta"
        else:
            os.makedirs(bundle_output_path)
    
    archivos_format= os.listdir(bundle_path)
    archivos = [name.rstrip('.bundles').rstrip('.bundlesdata').rstrip('.trk').rstrip('.tck') for name in archivos_format]
    archivos = sorted(list(set(archivos)))
    
    assert  format_change in ["tck2trk","bundles2tck","tck2bundles","trk2tck"],"Cambio de formato especificado no existe"
    
    for arch in archivos:
        if   format_change == "bundles2tck":
            
            
            
            streamlines        = bt.read_bundle(bundle_path+arch+".bundles")
            streamlines        = np.array(streamlines)
            streamlines_tck    = transform_streamlines(streamlines, T_inv)
            
            centroids_tractogram_file = Tractogram(streamlines = streamlines_tck)
            centroids_tractogram_file.affine_to_rasmm = np.eye(4)
            centroids_tck = TF.TckFile(centroids_tractogram_file,header = {'timestamp':0})
            centroids_tck.save(bundle_output_path+arch+".tck")
            
        elif format_change == "tck2bundles":
        
            #bbox_valid_check, si es True entrega un error si la posicion del tracto esta fuera del header especificado.
            #Si se quiere asegurar que esten "alineados" especifcar el header de una imagen del mismo sujeto y bbox_valid_check =True
            fibs_tck = load_tractogram(bundle_path+arch+".tck",header,bbox_valid_check=False).streamlines 
            fibs_aff = apply_aff_bundle_parallel(fibs_tck,'',T,nthreads)
            bt.write_bundle(bundle_output_path+arch+".bundles",fibs_aff)
            
        
        elif format_change == "trk2tck":
            
            fibs_trk = load_tractogram(bundle_path+arch+".trk",'same',bbox_valid_check=False)
            save_tractogram(fibs_trk, bundle_output_path+arch+".tck",bbox_valid_check=False)
        
        elif format_change == "tck2trk":
            
            fibs_tck = load_tractogram(bundle_path+arch+".tck",header,bbox_valid_check=False)
            save_tractogram(fibs_tck, bundle_output_path+arch+".trk",bbox_valid_check=False)
    


#Ejemplos de uso
#file_fchange("test_data/fibras_bundles/41.bundles","test_data/41.tck",format_change="bundles2tck")
#file_fchange("test_data/fibras_tck/41.tck","test_data/41.bundles",format_change="tck2bundles")
#folder_fchange("test_data/fibras_tck/", "test_data/fibas_bundles_prueba/", "tck2bundles",crear_carpeta=True)
