# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 19:33:42 2025

@author: snava
"""

import numpy as np
import os
from joblib import Parallel, delayed
import read_write_bundle as bt
import nibabel as nb
import nibabel.streamlines.tck as TF
from nibabel.streamlines.tractogram import Tractogram
from dipy.tracking.streamline import transform_streamlines
from dipy.io.streamline import load_tractogram, save_tractogram
import argparse

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
        
        
def verificar_archivos_sujetos(folder):
    """
    Verifica que los archivos diff.nii.gz, diff.bval y diff.bvec estén presentes en cada carpeta de sujetos.
    """
    sujetos = os.listdir(folder)
    sujetos = sorted(sujetos)
    sujetos = [os.path.join(folder, suj) for suj in sujetos if os.path.isdir(os.path.join(folder, suj))]
    flag = False
    for suj in sujetos:
        archivos_requeridos = ["diff.nii.gz", "diff.bval", "diff.bvec","T1w_acpc_dc_restore_brain.nii.gz"]
        faltantes = [archivo for archivo in archivos_requeridos if not os.path.exists(os.path.join(suj, archivo))]
	
        if faltantes:
            flag = True 
            print(f"Error: En la carpeta {suj} faltan los siguientes archivos: {', '.join(faltantes)}")
            

    if flag:
        return False
    print("Todos los archivos requeridos están presentes en las carpetas de los sujetos.")
    return True

T_inv=np.linalg.inv(np.load('tck2bundles.npy'))
T=np.load('tck2bundles.npy')

# Funciones auxiliares (apply_aff_bundle_parallel, apply_aff_fiber, apply_aff_point, file_fchange)
# ... (Mantén estas funciones igual que en tu código original)

def main():
    # Configurar el parser de argumentos
    parser = argparse.ArgumentParser(
        description="""Pipeline para sub-dividir los fascículos del atlas DWM creado por Sebastian Navarrete sebastian.navarrete@biomedica.udec.cl \n"""
    )

    # Añadir argumentos
    parser.add_argument(
        '-folder', 
        type=str, 
        required=True, 
        help="Ruta a la carpeta donde se encuentran los datos de sujetos a procesar."
    )

    args = parser.parse_args()

    # Validar que la carpeta exista
    if not os.path.isdir(args.folder):
        print(f"Error: La carpeta '{args.folder}' no existe.")
        return


    # Listar sujetos
    Sujetos = os.listdir(args.folder)
    Sujetos = sorted(Sujetos)
    Sujetos = [args.folder + "/" + suj + "/" for suj in Sujetos]
    print("Sujetos a procesar:", Sujetos)
    
    
    
    assert os.path.exists("MNI152_T1_3mm_brain.nii.gz") , \
        "Falta archivo 'MNI152_T1_3mm_brain.nii.gz' imagen T1w en el espacio MNI152, este se utiliza para el procesamiento de registro de imagenes."
        

    if not verificar_archivos_sujetos(args.folder):
        return
    
    # Verificar si se desea realizar segmentación de fascículos

    assert os.path.exists("main_index") and os.path.exists("atlas_faisceaux_UDD/") and os.path.exists("atlas_faisceaux_UDD.txt"), \
                "Archivos necesarios para segmentación de fibras largas no encontrados."

    # Procesar cada sujeto
    for suj in Sujetos:
        print(f"Procesando sujeto: {suj}")
        os.chdir(suj)

        print("--Inicio de proceso de segmentación de fascículos--")
        T1w = nb.as_closest_canonical(nb.load("T1w_acpc_dc_restore_brain_diff_space.nii.gz"))
        T1w_header = T1w.header
        fibs_tck = load_tractogram("tractography_prob_sift_3M_21p.tck", T1w_header, bbox_valid_check=False).streamlines
        tractography = np.asarray(fibs_tck)


        print("Calculo de Fascículos de fibras largas")
        os.makedirs("DWM", exist_ok=True)
        os.makedirs("DWM/results_folder_sub_bundles_DWI", exist_ok=True)
        os.chdir("..")
        os.chdir("..")
        os.system("./main_index 21 "+ suj+"tractography_prob_sift_3M_21p_MNI.bundles subject atlas_faisceaux_UDD/ atlas_faisceaux_UDD.txt "+suj+"DWM/results_folder_sub_bundles "+suj+"DWM/indices_folder_sub_bundles")
        DWM_fasc_index = os.listdir(suj + "DWM/indices_folder_sub_bundles")
        for fasc_ind in DWM_fasc_index:
            nombre_fasc = fasc_ind.rstrip(".txt")
            indices = open(suj + "DWM/indices_folder_sub_bundles/" + fasc_ind).readlines()
            indices = [int(ind.rstrip("\n")[:-7]) for ind in indices]
            bundle = tractography[indices]
            centroids_tractogram_file = Tractogram(streamlines=bundle)
            centroids_tractogram_file.affine_to_rasmm = np.eye(4)
            centroids_tck = TF.TckFile(centroids_tractogram_file, header={'timestamp': 0})
            centroids_tck.save(suj + "DWM/results_folder_sub_bundles_DWI/" + nombre_fasc[:] + ".tck")
        os.chdir(suj)

        os.chdir("..")
        os.chdir("..")


if __name__ == "__main__":
    main()
