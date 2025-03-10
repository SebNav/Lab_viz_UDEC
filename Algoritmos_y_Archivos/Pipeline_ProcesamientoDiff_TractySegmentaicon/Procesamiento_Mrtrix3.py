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
        description="""Pipeline de preprocesamiento y calculo de tractografía creado por Sebastian Navarrete sebastian.navarrete@biomedica.udec.cl basado en el software Mrtrix3. 
         Este pipeline lleva a cabo el procesamiento de imagenes de difusión, calculo de modelo de difusion (DTI o CSD), calculo de Tractograma (Det o Prob), filtrado de tractograma, segmentación de fasciculos de fibras cortas y largas, ademas del registro de imagenes.
         Para que este pipeline funcione correctamente se necesitan tener instalada en el mismo ambiente el software Mrtrix3, fsl y ants, Ademas de una serie de archivos, para mayor instalacion en la instalacion y archivos dirigase a https://github.com/SebNav/Lab_viz_UDEC/tree/main/Algoritmos_y_Archivos/Pipeline_Calculo_Tractograma \n"""
    )

    # Añadir argumentos
    parser.add_argument(
        '-folder', 
        type=str, 
        required=True, 
        help="Ruta a la carpeta donde se encuentran los datos de sujetos a procesar."
    )
    parser.add_argument(
        '-segmentacion', 
        nargs='+', 
        choices=['SWM', 'DWM'], 
        default=None, 
        help="Tipo de segmentación de fascículos: SWM (fibras cortas), DWM (fibras largas), o ambas."
    )
    parser.add_argument(
        '-diff_model', 
        type=str, 
        choices=['CSD'], 
        default='CSD', 
        help="Modelo de difusión a utilizar (por defecto: CSD)."
    )
    parser.add_argument(
        '-flag_tracto', 
        type=str, 
        choices=['prob'], 
        default='prob', 
        help="Tipo de algoritmo de tractografía a utilizar (por defecto: prob)."
    )
    parser.add_argument(
        '-phase_encoding', 
        type=str, 
        default='j-', 
        help="Dirección de phase encoding (por defecto: j-)."
    )
    # Parsear los argumentos
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
    if args.segmentacion is not None:
        if "SWM" in args.segmentacion:
            assert os.path.exists("main_index") and os.path.exists("AtlasRo") and os.path.exists("AtlasRo_estables.txt"), \
                "Archivos necesarios para segmentación de fibras cortas no encontrados."
        if "DWM" in args.segmentacion:
            assert os.path.exists("main_index") and os.path.exists("atlas_faisceaux") and os.path.exists("atlas_faisceaux.txt"), \
                "Archivos necesarios para segmentación de fibras largas no encontrados."

    # Procesar cada sujeto
    for suj in Sujetos:
        print(f"Procesando sujeto: {suj}")
        os.chdir(suj)

        if os.path.exists("tractography_prob_sift_3M_21p_MNI.tck"):
            print("Archivos ya calculados para este sujeto.")
            os.chdir("..")
            os.chdir("..")
            continue

        # Preprocesamiento de datos de difusión
        os.system("mrconvert diff.nii.gz dwi.mif -fslgrad diff.bvec diff.bval -force")
        print("--DWI Denoise--")
        os.system("dwidenoise dwi.mif dwi_den.mif -force")
        print("--DWI deGibbs--")
        os.system("mrdegibbs dwi_den.mif dwi_den_unr.mif -force")
        print("--DWI dwifslpreproc--")
        os.system(f"dwifslpreproc dwi_den_unr.mif dwi_preproc_unbiased.mif -rpe_none -pe_dir {args.phase_encoding}")
        print("--DWI mask--")
        os.system("dwi2mask dwi_preproc_unbiased.mif dwi_mask.mif")

        # Preprocesamiento de imagen T1w
        print("--Alineación T1w con DWI--")
        os.system("mrconvert T1w_acpc_dc_restore_brain.nii.gz T1w_acpc_dc_restore_brain.mif -force")
        os.system("dwiextract dwi_preproc_unbiased.mif - -bzero -force | mrmath - mean bzero.nii.gz -axis 3 -force")
        os.system("flirt -in bzero.nii.gz -ref T1w_acpc_dc_restore_brain.nii.gz -omat diff2struct_fsl.mat -dof 6")
        os.system("transformconvert diff2struct_fsl.mat bzero.nii.gz T1w_acpc_dc_restore_brain.mif flirt_import diff2struct_mrtrix.txt -force")
        os.system("mrtransform T1w_acpc_dc_restore_brain.mif -linear diff2struct_mrtrix.txt -inverse T1w_acpc_dc_restore_brain_diff_space.nii.gz -force")
        print("--Calculo de 5TT--")
        os.system("5ttgen fsl T1w_acpc_dc_restore_brain_diff_space.nii.gz 5TT.nii.gz -nthreads 24 -premasked -force")

        # Calculo de modelo de difusión
        if args.diff_model == "CSD":
            print("--Inicio de CSD, dhollander; msmt_csd--")
            os.system("dwi2response dhollander dwi_preproc_unbiased.mif response_wm.txt response_gm.txt response_csf.txt -voxels RF_voxels.mif -nthreads 24 -force")
            os.system("dwi2fod msmt_csd dwi.mif response_wm.txt wmfod.mif response_gm.txt gm.mif response_csf.txt csf.mif -mask dwi_mask.mif")

        # Calculo de tractografía
        if args.flag_tracto == "prob":
            print("--Inicio de calculo de tractografía Probabilistica--")
            os.system("tckgen -algorithm iFOD2 wmfod.mif output.tck -act 5TT.nii.gz -backtrack -crop_at_gmwmi -seed_image dwi_mask.mif -maxlength 250 -minlength 40 -select 10M -cutoff 0.06")
            os.system("tcksift output.tck wmfod.mif -act 5TT.nii.gz -term_number 3M output_sift.tck")
            os.system("tckresample output_sift.tck tractography_prob_sift_3M_21p.tck -num_points 21")
            os.system("rm output.tck")
            os.system("rm output_sift.tck")

        # Tractografía a espacio MNI
        print("--Registro de tractografía al espacio MNI152--")
        os.system("flirt -ref ../../MNI152_T1_3mm_brain.nii.gz -in T1w_acpc_dc_restore_brain_diff_space.nii.gz -omat T12MNI_affine.mat")
        os.system("fnirt --ref=../../MNI152_T1_3mm_brain.nii.gz --in=T1w_acpc_dc_restore_brain_diff_space.nii.gz --aff=T12MNI_affine.mat --cout=warps_T12MNI")
        os.system("invwarp --ref=T1w_acpc_dc_restore_brain_diff_space.nii.gz --warp=warps_T12MNI --out=warps_MNI2T1")
        os.system("warpinit ../../MNI152_T1_3mm_brain.nii.gz inv_identity_warp_no.nii.gz -force")
        os.system("applywarp --ref=T1w_acpc_dc_restore_brain_diff_space.nii.gz --in=inv_identity_warp_no.nii.gz --warp=warps_MNI2T1.nii.gz --out=mrtrix_warp_MNI2dwi.nii.gz")
        os.system("tcktransform tractography_prob_sift_3M_21p.tck mrtrix_warp_MNI2dwi.nii.gz tractography_prob_sift_3M_21p_MNI.tck -force")
        os.system("warpinvert mrtrix_warp_MNI2dwi.nii.gz mrtrix_warp_dwi2MNI.nii.gz -force")

        print("--Transformacion de formato de tractografia (.tck a .bundles)--")
        file_fchange("tractography_prob_sift_3M_21p_MNI.tck", "tractography_prob_sift_3M_21p_MNI.bundles", "tck2bundles",header="../../MNI152_T1_1mm.nii.gz")

        # Segmentación de fascículos (si se especifica)
        if  args.segmentacion is not None:
            print("--Inicio de proceso de segmentación de fascículos--")
            T1w = nb.as_closest_canonical(nb.load("T1w_acpc_dc_restore_brain_diff_space.nii.gz"))
            T1w_header = T1w.header
            fibs_tck = load_tractogram("tractography_prob_sift_3M_21p.tck", T1w_header, bbox_valid_check=False).streamlines
            tractography = np.asarray(fibs_tck)

            if "SWM" in args.segmentacion:
                print("Calculo de Fascículos de fibras cortas")
                os.makedirs("SWM", exist_ok=True)
                os.makedirs("SWM/results_folder_DWI", exist_ok=True)
                os.chdir("..")
                os.chdir("..")
                os.system("./main_index 21 "+ suj+"tractography_prob_sift_3M_21p_MNI.bundles subject AtlasRo/ AtlasRo_estables.txt "+ suj+"SWM/results_folder "+ suj+"SWM/indices_folder")
                SWM_fasc_index = os.listdir(suj + "SWM/indices_folder")
                for fasc_ind in SWM_fasc_index:
                    nombre_fasc = fasc_ind.rstrip(".txt")
                    indices = open(suj + "SWM/indices_folder/" + fasc_ind).readlines()
                    indices = [int(ind.rstrip("\n")[:-7]) for ind in indices]
                    bundle = tractography[indices]
                    centroids_tractogram_file = Tractogram(streamlines=bundle)
                    centroids_tractogram_file.affine_to_rasmm = np.eye(4)
                    centroids_tck = TF.TckFile(centroids_tractogram_file, header={'timestamp': 0})
                    centroids_tck.save(suj + "SWM/results_folder_DWI/" + nombre_fasc + ".tck")
                os.chdir(suj)

            if "DWM" in args.segmentacion:
                print("Calculo de Fascículos de fibras largas")
                os.makedirs("DWM", exist_ok=True)
                os.makedirs("DWM/results_folder_DWI", exist_ok=True)
                os.chdir("..")
                os.chdir("..")
                os.system("./main_index 21 "+ suj+"tractography_prob_sift_3M_21p_MNI.bundles subject atlas_faisceaux/ atlas_faisceaux.txt "+suj+"DWM/results_folder "+suj+"DWM/indices_folder")
                DWM_fasc_index = os.listdir(suj + "DWM/indices_folder")
                for fasc_ind in DWM_fasc_index:
                    nombre_fasc = fasc_ind.rstrip(".txt")
                    indices = open(suj + "DWM/indices_folder/" + fasc_ind).readlines()
                    indices = [int(ind.rstrip("\n")[:-7]) for ind in indices]
                    bundle = tractography[indices]
                    centroids_tractogram_file = Tractogram(streamlines=bundle)
                    centroids_tractogram_file.affine_to_rasmm = np.eye(4)
                    centroids_tck = TF.TckFile(centroids_tractogram_file, header={'timestamp': 0})
                    centroids_tck.save(suj + "DWM/results_folder_DWI/" + nombre_fasc[:-4] + ".tck")
                os.chdir(suj)

        os.chdir("..")
        os.chdir("..")


if __name__ == "__main__":
    main()
