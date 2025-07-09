"""
Created on Thu Jan 30 13:44:52 2025

@author: natal
"""

from time import time
import numpy as np
import os
import read_write_bundle as bt
import nibabel as nb
from nibabel.streamlines.tractogram import Tractogram
from dipy.tracking.streamline import transform_streamlines
from dipy.io.streamline import load_tractogram, save_tractogram
import argparse
from distutils.util import strtobool
from tqdm import tqdm
from dipy.tracking.utils import length
from dipy.tracking.streamline import set_number_of_points  
import pandas as pd


def main():
    # Configurar el parser de argumentos
    parser = argparse.ArgumentParser(
        description="""Pipeline de tractometria creado por Sebastian Navarrete sebastian.navarrete@biomedica.udec.cl. 
         Este pipeline lleva a cabo la extracion de los valores promedios de las metricas de difusion (FA, ADC,MD y/o RD) de los fascículos previamente segmentados y posterior guardado en un excel."""
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
        choices=['SWM', 'DWM','DWM_sub'], 
        default=None, 
        help="Tipo de segmentación de fascículos analizar: SWM (fibras cortas), DWM (fibras largas), o ambas."
    )
    parser.add_argument(
        '-difusion_metric', 
        type=str, 
        choices=['FA','ADC','MD','RD'], 
        default='FA', 
        help="Metrica de difusión a utilizar para el calculo ."
    )
    parser.add_argument(
        '-get_image', 
        type=lambda x: bool(strtobool(x)),
        default='True', 
        help="Si las imagenes de FA, ADC, MD o RD no fueron previamente calculadas, se calcularan mediante el software Mrtrix3."
    )
    parser.add_argument(
        '-precise', 
        type=lambda x: bool(strtobool(x)),
        default='False', 
        help="""Método de remuestreo de fibras:
    - False (default): Todas las fibras usan el mismo número de puntos (basado en la fibra más larga / (tamaño_vóxel*0.75)).
      * Ventaja: Rápido.
      * Limitación: Fibras cortas pueden tener puntos redundantes en el mismo vóxel.
    - True: Cada fibra se remuestrea con puntos proporcionales a su longitud.
      * Ventaja: Precisión mejorada.
      * Limitación: Mayor tiempo de cálculo."""
    )
    
# =============================================================================
#     args = type('', (), {})()  # Create empty object
#     args.folder = "Controles"  # Replace with your folder path
#     args.segmentacion =  'DWM'  # Or just one of them
#     args.difusion_metric = 'FA'
#     args.get_image = False
#     args.precise = True
# =============================================================================
    
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
	
    if args.precise:
    	  precise_str = "precise"
    else:
   	  precise_str = "fast"
  
    
    #Se cargan los nombres de los fasciculos del cual se obtendra las metricas de difusión
    
    
    
   

    
    

    if args.get_image:
        print("Se calculan las imagenes FA, MD, RD y ADC de los sujetos")
        for suj in Sujetos:
            print(f"Procesando sujeto: {suj}")
            os.chdir(suj)
            os.system("dwi2tensor dwi_preproc_unbiased.mif tensor.mif -force -mask dwi_mask.mif")
            os.system("tensor2metric tensor.mif -fa FA.nii.gz -ad AD.nii.gz -rd RD.nii.gz -adc ADC.nii.gz -force")
            
            os.chdir("..")
            os.chdir("..")
    
    
    
    if args.difusion_metric is not None:
        print("Se calculan los valores por sujeto y por fasciculo- FA, MD, RD y ADC de los sujetos")
        
        for segmentation in args.segmentacion:
            if 'SWM' == segmentation:
                Fasc_cortos = open("AtlasRo_estables.txt","r").readlines()
                Fasc_cortos =  [f.split("\t")[0].rstrip("\n") for f in Fasc_cortos]
                bundles = Fasc_cortos
                b_path = "SWM/results_folder_DWI/"                 
                ant = ""
            if 'DWM' == segmentation:
                Fasc_largos = open("atlas_faisceaux.txt","r").readlines()
                Fasc_largos =  [f.split("\t")[0][6:-4].rstrip("\n") for f in Fasc_largos]
                bundles = Fasc_largos
                b_path = "SWM/results_folder_DWI/"
                ant = "atlas_"
                
            if 'DWM' == segmentation:
                Fasc_largos_sub = open("atlas__faisceaux_UDD.txt","r").readlines()
                Fasc_largos_sub =  [f.split("\t")[0].rstrip("\n") for f in Fasc_largos]
                bundles = Fasc_largos_sub
                b_path = "DWM_sub/results_folder_sub_bundles_DWI"   
                ant = "atlas_"
    
            all_data_FA  =  []
            all_data_RD  =  []
            all_data_ADC = []
            all_data_AD  =  []
            all_data_len =  []
            
            for suj in Sujetos:
            
                #tipo = suj[10:-1]
                print(suj)
                
                s_path = suj+"/"
                
                data_FA = []
                data_RD = []
                data_ADC = []
                data_AD = []
                data_len = []
                
                ornt    = nb.orientations.axcodes2ornt(('R', 'A', 'S')) #Orientacion de fibras .bundles
                
                FA        = nb.as_closest_canonical(nb.load(s_path+"FA.nii.gz"))
                FA_voxdim = FA.header["pixdim"][1:4]
                FA_affine = FA.affine
                FA_traslacion    = FA.affine
                FA_traslacion    = FA_traslacion[:3,3]
                FA  = nb.apply_orientation(FA.dataobj, ornt)  
                   
                RD = nb.as_closest_canonical(nb.load(s_path+"RD.nii.gz"))
                RD_voxdim = RD.header["pixdim"][1:4]
                RD_affine = RD.affine
                RD  = nb.apply_orientation(RD.dataobj, ornt)  
                   
                
                ADC = nb.as_closest_canonical(nb.load(s_path+"ADC.nii.gz"))
                #ADC_voxdim = ADC.header["pixdim"][1:4]
                ADC_affine = ADC.affine
                ADC_traslacion    = ADC_affine
                ADC_traslacion    = ADC_traslacion[:3,3]
                ADC  = nb.apply_orientation(ADC.dataobj, ornt)  
                
                AD = nb.as_closest_canonical(nb.load(s_path+"AD.nii.gz"))
                AD_voxdim = AD.header["pixdim"][1:4]
                AD_affine = AD.affine
                AD  = nb.apply_orientation(AD.dataobj, ornt)  
                
                T1w        = nb.as_closest_canonical(nb.load(s_path+"T1w_acpc_dc_restore_brain_diff_space.nii.gz"))
                T1w_traslacion    = T1w.affine
                T1w_traslacion    = T1w_traslacion[:3,3]
                T1w_header = T1w.header 
                
                rasmm_to_voxel_affine   = np.linalg.inv(FA_affine)
                
                for bun in tqdm(bundles[:]):
                    
                    bun = ant+bun
                    bun_path = suj+b_path
                    
                    if os.path.exists(bun_path+bun+".tck"):
                
                        
                        bundle = load_tractogram(bun_path+bun+".tck",T1w_header,bbox_valid_check=False).streamlines
                        
                        
                        largo = list(length(bundle)) 
                        
                        if args.precise:
                            
                            bundle_values_FA  = []
                            bundle_values_RD  = []
                            bundle_values_ADC = []
                            bundle_values_AD  = []
                            
                            
                            for fibra, l in zip(bundle,largo): # Recorre cada fibra y su largo
                                cant_pto = l/(FA_voxdim[0]*0.75)         # Se eligen la cantidad de puntos personalizada a remuestrar
                                fibra = set_number_of_points(fibra, int(cant_pto))
                                
                                streamlines_voxel   = transform_streamlines(fibra, rasmm_to_voxel_affine)
                                streamlines_voxel   = np.round(np.array(streamlines_voxel)).astype(int)
                                
                                
                                bundle_values_FA.extend(FA[streamlines_voxel[:,0],streamlines_voxel[:,1],streamlines_voxel[:,2]])
                                bundle_values_RD.extend(RD[streamlines_voxel[:,0],streamlines_voxel[:,1],streamlines_voxel[:,2]])
                                bundle_values_ADC.extend(ADC[streamlines_voxel[:,0],streamlines_voxel[:,1],streamlines_voxel[:,2]])
                                bundle_values_AD.extend(AD[streamlines_voxel[:,0],streamlines_voxel[:,1],streamlines_voxel[:,2]])
                                

                            
                            data_FA.append(np.mean(bundle_values_FA))
                            data_RD.append(np.mean(bundle_values_RD))
                            data_ADC.append(np.mean(bundle_values_ADC))
                            data_AD.append(np.mean(bundle_values_AD))
                            data_len.append(len(bun))
                        
                        
                        
                        else:
                            
                            cant_pto = max(largo)/(FA_voxdim[0]*0.75)
                            bundle = set_number_of_points(bundle,int(cant_pto))
                            streamlines_voxel   = transform_streamlines(bundle, rasmm_to_voxel_affine)
                            streamlines_voxel   = np.round(np.array(streamlines_voxel)).astype(int)
                            
    
                            
                            data_FA.append(np.mean(FA[streamlines_voxel[:,:,0],streamlines_voxel[:,:,1],streamlines_voxel[:,:,2]]))
                            data_RD.append(np.mean(RD[streamlines_voxel[:,:,0],streamlines_voxel[:,:,1],streamlines_voxel[:,:,2]]))
                            data_ADC.append(np.mean(ADC[streamlines_voxel[:,:,0],streamlines_voxel[:,:,1],streamlines_voxel[:,:,2]]))
                            data_AD.append(np.mean(AD[streamlines_voxel[:,:,0],streamlines_voxel[:,:,1],streamlines_voxel[:,:,2]]))
                            data_len.append(len(bun))

                    else:
                        print(bun,suj,"no existe")
                        data_FA.append(0)
                        data_RD.append(0)
                        data_ADC.append(0)
                        data_AD.append(0)
                        data_len.append(0)
                        
                all_data_FA.append(data_FA)
                all_data_RD.append(data_RD)
                all_data_ADC.append(data_ADC)
                all_data_AD.append(data_AD)
                all_data_len.append(data_len)

            index = [suj[10:-1] for suj in Sujetos]

            df_FA = pd.DataFrame(all_data_FA,columns=bundles[:],index = index[:])
            df_FA.to_excel("FA_bundles_"+b_path+"_Controles_"+precise_str+".xlsx")

            df_RD = pd.DataFrame(all_data_RD,columns=bundles[:],index = index[:])
            df_RD.to_excel("RD_bundles_"+b_path+"_Controles_"+precise_str+".xlsx")

            df_ADC = pd.DataFrame(all_data_ADC,columns=bundles[:],index = index[:])
            df_ADC.to_excel("ADC_bundles_"+b_path+"_Controles_"+precise_str+".xlsx")

            df_AD = pd.DataFrame(all_data_AD,columns=bundles[:],index = index[:])
            df_AD.to_excel("AD_bundles_"+b_path+"_Controles_"+precise_str+".xlsx")
            
            df_len = pd.DataFrame(all_data_len,columns=bundles[:],index = index[:])
            df_AD.to_excel("cantidad_fib_bundles_"+b_path+"_Controles_"+precise_str+".xlsx")
                
if __name__ == "__main__":
    ini   = time()
    main()
    fin   = time()
    
    print("Tiempo total: ",fin-ini, "segundos.")
