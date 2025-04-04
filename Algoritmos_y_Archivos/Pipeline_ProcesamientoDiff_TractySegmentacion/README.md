# Pipeline Procesamiento de datos de difusión, calculo de tractogramas y segmentación de fascículos (SWM/DWM)

Pipeline de preprocesamiento y cálculo de tractografía creado por Sebastián Navarrete (sebastian.navarrete@biomedica.udec.cl), basado en el software MRtrix3. Este pipeline lleva a cabo el procesamiento de imágenes de difusión (Degibbs, Denoise, Masking, 5TT), cálculo de modelos de difusión (DTI o CSD), cálculo de tractogramas (Det o Prob), filtrado de tractogramas (SIFT), segmentación de fascículos de fibras cortas y largas (opcional), además del registro de imágenes. Para que este pipeline funcione correctamente, se necesita tener instalado en el mismo entorno el software MRtrix3, FSL y ANTs. Además, es necesario tener en la misma carpeta todos los archivos e imágenes necesarios para el proceso de registro, segmentación, etc.

Para que el pipeline funcione de forma correcta, se debe tener una carpeta con los datos por sujeto, separados por carpetas, como se muestra en el diagrama de abajo. Para cada sujeto, se necesitan las imágenes de difusión en formato .nii.gz, con los archivos bval y bvec correspondientes, además de una imagen estructural T1w con registro ACPC. Las imágenes de difusión deben tener el nombre "diff".

```
Folder  
   |  
   |--> Sujeto01  
   |	|  
   |	|--> diff.nii.gz  
   |	|--> diff.bval  
   |	|--> diff.bvec  
   |	L--> T1w_acpc_dc_restore_brain.nii.gz  
   |  
   |--> Sujeto02  
   |	|  
   |	|--> diff.nii.gz  
   |	|--> diff.bval  
   |	|--> diff.bvec  
   |	|--> T1w_acpc_dc_restore_brain.nii.gz  
   .  
   .  
   .  
   |--> SujetoXX  
   	|  
   	|--> diff.nii.gz  
   	|--> diff.bval  
   	|--> diff.bvec  
   	|--> T1w_acpc_dc_restore_brain.nii.gz  
```


## Modo de uso

El pipeline está escrito en Python y, por medio de la biblioteca OS, ejecuta diversos comandos de MRtrix3.

Para ejecutar el pipeline, se debe correr el siguiente comando en la terminal, desde la carpeta donde se encuentran todos los archivos necesarios y la carpeta con la información del sujet:

```
python3 Procesamiento_Mrtrix3.py -folder Nombre_Carpeta_sujetos
```

Si se desea que el pipeline ademas realice el proceso de segmentación y registro de fascículos de fibras cortas (SWM) y/o largas (DWM) al espacio MNI152, se debe ejecutar de la siguiente manera. Aquí se puede especificar uno (DWM o SWM) o ambos tipos de fascículos para segmentar:

```
python3 Procesamiento_Mrtrix3.py -folder Nombre_Carpeta_sujetos -segmentacion SWM DWM
```


**Por ahora el pipeline permite solamente el metodo CSD y el calculo probabilistico de tractografías**
