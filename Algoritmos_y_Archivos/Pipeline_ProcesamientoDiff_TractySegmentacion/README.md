# Pipeline Procesamiento de datos de difusión, calculo de tractogramas y segmentación de fascículos (SWM/DWM)

Pipeline de preprocesamiento y calculo de tractografía creado por Sebastian Navarrete sebastian.navarrete@biomedica.udec.cl basado en el software Mrtrix3.  Este pipeline lleva a cabo el procesamiento de imagenes de difusión (Degibbs, DeNoise,Masking,5TT), calculo de modelo de difusion (DTI o CSD), calculo de Tractograma (Det o Prob), filtrado de tractograma (SIFT), segmentación de fasciculos de fibras cortas y largas (Opcional), ademas del registro de imagenes. Para que este pipeline funcione correctamente se necesitan tener instalada en el mismo ambiente el software Mrtrix3, fsl y ants. Ademas de tener en la misma carpeta todas los archivos e imagenes necesarios para el proceso de registro, sementacion , etc.


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
