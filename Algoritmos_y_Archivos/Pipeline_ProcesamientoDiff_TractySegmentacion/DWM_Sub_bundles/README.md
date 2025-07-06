# DWM atlas Sub-bundles

Para sub-dividir los fascículos guardan la carpeta "atlas_faisceaux_udd" y el archivo de texto atlas_faisceaux_udd.txt en la carpeta donde tienes tienes guardados todos los otros archivos en la misma carpeta donde esta el atlas normal. Guarda en esa misma carpeta el codigo "DWM_sub_bundles.py" y correlo el archivo en la terminal de la carpeta dandole el nombre de la carpeta donde estan todos tos sujetos procesados.



```
python3 Procesamiento_Mrtrix3.py -folder Nombre_Carpeta_sujetos
```

Si se desea que el pipeline ademas realice el proceso de segmentación y registro de fascículos de fibras cortas (SWM) y/o largas (DWM) al espacio MNI152, se debe ejecutar de la siguiente manera. Aquí se puede especificar uno (DWM o SWM) o ambos tipos de fascículos para segmentar:

```
python3 Procesamiento_Mrtrix3.py -folder Nombre_Carpeta_sujetos -segmentacion SWM DWM
```


**Por ahora el pipeline permite solamente el metodo CSD y el calculo probabilistico de tractografías**
