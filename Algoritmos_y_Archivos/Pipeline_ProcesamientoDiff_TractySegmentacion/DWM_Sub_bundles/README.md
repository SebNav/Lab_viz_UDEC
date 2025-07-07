# DWM Atlas Sub-Bundles

Este Pipeline permite subdividir los fascículos del atlas DWM en sub-bundles.

## Instrucciones

1. Guarda la carpeta `atlas_faisceaux_udd` y el archivo de texto `atlas_faisceaux_udd.txt` en el mismo directorio donde se encuentra el atlas original y todos los demás archivos necesarios.

2. Coloca también el archivo `DWM_sub_bundles.py` en ese mismo directorio.

3. Abre una terminal en esa carpeta y ejecuta el siguiente comando, reemplazando `Nombre_Carpeta_sujetos` por el nombre de la carpeta que contiene todos los sujetos procesados:


```
python3 DWM_sub_bundles.py -folder Nombre_Carpeta_sujetos
```

## Salida

El script creará una nueva carpeta en la ruta `DWM` de cada sujeto. Esta carpeta contendrá las subdivisiones del atlas en:

- **Espacio MNI**: `result_folders_sub_bundles`
- **Espacio de difusión**: `result_folders_sub_bundles_DWI`



