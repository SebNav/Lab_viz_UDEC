# Conversión de Formato de Fibras (TCK, TRK, Bundles)

Este repositorio contiene un script en Python para realizar conversiones entre distintos formatos de archivos de fibras. Los formatos soportados incluyen:

- `tck2trk`  
- `trk2tck`  
- `tck2bundles`  
- `bundles2tck`  

> **Nota:** Para convertir archivos del formato `bundles` a `trk`, primero se debe convertir a `tck`.

## Funciones principales

El script ofrece dos funciones principales:

- `folder_fchange`: pensada para cambiar el formato de múltiples fascículos, por ejemplo, los resultantes de la segmentación de una tractografía.
- `file_fchange`: útil para convertir una única tractografía de un formato a otro.

## Parámetros de entrada

Ambas funciones requieren los siguientes parámetros:

- `bundle_path`: ruta al archivo (en `file_fchange`) o a la carpeta que contiene los archivos a convertir (en `folder_fchange`).
- `bundle_output_path`: ruta (y nombre) del archivo de salida o carpeta donde se guardarán los archivos convertidos. Si la carpeta no existe, se debe activar la opción `crear_carpeta=True`.
- `format_change`: tipo de conversión a realizar, por ejemplo `"tck2trk"`.
- `nthreads`: número de núcleos de CPU que se usarán para la conversión.
- `crear_carpeta`: (solo para `folder_fchange`) define si se debe crear la carpeta de salida en caso de que no exista.

## Ejemplo de uso

```python
from script_name import folder_fchange, file_fchange

# Para convertir una carpeta completa de archivos .tck a .trk
folder_fchange(
    bundle_path="ruta/a/carpeta_tck",
    bundle_output_path="ruta/a/carpeta_trk",
    format_change="tck2trk",
    nthreads=4,
    crear_carpeta=True
)

# Para convertir un solo archivo .trk a .tck
file_fchange(
    bundle_path="ruta/al/archivo.trk",
    bundle_output_path="ruta/salida/archivo.tck",
    format_change="trk2tck",
    nthreads=2
)

