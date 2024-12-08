# Algoritmo de segmentación de fibras de fascículos para datos de tractografías en el espacio MNI152

## Descripción

Método de segmentación de fibras basado en la distancia euclidiana máxima de las fibras de un tractograma con respecto a las fibras del atlas de fascículos a utilizar para la segmentación. Este método fue propuesto originalmente por Labra [N. Labra et al., 2017] y posteriormente paralelizado por Vázquez [Vázquez et al., 2019].

En este repositorio se tienen los archivos necesarios para realizar la segmentacion de fascículos de fibras largas (DWM atlas [Guevara et al., 2012] y fasciculos de fibras cortas (SWM atlas [Godoy et al., 2021]). 


| Carpeta/Archivo         | Descripción                                                                                 |
|-------------------------|---------------------------------------------------------------------------------------------|
| AtlasRo/                | Carpeta con los fascículos de fibras cortas SWM (.bundles)                                   |
| AtlasRo_tck/            | Carpeta con los fascículos de fibras cortas SWM (.tck)                                       |
| atlas_faisceaux/        | Carpeta con los fascículos de fibras largas DWM (.bundles)                                   |
| atlas_faisceaux_tck/    | Carpeta con los fascículos de fibras largas DWM (.tck)                                       |
| AtlasRo.txt             | Archivo de texto que especifica la distancia Euclideana máxima utilizada para la segmentación de todos los fascículos de fibras cortas |
| AtlasRo_estables.txt    | Archivo de texto que especifica la distancia Euclideana máxima utilizada para la segmentación de los 208 fascículos de fibras cortas más estables |
| atlas_faisceaux.txt     | Archivo de texto que especifica la distancia Euclideana máxima utilizada para la segmentación de los fascículos de fibras largas |
| main_index              | Ejecutable de código escrito en C++ que genera la segmentación de datos                      |



## Modo de Uso

> [!IMPORTANT]  
> **Para que la segmentación funcione correctamente el tractograma a segmentar debe estar en el espacio de referencia MNI152 y las fibras remuestradas a 21 puntos equidistantes (en el formato .bundles).**

El siguiente enlace contiene archivos de ejemplo para probar el funcionamiento de la segmentación en su dispositivo:

**https://drive.google.com/drive/folders/1sfcVBNMhlxCPAVE_BoPI9hZH3tKV8LYW?usp=drive_link**

Aplicación:

1.- Abre una terminal en la carpeta donde se tiene el ejecutable.

2.- Correr el ejecutable especificando la carpeta de fascículos y archivo de distancias dependiendo de los fascículos a segmentar:

2.1. Codigo para Segmentar fascículos de fibras cortas (estables) utilizando la carpeta de ejemplo:

```console
./main_index 21 test_data/subject_x_21p_MNI.bundles subject  AtlasRo/ AtlasRo_estables.txt test_data/SWM/results_folder test_data/SWM/indices_folder
```

2.2. Codigo para Segmentar fascículos de fibras largas utilizando la carpeta de ejemplo:

```console
./main_index 21 test_data/subject_x_21p_MNI.bundles subject atlas_faisceaux/  atlas_faisceaux.txt test_data/DWM/results_folder test_data/DWM/indices_folder

```

## Resultados

- result_folder: Se almacenarán las fibras de los fascículos segmentados con el prefijo especificado (en este caso, "subject").
- indices_folder: Se generará un archivo de texto por cada fascículo segmentado, con el mismo nombre del fascículo. Este archivo contiene los índices de las fibras segmentadas y permite rastrearlas entre tractogramas, facilitando el cambio de espacio de referencia o de formato de manera rápida.


## Citation