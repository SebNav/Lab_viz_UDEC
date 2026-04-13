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



## Segmentacion de fascículos usando Docker Image

## Requisitos

- [Docker](https://docs.docker.com/get-docker/) instalado y funcionando en su equipo.
- Archivo de tractografía en **espacio MNI152**, remuestreado a **21 puntos equidistantes**, en formato `.tck` o `.bundles`.

---

## 1. Descarga de la imagen

Descargue el archivo de imagen Docker (`.tar`) desde el siguiente enlace:

> **[Google Drive — Imagen Docker](https://drive.google.com/file/d/1F-_DdahEZs4TJQ3TMFwCSnHaWzPsjC3G/view?usp=sharing)**


---

## 2. Carga de la imagen

Abra una terminal en la carpeta donde descargó el archivo `.tar` y ejecute:

```bash
docker load -i tract_segment.tar.gz
```

Verifique que la imagen se cargó correctamente:

```bash
docker images
```

Debería aparecer una entrada llamada **`tract_segment`** en la lista.

---

## 3. Ejecutar la segmentación

La imagen requiere que sus datos sean **montados** en el contenedor en tiempo de ejecución mediante la opción `-v` — los archivos de tractografía no están incluidos dentro de la imagen.

### Sintaxis general

```bash
docker run --rm \
  -v /ruta/a/sus/datos:/data \
  tract_segment \
  -i /data/<archivo_tractografia> \
  -o /data/<carpeta_salida> \
  -s <tipo_segmentacion> \
  [-o_f <formato_salida>]
```

| Argumento | Requerido | Descripción |
|-----------|-----------|-------------|
| `-i` / `--input` | Sí | Ruta al archivo de tractografía de entrada (`.tck` o `.bundles`) **dentro del contenedor** (es decir, bajo `/data/`) |
| `-o` / `--output` | Sí | Ruta a la carpeta de salida **dentro del contenedor** donde se guardarán los resultados |
| `-s` / `--segmentation` | Sí | Tipo de segmentación: `DWM`, `SWM` o `SWM_e` |
| `-o_f` / `--output_format` | No | Formato de salida de los fascículos segmentados: `tck` o `bundles` (por defecto: `bundles`) |

> **Nota:** Reemplace `/ruta/a/sus/datos` con la **ruta absoluta** en su máquina local donde se encuentra el archivo de tractografía. Los resultados se guardarán en esa misma carpeta.

---

### Ejemplos

#### Segmentación de materia blanca profunda (DWM) — salida en `.bundles`

```bash
docker run --rm \
  -v /home/usuario/mis_datos:/data \
  tract_segment:latest \
  -i /data/sujeto_21p_MNI.bundles \
  -o /data/salida \
  -s DWM
```

#### Segmentación de materia blanca superficial — todos los fascículos (SWM) — salida en `.tck`

```bash
docker run --rm \
  -v /home/usuario/mis_datos:/data \
  tract_segment:latest \
  -i /data/sujeto_21p_MNI.tck \
  -o /data/salida \
  -s SWM \
  -o_f tck
```

#### Segmentación de materia blanca superficial — solo fascículos estables (SWM_e) — salida en `.bundles`

```bash
docker run --rm \
  -v /home/usuario/mis_datos:/data \
  tract_segment:latest \
  -i /data/sujeto_21p_MNI.tck \
  -o /data/salida \
  -s SWM_e
```

---

## 4. Estructura de la salida

Al finalizar la segmentación, la carpeta de salida contendrá:

```
salida/
├── results_folder_<tipo_segmentacion>/    # Fascículos segmentados (.bundles o .tck)
└── indices_folder_<tipo_segmentacion>/    # Un .txt por fascículo con los índices de fibras
```

Si se utilizó `-o_f tck`, la carpeta de resultados en `.bundles` se elimina automáticamente y es reemplazada por:

```
salida/
├── results_folder_<tipo_segmentacion>_tck/   # Fascículos segmentados (.tck)
└── indices_folder_<tipo_segmentacion>/        # Índices de fibras (.txt)
```

- **results_folder**: contiene un archivo por cada fascículo segmentado, nombrado según el fascículo del atlas.  
- **indices_folder**: contiene un `.txt` por fascículo con los índices de las fibras coincidentes en el tractograma original, útil para rastrear fibras entre formatos o espacios de referencia.

---

## Tipos de segmentación

| Tipo | Atlas | Fascículos |
|------|-------|------------|
| `DWM` | Atlas DWM ([Guevara et al., 2012](10.1016/j.neuroimage.2012.02.071)) | Fascículos de materia blanca profunda / largo alcance |
| `SWM` | Atlas SWM ([Roman et al., 2022](https://doi.org/10.1016/j.neuroimage.2022.119550)) | Todos los fascículos de materia blanca superficial |
| `SWM_e` | Atlas SWM (subconjunto estable) | Los 209 fascículos superficiales más reproducibles |

---



