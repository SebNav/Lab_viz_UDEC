# Pipeline DWI: preprocesamiento, tractografía y segmentación de fascículos

Pipeline del Laboratorio de Visualización de la Universidad de Concepción, creado por Sebastián Navarrete (sebastian.navarrete@biomedica.udec.cl) y basado en MRtrix3, FSL y ANTs.

## Qué hace

A partir de las imágenes de difusión y del T1w de cada sujeto, el pipeline:

1. Preprocesa la DWI: denoise, corrección de distorsión por susceptibilidad, eddy (movimiento y corrientes inducidas), corrección de bias y máscara cerebral.
2. Calcula los mapas DTI (FA, MD, AD, RD).
3. Estima las orientaciones de fibras (msmt-CSD) y genera una tractografía probabilística de 3 millones de fibras (iFOD2 + ACT + SIFT).
4. Lleva la tractografía al espacio MNI152 y segmenta los fascículos de fibras largas (DWM, 36) y cortas (SWM, 209) con atlas.
5. Extrae las métricas DTI de cada fascículo en espacio nativo (`metricas/bundle_metrics_<DWM|SWM>.csv`).
6. Genera un reporte de QC en PDF por sujeto.

Funciona con dos tipos de adquisición, que solo difieren en la corrección de la distorsión por susceptibilidad: **investigativa** (multishell con una b0 adquirida con la fase invertida, AP/PA) y **clínica** (single-shell o multishell, solo AP o solo PA). El tipo se detecta automáticamente para cada sujeto.

## Cambios respecto a la versión anterior

- La distorsión por susceptibilidad, que la versión 1 no corregía, ahora se corrige con topup: en adquisiciones de investigación se usa la b0 adquirida con la fase invertida, y en adquisiciones clínicas (solamente 1 b-values y menos dirreciones), que no la tienen, una b0 sintética sin distorsión generada a partir del T1w con [Synb0-DisCo](https://github.com/MASILab/Synb0-DISCO).
- La máscara cerebral puede obtenerse con dos métodos: [HD-BET](https://github.com/MIC-DKFZ/HD-BET), más confiable y con mejor desempeño, que se usa si está instalado (**se recomienda encarecidamente instalarlo**), o `dwi2mask` de MRtrix en caso contrario.
- El registro al espacio MNI152 pasó de FSL (`flirt` + `fnirt`) a ANTs SyN sobre la plantilla de 1 mm, que es más robusto y se ajusta mejor a la anatomía del sujeto.
- Se agregaron un reporte de QC para identificar rápidamente problemas en cada sujeto
- Se añadio el proceso de extracción de las métricas DTI por fascículo 
- Se añadio proceso opcional de limpieza que deja solo los archivos necesarios.
  
## Requisitos

| Software | Uso |
|---|---|
| [MRtrix3](https://www.mrtrix.org/download/) ≥ 3.0 | Procesamiento general |
| [FSL](https://fsl.fmrib.ox.ac.uk/fsl/docs/#/install/index) ≥ 6.0.6 | topup, eddy, registro DWI → T1w |
| [ANTs](https://github.com/ANTsX/ANTs/releases) | Bias field y registro a MNI |
| Python 3 | `nibabel`, `numpy`, `scipy`, `pandas`, `dipy`, `joblib`, `matplotlib` |
| [HD-BET](https://github.com/MIC-DKFZ/HD-BET) (recomendado) | Máscara cerebral |
| [Synb0-DisCo](https://github.com/MASILab/Synb0-DISCO) (datos clínicos) | b0 sintética; requiere Docker y una licencia de FreeSurfer |

## Instalación

El pipeline usa los atlas, `main_index` y las plantillas MNI de la carpeta `Algoritmos_y_Archivos/Pipeline_ProcesamientoDiff_TractySegmentacion/` del repositorio [Lab_viz_UDEC](https://github.com/SebNav/Lab_viz_UDEC).

- **Si ya tienes `Algoritmos_y_Archivos`:** descarga `Procesamiento_DWI.py` y `qc_dwi.py` y guárdalos juntos a los archivos de `Pipeline_ProcesamientoDiff_TractySegmentacion/` (donde está `Procesamiento_Mrtrix3.py`).
- **Si no lo tienes:**

```
git clone https://github.com/SebNav/Lab_viz_UDEC.git
cd Lab_viz_UDEC/Algoritmos_y_Archivos
git clone https://github.com/SebNav/Pipeline_DWI.git
```

La instalación de HD-BET y Synb0-DisCo se explica en [INSTALACION.md](INSTALACION.md).

## Datos de entrada

Una carpeta con un subdirectorio por sujeto, cada uno con el T1w procesada (`T1w_acpc_dc_restore_brain.nii.gz` y, opcionalmente, `T1w_acpc_dc_restore.nii.gz`) y la DWI en uno de estos formatos:

- **Salida de `dcm2niix`** (`.nii.gz/.json/.bval/.bvec`, en cualquier subcarpeta): los runs con b > 0 se concatenan y la serie de solo b0 con la fase opuesta se usa como b0 inversa. La dirección de fase se lee de los `.json`.
- **`diff.nii.gz/.bval/.bvec`** (formato de la versión 1): la dirección de fase sedebe indicar con al usar el script `-pe_dir` o se debe guardan el `acq.json` original con la variable (`{"PhaseEncodingDirection": "j- (este valor cambia según la acquisicion puede ser -j,j,-i,i-k,k"}`).

## Uso

Validar los datos (muestra, para cada sujeto, las series encontradas y la corrección de distorsión que se usará):

```
python3 Procesamiento_DWI.py -folder Carpeta_sujetos -check
```

Procesar:

```
python3 Procesamiento_DWI.py -folder Carpeta_sujetos -dti -segmentacion SWM DWM -cleanup
python3 Procesamiento_DWI.py -folder Carpeta_clinicos -pe_dir j- -dti -segmentacion SWM DWM
```

| Argumento | Descripción |
|---|---|
| `-folder` | Carpeta con los sujetos (obligatorio) |
| `-sujetos S1 S2` | Procesa solo esos sujetos |
| `-check` | Solo valida los datos |
| `-dti` | Mapas DTI y, con `-segmentacion`, métricas por fascículo |
| `-segmentacion SWM DWM` | Segmenta fascículos cortos, largos o ambos |
| `-sdc auto\|topup\|synb0\|none` | Corrección de distorsión; `auto` (por defecto) elige según los datos |
| `-pe_dir`, `-readout` | Dirección de fase y readout si no están en los `.json` |
| `-mask auto\|hdbet\|dwi2mask` | Método de la máscara; `auto` usa HD-BET si está instalado |
| `-cleanup` | Al terminar, borra los archivos intermedios |
| `-no_qc`, `-degibbs`, `-nthreads N` | Sin reporte de QC, con `mrdegibbs`, número de hilos |

Si un sujeto falla, el error queda en su log y se sigue con el siguiente. Al volver a ejecutar el mismo comando, el pipeline retoma desde el último paso terminado. El reporte de QC se puede regenerar con `python3 qc_dwi.py Carpeta_sujetos`.

## Salidas

Los resultados de cada sujeto quedan en `Sujeto/proc_dwi_v2/`:

| Archivo | Contenido |
|---|---|
| `QC/QC_<sujeto>.pdf` | Reporte de QC: corrección de distorsión, registros, tractografía y calidad de cada fascículo (concordancia con el atlas, simetría, cobertura del campo de visión), con métricas marcadas OK / REVISAR |
| `metricas/bundle_metrics_DWM.csv`, `metricas/bundle_metrics_SWM.csv` | FA, AD, RD y MD por fascículo: número de fibras, volumen, media, desviación, mediana e IQR |
| `DWM/`, `SWM/` | Fascículos segmentados en espacio MNI (`results_folder/`) y nativo (`results_folder_DWI/`) |
| `dwi_preproc_unbiased.nii.gz`, `dwi_mask.nii.gz`, `dti_*.nii.gz` | DWI preprocesada, máscara y mapas DTI |
| `tractography_prob_sift_3M.tck`, `tractography_prob_sift_3M_21p_MNI.tck` | Tractografía en espacio nativo y en MNI |
| `metodos.json`, `pipeline_v2.log` | Método de corrección usado, parámetros, versiones del software y registro de cada comando |
