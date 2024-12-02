## MRtrix3 

![Alt text](https://www.mrtrix.org/images/frontpage/mrview.jpg)


> [!CAUTION]
> Si bien MRtrix3 puede instalarse en computadoras con Windows o Mac, en estos sistemas operativos diversas funciones del software no funcionan correctamente o son muy difíciles de hacer funcionar. Se recomienda utilizar Ubuntu 22.04 para el uso de MRtrix3.

> [!WARNING]  
> Algunos de los scripts de MRtrix3 dependen de funciones de otros software, por lo que se recomienda instalar FSL y ANTS para el correcto funcionamiento de MRtrix3. Además, se recomienda utilizar un computador con Python 3.11.8 (o una versión inferior) debido a un error entre las últimas versiones de MRtrix3 y Python 3.12.

MRtrix3 es un software que permite procesar imágenes de difusión para realizar tractografía, analizar la materia blanca, medir anisotropía, registrar imágenes y crear mapas de conectividad cerebral.

La forma más fácil de instalar MRtrix3 es mediante Anaconda/Miniconda. Con Conda instalado, simplemente abre un terminal y ejecuta el siguiente código:

``` console
conda install -c mrtrix3 mrtrix3
```

Una vez instalado, para verificar que se haya instalado correctamente, se puede abrir el visualizador de MRtrix3 con el siguiente comando:

``` console
Mrview
```

## FSL

## Ants

## Freesurfer


## BrainVISA / Anatomist

![Alt text](https://brainvisa.info/web/_static/images/control_window0.png)

>[!NOTE]
> Este software funciona correctamente en Ubuntu y Windows. La instalación en Mac es más complicada, especialmente si se utiliza un computador con chip Mx.

El software se utiliza principalmente para la visualización de imágenes **.nii.gz** y datos de tractografía en formato **.bundles.** En la página oficial, existe un tutorial para la descarga e instalación de las nuevas versiones de Brainvisa/Anatomist. Sin embargo, estas versiones son difíciles de instalar o no funcionan con los archivos que utilizamos.

Para facilitar la instalación de Brainvisa/Anatomist, se recomienda descargar el ejecutable (.bin) y seguir el proceso de instalación de versiones anteriores (4.6.1). Aunque estas versiones fueron desarrolladas para Ubuntu 16.04, funcionan correctamente en versiones más recientes.

https://brainvisa.info/web/download-4.6.html


