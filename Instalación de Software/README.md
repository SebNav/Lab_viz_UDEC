## MRtrix3 

![Alt text](https://www.mrtrix.org/images/frontpage/mrview.jpg)


> [!CAUTION]
> Si bien MRtrix3 puede instalarse en computadoras con Windows o Mac, en estos sistemas operativos diversas funciones del software no funcionan correctamente o son muy difíciles de hacer funcionar. Se recomienda utilizar Ubuntu 22.04 para el uso de MRtrix3.

> [!WARNING]  
> Algunos de los scripts de MRtrix3 dependen de funciones de otros software, por lo que se recomienda instalar FSL y ANTS para el correcto funcionamiento de MRtrix3. Además, se recomienda utilizar un computador con Python 3.11.8 (o una versión inferior) debido a un error entre las últimas versiones de MRtrix3 y Python 3.12. A la fecha 10/06/2025 sigue habiando problemas con python 3.12.X y Mrtirx3 al usar stripts de fsl como '5ttgen fsl'(T_T).

MRtrix3 es un software que permite procesar imágenes de difusión para realizar tractografía, analizar la materia blanca, medir anisotropía, registrar imágenes y crear mapas de conectividad cerebral.

La forma más fácil de instalar MRtrix3 es mediante Anaconda/Miniconda. Con Conda instalado, simplemente abre un terminal y ejecuta el siguiente código:

``` console
conda install -c mrtrix3 mrtrix3
```

Una vez instalado, para verificar que se haya instalado correctamente, se puede abrir el visualizador de MRtrix3 con el siguiente comando:

``` console
mrview
```

## FSL

![Alt text](https://s3.us-east-2.amazonaws.com/brainder/2015/fsl-rpi/screenshot_debian_lxde_rpi2.png)

> [!CAUTION]
> Esta guía describe el proceso de instalación de FSL en Ubuntu. Para otros sistemas operativos, sigue las instrucciones oficiales en https://fsl.fmrib.ox.ac.uk/fsl/docs/install/linux.html.

### Instalación

I) Abre una terminal y ejecuta el siguiente comando. Este descargará e instalará FSL automáticamente (el proceso tarda entre 10 y 15 minutos según la velocidad de conexión):

```console
curl -Ls https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/releases/getfsl.sh | sh -s
```

Al finalizar, deberías ver el mensaje: **`FSL successfully installed`**

II) Cierra y vuelve a abrir la terminal para que los cambios tengan efecto.

### Verificación

1.- Escribe `echo $FSLDIR` en la terminal. Debería imprimir la ruta donde se instaló FSL, por ejemplo: `/home/labimagenes/fsl`.

2.- Ejecuta `fslmaths` en la terminal. Debería mostrar el texto de ayuda del comando.

3.- Abre el GUI de FSL escribiendo en la terminal `fsl &`.

4.- Abre el visualizador FSLeyes escribiendo `fsleyes -std &`. Esto debería abrir FSLeyes con la plantilla MNI152 T1.

### Solución de problemas

Si la instalación no terminó con el mensaje de éxito, guarda todo el texto de la terminal y busca en tu carpeta de inicio un archivo de log llamado `fsl_installation_<fecha>.log`.

Si alguno de los pasos de verificación falla, es posible que sea necesario añadir manualmente la ubicación de FSL al archivo `.bashrc`. Para ello:

1.- Abre `.bashrc` con un editor de texto:

```console
nano ~/.bashrc
```

2.- Añade al final del archivo las siguientes líneas, reemplazando la ruta con la ubicación real de tu instalación (por ejemplo, `/home/nombre_de_usuario/fsl`):

```console
FSLDIR=/home/nombre_de_usuario/fsl
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
```

3.- Guarda los cambios y reinicia la terminal.

4.- Vuelve a intentar los pasos de verificación.

Para otros problemas, consulta la página oficial de Troubleshooting de FSL: https://fsl.fmrib.ox.ac.uk/fsl/docs/install/troubleshooting.html


## ANTs / ANTsX

[ANTs/ANTsX](https://github.com/ANTsX/ANTs?tab=readme-ov-file) (Advanced Normalization Tools ecosystem) es un conjunto de herramientas de procesamiento de imágenes biomédicas diseñado para el registro (alineación), segmentación y análisis estadístico de imágenes cerebrales.

Existen diversos métodos para instalar ANTs, incluyendo la compilación desde el código fuente (proceso que puede tardar más de una hora y requiere configurar múltiples dependencias de C++) o el uso de contenedores Docker/Singularity.

Sin embargo, para el setup que utilizamos en el laboratorio, la forma más rápida y sencilla es mediante Conda.

1.- Abre una terminal y ejecuta el siguiente comando:

``` console
conda install -c conda-forge ants
```

2.- Verificación de la instalación:

``` console
antsRegistration --version
```


## Freesurfer

![Alt text](https://andysbrainbook.readthedocs.io/en/stable/_images/06_Freeview_Example.png)

Para instalar FreeSurfer, se recomienda seguir la guía oficial disponible en el siguiente enlace:  
https://drive.google.com/file/d/1uNwv29fCeuMHrmTyXw94ZSuroNsPOxu-/view?pli=1

Comience desde la página 27, en el apartado:  
"To use the wget command from the terminal, cut and paste (or type in) the following commands into the terminal window:".

(Las páginas anteriores describen el proceso de instalación en una máquina virtual.)

## BrainVISA / Anatomist

![Alt text](https://brainvisa.info/web/_static/images/control_window0.png)

> [!NOTE]
> Linux es el único sistema operativo soportado de forma nativa. Los usuarios de Windows deben habilitar WSL2, y los usuarios de Mac deben instalar una máquina virtual de Linux.

El software se utiliza principalmente para la visualización de imágenes **.nii.gz** y datos de tractografía en formato **.bundles.**

Anteriormente, las versiones nuevas de BrainVISA eran difíciles de instalar o no funcionaban correctamente, por lo que se recomendaba instalar versiones antiguas como la 4.6.1 (disponibles en https://brainvisa.info/web/download-4.6.html). Sin embargo, a partir de la versión 6.0, el proceso de instalación es sencillo y funciona correctamente en Ubuntu 22.04 sin necesidad de máquina virtual.

La versión actual se distribuye mediante el gestor de paquetes **Pixi** a través del canal `neuro-forge`. Para instalarla, sigue los siguientes pasos:

### Instalación

I) Instala el gestor de paquetes **Pixi**:

```console
curl -fsSL https://pixi.sh/install.sh | sh
```

Cierra y vuelve a abrir la terminal para que los cambios tengan efecto.

II) Crea un directorio de trabajo e inicializa el entorno:

```console
mkdir ~/brainvisa
cd ~/brainvisa
pixi init -c https://brainvisa.info/neuro-forge -c conda-forge
```

III) Instala BrainVISA dentro del entorno:

```console
pixi add brainvisa
```

IV) Entra al entorno e inicia la aplicación:

```console
pixi shell
brainvisa
```

V) (Recomendado) Ejecuta el setup inicial para actualizar la base de datos compartida de BrainVISA:

```console
pixi run brainvisa -b --setup
```

Para más información o versiones anteriores, visita la página oficial de descarga: https://brainvisa.info/web/download.html


