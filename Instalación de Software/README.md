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

![Alt text](https://s3.us-east-2.amazonaws.com/brainder/2015/fsl-rpi/screenshot_debian_lxde_rpi2.png)

> [!CAUTION]
> Esta guia describe el proceso de instalación de FSL en Ubuntu, si se quiere instalar FSL en otro sistema operativo siga el proceso descrito en https://fsl.fmrib.ox.ac.uk/fsl/docs/#/install/linux.

I) Descargar el instalador de Python **fslinstaller.py**.

II) Abrir una terminal y ejecutar el script usando Python:

```console
python ~/Downloads/fslinstaller.py
```

Si el comando anterior no funciona, pruebe con python
3:
```console
python3 ~/Downloads/fslinstaller.py
```

III) Verificar la instalación:

1.- Escribe ```echo $FSLDIR``` en la terminal. Esto debería imprimir en pantalla la ubicación donde FSL fue instalado, por ejemplo: ```/home/labimagenes/fsl```.

2.- Abre el GUI de fsl escribiendo en el terminal ```fsl```.

3.- Abre el GUI de FSLeyes, escribiendo en el terminal ```fsleyes -std &```, esto deberia abrir FSL con una template MNI152 T1.

IV) Solución de problemas:
Si alguno de los pasos de verificación no funcionó, es posible que sea necesario añadir la ubicación de FSL en el archivo .bashrc. Para ello, sigue estos pasos:

1.- Abre una terminal y edita el archivo .bashrc escribiendo:

```console
open .bashrc
```

2.- Añade al final del archivo la siguiente línea, reemplazando con la ruta donde está instalado FSL (por ejemplo, /home/nombre_de_usuario/fsl):

``` console
export FSLDIR=/home/labimagenes/fsl
```

3.- Guarda los cambios y reinicia la terminal.

4.- Vuelve a intentar los pasos de verificación.

Si presentan otro tipo de problemas revisa el apartado de Troubleshooting oficial de FSL https://fsl.fmrib.ox.ac.uk/fsl/docs/#/install/troubleshooting


## Ants

## Freesurfer

![Alt text](https://www.researchgate.net/profile/Efren-Murillo-Zamora/publication/374489831/figure/fig1/AS:11431281196265744@1696597965800/MRI-volumetry-of-the-brain-with-FreeSurfer-software-A-a-3D-image-of-the-pial-area-The.png)

Para instalar Freesurfer se recomienda seguir la guia oficial https://drive.google.com/file/d/1uNwv29fCeuMHrmTyXw94ZSuroNsPOxu-/view?pli=1 desde la pagina 27 en apartado "To use the wget command from the terminal, cut and paste (or type in) the following commands into the terminal window:" en adelante (Las paginas anteriores describen el proceso de instalación de una maquina virtual)

## BrainVISA / Anatomist

![Alt text](https://brainvisa.info/web/_static/images/control_window0.png)

>[!NOTE]
> Este software funciona correctamente en Ubuntu y Windows. La instalación en Mac es más complicada, especialmente si se utiliza un computador con chip Mx.

El software se utiliza principalmente para la visualización de imágenes **.nii.gz** y datos de tractografía en formato **.bundles.** En la página oficial, existe un tutorial para la descarga e instalación de las nuevas versiones de Brainvisa/Anatomist. Sin embargo, estas versiones son difíciles de instalar o no funcionan con los archivos que utilizamos.

Para facilitar la instalación de Brainvisa/Anatomist, se recomienda descargar el ejecutable (.bin) y seguir el proceso de instalación de versiones anteriores (4.6.1). Aunque estas versiones fueron desarrolladas para Ubuntu 16.04, funcionan correctamente en versiones más recientes.

https://brainvisa.info/web/download-4.6.html


