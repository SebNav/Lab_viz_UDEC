# Instalación de HD-BET y Synb0-DisCo

Ambos son opcionales, pero mejoran el resultado: HD-BET da una máscara cerebral más precisa que `dwi2mask`, y Synb0-DisCo permite corregir la distorsión por susceptibilidad en datos clínicos. Para comprobar si el pipeline los encuentra, se usa `python3 Procesamiento_DWI.py -folder Carpeta_sujetos -check`.

## HD-BET (recomendado)

[HD-BET](https://github.com/MIC-DKFZ/HD-BET) es una red neuronal para extraer el cerebro. Se instala con `pip` en el mismo entorno de Python del pipeline:

```
pip install hd-bet
hd-bet -h
```

Los pesos de la red se descargan solos la primera vez que se usa. Con una GPU NVIDIA tarda ~40 s por sujeto; sin GPU funciona igual, pero más lento. El pipeline lo usa automáticamente si encuentra el comando `hd-bet`.

## Synb0-DisCo (datos clínicos)

[Synb0-DisCo](https://github.com/MASILab/Synb0-DISCO) genera, a partir del T1w, una b0 sin distorsión que reemplaza a la b0 con fase invertida que no tienen las adquisiciones clínicas. Corre dentro de un contenedor y necesita tres cosas:

**1. Docker** ([instalación en Ubuntu](https://docs.docker.com/engine/install/ubuntu/)), con permiso para usarlo sin `sudo`:

```
sudo usermod -aG docker $USER    # luego cerrar sesión y volver a entrar
docker run hello-world
```

**2. La imagen de Synb0-DisCo** (~12 GB):

```
docker pull leonyichencai/synb0-disco:v3.1
```

**3. Una licencia de FreeSurfer**, gratuita, que se obtiene en el [registro de FreeSurfer](https://surfer.nmr.mgh.harvard.edu/registration.html). El pipeline la busca en este orden: `-fs_license`, la variable `FS_LICENSE`, `$FREESURFER_HOME/license.txt` y `~/license.txt`. Lo más simple es guardarla como `~/license.txt` o agregar a `~/.bashrc`:

```
export FS_LICENSE=/ruta/a/license.txt
```

Synb0-DisCo no requiere GPU y tarda ~13 min por sujeto. Si falta alguno de los tres requisitos, los sujetos clínicos se procesan sin corrección de susceptibilidad (solo eddy), y `-check` y el reporte de QC lo indican.

**En un cluster sin Docker**, se puede usar [Apptainer](https://apptainer.org/docs/user/latest/):

```
apptainer pull synb0.sif docker://leonyichencai/synb0-disco:v3.1
python3 Procesamiento_DWI.py -folder Carpeta_clinicos -container apptainer -synb0_sif synb0.sif ...
```

Si se usa Synb0-DisCo, hay que citar: Schilling et al., *Magn Reson Imaging* 2019;64:62-70 y *PLoS One* 2020;15(7):e0236418.
