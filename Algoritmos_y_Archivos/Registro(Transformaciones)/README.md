# Registro de Imagenes y Tractogramas

![Alt text](https://3dqlab.stanford.edu/wp-content/uploads/2023/04/registered-final.png)

El registro en neuroimagen es el proceso de alinear diferentes conjuntos de datos, como imágenes cerebrales o tractogramas, en un mismo espacio de referencia para permitir su comparación y análisis conjunto.

Registro de imágenes: Consiste en alinear imágenes de distintas modalidades (RM estructural, funcional o dMRI) o de diferentes sujetos mediante transformaciones rígidas, afines o no lineales (warp).

Registro de tractogramas: Es el proceso mediante el cual se alinea la posición espacial de las fibras provenientes de un tractograma de un espacio a otro. Usualmente, las fibras se alinean a un espacio de referencia, como el MNI152, con el objetivo de realizar procedimientos específicos o llevar a cabo análisis detallados en dicho espacio.

La imagen que se somete a la transformación para ser alineada a otro espacio se conoce comúnmente como moving image o input image. Por otro lado, la imagen en el espacio deseado que utilizamos para calcular la transformada se denomina imagen de referencia, template o fixed image.

## Transformaciones afines 

![](https://github.com/SebNav/Lab_viz_UDEC/blob/main/Algoritmos_y_Archivos/Registro(Transformaciones)/Affine_transform.gif)

El registro mediante transformaciones afines consiste en una transformación lineal seguida de una traslación, que se aplica a cada vóxel de la imagen.
La transformación afín se representa mediante una matriz de 4x4, en la cual los 12 componentes corresponden a diferentes acciones: escalado, rotación, cizallamiento y traslación.
Existen distintos tipos de transformaciones para imágenes 3D según el número de grados de libertad (Degrees of Freedom, DOF) utilizados: 12 (afín), 9 (tradicional), 7 (escalado global) y 6 (cuerpo rígido).

![Alt text](https://community.mrtrix.org/uploads/default/original/2X/a/a589b481ff0e5c763d9740824f1787487c04276e.png)


## Transformaciones no-lineales (Warps)


Los warps pueden almacenar la información de las transformaciones de dos formas distintas:

**Deformation fields:** Se define como una imagen en la que cada vóxel define la posición correspondiente en la otra imagen (en coordenadas del espacio del escáner).

**Displacement fields:** Un campo de desplazamiento almacena los desplazamientos (en mm) a la otra imagen desde la posición de cada vóxel (en el espacio del escáner). 

Es relevante manejar esta información ya que distintos Software pueden preferir generar y/o usar un tipo especifico de warps, generando errores o resultados no esperados si se ultiliza uno distinto. Si es que se tiene un tipo de warp y se necesita el otro existen funciones como warpconvert de Mrtrix3 que permite modificar el tipo de warp.

## Registro de Tractogramas

## Ejemplos

> [!IMPORTANT]
> Los ejemplos descritos muestran el proceso de registro utilizando funciones de los Software FSL y Mrtrix3. Todos los archivos utilizados en los ejemplos se encuentran en la carpeta Test_data

### Registro de imagen estructural T1w con imagenes DWI

![Alt text](https://github.com/SebNav/Lab_viz_UDEC/blob/main/Algoritmos_y_Archivos/Registro(Transformaciones)/T1w_dwi_overlay2.png)

Para este caso, tenemos una imagen estructural T1w que necesitamos registrar al espacio de las imágenes de difusión. Este es un proceso común, ya que se requiere la imagen T1w en el espacio de difusión para realizar diversos procedimientos. Por ejemplo, se utiliza para mejorar la calidad de la tractografía aplicando ACT (Anatomically Constrained Tractography) mediante la imagen 5TT (5 Tissue Types), que se genera a partir de la T1w. Además, tener una imagen del cerebro de buena calidad en el espacio de difusión permite crear mejores máscaras, transformaciones, entre otros elementos esenciales.

1.- **Remover Craneo:** La imagen estructural T1w presenta el cerebro, creaneo y otros tejidos de la cabeza. para mejorar la calidad del registro es necesario remover todo lo que no sea cerebro, para ello existen muchos funciones de distintos programas aqui una lista de distintos comandos que existem:

- [Freesurfer mri_synthstrip](https://surfer.nmr.mgh.harvard.edu/docs/synthstrip/) (Personalmente uso este debido a la calidad de sus resultados pero ojo es un modelo DL puede generar malos resultados en la segmentación de infantes)
- [fsl bet](https://web.mit.edu/fsl_v5.0.10/fsl/doc/wiki/BET(2f)UserGuide.html)
- [antsBrainExtraction.sh](https://github.com/ANTsX/ANTs/blob/master/Scripts/antsBrainExtraction.sh)

```console
mri_synthstrip -i T1w_acpc_dc.nii.gz -o T1w_acpc_dc_brain.nii.gz
```

![Alt text](https://github.com/SebNav/Lab_viz_UDEC/blob/main/Algoritmos_y_Archivos/Registro(Transformaciones)/Brain_striping.png)


2. **Promediado de b-zero's (Opcional):** Debido a que usaremos como referencia una imagen de difusión, la cual presenta múltiples volúmenes, se recomienda extraer y promediar las imágenes b-zero para obtener un promedio. Esto se debe a que este tipo de imágenes está compuesto por múltiples volúmenes calculados en diferentes momentos del tiempo, lo que genera pequeñas diferencias en la posición de los volúmenes debido al movimiento del paciente durante la adquisición. Si no se aplica este paso, la función de FSL que utilizaremos para la transformación tomará únicamente el primer volumen como referencia para la transformación.

- Los siguientes comandos de Mrtrix3 extraen y promedios estos volumenes, para este paso se necesita los valores y las direcciones de los b-values.

```console
dwiextract dwi_preproc_unbiased.mif - -bzero -force | mrmath - mean bzero.nii.gz -axis 3 -force"
```


3. **Calculo de Transformación Afín:** Se utiliza el software FSL para calcular la transformada y MRtrix3 para aplicarla. Para garantizar que el resultado de la transformada sea de buena calidad, se recomienda utilizar como referencia la imagen con mayor resolución. Esto no supone un problema si es necesario dejar la imagen que se desea mover como referencia, ya que en ese caso simplemente se aplica la transformada inversa.

En este caso, como contamos con una imagen T1w (alta resolución) y una imagen DWI (menor resolución), utilizaremos la T1w como referencia. Primero, utilizaremos el comando flirt de FSL para calcular la matriz afín. Luego, emplearemos el comando transformconvert de MRtrix3 para convertir el formato de la matriz a uno compatible con MRtrix3. Finalmente, aplicaremos la transformada inversa de la matriz afín calculada para obtener nuestra imagen T1w en el espacio de la DWI.


```console
flirt -in mean_bzero.nii.gz -ref T1w_acpc_dc_brain.nii.gz -omat struct2dwi.mat -dof 12 
transformconvert struct2dwi.mat mean_bzero.nii.gz T1w_acpc_dc_brain.nii.gz flirt_import struct2dwi_mrtrix.txt -force
mrtransform T1w_acpc_dc_brain.nii.gz -linear struct2dwi_mrtrix.txt -inverse T1w_acpc_dc_brain_dwi_space.nii.gz -force
```

- ref: Imagen en el espacio de referencia
- in: Imagen de entrada
- omat: Nombre del archivo con la transfomación afín
- dof: Cantidad de Degrees of Freedom que se utilizaran para calcular la transformación
- out: Imagen transformada al espacio deseado.

Como podemos observar en la figura, la imagen T1w (en gris) está correctamente alineada con la imagen de difusión (en rojo). Las manchas rojas que se observan debajo del cerebro corresponden a los ojos y otros tejidos del sujeto. Si estos tejidos causan problemas en el registro, se podrían eliminar antes de aplicar la transformada, como se realizó con la imagen T1w.

Si la imagen de difusión utilizada no fue corregida por distorsiones EPI (EPI distortion correction) durante el preprocesamiento, es normal que la parte anterior del cerebro no quede perfectamente alineada. Esto ocurre porque, al calcular la imagen DWI, se utiliza un barrido de difusión en una dirección, lo que genera artefactos de achatamiento o estiramiento en la parte anterior del cerebro. Este efecto depende de la codificación de fase (phase encoding) utilizada para obtener la DWI, que usualmente es A-P (anterior-posterior) o P-A (posterior-anterior).

![Alt text]()

### Transformación Lineal de Imagenes


### Transformación No-Linear de Imagenes


### Transformación de Tractogramas
