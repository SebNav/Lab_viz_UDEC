# Registro de Imagenes y Tractogramas

![Alt text](https://3dqlab.stanford.edu/wp-content/uploads/2023/04/registered-final.png)

El registro en neuroimagen es el proceso de alinear diferentes conjuntos de datos, como imágenes cerebrales o tractogramas, en un mismo espacio de referencia para permitir su comparación y análisis conjunto.

Registro de imágenes: Consiste en alinear imágenes de distintas modalidades (RM estructural, funcional o dMRI) o sujetos mediante transformaciones rígidas, afines o no lineales (warp).

Registro de tractogramas: Proceso en el cual se alinea la posición espacial de las fibras provenientes de un tractograma de un espacio a otro. Usualmente, las fibras se alinean a un espacio de referencia, como el MNI152, con el objetivo de realizar algún procedimiento en este espacio o llevar a cabo algún análisis.


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
> Los ejemplos descritos muestran el proceso de registro utilizando funciones de los Software FSL y Mrtrix3.

### Registro de imagen estructural T1w con imagenes DWI

![Alt text](https://github.com/SebNav/Lab_viz_UDEC/blob/main/Algoritmos_y_Archivos/Registro(Transformaciones)/T1w_dwi_overlay2.png)

Para este caso, tenemos una imagen estructural T1w que necesitamos registrar al espacio de las imágenes de difusión. Este es un proceso común, ya que se requiere la imagen T1w en el espacio de difusión para realizar diversos procedimientos. Por ejemplo, se utiliza para mejorar la calidad de la tractografía aplicando ACT (Anatomically Constrained Tractography) mediante la imagen 5TT (5 Tissue Types), que se genera a partir de la T1w. Además, tener una imagen del cerebro de buena calidad en el espacio de difusión permite crear mejores máscaras, transformaciones, entre otros elementos esenciales.

1.-


### Transformación Lineal de Imagenes


### Transformación No-Linear de Imagenes


### Transformación de Tractogramas
