# Registro de Imagenes y Tractogramas

![Alt text](https://3dqlab.stanford.edu/wp-content/uploads/2023/04/registered-final.png)

El registro en neuroimagen es el proceso de alinear diferentes conjuntos de datos, como imágenes cerebrales o tractogramas, en un mismo espacio de referencia para permitir su comparación y análisis conjunto.

Registro de imágenes: Consiste en alinear imágenes de distintas modalidades (RM estructural, funcional o dMRI) o sujetos mediante transformaciones rígidas, afines o no lineales (warp).

Registro de tractogramas: Proceso en el cual se alinea la posición espacial de las fibras provenientes de un tractograma de un espacio a otro. Usualmente, las fibras se alinean a un espacio de referencia, como el MNI152, con el objetivo de realizar algún procedimiento en este espacio o llevar a cabo algún análisis.


## Transformaciones afines 

![](https://github.com/SebNav/Lab_viz_UDEC/blob/main/Algoritmos_y_Archivos/Registro(Transformaciones)/Affine_transform.gif)


## Transformaciones no-lineales (Warps)

Los warps pueden almacenar la información de las transformaciones de dos formas distintas:

**Deformation fields:** Se define como una imagen en la que cada vóxel define la posición correspondiente en la otra imagen (en coordenadas del espacio del escáner).

**Displacement fields:** Un campo de desplazamiento almacena los desplazamientos (en mm) a la otra imagen desde la posición de cada vóxel (en el espacio del escáner). 

Es relevante manejar esta información ya que distintos Software pueden preferir generar y/o usar un tipo especifico de warps, generando errores o resultados no esperados si se ultiliza uno distinto. Si es que se tiene un tipo de warp y se necesita el otro existen funciones como warpconvert de Mrtrix3 que permite modificar el tipo de warp.

## Registro de Tractogramas

## Ejemplos

> [!IMPORTANT]
> Los ejemplos descritos muestran el proceso de registro utilizando funciones de los Software FSL y Mrtrix3.

> [!CAUTION]
> Si la imagen de referencia no presenta cráneo y la imagen que se va a mover sí lo tiene, se recomienda eliminar el cráneo de la imagen que se va a mover para mejorar la calidad del registro. Si la imagen de referencia presenta cráneo, se recomienda que la imagen que se va a mover también lo tenga.




### Transformación Lineal de Imagenes

### Transformación No-Linear de Imagenes


### Transformación de Tractogramas