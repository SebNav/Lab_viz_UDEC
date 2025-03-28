# Pipeline de obtención de valores promedios de metricas de difusion (FA,ADC,RD,AD) (SWM/DWM)

Pipeline automatico para la obtencion de los valores de promedio para las metricas de difusión FA, ADC, RD, AD para los fascículos de fibras largas(DWM) y/o fibras cortas(SWM). El algoritmo como salida un archivo excel (para cada metrica un excel distinto), donde se tiene los sujetos (de la carpeta que se pase para el procesamiento) con los valores promedio para cada fascículos.

```
Folder  
   |  
   |--> Sujeto01  
   |	|  
   |	|--> diff.nii.gz  
   |	|--> diff.bval  
   |	|--> diff.bvec  
   |	L--> T1w_acpc_dc_restore_brain.nii.gz  
   |  
   |--> Sujeto02  
   |	|  
   |	|--> diff.nii.gz  
   |	|--> diff.bval  
   |	|--> diff.bvec  
   |	|--> T1w_acpc_dc_restore_brain.nii.gz  
   .  
   .  
   .  
   |--> SujetoXX  
   	|  
   	|--> diff.nii.gz  
   	|--> diff.bval  
   	|--> diff.bvec  
   	|--> T1w_acpc_dc_restore_brain.nii.gz  
```




## Uso del programa

Ejecuta el script desde la terminal con:

```bash
python3 Extraccion_metricas_v2.py [OPCIONES]
```

### Opciones disponibles:

#### Obligatorias:
- `-folder RUTA`  
  Ruta a la carpeta con los datos de los sujetos a procesar (requerido)
  
- `-segmentacion {SWM,DWM} [SWM,DWM...]`  
  Se especifica de cual fascículos se quiere obtener las metricas de difusión:  
  - `SWM` para fibras cortas (Superficial White Matter)  
  - `DWM` para fibras largas (Deep White Matter)  
  Puedes especificar uno o ambos tipos

#### Opcionales:

- `-difusion_metric {FA,ADC,MD,RD}` #Este paso esta en proceso 
  Métrica de difusión para los cálculos:  
  - `FA`: Anisotropía Fraccional  
  - `ADC`: Coeficiente de Difusión Aparente  
  - `MD`: Difusividad Media  
  - `RD`: Difusividad Radial  

- `-get_image {True,False,y,n,yes,no}{Default:True}`  
    Calculas las imagenes (FA,ADC,RD,AD) del cual se extraeran los valores:  
  - `True`: Calcula estas imagenes utilizando Mrtrix3  
  - `False`: No se calcula, estas imagenes ya existen

- `-precise {True,False}`  
  Método de remuestreo de fibras:  
  
  - `False` (default): Todas las fibras usan el mismo número de puntos (basado en la fibra más larga / (tamaño_vóxel*0.75)).
      * Ventaja: Rápido.
      * Limitación: Fibras cortas pueden tener puntos redundantes en el mismo vóxel, lo que influye muy levemente en el valor promedio del fascículo.
  - `True`: Cada fibra se remuestrea con puntos proporcionales a su longitud.
      * Ventaja: Precisión mejorada.
      * Limitación: MAYOR tiempo de cálculo."""

### Ejemplos de uso:

1. Procesamiento básico:
```bash
python nombre_del_script.py -folder /datos/pacientes -segmentacion SWM
```

2. Procesamiento completo con todas las opciones:
```bash
python nombre_del_script.py -folder /datos/pacientes -segmentacion SWM DWM -difusion_metric FA -get_image True -precise True
```

## Detalles técnicos

### Diferencias en el metodo de remuestreo (precise)
