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
python nombre_del_script.py [OPCIONES]
```

### Opciones disponibles:

#### Obligatorias:
- `-folder RUTA`  
  Ruta a la carpeta con los datos de los sujetos a procesar (requerido)

#### Opcionales:
- `-h`, `--help`  
  Muestra este mensaje de ayuda y sale

- `-segmentacion {SWM,DWM} [SWM,DWM...]`  
  Tipo de segmentación de fascículos:  
  - `SWM` para fibras cortas (Superficial White Matter)  
  - `DWM` para fibras largas (Deep White Matter)  
  Puedes especificar uno o ambos tipos

- `-difusion_metric {FA,ADC,MD,RD}`  
  Métrica de difusión para los cálculos:  
  - `FA`: Anisotropía Fraccional  
  - `ADC`: Coeficiente de Difusión Aparente  
  - `MD`: Difusividad Media  
  - `RD`: Difusividad Radial  

- `-get_image {True,False}`  
  Si las imágenes de métricas no existen:  
  - `True`: Calcula usando Mrtrix3  
  - `False`: Usa las existentes (default)  

- `-precise {True,False}`  
  Método de remuestreo de fibras:  
  - `False` (default): Método rápido (mismo número de puntos para todas las fibras)  
  - `True`: Método preciso (puntos proporcionales a la longitud de cada fibra)  

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

### Métricas disponibles:
- **FA (Anisotropía Fraccional)**: Evalúa la integridad de la materia blanca
- **ADC/MD/RD**: Proporcionan información complementaria sobre la microestructura

### Métodos de segmentación:
- **SWM**: Optimizado para fibras cortas (<4cm)
- **DWM**: Especializado en tractos largos profundos

### Rendimiento:
- El modo `-precise False` es más rápido pero menos preciso
- El modo `-precise True` ofrece mayor precisión a costa de tiempo de cálculo
```

Este formato:
1. Es claro y fácil de leer directamente en GitHub
2. Usa formato Markdown estándar
3. Incluye ejemplos prácticos de uso
4. Explica las opciones técnicas de manera accesible
5. Destaca las opciones obligatorias vs opcionales

Puedes copiarlo directamente a tu README.md y personalizar los nombres de script o añadir cualquier detalle específico de tu implementación.
