# Algoritmo de segmentación de fibras de fascículos para datos de tractografías en el espacio MNI152

## Descripción

Método de segmentación de fibras basado en la distancia euclidiana máxima de las fibras de un tractograma con respecto a las fibras del atlas de fascículos a utilizar para la segmentación. Este método fue propuesto originalmente por Labra [N. Labra et al., 2017] y posteriormente paralelizado por Vázquez [Vázquez et al., 2019].

En este repositorio se tienen los archivos necesarios para realizar la segmentacion de fascículos de fibras largas (DWM atlas [Guevara et al., 2012] y fasciculos de fibras cortas (SWM atlas [Godoy et al., 2021]). 


* AtlasRo/             : Carpeta con los fascículos de fibras cortas SWM (.bundles)
* AtlasRo_tck/         : Carpeta con los fascículos de fibras cortas SWM (.tck)
* atlas_faisceaux/     : Carpeta con los fascículos de fibras largas DWM (.bundles)
* atlas_faisceaux_tck/ : Carpeta con los fascículos de fibras largas DWM (.tck)
* AtlasRo.txt          : Archivo de texto que especifica la distancia Euclideana maxima utilizada para la segmentacion **de todos** los fascículos de fibras cortas 
* AtlasRo_estables.txt : Archivo de texto que especifica la distancia Euclideana maxima utilizada para la segmentacion de los 208 fascículos de fibras cortas más estables
* atlas_faisceaux.txt  : Archivo de texto que especifica la distancia Euclideana maxima utilizada para la segmentacion de los fascículos de fibras largas
* main_index           : Ejecutable de codigo escrito en C++ que genera la segmentación de datos


## Modo de Uso

> [!IMPORTANT]  
> **Para que la segmentación funcione correctamente el tractograma a segmentar debe estar en el espacio de referencia MNI152 y en el formato .bundles.**




## Citation