
Este flujo de trabajo carga datos de expresión de pacientes con cáncer de mama del TCGA y ejecuta un análisis de expresión diferencial de muestras del subtipo triple-negativo (TNBC) y luminal A (Lum A).

Para ejecutar el flujo introducir en consola el comando
"$ bash tcga_dea.sh" y añadir un parámetro de entre los siguientes:

- 'All' en caso de querer ejecutar todo el flujo de trabajo desde la descarga de datos hasta el análisis.
$ bash tcga_dea.sh All

- 'Data' en caso de sólo querer cargar los datos.
$ bash tcga_dea.sh Data

- 'Analysis' si los datos ya están cargados y sólo se desea efectuar el análisis.
$ bash tcga_dea.sh Analysis

Los datos cargados se guardan en la carpeta Data.

Los resultados se almacenan en la carpeta Results. Tras ejecutar el análisis deberíamos encontrar:
- Dos archivos .csv correspondientes al análisis de expresión diferencial y al análisis de enriquecimiento funcional.
- Un volcano plot
- Un dotplot
- Un heatmap de la matriz de datos de expresión de los genes significativos obtenidos del análisis de expresión diferencial
- Una imagen de los clusters realizados
- Una curva de Kaplan-Meier con las probabilidades de supervivencia de los dos grupos obtenidos en el clustering.

Los scripts de R se encuentran en la carpeta modulesR.

Junto al script bash también encontramos un informe más detallado en formato html confeccionado en Rmarkdown.

Versión R utilizada: 3.5.1

Paquetes R requeridos:
- RTCGAToolbox
- readxl
- dplyr
- limma
- edgeR
- calibrate
- org.Hs.eg.db
- biomaRt
- clusterProfiler
- gplots
- factoextra
- pROC
- survival
- survminer

