# Análisis integrómico de asociaciones cis-mQTL en ictus isquémico

## Descripción

Este repositorio contiene los scripts utilizados para el análisis bioinformático y estadístico desarrollado en el Trabajo de Fin de Grado del Grado en Ciencias Biomédicas.

El estudio se centra en la identificación de asociaciones cis-mQTL y en la evaluación de posibles diferencias en la regulación epigenética entre individuos con ictus isquémico y controles, mediante la integración de datos genéticos y de metilación del ADN procedentes de la cohorte española GENERACIÓN 1-2.

---

## Estructura del repositorio

### generacion1_2/

Scripts relacionados con el procesamiento y preparación de datos de la cohorte española GENERACIÓN 1-2.

Incluye:
- preparación de datos,
- control de calidad,
- filtrado de variantes y CpGs.

#### Scripts
- `01_preparacion_datos_GENERACION1_2.Rmd`
- `02_control_calidad_y_filtrado.Rmd`

---

### mQTL/

Scripts empleados en el análisis cis-mQTL y en la interpretación funcional de los resultados.

Incluye:
- exploración de datos,
- análisis de heterogeneidad entre casos y controles,
- análisis de enriquecimiento funcional y pathways.

#### Scripts
- `03_data_exploration_mQTL.R`
- `04_heterogeneity_analysis.R`
- `05_enrichment_pathways.R`

---

## Software y paquetes utilizados

### Software
- R version 4.X

### Principales paquetes utilizados
- data.table
- dplyr
- ggplot2
- plyr
- clusterProfiler
- ReactomePA
- enrichR
- stringr

---

## Objetivo del estudio

El objetivo principal del trabajo fue analizar la relación entre variantes genéticas y niveles de metilación del ADN mediante asociaciones cis-mQTL, así como evaluar posibles diferencias entre individuos con ictus isquémico y controles, con el fin de explorar mecanismos epigenéticos potencialmente implicados en la enfermedad.

---

## Autora

Leyre Franquet Navío  
Grado en Ciencias Biomédicas  
Trabajo de Fin de Grado (TFG)

---

## Nota

Los datos originales utilizados en el estudio no se incluyen en este repositorio debido a restricciones de tamaño y confidencialidad. El repositorio contiene únicamente los scripts y recursos necesarios para documentar el flujo de análisis realizado.
