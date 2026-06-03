# Caracterización del Microentorno Inmune Tumoral en Tumores Sólidos Pediátricos utilizando Datos Transcriptómicos

Este repositorio contiene el pipeline bioinformático desarrollado para el Trabajo de Fin de Máster (TFM) en Ciencia de Datos. El proyecto implementa un flujo de trabajo reproducible en **R** para automatizar la importación, normalización y análisis de datos transcriptómicos masivos enfocados en la caracterización del microambiente tumoral en oncología pediátrica.

## 📄 Resumen del Proyecto

El tratamiento actual del cáncer infantil carece a menudo de un enfoque especializado, aplicando adaptaciones de terapias para adultos que resultan altamente tóxicas. Este estudio aborda dicha limitación mediante el análisis de datos de secuenciación masiva de ARN (*bulk RNA-seq*) procedentes del consorcio internacional **cBioPortal**, correspondientes a cuatro cohortes de tumores sólidos pediátricos:
* **Tumores cerebrales** (CPTAC/CHOP 2020)
* **Neuroblastoma** (TARGET 2018)
* **Tumor de Wilms** (TARGET 2018)
* **Tumor Rabdoide** (TARGET 2018)

Mediante el uso de herramientas computacionales de **deconvolución celular** (CIBERSORTx, ConsensusTME y ESTIMATE), se infiere la composición de las poblaciones inmunes y estromales. Posteriormente, se aplican técnicas de *clustering* no supervisado y modelos de supervivencia (Kaplan-Meier y regresión de Cox) para identificar fenotipos inmunológicos correlacionados con variables clínicas y determinar biomarcadores celulares con valor pronóstico directo.

## 🎯 Objetivos

* **Objetivo Principal:** Caracterizar y comparar de forma exhaustiva el microambiente inmune de los principales tumores sólidos pediátricos para identificar patrones diferenciales y vulnerabilidades terapéuticas.
* **Objetivos Secundarios:**
  * Recopilar y normalizar datos genómicos a la métrica homogénea TPM (*Transcripts Per Million*).
  * Estimar la composición celular mediante deconvolución bioinformática.
  * Agrupar de manera no supervisada las muestras según su perfil inmunológico.
  * Evaluar asociaciones clinicopatológicas y moleculares clave (ej. impacto de la amplificación de *MYCN* o la carga mutacional *TMB*).
  * Determinar el valor pronóstico de los fenotipos inmunitarios identificados.

## 🛠️ Tecnologías y Paquetes Utilizados

El entorno de desarrollo principal ha sido **R (versión 4.4.1)**. Los paquetes fundamentales empleados se dividen según su propósito:

* **Manipulación y Limpieza de Datos:** `tidyverse`, `dplyr`, `janitor`, `tibble`, `purrr`.
* **Análisis Exploratorio y Control de Calidad:** `FactoMineR` y `factoextra` (para Análisis de Componentes Principales - PCA), `sva` (evaluación de efectos de lote con ComBat).
* **Deconvolución y Microentorno Tumoral:** `ConsensusTME`, `CIBERSORT`, `estimate`.
* **Visualización y Reporte:** `pheatmap` (mapas de calor con agrupamiento jerárquico), `ggplot2`, `patchwork`, `flextable`, `officer`.
* **Análisis Clínico:** `survival` (modelos de curvas de Kaplan-Meier y regresión de Cox).
