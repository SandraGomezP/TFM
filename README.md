# Caracterización del Microentorno Inmune Tumoral en Tumores Sólidos Pediátricos utilizando Datos Transcriptómicos

Este repositorio contiene el pipeline bioinformático, los scripts de análisis estadístico y los flujos de trabajo desarrollados para mi **Trabajo de Fin de Máster (TFM) en Ciencia de Datos**. El objetivo principal del proyecto es caracterizar y comparar exhaustivamente la arquitectura celular del microambiente inmunitario en los principales tumores sólidos pediátricos a partir de datos ómicos públicos, permitiendo identificar patrones diferenciales y potenciales biomarcadores con valor pronóstico clínico.

## 👥 Autores y Dirección
* **Autora:** Sandra Gómez Peña
* **Directoras/Director de TF:** Rebeca Sanz Pamplona, Andrea Moreno Manuel, Samuel Paul Gallegos
* **Profesora Responsable de Asignatura (PRA):** Laia Subirats Maté
* **Institución:** Máster Universitario en Ciencia de Datos

---

## 🎯 Objetivos del Proyecto
* **General:** Caracterizar el microambiente inmune de tumores sólidos pediátricos (neuroblastoma, tumores cerebrales, tumor de Wilms y tumor rabdoide) mediante análisis bioinformático de datos transcriptómicos masivos (*bulk RNA-seq*).
* **Secundarios:**
  * Recopilar, curar y normalizar expresiones génicas procedentes de consorcios internacionales.
  * Estimar la composición celular infiltrante aplicando algoritmos de deconvolución computacional.
  * Identificar fenotipos inmunológicos comunes y diferenciales mediante clustering no supervisado.
  * Evaluar asociaciones biológicas entre perfiles celulares, alteraciones moleculares (p. ej., amplificación de *MYCN*, carga mutacional - TMB) y variables clínicas.
  * Determinar el impacto pronóstico de las poblaciones celulares utilizando modelos de supervivencia.

---

## 📊 Origen de los Datos
Los datos clínicos y transcriptómicos (*bulk RNA-seq*) fueron recuperados de manera anonimizada mediante la plataforma de acceso abierto **cBioPortal for Cancer Genomics**. Se integraron y armonizaron las siguientes 4 cohortes pediátricas:
1. **Pediatric Brain Cancer (CPTAC/CHOP, Cell 2020)** – `brain_cptac_2020` (188 muestras)
2. **Pediatric Neuroblastoma (TARGET, 2018)** – `nbl_target_2018_pub` (143 muestras)
3. **Pediatric Wilms' Tumor (TARGET, 2018)** – `wt_target_2018_pub` (130 muestras)
4. **Pediatric Rhabdoid Tumor (TARGET, 2018)** – `rt_target_2018_pub` (43 muestras)

*Nota metodológica:* Los niveles de expresión originales expresados en FPKM/RPKM fueron transformados homogéneamente a Transcritos por Millón (**TPM**) para garantizar la comparabilidad biológica e inter-muestra.

---

## 🛠️ Metodología y Pipeline Bioinformático
El flujo de trabajo (*pipeline*) automatizado e *in silico* se estructura en cuatro fases principales ejecutadas de forma iterativa bajo un entorno reproducible:

1. **Control de Calidad y Normalización:** Análisis exploratorio multivariante (PCA), estabilización de la varianza mediante `log2(x+1)` y mitigación de sesgos técnicos/efectos de lote (*batch effects*) con el algoritmo **ComBat**.
2. **Deconvolución Computacional (Cellular Estimation):** Inferencia cuantitativa de la abundancia del infiltrado celular inmune y estromal combinando tres aproximaciones complementarias:
   * **CIBERSORTx:** Resolución fina de hasta 22 subtipos de células inmunes (matriz LM22) mediante modelos de regresión lineal (SVR).
   * **ConsensusTME:** Robustez metodológica mediante puntuaciones de enriquecimiento por consenso y ssGSEA.
   * **ESTIMATE:** Evaluación global de la pureza tumoral e infiltración estromal/inmune.
3. **Clustering No Supervisado:** Agrupamiento jerárquico y algoritmos *K-means* para segregar a los pacientes en fenotipos inmunológicos diferenciados evaluando la estabilidad por Coeficiente de Silueta.
4. **Validación Clínica y Modelización de Supervivencia:** Contraste de hipótesis ajustados por tasa de falso descubrimiento (FDR), estimaciones de curvas de **Kaplan-Meier** y modelos de riesgos proporcionales multivariantes de **Cox**.

---

## 💻 Tecnologías y Paquetes Utilizados
El entorno técnico de desarrollo se ha implementado íntegramente en el lenguaje de programación **R (v4.4.1)**. Se emplearon de forma intensiva las siguientes librerías:

* **Manipulación y Limpieza de Datos:** `tidyverse`, `dplyr`, `janitor`, `tibble`, `purrr`.
* **Análisis Bioinformático y Deconvolución:** `ConsensusTME`, `CIBERSORT`, `estimate`, `sva` (algoritmo ComBat).
* **Análisis Estadístico y Supervivencia:** `survival`, `FactoMineR`, `factoextra`, `arsenal`.
* **Visualización y Reporte:** `ggplot2`, `pheatmap`, `patchwork`, `flextable`, `officer`.

