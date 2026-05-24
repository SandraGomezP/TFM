library(summarytools)
# library(readxl)
library(dplyr)
library(janitor)
library(tibble)
library(FactoMineR)
library(factoextra)
library(sva) # para combat
library(tidyverse)
library(purrr)
library(arsenal)
library(ConsensusTME)
library(CIBERSORT)
library(pheatmap)
library(patchwork)
library(flextable)
library(officer)
library(estimate)
library(survival)


setwd("Datasets")

## FUNCIONES

resumir_df <- function(clin_vars, df, nombreWord){
  tb_result <- tableby(
    formula = as.formula(paste("~", paste(clin_vars, collapse = " + "))),
    data = df,
    cat.stats = c("countpct"),  # cuenta y porcentaje
    numeric.stats = c("max", "min", "median","q1q3")
  )
  write2word(
    object = tb_result,
    file = nombreWord,
    title = "Tabla descriptiva de variables clínicas"
  )
}

norm_to_tpm <- function(df) { # Se normaliza igual para FPKM y RPKM
  (df / colSums(df, na.rm = TRUE)) * 1e6
}

preparar_densidad <- function(df, nombre) {
  df %>%
    select(where(is.numeric)) %>%
    select(-any_of(c("Entrez_Gene_Id", "Hugo_Symbol"))) %>%
    pivot_longer(everything(), names_to = "Sample", values_to = "Expression") %>%
    mutate(Dataset = nombre) %>%
    filter(Expression > 0) # Filtrar ceros para que no sesguen el logaritmo si hay demasiados
}

# preprocesar_pca <- function(df) {
#   df_procesado <- df %>%
#     # Filtro de vacíos y eliminación de Entrez (usamos any_of para que no falle si no existe)
#     filter(Hugo_Symbol != "" & !is.na(Hugo_Symbol)) %>%
#     select(-any_of("Entrez_Gene_Id")) %>%
#     distinct() %>%
#     # Transposición
#     t() %>%
#     as.data.frame() %>%
#     # Ponemos la primera fila como nombres de columnas y elimina la fila
#     row_to_names(row_number = 1) %>% 
#     .[, colSums(. != 0, na.rm = TRUE) > 0] %>%
#     setNames(make.unique(colnames(.))) %>% 
#     # Conversión a numérico y limpieza de ceros
#     mutate(across(everything(), as.numeric)) %>%
#     # Ordenar por muestras
#     arrange(rownames(.)) %>%
#   
#   return(df_procesado)
# }

preprocesar_pca <- function(df) {
  df_procesado <- df %>%
    # Filtro de vacíos y eliminación de Entrez (usamos any_of para que no falle si no existe)
    filter(Hugo_Symbol != "" & !is.na(Hugo_Symbol)) %>%
    select(-any_of("Entrez_Gene_Id")) %>%
    distinct(Hugo_Symbol, .keep_all = TRUE) %>%
    column_to_rownames("Hugo_Symbol") %>%
    t() %>%
    as.data.frame() %>%
    mutate(across(everything(), as.numeric)) %>%
    # Eliminar columnas que solo tienen 0
    select(where(~ sum(. != 0, na.rm = TRUE) > 0)) %>%
    setNames(make.unique(colnames(.))) %>% 
    # Ordenar por muestras
    arrange(rownames(.)) %>%
    
    return(df_procesado)
}

plot_pca_all <- function(pca_deconvolution_all, title){
  df_coords <- as.data.frame(pca_deconvolution_all$ind$coord)
  colnames(df_coords)[1:2] <- c("PC1", "PC2")
  df_coords$Dataset <- sub("_.*", "", rownames(df_coords))
  
  ggplot(df_coords, aes(x = PC1, y = PC2, color = Dataset)) +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    labs(title = paste("PCA", title),
         x = paste0("PC1 (", round(pca_deconvolution_all$eig[1, 2], 2), "%)"),
         y = paste0("PC2 (", round(pca_deconvolution_all$eig[2, 2], 2), "%)"))
}

plot_deconvolution_heatmap <- function(df, title){ # df_clin
  mat_heatmap <- df %>%
    select(-any_of(c("RMSE", "P-value", "Correlation", "P.value", "Immune_Score"))) %>%
    arrange(rownames(.)) %>%
    t(.) %>%
    .[apply(., 1, var, na.rm = TRUE) > 0, ] # Para quedarnos solo con las celulas inmunes con expresión

  
  annotation_data <- data.frame(
    Dataset = sub("_.*", "", colnames(mat_heatmap))
  )
  rownames(annotation_data) <- colnames(mat_heatmap)
  ann_colors <- list(Dataset = c(Brain = "#F8766D", Neuroblastoma = "#7CAE00", Wilms = "#C77CFF", Rhabdoid = "#00BFC4"))
  
  pheatmap(mat_heatmap, 
           scale = "column",
           annotation_col = annotation_data,
           annotation_colors = ann_colors,
           breaks = seq(-4.8, 4.8, length.out = 101), # -4.8, 4.8
           show_colnames = FALSE,
           main = title,
           color = colorRampPalette(c("blue", "white", "red"))(100))
}

estimate <- function(df_rnaseq_log2, dataset){
  df_estimate <- df_rnaseq_log2 %>%
    select(-any_of("Entrez_Gene_Id")) %>%
    filter(!is.na(Hugo_Symbol)) %>%
    group_by(Hugo_Symbol) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
    column_to_rownames(var = "Hugo_Symbol")
  
  write.table(df_estimate, file = paste0(dataset, "_estimate.txt"), sep = "\t", quote = FALSE, col.names = NA)
  
  filterCommonGenes(input.f = paste0(dataset, "_estimate.txt"), 
                    output.f = paste0(dataset, "_genes_filtered.gct"), 
                    id = "GeneSymbol")
  
  estimateScore(input.ds = paste0(dataset, "_genes_filtered.gct"), 
                output.ds = paste0("estimate_", dataset, ".gct"))
  
  scores_df <- read.table(paste0("estimate_", dataset, ".gct"), skip = 2, header = TRUE, check.names = FALSE) %>%
    t() %>%
    as.data.frame() %>%
    `colnames<-`(.[which(rownames(.) == "NAME"), ]) %>%
    filter(!(rownames(.) %in% c("NAME", "Description"))) %>%
    mutate(across(everything(), ~ as.numeric(as.character(.))))
  
  return(scores_df)
}

correlate <- function(df_data, df_clin, clin_vars, feature_name = "Feature",
                      exclude_cols = c("RMSE", "P-value", "Correlation", "P.value", "Immune_Score", "NAME", "Description")) {
  
  # Función interna para limpiar IDs
  clean_ids <- function(x) {
    x <- gsub("^[^_]+_", "", x)
    x <- gsub("-", ".", x)
    x <- gsub("^X(?=[0-9])", "", x, perl = TRUE)
    return(toupper(x))
  }
  
  # Limpieza de Identificadores (rownames)
  rownames(df_data) <- clean_ids(rownames(df_data))
  df_clin$Sample.ID <- clean_ids(df_clin$Sample.ID)
  rownames(df_clin) <- df_clin$Sample.ID
  
  # Alinear muestras comunes
  common_samples <- intersect(rownames(df_data), rownames(df_clin))
  df_data <- df_data[common_samples, , drop = FALSE]
  df_clin <- df_clin[common_samples, , drop = FALSE]
  
  # Identificar variables a testear (excluyendo metadatos)
  features <- setdiff(colnames(df_data), exclude_cols)
  results_list <- list()
  
  for (feat in features) {
    # Asegurar que el dato sea numérico
    feat_data <- as.numeric(as.character(df_data[[feat]]))
    
    for (var in clin_vars) {
      clin_data <- df_clin[[var]]
      
      # Limpieza de NAs y alineación
      blacklist <- c("Unavailable", "Not Reported", "UNKNOWN", "Unknown")
      valid_idx <- !is.na(clin_data) & !is.na(feat_data) & !(as.character(clin_data) %in% blacklist)
      
      if (sum(valid_idx) < 3) next
      
      x <- clin_data[valid_idx]
      y <- feat_data[valid_idx]
      
      # Verificación de varianza en el feature
      if (sd(y, na.rm = TRUE) == 0 || is.na(sd(y, na.rm = TRUE))) next
      
      # CASO 1: Variable Clínica Numérica
      if (is.numeric(x)) {
        if (sd(x, na.rm = TRUE) == 0) next
        
        norm_x <- shapiro.test(x)$p.value > 0.05
        norm_y <- shapiro.test(y)$p.value > 0.05
        method <- if (norm_x && norm_y) "pearson" else "spearman"
        res <- cor.test(x, y, method = method, exact = FALSE)
        
        results_list[[paste(feat, var, sep = "_")]] <- data.frame(
          Feature = feat, Variable = var, Type = "Numeric",
          Method = method, P_val = res$p.value, Estimate = as.numeric(res$estimate)
        )
        
        # CASO 2: Variable Clínica Categórica
      } else if (is.factor(x) || is.character(x)) {
        x <- as.factor(x)
        levels_present <- droplevels(x)
        if (length(levels(levels_present)) < 2) next
        
        norm_y <- shapiro.test(y)$p.value > 0.05
        
        if (length(levels(levels_present)) == 2) {
          res <- if (norm_y) t.test(y ~ x) else wilcox.test(y ~ x)
          method <- if (norm_y) "t-test" else "Wilcoxon"
          p_val <- res$p.value
          
        } else {
          res <- if (norm_y) aov(y ~ x) else kruskal.test(y ~ x)
          method <- if (norm_y) "ANOVA" else "Kruskal-Wallis"
          
          # Extracción segura para evitar el error de filas
          if (norm_y) {
            sum_res <- summary(res)[[1]]
            p_val <- if(!is.null(sum_res[["Pr(>F)"]])) sum_res[["Pr(>F)"]][1] else NA
          } else {
            p_val <- res$p.value
          }
        }
        
        # Solo guardar si p_val es válido
        if(!is.na(p_val)) {
          results_list[[paste(feat, var, sep = "_")]] <- data.frame(
            Feature = feat, Variable = var, Type = "Categorical",
            Method = method, P_val = p_val, Estimate = NA
          )
        }
      }
    }
  }
  
  if (length(results_list) == 0) return(NULL)
  
  # Consolidar y ajustar nombres
  final_df <- do.call(rbind, results_list) %>%
    mutate(FDR = p.adjust(P_val, method = "fdr"))
  
  # MODIFICACIÓN MÍNIMA: Usar el argumento feature_name
  colnames(final_df)[1] <- feature_name
  rownames(final_df) <- NULL
  
  return(final_df)
}

plot_correlations <- function(deconvolution_df, df_clin, clin_vars, correlations_df, feature_name = "Feature") {
  clean_ids <- function(x) {
    x <- gsub("^[^_]+_", "", x)
    x <- gsub("-", ".", x)
    x <- gsub("^X(?=[0-9])", "", x, perl = TRUE)
    return(toupper(x))
  }
  
  rownames(deconvolution_df) <- clean_ids(rownames(deconvolution_df))
  df_clin$Sample.ID <- clean_ids(df_clin$Sample.ID)
  rownames(df_clin) <- df_clin$Sample.ID
  
  common_samples <- intersect(rownames(deconvolution_df), rownames(df_clin))
  deconvolution_df <- deconvolution_df[common_samples, , drop = FALSE]
  df_clin <- df_clin[common_samples, , drop = FALSE]
  
  exclude_cols <- c("RMSE", "P-value", "Correlation", "P.value", "Immune_Score")
  cell_types <- setdiff(colnames(deconvolution_df), exclude_cols)
  
  # Lista para guardar los "mega-gráficos" finales
  final_plots <- list()
  
  # Iterar por cada variable clínica
  for (var in clin_vars) {
    plot_list <- list()
    clin_data <- df_clin[[var]]
    
    # Determinar tipo de variable
    is_numeric_var <- is.numeric(clin_data)
    
    for (cellType in cell_types) {
      cell_data <- deconvolution_df[[cellType]]
      
      p_val_row <- correlations_df[correlations_df[[feature_name]] == cellType & 
                                     correlations_df$Variable == var, ]
      
      # Formatear el p-valor (si no existe, la variable no tiene varianza, así que no mostramos nada)
      p_label <- if(nrow(p_val_row) > 0) {
        sprintf("p = %.2e", p_val_row$P_val[1])
      } else {
        next
      }
      
      # Crear DF temporal para el gráfico
      df_temp <- data.frame(Clinical = clin_data, Feature = cell_data)
      df_temp <- df_temp[!is.na(df_temp$Clinical) & !is.na(df_temp$Feature), ]
      
      blacklist <- c("Unavailable", "Not Reported", "UNKNOWN", "Unknown")
      df_temp <- df_temp[!(df_temp$Clinical %in% blacklist), ]
      
      if (nrow(df_temp) < 3) next
      
      # Construir el gráfico base
      p <- ggplot(df_temp, aes(x = Clinical, y = Feature)) +
        labs(title = cellType, subtitle = p_label, x = NULL, y = NULL) +
        theme_bw(base_size = 10) +
        theme(
          axis.text.x = element_text(angle = 10, vjust = 1, hjust = 1, size = rel(0.8)),
          axis.text.y = element_text(size = rel(0.8)),
          panel.grid.minor = element_blank(),
          plot.margin = margin(5, 5, 5, 5))
      
      if (is_numeric_var) {
        # Scatter plot + Línea de regresión lineal
        p <- p + 
          geom_point(alpha = 0.5, color = "steelblue") +
          geom_smooth(method = "lm", color = "darkred", fill = "lightcoral")
      } else {
        # Boxplot + Puntos para variables categóricas
        df_temp$Clinical <- as.factor(df_temp$Clinical)
        p <- p + 
          geom_boxplot(aes(fill = Clinical), outlier.shape = NA, alpha = 0.7) +
          geom_jitter(width = 0.2, alpha = 0.4) +
          theme(legend.position = "none") # Quitar leyenda para ahorrar espacio
      }
      
      plot_list[[cellType]] <- p
    }
    
    # Organizar subgráficos (Facet-like) usando patchwork
    if (length(plot_list) > 0) {
      combined_plot <- wrap_plots(plot_list, ncol = 4) + 
        plot_annotation(
          title = paste("Correlation with:", var),
          theme = theme(plot.title = element_text(size = 16, face = "bold"))
        )
      
      final_plots[[var]] <- combined_plot
    }
  }
  
  return(final_plots)
}

tabla_correlations <- function(correlations_df, nombreWord, columnName = "Feature"){
  tabla_matriz <- correlations_df %>%
    select(columnName, Variable, P_val) %>%
    pivot_wider(names_from = Variable, values_from = P_val) %>%
    mutate(across(-columnName, ~ formatC(.x, format = "e", digits = 2)))
  
  # 2. Crear la flextable
  ft <- flextable(tabla_matriz) %>%
    set_header_labels(columnName = "Tipo Celular") %>%
    theme_booktabs() %>% 
    autofit() %>%
    bold(part = "header") %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "all")
  
  # 3. Exportar a Word
  doc <- read_docx() 

  doc <- officer::body_add_par(doc, "Matriz de P-valores", style = "heading 1")
  doc <- flextable::body_add_flextable(doc, value = ft)
  
  print(doc, target = paste(nombreWord, ".docx"))
}

plot_significant_correlations <- function(deconvolution_df, df_clin, clin_vars, correlations_df, feature_name = "Feature") {

  sig_correlations <- correlations_df %>%
    filter(P_val < 0.05)

  # Identificar qué variables clínicas sobrevivieron al filtro
  sig_vars <- intersect(clin_vars, unique(sig_correlations$Variable))
  
  plots <- plot_correlations(
    deconvolution_df = deconvolution_df,
    df_clin = df_clin,
    clin_vars = sig_vars,
    correlations_df = sig_correlations,
    feature_name = feature_name
  )
  
  return(plots)
}


## IMPORTAR DATASETS

brain <- read.delim("brain_cptac_2020_clinical_data.tsv", header=TRUE)
brain_rnaseq_log2 <- read.delim("brain_data_mrna_seq_v2_rsem.txt") %>%
  mutate(across(where(is.numeric), ~ pmax(0, .x)))
brain_rnaseq <- brain_rnaseq_log2 %>%
  mutate(across(where(is.numeric), ~ 2^.x - 1)) %>%
  mutate(across(-c(Hugo_Symbol), ~ ( . / sum(. , na.rm = TRUE)) * 1e6))
brain_rnaseq_log2 <- brain_rnaseq %>%
  mutate(across(where(is.numeric), ~ log2(.x + 1)))

neuroblastoma <- read.delim("nbl_target_2018_pub_clinical_data.tsv", header=TRUE)
neuroblastoma_rnaseq <- read.delim("nbl_data_mrna_seq_rpkm.txt") %>%
  mutate(across(-c(Hugo_Symbol, Entrez_Gene_Id), ~ ( . / sum(. , na.rm = TRUE)) * 1e6)) # Normalizar a TPM
neuroblastoma_rnaseq_log2 <- neuroblastoma_rnaseq %>%
  mutate(across(where(is.numeric) & -Entrez_Gene_Id, ~ log2(.x + 1)))
  # mutate_at(vars(-(1:2)), ~ log2(as.numeric(as.character(.)) + 1))

wilms <- read.delim("wt_target_2018_pub_clinical_data.tsv", header=TRUE)
wilms_rnaseq <- read.delim("wt_data_mrna_seq_rpkm.txt") %>%
  mutate(across(-c(Hugo_Symbol, Entrez_Gene_Id), ~ ( . / sum(. , na.rm = TRUE)) * 1e6))
wilms_rnaseq_log2 <- wilms_rnaseq %>%
  mutate(across(where(is.numeric) & -Entrez_Gene_Id, ~ log2(.x + 1)))
  
rhabdoid <- read.delim("rt_target_2018_pub_clinical_data.tsv", header=TRUE)
rhabdoid_rnaseq <- read.delim("rt_data_mrna_seq_rpkm.txt") %>%
  mutate(across(-c(Hugo_Symbol, Entrez_Gene_Id), ~ ( . / sum(. , na.rm = TRUE)) * 1e6))
rhabdoid_rnaseq_log2 <- rhabdoid_rnaseq %>%
  mutate(across(where(is.numeric) & -Entrez_Gene_Id, ~ log2(.x + 1)))


## DESCRIPTIVA

variables_originales <- list(
  Brain_CPTAC = colnames(brain),
  Neuroblastoma_TARGET = colnames(neuroblastoma),
  Wilms_TARGET = colnames(wilms),
  Rhabdoid_TARGET = colnames(rhabdoid)
)

# Para verlas en consola de forma limpia:
lapply(variables_originales, function(x) paste(x, collapse = ", "))

brain_clin <- brain %>%
  filter(Sample.ID %in% (colnames(brain_rnaseq) %>% str_remove("^X") %>% str_replace_all("\\.", "-"))) 
  # mutate(Sample.ID = paste0("Brain_X", gsub("-", ".", Sample.ID)))
# brain_clin_vars <- c("AGE", "Sex", "Race", "Mutation.Count", "Cancer.Type")
# resumir_df(brain_clin_vars, brain_clin, "brain")
brain_clin_vars <- brain_clin %>%
  select(-Study.ID, -Patient.ID, -Sample.ID, -External.Patient.ID) %>%
  colnames()

neuroblastoma_clin <- neuroblastoma %>%
  filter(Sample.ID %in% (colnames(neuroblastoma_rnaseq) %>% str_remove("^X") %>% str_replace_all("\\.", "-"))) %>%
  mutate(Sample.ID = paste0("Brain_X", gsub("-", ".", Sample.ID)))
# neuroblastoma_clin_vars <- c("Diagnosis.Age", "Sex", "Race.Category", "Mutation.Count", "Cancer.Type", "Fraction.Genome.Altered", "TMB..nonsynonymous.")
# resumir_df(neuroblastoma_clin_vars, neuroblastoma_clin, "neuroblastoma")
neuroblastoma_clin_vars <- neuroblastoma_clin %>%
  select(-Study.ID, -Patient.ID, -Sample.ID) %>%
  colnames()

wilms_clin <- wilms %>%
  filter(Sample.ID %in% (colnames(wilms_rnaseq) %>% str_remove("^X") %>% str_replace_all("\\.", "-")))
# wilms_clin_vars <- c("Diagnosis.Age", "Sex", "Race.Category", "Mutation.Count", "Cancer.Type", "Fraction.Genome.Altered", "TMB..nonsynonymous.")
# resumir_df(wilms_clin_vars, wilms_clin, "wilms")
wilms_clin_vars <- wilms_clin %>%
  select(-Study.ID, -Patient.ID, -Sample.ID) %>%
  colnames()

rhabdoid_clin <- rhabdoid %>%
  filter(Sample.ID %in% (colnames(rhabdoid_rnaseq) %>% str_remove("^X") %>% str_replace_all("\\.", "-")))
# rhabdoid_clin_vars <- c("Diagnosis.Age", "Sex", "Race.Category", "Mutation.Count", "Cancer.Type", "TMB..nonsynonymous.")
# resumir_df(rhabdoid_clin_vars, rhabdoid_clin, "rhabdoid")
rhabdoid_clin_vars <- rhabdoid_clin %>%
  select(-Study.ID, -Patient.ID, -Sample.ID) %>%
  colnames()

all_clin <- 


## DENSITY PLOT

lista_rnaseq_log2 <- list(
  Brain = brain_rnaseq_log2,
  Neuroblastoma = neuroblastoma_rnaseq_log2,
  Wilms = wilms_rnaseq_log2,
  Rhabdoid = rhabdoid_rnaseq_log2
)

all_density_log2 <- map2_dfr(lista_rnaseq_log2, names(lista_rnaseq_log2), preparar_densidad) # Aplicamos la función a toda la lista

ggplot(all_density_log2, aes(x = Expression, fill = Dataset)) +
  geom_density(alpha = 0.5) + # Alpha para ver dónde se solapan
  scale_fill_manual(values = c("Brain" = "#F8766D", 
                               "Neuroblastoma" = "#7CAE00", 
                               "Wilms" = "#C77CFF", 
                               "Rhabdoid" = "#00BFC4")) +
  labs(
    title = "Comparación de Distribución de Expresión",
    x = "log2(Expresión + 1)",
    y = "Densidad"
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom")


## PCA

brain_pca <- preprocesar_pca(brain_rnaseq)
pca_brain <- PCA(brain_pca, graph = FALSE)
fviz_pca_ind(pca_brain, label = "none", title = "PCA Brain Tumor cBioPortal", col.ind = "#F8766D") +
  theme(text = element_text(size = 15))

neuroblastoma_pca <- preprocesar_pca(neuroblastoma_rnaseq)
  # mutate(across(where(is.numeric), ~ replace_na(., 0)))
pca_neuroblastoma <- PCA(neuroblastoma_pca, graph = FALSE)
fviz_pca_ind(pca_neuroblastoma, label = "none", title = "PCA neuroblastoma Tumor cBioPortal", col.ind = "#7CAE00") +
  theme(text = element_text(size = 15))

wilms_pca <- preprocesar_pca(wilms_rnaseq)
pca_wilms <- PCA(wilms_pca, graph = FALSE)
fviz_pca_ind(pca_wilms, label = "none", title = "PCA Wilm's Tumor cBioPortal", col.ind = "#C77CFF") +
  theme(text = element_text(size = 15))

rhabdoid_pca <- preprocesar_pca(rhabdoid_rnaseq)
pca_rhabdoid <- PCA(rhabdoid_pca, graph = FALSE)
fviz_pca_ind(pca_rhabdoid, label = "none", title = "PCA Rhabdoid's Tumor cBioPortal", col.ind = "#00BFC4") +
  theme(text = element_text(size = 15))


sets <- map(list(brain_pca,neuroblastoma_pca,wilms_pca,rhabdoid_pca), ~ .x %>% 
              as.data.frame() %>% 
              tibble::rownames_to_column("SampleID"))
names(sets) <- c("Brain", "Neuroblastoma", "Wilms", "Rhabdoid")
all_pca <- bind_rows(sets, .id = "DatasetName") %>%
  mutate(UniqueID = paste(DatasetName, SampleID, sep = "_")) %>% # Creamos un ID único combinando el origen y el ID de muestra
  select(, -c(SampleID, DatasetName)) %>%
  tibble::column_to_rownames("UniqueID") 
  # mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

pca_all_precombat <- PCA(all_pca, graph = FALSE)
fviz_pca_ind(pca_all_precombat, label = "none", 
             habillage = as.factor(gsub("_.*", "", rownames(pca_all_precombat$ind$coord))),
             pointshape = 19,   # El código 19 es el punto sólido redondo
             title = "PCA Conjunto") +
  theme(text = element_text(size = 15))

# Aplicar ComBat para reducir el efecto lote
matrix_for_combat <- all_pca %>% 
  t() %>% 
  as.matrix() %>%
  .[apply(., 1, var) > 0, ]
# all_deconvolution <- ComBat(dat = matrix_for_combat, batch = as.factor(all_pca$DatasetName)) %>%
#   t() %>%
#   as.data.frame() %>%
#   .[, apply(., 2, var, na.rm = TRUE) > 0]
# 
# pca_all_postcombat <- PCA(all_deconvolution, graph = FALSE, scale.unit = TRUE)
# plot_pca_all(pca_all_postcombat, "Conjunto Post-ComBat")


## PCA con expresión en log2

brain_pca_log2 <- preprocesar_pca(brain_rnaseq_log2)
pca_brain_log2 <- PCA(brain_pca_log2, graph = FALSE)
fviz_pca_ind(pca_brain_log2, label = "none", title = "PCA Brain Tumor log2 Expression", col.ind = "#F8766D") +
  theme(text = element_text(size = 15))

neuroblastoma_pca_log2 <- preprocesar_pca(neuroblastoma_rnaseq_log2)
  # select(where(~ !is.numeric(.x) || (any(.x != 0, na.rm = TRUE) & !all(is.na(.x)))))
  # mutate(across(where(is.numeric), ~ replace_na(., 0)))
pca_neuroblastoma_log2 <- PCA(neuroblastoma_pca_log2, graph = FALSE)
fviz_pca_ind(pca_neuroblastoma_log2, label = "none", title = "PCA neuroblastoma Tumor log2 Expression", col.ind = "#7CAE00") +
  theme(text = element_text(size = 15))

wilms_pca_log2 <- preprocesar_pca(wilms_rnaseq_log2)
pca_wilms_log2 <- PCA(wilms_pca_log2, graph = FALSE)
fviz_pca_ind(pca_wilms_log2, label = "none", title = "PCA Wilm's Tumor log2 Expression", col.ind = "#C77CFF") +
  theme(text = element_text(size = 15))

rhabdoid_pca_log2 <- preprocesar_pca(rhabdoid_rnaseq_log2)
pca_rhabdoid_log2 <- PCA(rhabdoid_pca_log2, graph = FALSE)
fviz_pca_ind(pca_rhabdoid_log2, label = "none", title = "PCA Rhabdoid's Tumor log2 Expression", col.ind = "#00BFC4") +
  theme(text = element_text(size = 15))


sets <- map(list(brain_pca_log2,neuroblastoma_pca_log2,wilms_pca_log2,rhabdoid_pca_log2), ~ .x %>% 
              as.data.frame() %>% 
              tibble::rownames_to_column("SampleID"))
names(sets) <- c("Brain", "Neuroblastoma", "Wilms", "Rhabdoid")
all_pca_log2 <- bind_rows(sets, .id = "DatasetName") %>%
  mutate(UniqueID = paste(DatasetName, SampleID, sep = "_")) %>% # Creamos un ID único combinando el origen y el ID de muestra
  select(, -c(SampleID, DatasetName)) %>%
  tibble::column_to_rownames("UniqueID") 
# mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

pca_all_precombat_log2 <- PCA(all_pca_log2, graph = FALSE)
fviz_pca_ind(pca_all_precombat_log2, label = "none", 
             habillage = as.factor(gsub("_.*", "", rownames(pca_all_precombat_log2$ind$coord))),
             pointshape = 19,   # El código 19 es el punto sólido redondo
             title = "PCA Conjunto Expresión en log2") +
  theme(text = element_text(size = 15))


## ESTIMATE

scores_brain <- estimate(brain_rnaseq_log2, "brain") 
scores_neuroblastoma <- estimate(neuroblastoma_rnaseq_log2, "neuroblastoma") 
scores_wilms <- estimate(wilms_rnaseq_log2, "wilms") 
scores_rhabdoid <- estimate(rhabdoid_rnaseq_log2, "rhabdoid") 

correlations_brain_estimate <- correlate(scores_brain, brain_clin, brain_clin_vars, "ScoreType")
correlations_neuroblastoma_estimate <- correlate(scores_neuroblastoma, neuroblastoma_clin, neuroblastoma_clin_vars, "ScoreType")
correlations_wilms_estimate <- correlate(scores_wilms, wilms_clin, wilms_clin_vars, "ScoreType")
correlations_rhabdoid_estimate <- correlate(scores_rhabdoid, rhabdoid_clin, rhabdoid_clin_vars, "ScoreType")

# tabla_correlations(correlations_brain_estimate, "BrainEstimate", "ScoreType")
# tabla_correlations(correlations_neuroblastoma_estimate, "NeuroblastomaEstimate", "ScoreType")
# tabla_correlations(correlations_wilms_estimate, "WilmsEstimate", "ScoreType")
# tabla_correlations(correlations_rhabdoid_estimate, "RhabdoidEstimate", "ScoreType")

## DECONVOLUCIÓN pre-combat

ConsensusTME::cancerAll
sig_matrix <- as.matrix(system.file("extdata", "LM22.txt", package = "CIBERSORT"))
# options(parallelly.fork.enable = TRUE)

brain_deconvolution_preCombat <- as.data.frame(t(brain_pca)) %>%
  rename_with(~ paste0("Brain_", .x))
  
consensus_brain_preCombat <- as.data.frame(t(consensusTMEAnalysis(bulkExp=as.matrix(brain_deconvolution_preCombat), 
                                                                  cancerType = "Unfiltered", statMethod = "gsva")))
# cibersort_brain_preCombat <- as.data.frame(cibersort(sig_matrix, as.matrix(brain_deconvolution_preCombat), QN=FALSE))

neuroblastoma_deconvolution_preCombat <- as.data.frame(t(neuroblastoma_pca)) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  rename_with(~ paste0("Neuroblastoma_", .x))
consensus_neuroblastoma_preCombat <- as.data.frame(t(consensusTMEAnalysis(bulkExp=as.matrix(neuroblastoma_deconvolution_preCombat), 
                                                                          cancerType = "Unfiltered", statMethod = "gsva")))
# cibersort_neuroblastoma_preCombat <- as.data.frame(cibersort(sig_matrix, as.matrix(neuroblastoma_deconvolution_preCombat), QN=FALSE))

wilms_deconvolution_preCombat <- as.data.frame(t(wilms_pca)) %>%
  rename_with(~ paste0("Wilms_", .x))
consensus_wilms_preCombat <- as.data.frame(t(consensusTMEAnalysis(bulkExp=as.matrix(wilms_deconvolution_preCombat), 
                                                                  cancerType = "Unfiltered", statMethod = "gsva")))
# cibersort_wilms_preCombat <- as.data.frame(cibersort(sig_matrix, as.matrix(wilms_deconvolution_preCombat), QN=FALSE))

rhabdoid_deconvolution_preCombat <- as.data.frame(t(rhabdoid_pca)) %>%
  rename_with(~ paste0("Rhabdoid_", .x))
consensus_rhabdoid_preCombat <- as.data.frame(t(consensusTMEAnalysis(bulkExp=as.matrix(rhabdoid_deconvolution_preCombat), 
                                                                     cancerType = "Unfiltered", statMethod = "gsva")))
# cibersort_rhabdoid_preCombat <- as.data.frame(cibersort(sig_matrix, as.matrix(rhabdoid_deconvolution_preCombat), QN=FALSE))

# PCA para verificar que no hay efecto lote
consensus_all_preCombat <- bind_rows(consensus_brain_preCombat, consensus_neuroblastoma_preCombat, 
                                     consensus_wilms_preCombat, consensus_rhabdoid_preCombat) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) # Reemplazar NA por 0
pca_consensus_all_preCombat <- PCA(consensus_all_preCombat, graph = FALSE, scale.unit = TRUE)
plot_pca_all(pca_consensus_all_preCombat, "Deconvolución con ConsensusTME")

# cibersort_all_preCombat <- bind_rows(cibersort_brain_preCombat, cibersort_neuroblastoma_preCombat, 
#                                      cibersort_wilms_preCombat, cibersort_rhabdoid_preCombat) %>%
#   mutate(across(everything(), ~ replace_na(.x, 0)))
# pca_cibersort_all_preCombat <- PCA(cibersort_all_preCombat, graph = FALSE, scale.unit = TRUE)
# plot_pca_all(pca_cibersort_all_preCombat, "Deconvolución con CIBERSORT")

# plot_deconvolution_heatmap(consensus_brain_preCombat, "Deconvolución Brain Consensus")
# plot_deconvolution_heatmap(consensus_neuroblastoma_preCombat, "Deconvolución Neuroblastoma Consensus")
# plot_deconvolution_heatmap(consensus_wilms_preCombat, "Deconvolución Wilms Consensus")
# plot_deconvolution_heatmap(consensus_rhabdoid_preCombat, "Deconvolución Rhabdoid Consensus")
plot_deconvolution_heatmap(consensus_all_preCombat, "Deconvolución All Consensus")
# plot_deconvolution_heatmap(cibersort_brain_preCombat, "Deconvolución Brain CIBERSORT")
# plot_deconvolution_heatmap(cibersort_neuroblastoma_preCombat, "Deconvolución Neuroblastoma CIBERSORT")
# plot_deconvolution_heatmap(cibersort_wilms_preCombat, "Deconvolución Wilms CIBERSORT")
# plot_deconvolution_heatmap(cibersort_rhabdoid_preCombat, "Deconvolución Rhabdoid CIBERSORT")
# plot_deconvolution_heatmap(cibersort_all_preCombat, "Deconvolución All CIBERSORT")


## DECONVOLUCIÓN post-combat

# brain_deconvolution_postCombat <- all_pca[grep("^Brain", rownames(all_pca), value = TRUE), ] %>%
#   select(-DatasetName) %>%
#   t() %>%
#   as.data.frame()
# consensus_brain_postCombat <- as.data.frame(t(consensusTMEAnalysis(bulkExp=as.matrix(brain_deconvolution_postCombat), 
#                                                                    cancerType = "Unfiltered")))
# cibersort_brain_postCombat <- as.data.frame(cibersort(sig_matrix, as.matrix(brain_deconvolution_postCombat), QN=FALSE))
# 
# neuroblastoma_deconvolution_postCombat <- all_pca[grep("^Neuroblastoma", rownames(all_pca), value = TRUE), ] %>%
#   select(-DatasetName) %>%
#   t() %>%
#   as.data.frame()
# consensus_neuroblastoma_postCombat <- as.data.frame(t(consensusTMEAnalysis(bulkExp=as.matrix(neuroblastoma_deconvolution_postCombat), 
#                                                                            cancerType = "Unfiltered")))
# cibersort_neuroblastoma_postCombat <- as.data.frame(cibersort(sig_matrix, as.matrix(neuroblastoma_deconvolution_postCombat), QN=FALSE))
# 
# wilms_deconvolution_postCombat <- all_pca[grep("^Wilms", rownames(all_pca), value = TRUE), ] %>%
#   select(-DatasetName) %>%
#   t() %>%
#   as.data.frame()
# consensus_wilms_postCombat <- as.data.frame(t(consensusTMEAnalysis(bulkExp=as.matrix(wilms_deconvolution_postCombat), 
#                                                         cancerType = "Unfiltered")))
# cibersort_wilms_postCombat <- as.data.frame(cibersort(sig_matrix, as.matrix(wilms_deconvolution_postCombat), QN=FALSE))
# 
# rhabdoid_deconvolution_postCombat <- all_pca[grep("^Rhabdoid", rownames(all_pca), value = TRUE), ] %>%
#   select(-DatasetName) %>%
#   t() %>%
#   as.data.frame()
# consensus_rhabdoid_postCombat <- as.data.frame(t(consensusTMEAnalysis(bulkExp=as.matrix(rhabdoid_deconvolution_postCombat), 
#                                                                       cancerType = "Unfiltered")))
# cibersort_rhabdoid_postCombat <- as.data.frame(cibersort(sig_matrix, as.matrix(rhabdoid_deconvolution_postCombat), QN=FALSE))

# PCA para verificar que no hay efecto lote
# consensus_all_postCombat <- bind_rows(consensus_brain_postCombat, consensus_neuroblastoma_postCombat, 
#                            consensus_wilms_postCombat, consensus_rhabdoid_postCombat) %>%
#   mutate(across(everything(), ~ replace_na(.x, 0))) # Reemplazar NA por 0
# pca_consensus_all_postCombat <- PCA(consensus_all_postCombat, graph = FALSE, scale.unit = TRUE)
# plot_pca_all(pca_consensus_all_postCombat, "Deconvolución con ConsensusTME")
# 
# cibersort_all_postCombat <- bind_rows(cibersort_brain_postCombat, cibersort_neuroblastoma_postCombat, 
#                                       cibersort_wilms_postCombat, cibersort_rhabdoid_postCombat) %>%
#   mutate(across(everything(), ~ replace_na(.x, 0)))
# pca_cibersort_all_postCombat <- PCA(cibersort_all_postCombat, graph = FALSE, scale.unit = TRUE)
# plot_pca_all(pca_cibersort_all_postCombat, "Deconvolución con CIBERSORT")

# plot_deconvolution_heatmap(consensus_brain_postCombat, "Deconvolución Brain Consensus")
# plot_deconvolution_heatmap(consensus_neuroblastoma_postCombat, "Deconvolución Neuroblastoma Consensus")
# plot_deconvolution_heatmap(consensus_wilms_postCombat, "Deconvolución Wilms Consensus")
# plot_deconvolution_heatmap(consensus_rhabdoid_postCombat, "Deconvolución Rhabdoid Consensus")
# plot_deconvolution_heatmap(consensus_all_postCombat, "Deconvolución All Consensus")
# plot_deconvolution_heatmap(cibersort_brain_postCombat, "Deconvolución Brain CIBERSORT")
# plot_deconvolution_heatmap(cibersort_neuroblastoma_postCombat, "Deconvolución Neuroblastoma CIBERSORT")
# plot_deconvolution_heatmap(cibersort_wilms_postCombat, "Deconvolución Wilms CIBERSORT")
# plot_deconvolution_heatmap(cibersort_rhabdoid_postCombat, "Deconvolución Rhabdoid CIBERSORT")
# plot_deconvolution_heatmap(cibersort_all_postCombat, "Deconvolución All CIBERSORT")


## ESTUDIO DE GENES: B4GALNT1, CNTRF, ERBB2, NTRK1

# compararExpresionGenes <- function(genes_interes, df_pre_post){
#   df_expresion <- df_pre_post %>%
#     rownames_to_column(var = "Sample") %>%
#     mutate(Dataset = str_extract(Sample, "^[^_]+")) %>% # Extraemos el Dataset del prefijo
#     select(Dataset, all_of(genes_interes)) %>% # Seleccionamos Dataset y los 4 genes
#     pivot_longer(cols = all_of(genes_interes),  # Pasamos a formato largo para que los genes estén en el eje X
#                  names_to = "Gene", 
#                  values_to = "Expression") %>%
#     mutate(Expression_Log2 = log2(Expression + 1))
#   
#   lista_datasets <- unique(df_expresion$Dataset)
#   
#   for (gen in genes_interes) {
#     df_subset <- df_expresion %>% 
#       filter(Gene == gen)
#     
#     p <- ggplot(df_subset, aes(x = Dataset, y = Expression_Log2, fill = Dataset)) +
#       geom_boxplot(outlier.shape = NA, alpha = 0.7) +
#       geom_jitter(width = 0.2, size = 1.2, alpha = 0.4) +
#       theme_minimal() +
#       labs(title = gen,
#            x = "Dataset",
#            y = "Expression (log2)") +
#       theme(legend.position = "none",
#             axis.text.x = element_text(face = "bold")) +
#       scale_fill_brewer(palette = "Set2")  # Para que los genes siempre tengan el mismo color
#     
#     print(p)
#   }
# }

# compararExpresionGenes(c("B4GALNT1", "CNTFR", "ERBB2", "NTRK1"), all_pca)
# compararExpresionGenes(c("B4GALNT1", "CNTFR", "ERBB2", "NTRK1"), all_deconvolution)


## CORRELACIONES CON VARIABLES CLÍNICAS

# correlations_brain_cibersort_preCombat <- correlate(cibersort_brain_preCombat, brain_clin, brain_clin_vars, "CellType")
# correlations_brain_consensus_preCombat <- correlate(consensus_brain_preCombat, brain_clin, brain_clin_vars)
# correlations_neuroblastoma_cibersort_preCombat <- correlate(cibersort_neuroblastoma_preCombat, neuroblastoma_clin, neuroblastoma_clin_vars, "CellType")
# correlations_neuroblastoma_consensus_preCombat <- correlate(consensus_neuroblastoma_preCombat, neuroblastoma_clin, neuroblastoma_clin_vars, "CellType")
# correlations_wilms_cibersort_preCombat <- correlate(cibersort_wilms_preCombat, wilms_clin, wilms_clin_vars, "CellType")
# correlations_wilms_consensus_preCombat <- correlate(consensus_wilms_preCombat, wilms_clin, wilms_clin_vars, "CellType")
# correlations_rhabdoid_cibersort_preCombat <- correlate(cibersort_rhabdoid_preCombat, rhabdoid_clin, rhabdoid_clin_vars, "CellType")
# correlations_rhabdoid_consensus_preCombat <- correlate(consensus_rhabdoid_preCombat, rhabdoid_clin, rhabdoid_clin_vars, "CellType")

# tabla_correlations(correlations_brain_cibersort_preCombat, "BrainCibersort_preCombat")
# tabla_correlations(correlations_neuroblastoma_cibersort_preCombat, "NeuroblastomaCibersort_preCombat")
# tabla_correlations(correlations_wilms_cibersort_preCombat, "WilmsCibersort_preCombat")
# tabla_correlations(correlations_rhabdoid_cibersort_preCombat, "RhabdoidCibersort_preCombat")
# tabla_correlations(correlations_brain_consensus_preCombat, "BrainConsensus_preCombat")
tabla_correlations(correlations_neuroblastoma_consensus_preCombat, "NeuroblastomaConsensus_preCombat", "CellType")
# tabla_correlations(correlations_wilms_consensus_preCombat, "WilmsConsensus_preCombat", "CellType")
# tabla_correlations(correlations_rhabdoid_consensus_preCombat, "RhabdoidConsensus_preCombat", "CellType")

# plot_significant_correlations(cibersort_brain_preCombat, brain_clin, brain_clin_vars, correlations_brain_cibersort_preCombat)
plot_significant_correlations(consensus_brain_preCombat, brain_clin, brain_clin_vars, correlations_brain_consensus_preCombat)
# plot_significant_correlations(cibersort_neuroblastoma_preCombat, neuroblastoma_clin, neuroblastoma_clin_vars, correlations_neuroblastoma_cibersort_preCombat)
plot_significant_correlations(consensus_neuroblastoma_preCombat, neuroblastoma_clin, neuroblastoma_clin_vars, correlations_neuroblastoma_consensus_preCombat, "CellType")
# plot_significant_correlations(cibersort_wilms_preCombat, wilms_clin, wilms_clin_vars, correlations_wilms_cibersort_preCombat)
# plot_significant_correlations(consensus_wilms_preCombat, wilms_clin, wilms_clin_vars, correlations_wilms_consensus_preCombat, "CellType")
# plot_significant_correlations(cibersort_rhabdoid_preCombat, rhabdoid_clin, rhabdoid_clin_vars, correlations_rhabdoid_cibersort_preCombat)
# plot_significant_correlations(consensus_rhabdoid_preCombat, rhabdoid_clin, rhabdoid_clin_vars, correlations_rhabdoid_consensus_preCombat, "CellType")

# correlations_brain_cibersort_postCombat <- correlate(cibersort_brain_postCombat, brain_clin, brain_clin_vars, "CellType")
# correlations_brain_consensus_postCombat <- correlate(consensus_brain_postCombat, brain_clin, brain_clin_vars, "CellType")
# correlations_neuroblastoma_cibersort_postCombat <- correlate(cibersort_neuroblastoma_postCombat, neuroblastoma_clin, neuroblastoma_clin_vars, "CellType")
# correlations_neuroblastoma_consensus_postCombat <- correlate(consensus_neuroblastoma_postCombat, neuroblastoma_clin, neuroblastoma_clin_vars, "CellType")
# correlations_wilms_cibersort_postCombat <- correlate(cibersort_wilms_postCombat, wilms_clin, wilms_clin_vars, "CellType")
# correlations_wilms_consensus_postCombat <- correlate(consensus_wilms_postCombat, wilms_clin, wilms_clin_vars, "CellType")
# correlations_rhabdoid_cibersort_postCombat <- correlate(cibersort_rhabdoid_postCombat, rhabdoid_clin, rhabdoid_clin_vars, "CellType")
# correlations_rhabdoid_consensus_postCombat <- correlate(consensus_rhabdoid_postCombat, rhabdoid_clin, rhabdoid_clin_vars, "CellType")

# tabla_correlations(correlations_brain_cibersort_postCombat, "BrainCibersort_postCombat")
# tabla_correlations(correlations_neuroblastoma_cibersort_postCombat, "NeuroblastomaCibersort_postCombat")
# tabla_correlations(correlations_wilms_cibersort_postCombat, "WilmsCibersort_postCombat")
# tabla_correlations(correlations_rhabdoid_cibersort_postCombat, "RhabdoidCibersort_postCombat")
# tabla_correlations(correlations_brain_consensus_postCombat, "BrainConsensus_postCombat")
# tabla_correlations(correlations_neuroblastoma_consensus_postCombat, "NeuroblastomaConsensus_postCombat")
# tabla_correlations(correlations_wilms_consensus_postCombat, "WilmsConsensus_postCombat")
# tabla_correlations(correlations_rhabdoid_consensus_postCombat, "RhabdoidConsensus_postCombat")

# plot_correlations(cibersort_brain, brain_clin, brain_clin_vars, correlations_brain_cibersort)
# plot_correlations(consensus_brain_preCombat, brain_clin, brain_clin_vars, correlations_brain_consensus_preCombat)
# plot_correlations(cibersort_neuroblastoma, neuroblastoma_clin, neuroblastoma_clin_vars, correlations_neuroblastoma_cibersort)
# plot_correlations(consensus_neuroblastoma_preCombat, neuroblastoma_clin, neuroblastoma_clin_vars, correlations_neuroblastoma_consensus_preCombat, feature_name = "CellType")
# plot_correlations(cibersort_wilms, wilms_clin, wilms_clin_vars, correlations_wilms_cibersort)
# plot_correlations(consensus_wilms_preCombat, wilms_clin, wilms_clin_vars, correlations_wilms_consensus_preCombat, feature_name = "CellType")
# plot_correlations(cibersort_rhabdoid, rhabdoid_clin, rhabdoid_clin_vars, correlations_rhabdoid_cibersort)
# plot_correlations(consensus_rhabdoid_preCombat, rhabdoid_clin, rhabdoid_clin_vars, correlations_rhabdoid_consensus_preCombat, feature_name = "CellType")

# plot_significant_correlations(cibersort_brain_postCombat, brain_clin, brain_clin_vars, correlations_brain_cibersort_postCombat)
# plot_significant_correlations(consensus_brain_postCombat, brain_clin, brain_clin_vars, correlations_brain_consensus_postCombat)
# plot_significant_correlations(cibersort_neuroblastoma_postCombat, neuroblastoma_clin, neuroblastoma_clin_vars, correlations_neuroblastoma_cibersort_postCombat)
# plot_significant_correlations(consensus_neuroblastoma_postCombat, neuroblastoma_clin, neuroblastoma_clin_vars, correlations_neuroblastoma_consensus_postCombat)
# plot_significant_correlations(cibersort_wilm_postCombats, wilms_clin, wilms_clin_vars, correlations_wilms_cibersort_postCombat)
# plot_significant_correlations(consensus_wilms_postCombat, wilms_clin, wilms_clin_vars, correlations_wilms_consensus_postCombat)
# plot_significant_correlations(cibersort_rhabdoid_postCombat, rhabdoid_clin, rhabdoid_clin_vars, correlations_rhabdoid_cibersort_postCombat)
# plot_significant_correlations(consensus_rhabdoid_postCombat, rhabdoid_clin, rhabdoid_clin_vars, correlations_rhabdoid_consensus_postCombat)





# plot_single_correlation <- function(deconvolution_df, df_clin, cell_type, clin_var, correlations_df) {
#   
#   all_plots_list <- plot_correlations(
#     deconvolution_df = deconvolution_df, 
#     df_clin = df_clin, 
#     clin_vars = clin_var, 
#     correlations_df = correlations_df,
#     feature_name = "CellType"
#   )
#   
#   # Extract the patchwork object for that clinical variable
#   combined_obj <- all_plots_list[[clin_var]]
#   
#   # if (is.null(combined_obj)) {
#   #   stop("No plots were generated for this clinical variable.")
#   # }
#   
#   # Because wrap_plots was used, the sub-plots are stored internally:
#   all_sub_plots <- c(list(combined_obj$patches$plots[[1]]), combined_obj$patches$plots)
#   
#   target_plot <- NULL
#   
#   # We iterate through the patches to find the one with the right title
#   potential_plots <-  combined_obj$patches$plots
#   # Add the "active" plot which is the last one added/main one
#   potential_plots[[length(potential_plots)+1]] <- combined_obj
#   
#   for(p in potential_plots) {
#     if(!is.null(p$labels$title) && p$labels$title == cell_type) {
#       target_plot <- p
#       break
#     }
#   }
#   
#   if (is.null(target_plot)) {
#     stop(paste("Could not find cell type:", cell_type, "in the plots."))
#   }
# 
#   target_plot <- target_plot + 
#     labs(caption = paste("Dataset:", sub("_.*", "", rownames(deconvolution_df)[1])),
#          x = clin_var, # Add X axis label back since it's a single plot now
#          y = "Proportion") + 
#     theme_bw(base_size = 25)
#   
#   return(target_plot)
# }


plot_single_correlation <- function(deconvolution_df, df_clin, cell_type, clin_var, correlations_df) {
  
  all_plots_list <- plot_correlations(
    deconvolution_df = deconvolution_df, 
    df_clin = df_clin, 
    clin_vars = clin_var, 
    correlations_df = correlations_df,
    feature_name = "Feature" # Si da errror, esto es Feature o CellType
  )
  
  combined_obj <- all_plots_list[[clin_var]]
  print(combined_obj)
  
  # Extraer todos los sub-plots de forma más robusta
  # Patchwork guarda el gráfico "base" y luego los "patches"
  potential_plots <- c(list(combined_obj), combined_obj$patches$plots)
  # print(potential_plots)
  
  target_plot <- NULL
  found_titles <- c() # Para debugging
  
  # Normalizar el nombre que buscamos (quitar puntos/guiones y pasar a minúsculas)
  clean_target <- tolower(gsub("[._ ]", "", cell_type))
  
  for(p in potential_plots) {
    current_title <- p$labels$title
    if(!is.null(current_title)) {
      found_titles <- c(found_titles, current_title)
      
      # Normalizar el título encontrado para la comparación
      clean_title <- tolower(gsub("[._ ]", "", current_title))
      
      if(clean_title == clean_target) {
        target_plot <- p
        break
      }
    }
  }
  # cat("Títulos encontrados en el objeto:", paste(found_titles, collapse = ", "), "\n")
  if (is.null(target_plot)) {
    cat("Títulos encontrados en el objeto:", paste(found_titles, collapse = ", "), "\n")
    stop(paste("No se encontró el tipo celular:", cell_type))
  }
  
  # Personalización final
  dataset_name <- sub("_.*", "", rownames(deconvolution_df)[1])
  
  target_plot <- target_plot + 
    labs(caption = paste("Dataset:", dataset_name),
         x = clin_var, 
         y = "Proportion",
         color = clin_var,  # Cambia el título de la leyenda de color
         fill = clin_var) + 
    theme_bw(base_size = 25)  
    # theme(
    #   axis.text.x = element_text(angle = 45, hjust = 1)
    # )
  
  return(target_plot)
}

plot_single_correlation(
  consensus_brain_preCombat, 
  brain_clin, 
  "Mast_cells", 
  "Chemotherapy", 
  correlations_brain_consensus_preCombat
)

plot_single_correlation(
  consensus_neuroblastoma_preCombat, 
  neuroblastoma_clin, 
  "Cytotoxic_cells", 
  "MYCN", 
  correlations_neuroblastoma_consensus_preCombat
)

plot_single_correlation(
  consensus_wilms_preCombat, 
  wilms_clin, 
  "NK_cells", 
  "Fraction.Genome.Altered", 
  correlations_wilms_consensus_preCombat
)

plot_single_correlation(
  consensus_rhabdoid_preCombat, 
  rhabdoid_clin, 
  "NK_cells", 
  "TMB..nonsynonymous.", 
  correlations_rhabdoid_consensus_preCombat
)

plot_single_correlation(
  consensus_wilms_preCombat, 
  wilms_clin, 
  "NK_cells", 
  "Fraction.Genome.Altered", 
  correlations_wilms_consensus_preCombat
)


cluster_deconvolution <- function(df, title){
  # 1. Preparación inicial de la matriz
  mat_heatmap <- df %>%
    select(-any_of(c("RMSE", "P-value", "Correlation", "P.value", "Immune_Score"))) %>%
    t() %>%
    .[apply(., 1, var, na.rm = TRUE) > 0, ]
  
  # 2. Primera ejecución "silenciosa" para extraer el clustering exacto
  # Usamos los mismos parámetros que quieres en el gráfico final
  res_temp <- pheatmap(mat_heatmap, 
                       scale = "column", 
                       silent = TRUE)
  
  clusters <- cutree(res_temp$tree_col, k = 4)
  clusters_named <- as.character(clusters)
  clusters_named[clusters_named == "1"] <- "Inflammatory"
  clusters_named[clusters_named == "2"] <- "T-cell Enriched"
  clusters_named[clusters_named == "3"] <- "Suppresor"
  clusters_named[clusters_named == "4"] <- "Citotoxic"
  
  # 3. Extraer los clusters del árbol que pheatmap generó
  # Esto garantiza que 'cutree' use el mismo dendrograma que ves en pantalla
  clusters <- cutree(res_temp$tree_col, k = 4)
  
  # 4. Crear el dataframe de anotación con los clusters extraídos
  annotation_data <- data.frame(
    Dataset = sub("_.*", "", colnames(mat_heatmap)),
    Cluster = factor(clusters_named, levels = c("Inflammatory", "T-cell Enriched", "Suppresor", "Citotoxic"))
  )
  rownames(annotation_data) <- colnames(mat_heatmap)
  
  # Colores para los datasets y los clusters
  ann_colors <- list(
    Dataset = c(Brain = "#F8766D", Neuroblastoma = "#7CAE00", Wilms = "#C77CFF", Rhabdoid = "#00BFC4"),
    Cluster = c("Inflammatory" = "#7570B3", "T-cell Enriched" = "#E7298A", "Suppresor" = "#66A61E", "Citotoxic" = "#E6AB02")
  )
  
  # 5. Generar el Heatmap FINAL con la anotación incorporada
  # pheatmap usará el orden del dendrograma ya calculado
  final_plot <- pheatmap(mat_heatmap,
                         scale = "column",
                         annotation_col = annotation_data,
                         annotation_colors = ann_colors,
                         breaks = seq(-4.8, 4.8, length.out = 101),
                         show_colnames = FALSE,
                         main = title,
                         cutree_cols = 4,
                         color = colorRampPalette(c("blue", "white", "red"))(100))
  
  # 6. Retornar el dataframe de clusters ordenado
  df_clusters <- data.frame(
    Sample = names(clusters),
    Cluster = annotation_data$Cluster
  ) %>% arrange(Cluster)
  
  return(df_clusters)
}

clusters_all <- cluster_deconvolution(consensus_all_preCombat, "Deconvolución All Consensus")

# plot_survival_clusters <- function(df_clin, cluster_all, cancer_type) {
#   
#   old_par <- par(cex = 1.5, font.main = 2) 
#   # Al final de la función restauramos la configuración original
#   on.exit(par(old_par))
#   
#   cluster_colors <- c("1" = "#BF0505", "2" = "#FF944B", "3" = "#DFDA00", "4" = "#663300")
#   
#   # Limpieza de IDs (Normalización para el merge)
#   clean_id <- function(x) {
#     x <- as.character(x)
#     x <- gsub(paste0("^", cancer_type, "_"), "", x, ignore.case = TRUE)
#     x <- gsub("^X", "", x)
#     x <- gsub("[^A-Za-z0-9]", "", x)
#     return(toupper(x))
#   }
# 
#   # Esta función extrae solo números y decimales de un string como "12.5 months"
#   clean_numeric <- function(x) {
#     x <- as.character(x)
#     # Elimina todo lo que no sea número o punto decimal
#     x_clean <- gsub("[^0-9.]", "", x)
#     return(as.numeric(x_clean))
#   }
#   
#   # Selección de columnas
#   if (tolower(cancer_type) == "brain") {
#     time_col <- "OS.Months"; status_col <- "OS.Status"
#   } else {
#     time_col <- "Overall.Survival..Months."; status_col <- "Overall.Survival.Status"
#   }
#   
#   # Filtrado y Match
#   curr_clusters <- cluster_all[grepl(paste0("^", cancer_type), cluster_all$Sample, ignore.case = TRUE), ]
#   
#   df_clin$ID_match <- clean_id(df_clin$Sample.ID)
#   curr_clusters$ID_match <- clean_id(curr_clusters$Sample)
#   
#   df_merged <- merge(df_clin, curr_clusters, by = "ID_match")
#   
#   # Aplicamos el limpiador de números a las columnas críticas
#   df_merged$time_clean <- clean_numeric(df_merged[[time_col]])
#   df_merged$status_clean <- clean_numeric(df_merged[[status_col]])
#   
#   # Quitamos filas donde el tiempo sea 0 o NA (survfit fallaría con ellas)
#   df_final <- df_merged[!is.na(df_merged$time_clean) & !is.na(df_merged$status_clean), ]
#   df_final <- df_final[df_final$time_clean > 0, ]
#   
#   # Ejecución del análisis
#   surv_obj <- Surv(df_final$time_clean, df_final$status_clean)
#   fit <- survfit(surv_obj ~ Cluster, data = df_final)
#   
#   # Log-rank test
#   sd <- survdiff(surv_obj ~ Cluster, data = df_final)
#   p_val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
#   p_txt <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 4)))
# 
#   plot(fit, col = 1:length(unique(df_final$Cluster)), lwd = 2,
#        main = toupper(cancer_type), xlab = "Months", ylab = "Survival Probability")
#   
#   legend("topright", legend = paste("Cluster", levels(as.factor(df_final$Cluster))), 
#          col = 1:length(unique(df_final$Cluster)), lty = 1, bty = "n", lwd = 2,
#          y.intersp = 0.2,   # <--- Reduce el interlineado (puedes probar con 0.6 o 0.8)
#          seg.len = 1.5)
#   
#   mtext(p_txt, side = 1, line = -2, adj = 0.95, font = 2, cex = 1.5)
# }

plot_survival_clusters <- function(df_clin, cluster_all, cancer_type) {
  
  # --- Configuración Visual ---
  old_par <- par(cex = 1.4, font.main = 2) 
  on.exit(par(old_par))
  
  # Paleta de colores personalizada
  cluster_colors <- c("Inflammatory" = "#7570B3", "T-cell Enriched" = "#E7298A", "Suppresor" = "#66A61E", "Citotoxic" = "#E6AB02")
  
  # 1. Limpieza de IDs y Datos
  clean_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", gsub("^X", "", gsub(paste0("^", cancer_type, "_"), "", as.character(x), ignore.case = TRUE))))
  clean_numeric <- function(x) as.numeric(gsub("[^0-9.]", "", as.character(x)))
  
  # 2. Selección de columnas según tipo de cáncer
  if (tolower(cancer_type) == "brain") {
    time_col <- "OS.Months"; status_col <- "OS.Status"
  } else {
    time_col <- "Overall.Survival..Months."; status_col <- "Overall.Survival.Status"
  }
  
  # 3. Match y Filtrado
  curr_clusters <- cluster_all[grepl(paste0("^", cancer_type), cluster_all$Sample, ignore.case = TRUE), ]
  df_clin$ID_match <- clean_id(df_clin$Sample.ID)
  curr_clusters$ID_match <- clean_id(curr_clusters$Sample)
  
  df_merged <- merge(df_clin, curr_clusters, by = "ID_match")
  df_merged$time_clean <- clean_numeric(df_merged[[time_col]])
  df_merged$status_clean <- clean_numeric(df_merged[[status_col]])
  
  df_final <- df_merged[!is.na(df_merged$time_clean) & df_merged$time_clean > 0, ]
  
  # Aseguramos que Cluster sea un factor para mantener el orden de los colores
  df_final$Cluster <- as.factor(df_final$Cluster)
  active_clusters <- levels(df_final$Cluster)
  cols_to_use <- cluster_colors[active_clusters]
  
  # 4. Análisis de Supervivencia
  fit <- survfit(Surv(time_clean, status_clean) ~ Cluster, data = df_final)
  sd <- survdiff(Surv(time_clean, status_clean) ~ Cluster, data = df_final)
  p_val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  p_txt <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 4)))
  
  # 5. Gráfico
  plot(fit, 
       col = cols_to_use, 
       lwd = 3, 
       main = toupper(cancer_type), 
       xlab = "OS (Months)", 
       ylab = "Survival Probability",
       cex.main = 1.7,                # Título muy grande
       cex.lab = 1.7,                   # Etiquetas de ejes grandes
       cex.axis = 1.7)
  
  par(xpd = TRUE)
  
  # Leyenda con colores vinculados y espaciado corregido
  legend(x = max(df_final$time_clean - 60, na.rm = TRUE), y = 1.25, # MIRAR POR QUE NO FUNCIONA TOPRIGHT
         legend = paste(active_clusters), 
         col = cols_to_use, 
         lty = 1, lwd = 3, bty = "n",
         y.intersp = 0.15,  # Espaciado compacto
         seg.len = 1,
         cex = 1.5)
  
  # P-valor más grande
  mtext(p_txt, side = 1, line = -2, adj = 0.95, font = 2, cex = 2)
}

plot_survival_clusters(brain_clin, clusters_all, "brain")
plot_survival_clusters(neuroblastoma_clin, clusters_all, "neuroblastoma")
plot_survival_clusters(wilms_clin, clusters_all, "wilms")
plot_survival_clusters(rhabdoid_clin, clusters_all, "rhabdoid")


# cox_analysis <- function(deconv_df, clinical_df, time_col, event_col, sample_id) {
#   
#   deconv_df <- deconv_df %>% 
#     as.data.frame() %>% 
#     rownames_to_column(var = sample_id)
#   
#   # 1. Unir dataframes por ID de muestra
#   # Asegúrate de que sample_id exista en ambos
#   combined_df <- inner_join(deconv_df, clinical_df, by = sample_id)
#   
#   if(is.character(combined_df[[event_col]])) {
#     combined_df[[event_col]] <- as.numeric(substr(combined_df[[event_col]], 1, 1))
#   }
#   
#   
#   # 2. Identificar qué columnas son los tipos celulares 
#   # (asumiendo que son todas las de deconv_df menos el ID)
#   cell_types <- setdiff(colnames(deconv_df), sample_id)
#   
#   results <- list()
#   
#   # 3. Bucle para ajustar Cox a cada tipo celular
#   for (cell in cell_types) {
#     # Creamos la fórmula dinámicamente
#     formula_str <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ `", cell, "`"))
#     
#     # Ajustar el modelo
#     fit <- tryCatch({
#       coxph(formula_str, data = combined_df)
#     }, error = function(e) return(NULL))
#     
#     if (!is.null(fit)) {
#       s <- summary(fit)
#       results[[cell]] <- data.frame(
#         CellType = cell,
#         HazardRatio = s$conf.int[1],
#         Lower_95 = s$conf.int[3],
#         Upper_95 = s$conf.int[4],
#         PValue = s$coefficients[5]
#       )
#     }
#   }
#   # 4. Consolidar y ajustar p-values por comparaciones múltiples (FDR)
#   final_table <- bind_rows(results) %>%
#     mutate(Adj_PValue = p.adjust(PValue, method = "fdr")) %>%
#     arrange(PValue)
#   
#   return(final_table)
# }




# cox_analysis <- function(deconv_df, clinical_df, time_col, event_col, sample_id) {
#   
#   # 1. Asegurar que deconv_df tenga el ID como columna
#   if (!(sample_id %in% colnames(deconv_df))) {
#     deconv_df <- deconv_df %>% as.data.frame() %>% rownames_to_column(var = sample_id)
#   }
#   
#   # 2. NORMALIZACIÓN AGRESIVA DE IDs
#   # Función interna para limpiar IDs: 
#   # - Cambia guiones por puntos
#   # - Elimina prefijos de texto seguidos de punto o guion bajo (ej. Brain_X o Neuroblastoma_)
#   clean_ids <- function(x) {
#     x <- gsub("-", ".", x) # Estandarizamos a puntos
#     x <- gsub("^[A-Za-z]+[._]+", "", x) # Borra "CualquierTexto_" o "CualquierTexto."
#     return(x)
#   }
#   
#   deconv_df[[sample_id]] <- clean_ids(deconv_df[[sample_id]])
#   clinical_df[[sample_id]] <- clean_ids(clinical_df[[sample_id]])
#   
#   # 3. UNIÓN
#   combined_df <- inner_join(deconv_df, clinical_df, by = sample_id)
#   
#   if (nrow(combined_df) == 0) {
#     stop("Error: No se pudieron alinear las muestras. Revisa los IDs manualmente.")
#   }
#   
#   # 4. Limpieza del Status (0:LIVING / 1:DECEASED)
#   if(is.character(combined_df[[event_col]])) {
#     combined_df[[event_col]] <- as.numeric(substr(combined_df[[event_col]], 1, 1))
#   }
#   
#   # 5. Bucle de Cox
#   cell_types <- setdiff(colnames(deconv_df), sample_id)
#   results <- list()
#   
#   for (cell in cell_types) {
#     formula_str <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ `", cell, "`"))
#     fit <- tryCatch({ coxph(formula_str, data = combined_df) }, error = function(e) NULL)
#     
#     if (!is.null(fit)) {
#       s <- summary(fit)
#       results[[cell]] <- data.frame(
#         CellType = cell,
#         HazardRatio = s$conf.int[1],
#         Lower_95 = s$conf.int[3],
#         Upper_95 = s$conf.int[4],
#         PValue = s$coefficients[5]
#       )
#     }
#   }
#   
#   if (length(results) == 0) return(NULL)
#   
#   bind_rows(results) %>%
#     mutate(Adj_PValue = p.adjust(PValue, method = "fdr")) %>%
#     arrange(PValue)
# }

# library(dplyr)
# 
# unify_sample_ids <- function(deconv_df, clinical_df, sample_id_col = "Sample.ID") {
#   
#   # 1. Preparar deconv_df (extraer rownames si es necesario)
#   if (!sample_id_col %in% colnames(deconv_df)) {
#     deconv_df <- as.data.frame(deconv_df) %>% 
#       tibble::rownames_to_column(var = sample_id_col)
#   }
#   
#   # 2. Función de limpieza universal
#   # Esta lógica elimina cualquier prefijo de texto que termine en _ o .
#   # y convierte todos los guiones en puntos para estandarizar.
#   clean_pattern <- function(ids) {
#     ids <- gsub("-", ".", ids)             # Todo a puntos
#     ids <- gsub("^[A-Za-z]+[_X.]+", "", ids) # Elimina "CualquierTexto_" o "CualquierTexto_X."
#     return(ids)
#   }
#   
#   # 3. Aplicar limpieza a ambos dataframes
#   deconv_df[[sample_id_col]]   <- clean_pattern(deconv_df[[sample_id_col]])
#   clinical_df[[sample_id_col]] <- clean_pattern(clinical_df[[sample_id_col]])
#   
#   # 4. Verificación de intersección
#   common <- intersect(deconv_df[[sample_id_col]], clinical_df[[sample_id_col]])
#   message(paste("Muestras unificadas con éxito:", length(common)))
#   
#   if(length(common) == 0) {
#     warning("¡Cuidado! No se encontraron coincidencias. Revisa el formato de los IDs.")
#   }
#   
#   return(list(deconv = deconv_df, clinical = clinical_df))
# }
 

cox_analysis <- function(deconv_df, clinical_df, time_col, event_col, sample_id = "Sample.ID") {
  
  # 1. Asegurar que deconv_df tenga el ID como columna (si viene con rownames)
  if (!(sample_id %in% colnames(deconv_df))) {
    deconv_df <- as.data.frame(deconv_df) %>% 
      rownames_to_column(var = sample_id)
  }
  
  # 2. LÓGICA DE UNIFICACIÓN UNIVERSAL
  # Función interna para limpiar IDs:
  # - Cambia guiones por puntos (estándar de R para nombres de columnas)
  # - Elimina prefijos como "Brain_X", "Neuroblastoma_", "Wilms_", "Rhabdoid_"
  clean_id <- function(x) {
    x <- gsub("^[A-Za-z]+_", "", x) # Quita "Neuroblastoma_", "Wilms_", etc.
    x <- gsub("^X", "", x)          # Quita la "X" que R añade a los números
    x <- gsub("-", ".", x)          # Convierte guiones en puntos
    return(x)
  }
  
  # Aplicamos la limpieza a ambos dataframes internamente
  deconv_clean <- deconv_df
  deconv_clean[[sample_id]] <- clean_id(deconv_clean[[sample_id]])
  
  clin_clean <- clinical_df
  clin_clean[[sample_id]] <- clean_id(clin_clean[[sample_id]])
  
  # 3. UNIÓN Y VERIFICACIÓN
  combined_df <- inner_join(deconv_clean, clin_clean, by = sample_id)
  
  n_match <- nrow(combined_df)
  if (n_match == 0) {
    stop("Error: No hay coincidencia de IDs. Revisa si el sample_id es correcto.")
  } else {
    message(paste("Éxito: Se han emparejado", n_match, "muestras."))
  }
  
  # 4. LIMPIEZA DE EVENTO (Soporta "0:LIVING", "1:DECEASED" o numérico)
  if(is.character(combined_df[[event_col]])) {
    combined_df[[event_col]] <- as.numeric(substr(combined_df[[event_col]], 1, 1))
  }
  
  # 5. EJECUCIÓN DEL MODELO POR TIPO CELULAR
  cell_types <- setdiff(colnames(deconv_df), sample_id)
  results <- list()
  
  for (cell in cell_types) {
    # Usamos backticks para manejar nombres de células con espacios o caracteres especiales
    formula_str <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ `", cell, "`"))
    
    fit <- tryCatch({
      coxph(formula_str, data = combined_df)
    }, error = function(e) return(NULL))
    
    if (!is.null(fit)) {
      s <- summary(fit)
      results[[cell]] <- data.frame(
        CellType = cell,
        HazardRatio = s$conf.int[1],
        Lower_95 = s$conf.int[3],
        Upper_95 = s$conf.int[4],
        PValue = s$coefficients[5]
      )
    }
  }
  
  # 6. RESULTADOS FINALES
  if (length(results) == 0) return(NULL)
  
  bind_rows(results) %>%
    mutate(Adj_PValue = p.adjust(PValue, method = "fdr")) %>%
    arrange(PValue)
}

plot_cox_forest <- function(cox_results, title = "Asociación de Tipos Celulares con Supervivencia") {
  
  # 1. Preparar datos: ordenar por HR y crear etiquetas de significancia
  plot_data <- cox_results %>%
    mutate(
      Significance = ifelse(Adj_PValue < 0.05, "Significativo (FDR < 0.05)", "No Significativo"),
      CellType = reorder(CellType, HazardRatio) # Ordenar para que sea legible
    )
  
  # 2. Construir el gráfico
  ggplot(plot_data, aes(x = HazardRatio, y = CellType)) +
    # Línea de referencia en HR = 1 (efecto nulo)
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
    # Intervalos de confianza
    geom_errorbarh(aes(xmin = Lower_95, xmax = Upper_95, color = Significance), height = 0.3) +
    # Punto del Hazard Ratio
    geom_point(aes(color = Significance), size = 3) +
    # Escala logarítmica (opcional pero recomendada para HR)
    scale_x_log10() +
    # Estética
    theme_minimal(base_size = 20) +
    labs(
      title = title,
      subtitle = "Hazard Ratio e Intervalos de Confianza (95%)",
      x = "Hazard Ratio (Escala Log)",
      y = "Tipo Celular",
      color = ""
    ) +
    theme(legend.position = "bottom")
}

# Supongamos que tus columnas se llaman así:
cox_brain <- cox_analysis(
  deconv_df = consensus_brain_preCombat,      # Tu df de deconvolución
  clinical_df = brain_clin,   # Tu df con supervivencia
  time_col = "OS.Months",       # Nombre de la columna de tiempo
  event_col = "OS.Status",    # Nombre de la columna de evento (0/1)
  sample_id = "Sample.ID"      # La columna común para unir ambos df
)
plot_cox_forest(cox_brain, title = "Cox en dataset de Brain")

cox_neuroblastoma <- cox_analysis(
  deconv_df = consensus_neuroblastoma_preCombat,      # Tu df de deconvolución
  clinical_df = neuroblastoma_clin,   # Tu df con supervivencia
  time_col = "Overall.Survival..Months.",       # Nombre de la columna de tiempo
  event_col = "Overall.Survival.Status",    # Nombre de la columna de evento (0/1)
  sample_id = "Sample.ID"      # La columna común para unir ambos df
)
plot_cox_forest(cox_neuroblastoma, title = "Cox en dataset de Neuroblastoma")

cox_wilms <- cox_analysis(
  deconv_df = consensus_wilms_preCombat,      # Tu df de deconvolución
  clinical_df = wilms_clin,   # Tu df con supervivencia
  time_col = "Overall.Survival..Months.",       # Nombre de la columna de tiempo
  event_col = "Overall.Survival.Status",    # Nombre de la columna de evento (0/1)
  sample_id = "Sample.ID"      # La columna común para unir ambos df
)
plot_cox_forest(cox_wilms, title = "Cox en dataset de Wilms")

cox_rhabdoid <- cox_analysis(
  deconv_df = consensus_rhabdoid_preCombat,      # Tu df de deconvolución
  clinical_df = rhabdoid_clin,   # Tu df con supervivencia
  time_col = "Overall.Survival..Months.",       # Nombre de la columna de tiempo
  event_col = "Overall.Survival.Status",    # Nombre de la columna de evento (0/1)
  sample_id = "Sample.ID"      # La columna común para unir ambos df
)
plot_cox_forest(cox_rhabdoid, title = "Cox en dataset de Rhabdoid")

# Ver los resultados más significativos
print(head(mi_tabla_supervivencia))





cox_pancancer_analysis <- function(list_of_datasets, sample_id = "Sample.ID") {
  
  # 1. Función interna para normalizar IDs (Puntos vs Guiones y Prefijos)
  clean_id <- function(x) {
    x <- gsub("^[A-Za-z]+_", "", x) # Quita "Neuroblastoma_", etc.
    x <- gsub("^X", "", x)          # Quita "X" inicial de R
    x <- gsub("-", ".", x)          # Estandariza a puntos
    return(x)
  }
  
  # 2. Procesamiento y Unión de Datasets
  # Asumimos que cada elemento de la lista tiene $dec (deconvolución) y $clin (clínica)
  combined_data <- lapply(names(list_of_datasets), function(name) {
    d <- list_of_datasets[[name]]
    
    # Asegurar que ID sea columna en deconv
    if (!(sample_id %in% colnames(d$dec))) {
      d$dec <- as.data.frame(d$dec) %>% rownames_to_column(var = sample_id)
    }
    
    # Limpiar IDs
    d$dec[[sample_id]] <- clean_id(d$dec[[sample_id]])
    d$clin[[sample_id]] <- clean_id(d$clin[[sample_id]])
    
    # Identificar columnas de supervivencia dinámicamente
    t_col <- grep("OS.*Month|Overall.*Month", colnames(d$clin), value = TRUE, ignore.case = TRUE)[1]
    e_col <- grep("OS.*Status|Overall.*Status", colnames(d$clin), value = TRUE, ignore.case = TRUE)[1]
    
    # Unir y añadir etiqueta de tumor
    inner_join(d$dec, d$clin, by = sample_id) %>%
      mutate(
        TumorType = name,
        OS_Months = as.numeric(.[[t_col]]),
        # Limpieza de evento (ej: "0:LIVING" -> 0)
        OS_Status = as.numeric(gsub(":[A-Za-z]+", "", .[[e_col]]))
      )
  }) %>% bind_rows()
  
  # 3. Ejecución del Modelo Cox Estratificado
  # (Compara pacientes dentro de su mismo tipo de tumor)
  cell_types <- setdiff(colnames(list_of_datasets[[1]]$dec), sample_id)
  results <- list()
  
  for (cell in cell_types) {
    # La clave Pan-Cancer: strata(TumorType)
    formula_str <- as.formula(paste0("Surv(OS_Months, OS_Status) ~ `", cell, "` + strata(TumorType)"))
    
    fit <- tryCatch({
      coxph(formula_str, data = combined_data)
    }, error = function(e) return(NULL))
    
    if (!is.null(fit)) {
      s <- summary(fit)
      results[[cell]] <- data.frame(
        CellType = cell,
        HazardRatio = s$conf.int[1],
        Lower_95 = s$conf.int[3],
        Upper_95 = s$conf.int[4],
        PValue = s$coefficients[5]
      )
    }
  }
  
  # 4. Formateo Final
  bind_rows(results) %>%
    mutate(Adj_PValue = p.adjust(PValue, method = "fdr")) %>%
    arrange(PValue)
}

mis_datos <- list(
  Neuroblastoma = list(dec = consensus_neuroblastoma_preCombat, clin = neuroblastoma_clin),
  Brain = list(dec = consensus_brain_preCombat, clin = brain_clin),
  Wilms = list(dec = consensus_wilms_preCombat, clin = wilms_clin),
  Rhabdoid = list(dec = consensus_rhabdoid_preCombat, clin = rhabdoid_clin)
)

# Ejecutar
mi_tabla_pancancer <- cox_pancancer_analysis(mis_datos)
plot_cox_forest(mi_tabla_pancancer, title = "Cox en dataset de All")



resumir_df(colnames(brain_clin), brain_clin, "TablaBrain.docx")
variables_limpias <- setdiff(colnames(neuroblastoma_clin), "Tumor.Disease.Anatomic.Site")
resumir_df(variables_limpias, neuroblastoma_clin, "TablaNeuroblastoma.docx")
resumir_df(colnames(wilms_clin), wilms_clin, "TablaWilms")
resumir_df(colnames(rhabdoid_clin), rhabdoid_clin, "TablaRhabdoid")


library(openxlsx)
write.xlsx(correlations_brain_consensus_preCombat, file = "CorrelationsBrain.xlsx")
write.xlsx(correlations_neuroblastoma_consensus_preCombat, file = "CorrelationsNeuroblastoma.xlsx")
write.xlsx(correlations_wilms_consensus_preCombat, file = "CorrelationsWilms.xlsx")
write.xlsx(correlations_rhabdoid_consensus_preCombat, file = "CorrelationsRhabdoid.xlsx")
