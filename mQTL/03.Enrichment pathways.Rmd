---
title: "Untitled"
author: "Leyre Franquet"
date: "2026-04-17"
output: html_document
---

```{r}
library(data.table)
library(plyr)
library(gridExtra)
library(grid)

results_cases <-  ldply(paste0("../results/osca_maf5/mqtl_description_stat_cases_chr",1:22,".txt.gz"), fread)  

results_controls <-  ldply(paste0("../results/osca_maf5/mqtl_description_stat_controls_chr",1:22,".txt.gz"), fread)  

# CASES
hetero_cases <- results_cases[
  (!results_cases$ProbeID %in% results_controls$ProbeID) &
    results_cases$P_het < 0.05 &
    results_cases$I2_het > 80,
]

# CONTROLS
hetero_controls <- results_controls[
  (!results_controls$ProbeID %in% results_cases$ProbeID) &
    results_controls$P_het < 0.05 &
    results_controls$I2_het > 80,
]

cat("Total filas hetero_cases:", nrow(hetero_cases), "\n")
cat("Total filas hetero_controls:", nrow(hetero_controls), "\n")


#CARGAR MANIFEST_450K
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- as.data.table(ann)

#total cpg analizadas
all_cpg <- unique(c(results_cases$ProbeID, results_controls$ProbeID))
cat("Número total de CpGs analizadas:", length(all_cpg), "\n")

#crear lista de CpGs
sig_cpg_cases <- unique(hetero_cases$ProbeID)
sig_cpg_controls <- unique(hetero_controls$ProbeID)

mapped_cases <- getMappedEntrezIDs(
  sig.cpg = sig_cpg_cases,
  all.cpg = all_cpg,
  array.type = "450K",
  anno = NULL,
  genomic.features = "ALL"
)

mapped_genes <- getMappedEntrezIDs(
  sig.cpg = all_cpg,
  all.cpg = all_cpg,
  array.type = "450K",
  anno = NULL,
  genomic.features = "ALL"
)

#filtrar la anotación para mis cpg
res_cases <- ann[ann$Name %in% sig_cpg_cases, ]
res_controls <- ann[ann$Name %in% sig_cpg_controls, ]

#dejar solo las columnas que necesito
res_cases_small <- res_cases[, c(
  "Name",
  "chr",
  "pos",
  "UCSC_RefGene_Name",
  "UCSC_RefGene_Group",
  "Relation_to_Island"
)]

res_controls_small <- res_controls[, c(
  "Name",
  "chr",
  "pos",
  "UCSC_RefGene_Name",
  "UCSC_RefGene_Group",
  "Relation_to_Island"
)]

res_cases_small <- as.data.table(res_cases_small)
res_controls_small <- as.data.table(res_controls_small)

setDT(hetero_cases)
setDT(hetero_controls)

hetero_cases_annot <- merge(
  hetero_cases,
  res_cases_small,
  by.x = "ProbeID",
  by.y = "Name",
  all.x = TRUE
)

hetero_controls_annot <- merge(
  hetero_controls,
  res_controls_small,
  by.x = "ProbeID",
  by.y = "Name",
  all.x = TRUE
)
```

#pathway enrichment
```{r}
BiocManager::install("missMethyl")
```

```{r}
library(missMethyl)
library(data.table)
library(dplyr)
library(ggplot2)


mapped_cases <- getMappedEntrezIDs(
  sig.cpg = sig_cpg_cases,
  all.cpg = all_cpg,
  array.type = "450K",
  anno = NULL,
  genomic.features = "ALL"
)

mapped_controls <- getMappedEntrezIDs(
  sig.cpg = sig_cpg_controls,
  all.cpg = all_cpg,
  array.type = "450K",
  anno = NULL,
  genomic.features = "ALL"
)

cat("\nElementos en mapped_cases:\n")
print(names(mapped_cases))

cat("\nElementos en mapped_controls:\n")
print(names(mapped_controls))
```

#ENRIQUECIMIENTO FUNCIONAL CON GO
```{r}
go_cases <- gometh(
  sig.cpg = sig_cpg_cases,
  all.cpg = all_cpg,
  collection = "GO",
  array.type = "450K",
  prior.prob = TRUE,
  anno = NULL
)

go_controls <- gometh(
  sig.cpg = sig_cpg_controls,
  all.cpg = all_cpg,
  collection = "GO",
  array.type = "450K",
  prior.prob = TRUE,
  anno = NULL
)

go_cases_top10 <- go_cases %>%
  arrange(FDR, P.DE) %>%
  slice_head(n = 10)

go_controls_top10 <- go_controls %>%
  arrange(FDR, P.DE) %>%
  slice_head(n = 10)
```

#over representation analysis (ORA)
```{r}
kegg_cases <- gometh(
  sig.cpg = sig_cpg_cases,
  all.cpg = all_cpg,
  collection = "KEGG",
  array.type = "450K",
  prior.prob = TRUE,
  anno = NULL
)

kegg_controls <- gometh(
  sig.cpg = sig_cpg_controls,
  all.cpg = all_cpg,
  collection = "KEGG",
  array.type = "450K",
  prior.prob = TRUE,
  anno = NULL
)

kegg_cases_top10 <- kegg_cases %>%
  arrange(FDR, P.DE) %>%
  slice_head(n = 10)

kegg_controls_top10 <- kegg_controls %>%
  arrange(FDR, P.DE) %>%
  slice_head(n = 10)
```

#REACTOME
```{r}
library(ReactomePA)
library(clusterProfiler)

genes_cases <- mapped_cases$sig.eg
genes_controls <- mapped_controls$sig.eg

# limpiar formato por si vienen raros
genes_cases <- unique(as.character(genes_cases))
genes_controls <- unique(as.character(genes_controls))

genes_cases <- genes_cases[!is.na(genes_cases) & genes_cases != ""]
genes_controls <- genes_controls[!is.na(genes_controls) & genes_controls != ""]

cat("\nGenes Entrez en casos:", length(genes_cases), "\n")
cat("Genes Entrez en controles:", length(genes_controls), "\n")

cat("\nPrimeros genes casos:\n")
print(head(genes_cases))

cat("\nPrimeros genes controles:\n")
print(head(genes_controls))

if (length(genes_cases) == 0) stop("No hay genes Entrez válidos en casos")
if (length(genes_controls) == 0) stop("No hay genes Entrez válidos en controles")

reactome_cases <- enrichPathway(
  gene = genes_cases,
  organism = "human",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  universe = unique(as.character(mapped_genes$universe)),
  readable = TRUE
)

reactome_controls <- enrichPathway(
  gene = genes_controls,
  organism = "human",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  universe = unique(as.character(mapped_genes$universe)),
  readable = TRUE
)

reactome_cases_top10 <- reactome_cases@result %>%
  arrange(p.adjust, pvalue) %>%
  slice_head(n = 10)

reactome_controls_top10 <- reactome_controls@result %>%
  arrange(p.adjust, pvalue) %>%
  slice_head(n = 10)

```

#pathway analysisi plot
```{r}
library(dplyr)
library(ggplot2)
library(stringr)

go_cases_plot <- go_cases_top10 %>%
  transmute(
    Database = "GO",
    Group = "Cases",
    Description = TERM,
    score = -log10(P.DE)
  )

go_controls_plot <- go_controls_top10 %>%
  transmute(
    Database = "GO",
    Group = "Controls",
    Description = TERM,
    score = -log10(P.DE)
  )

kegg_cases_plot <- kegg_cases_top10 %>%
  transmute(
    Database = "KEGG",
    Group = "Cases",
    Description = Description,
    score = -log10(P.DE)
  )

kegg_controls_plot <- kegg_controls_top10 %>%
  transmute(
    Database = "KEGG",
    Group = "Controls",
    Description = Description,
    score = -log10(P.DE)
  )

reactome_cases_plot <- reactome_cases_top10 %>%
  transmute(
    Database = "Reactome",
    Group = "Cases",
    Description = Description,
    score = -log10(pvalue)
  )

reactome_controls_plot <- reactome_controls_top10 %>%
  transmute(
    Database = "Reactome",
    Group = "Controls",
    Description = Description,
    score = -log10(pvalue)
  )

#unir
cases_pathways <- bind_rows(
  go_cases_plot,
  kegg_cases_plot,
  reactome_cases_plot
)

controls_pathways <- bind_rows(
  go_controls_plot,
  kegg_controls_plot,
  reactome_controls_plot
)

cases_pathways <- cases_pathways %>%
  mutate(Description = str_wrap(Description, width = 40))

controls_pathways <- controls_pathways %>%
  mutate(Description = str_wrap(Description, width = 40))

## crear plots
#casos
plot_cases <- ggplot(cases_pathways,
                     aes(x = score, y = reorder(Description, score), fill = Database)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~Database, scales = "free_y", ncol = 1) +
  labs(
    title = "Top pathways enriched in stroke cases (ranked by nominal p-value)",
    x = expression(-log[10]("nominal p-value")),
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )
plot_cases

#controles
plot_controls <- ggplot(controls_pathways,
                        aes(x = score, y = reorder(Description, score), fill = Database)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~Database, scales = "free_y", ncol = 1) +
  labs(
    title = "Top pathways enriched in controls (ranked by nominal p-value)",
    x = expression(-log[10]("nominal p-value")),
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )
plot_controls

#guardar plots
ggsave("pathway_plot_cases_nominalP.png", plot_cases, width = 10, height = 9, dpi = 300)
ggsave("pathway_plot_controls_nominalP.png", plot_controls, width = 10, height = 9, dpi = 300)
```


```{r}

```

#listas de genes para WebStalt
```{r}
genes_all <- unique(mapped_genes$gene)

```


distribución por región génica
```{r}
library(data.table)
library(ggplot2)

#distribución de cpg por regiçon génica
setDT(hetero_cases_annot)
setDT(hetero_controls_annot)
setDT(ann)

#una misma cpg puede estar varias veces en SNP distintos por lo que solo quieor contar cpg y no asociaciones
cases_unique_regions <- unique(
  hetero_cases_annot[, .(ProbeID, UCSC_RefGene_Group)]
)

controls_unique_regions <- unique(
  hetero_controls_annot[, .(ProbeID, UCSC_RefGene_Group)]
)

cat("CpGs únicas en casos para anotación funcional:",
    nrow(cases_unique_regions), "\n")
cat("CpGs únicas en controles para anotación funcional:",
    nrow(controls_unique_regions), "\n")

# classify_region <- function(x){
# 
#   if(is.na(x)) return("Other")
# 
#   if(grepl("TSS200|TSS1500|1stExon|5'UTR", x)) return("Promoter")
# 
#   if(grepl("Body", x)) return("Body")
# 
#   if(grepl("3'UTR", x)) return("3UTR")
# 
#   return("Other")
# }
# 
# cases_unique_regions[, Region_grouped := sapply(UCSC_RefGene_Group, classify_region)]
# controls_unique_regions[, Region_grouped := sapply(UCSC_RefGene_Group, classify_region)]

#group_counts_cases <- cases_unique_regions[, .N, by = Region_grouped]
#group_counts_controls <- controls_unique_regions[, .N, by = Region_grouped]

#group_counts_cases[, Percentage := round(100 * N / sum(N),2)]
#group_counts_controls[, Percentage := round(100 * N / sum(N),2)]

# convertir NA o vacío en Intergenic
cases_unique_regions[
  is.na(UCSC_RefGene_Group) | UCSC_RefGene_Group == "",
  UCSC_RefGene_Group := "Intergenic"
]

controls_unique_regions[
  is.na(UCSC_RefGene_Group) | UCSC_RefGene_Group == "",
  UCSC_RefGene_Group := "Intergenic"
]

# separar anotaciones múltiples
cases_regions_long <- cases_unique_regions[
  ,
  .(Region = unlist(strsplit(UCSC_RefGene_Group, ";"))),
  by = ProbeID
]

controls_regions_long <- controls_unique_regions[
  ,
  .(Region = unlist(strsplit(UCSC_RefGene_Group, ";"))),
  by = ProbeID
]


# Quitar duplicados por seguridad
cases_regions_long <- unique(cases_regions_long)
controls_regions_long <- unique(controls_regions_long)

#análisis de todas las cpg
all_cpg_annot <- ann[Name %in% all_cpg, .(Name, UCSC_RefGene_Group)]

all_cpg_annot[
  is.na(UCSC_RefGene_Group) | UCSC_RefGene_Group == "",
  UCSC_RefGene_Group := "Intergenic"
]

all_regions_long <- all_cpg_annot[
  ,
  .(Region = unlist(strsplit(UCSC_RefGene_Group, ";"))),
  by = Name
]

all_regions_long <- unique(all_regions_long)

# contar
region_counts_all <- all_regions_long[, .N, by = Region][order(-N)]
region_counts_cases <- cases_regions_long[, .N, by = Region][order(-N)]
region_counts_controls <- controls_regions_long[, .N, by = Region][order(-N)]

# porcentajes
region_counts_all[, Percentage := round(100 * N / sum(N), 2)]
region_counts_cases[, Percentage := round(100 * N / sum(N), 2)]
region_counts_controls[, Percentage := round(100 * N / sum(N), 2)]

cat("\nRegiones detalladas en casos:\n")
print(region_counts_cases)

cat("\nRegiones detalladas en controles:\n")
print(region_counts_controls)

# guardar tablas
fwrite(region_counts_cases, "region_counts_cases_detailed.csv")
fwrite(region_counts_controls, "region_counts_controls_detailed.csv")

# tabla resumen que integla las anteriores
table_all <- region_counts_all[, .(Region, All = paste0(N, " (", Percentage, "%)"))]
table_cases <- region_counts_cases[, .(Region, Cases = paste0(N, " (", Percentage, "%)"))]
table_controls <- region_counts_controls[, .(Region, Controls = paste0(N, " (", Percentage, "%)"))]

final_table <- Reduce(function(x, y) merge(x, y, by = "Region", all = TRUE),
                      list(table_all, table_cases, table_controls))

print(final_table)
fwrite(final_table, "functional_regions_comparison.csv")

# gráfico de barras
region_counts_all[, Group := "All CpGs"]
region_counts_cases[, Group := "Cases hetero"]
region_counts_controls[, Group := "Controls hetero"]

plot_data <- rbindlist(list(
  region_counts_all,
  region_counts_cases,
  region_counts_controls
))

plot_data <- plot_data[, .(Region, Percentage, Group)]

plot_data[, Region := factor(
  Region,
  levels = c("Body", "TSS1500", "5'UTR", "TSS200", "1stExon", "3'UTR", "Intergenic")
)]

p <- ggplot(plot_data, aes(x = Region, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab("Percentage of CpGs") +
  xlab("Genomic region") +
  scale_fill_manual(values = c("grey40", "#1b9e77", "#d95f02")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

print(p)

ggsave(
  "CpG_genomic_distribution.png",
  p,
  width = 8,
  height = 5,
  dpi = 300
)

# stacked plot complementario
plot_data_stacked <- rbindlist(list(
  region_counts_all,
  region_counts_cases,
  region_counts_controls
))

plot_data_stacked <- plot_data_stacked[, .(Group, Region, Percentage)]

plot_data_stacked[, Region := factor(
  Region,
  levels = c("Body", "TSS1500", "5'UTR", "TSS200", "1stExon", "3'UTR", "Intergenic")
)]

plot_data_stacked[, Group := factor(
  Group,
  levels = c("All CpGs", "Cases hetero", "Controls hetero")
)]

p_stacked <- ggplot(plot_data_stacked, aes(x = Group, y = Percentage, fill = Region)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_classic() +
  xlab("") +
  ylab("Percentage of CpGs") +
  theme(
    axis.text.x = element_text(size = 12),
    legend.title = element_blank()
  )

print(p_stacked)

ggsave(
  "CpG_genomic_distribution_stacked.png",
  p_stacked,
  width = 7,
  height = 5,
  dpi = 300
)

```

discusión de genes
```{r}
library(data.table)

setDT(hetero_cases_annot)
setDT(hetero_controls_annot)

# Función para limpiar genes anotados
extract_genes <- function(dt) {
  genes <- dt[!is.na(UCSC_RefGene_Name) & UCSC_RefGene_Name != "",
              .(ProbeID, UCSC_RefGene_Name)]
  
  genes[, Gene := strsplit(UCSC_RefGene_Name, ";")]
  genes <- genes[, .(Gene = unlist(Gene)), by = ProbeID]
  genes[, Gene := trimws(Gene)]
  genes <- genes[Gene != ""]
  
  return(unique(genes))
}

genes_cases_dt <- extract_genes(hetero_cases_annot)
genes_controls_dt <- extract_genes(hetero_controls_annot)

# Genes más repetidos según número de CpGs asociadas
top_genes_cases <- genes_cases_dt[, .(
  n_CpGs = uniqueN(ProbeID)
), by = Gene][order(-n_CpGs)]

top_genes_controls <- genes_controls_dt[, .(
  n_CpGs = uniqueN(ProbeID)
), by = Gene][order(-n_CpGs)]

# Ver top 20
head(top_genes_cases, 20)
head(top_genes_controls, 20)

# Genes presentes solo en casos o solo en controles
genes_cases_unique <- setdiff(top_genes_cases$Gene, top_genes_controls$Gene)
genes_controls_unique <- setdiff(top_genes_controls$Gene, top_genes_cases$Gene)

top_unique_cases <- top_genes_cases[Gene %in% genes_cases_unique][order(-n_CpGs)]
top_unique_controls <- top_genes_controls[Gene %in% genes_controls_unique][order(-n_CpGs)]

head(top_unique_cases, 20)
head(top_unique_controls, 20)

# Guardar resultados
fwrite(top_genes_cases, "top_genes_cases.csv")
fwrite(top_genes_controls, "top_genes_controls.csv")
fwrite(top_unique_cases, "top_unique_genes_cases.csv")
fwrite(top_unique_controls, "top_unique_genes_controls.csv")

##tabla para anexo
library(data.table)

setDT(hetero_cases_annot)
setDT(hetero_controls_annot)

# Función para extraer genes
extract_genes <- function(dt) {
  genes <- dt[!is.na(UCSC_RefGene_Name) & UCSC_RefGene_Name != "",
              .(ProbeID, UCSC_RefGene_Name)]
  
  genes[, Gene := strsplit(UCSC_RefGene_Name, ";")]
  genes <- genes[, .(Gene = unlist(Gene)), by = ProbeID]
  genes[, Gene := trimws(Gene)]
  genes <- genes[Gene != ""]
  
  return(unique(genes))
}

# Extraer genes
genes_cases_dt <- extract_genes(hetero_cases_annot)
genes_controls_dt <- extract_genes(hetero_controls_annot)

# Contar número de CpGs por gen y ordenar
top_genes_cases <- genes_cases_dt[, .(
  n_CpGs = uniqueN(ProbeID)
), by = Gene][order(-n_CpGs)][1:20]

top_genes_controls <- genes_controls_dt[, .(
  n_CpGs = uniqueN(ProbeID)
), by = Gene][order(-n_CpGs)][1:20]

# Renombrar columnas para el anexo
setnames(top_genes_cases,
         c("Gene", "n_CpGs"),
         c("Gen", "Número de CpGs asociadas"))

setnames(top_genes_controls,
         c("Gene", "n_CpGs"),
         c("Gen", "Número de CpGs asociadas"))

# Exportar
fwrite(top_genes_cases, "top20_genes_casos.csv", sep = ";")
fwrite(top_genes_controls, "top20_genes_controles.csv", sep = ";")
```

