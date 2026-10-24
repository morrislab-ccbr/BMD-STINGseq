# Supplementary/revision figures:
#   - hFOB ATAC-seq sample PCA (nf-core/atacseq DESeq2 PCA)
#   - expression of osteoblast lineage markers in human bone marrow scRNA-seq (GSE317068)
#     and across hFOB differentiation (Supplementary Figure 1)
#   - hMSC replication: hFOB vs hMSC effects and per-target hMSC RT-qPCR results

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
external_dir <- Sys.getenv("EXTERNAL_DIR", unset = "path/to/external_inputs")
setwd(project_dir)
system("mkdir -p revisions")

library(data.table)
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggrepel)
library(ggrastr)
library(colorspace)
library(ggpubr)
library(ggtext)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ==========================================================================
# ATAC-seq sample PCA with k-means clusters
# ==========================================================================
DF.0 = fread(paste0(external_dir, "/hFOB_ATAC_nfcore/results/multiqc/narrow_peak/multiqc_data/multiqc_mlib_deseq2_pca_1-plot.txt"))

DF.0$Label <- gsub("^DAY([0-9]+)_REP1$", "Day \\1, Rep A", DF.0$Sample)
DF.0$Label <- gsub("^DAY([0-9]+)_REP2$", "Day \\1, Rep B", DF.0$Label)
DF.0$Label <- gsub("^DAY([0-9]+)_REP3$", "Day \\1, Rep C", DF.0$Label)
DF.0$Label <- gsub("^DAY([0-9]+)_REP4$", "Day \\1, Rep D", DF.0$Label)

set.seed(123)
km <- kmeans(DF.0[, c("x", "y")], centers = 3)

DF.0$Cluster <- factor(km$cluster)

ggout = ggplot(DF.0, aes(x = x, y = y, colour = Cluster)) +
  theme_classic(base_size = 6) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = cbPalette) +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = Inf, color = "black", min.segment.length = 0, segment.size = 0.2) +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black", linewidth = 0.1),
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        axis.text = element_text(color = "black"))

ggsave(filename = "revisions/260721_ATAC-QC-PCA.pdf", height = 4, width = 4, device = "pdf", plot = ggout)

# ==========================================================================
# Osteoblast lineage markers (up to 1,000 cells per cluster)
# ==========================================================================
H.0 = readRDS("revisions/Chai2026/GSE317068_Human_bone_scRNAseq_seurat_object_20251120.Rds")
H_cell_data_subset <- H.0@meta.data %>%
  dplyr::mutate(cell = colnames(H.0)) %>%
  dplyr::group_by(Cluster_annotation) %>%
  dplyr::filter(row_number() <= 1000) %>%
  dplyr::ungroup()
H_subset_cells = colnames(H.0) %in% H_cell_data_subset$cell
H.1 = H.0[,H_subset_cells]

correct_order <- c(
  # Mesenchymal / stromal
  "0. MSCs CXCL12+GALNT17+",
  "1. MSCs CXCL12+STMN2+",
  "6. Pericytes",
  "7. SMCs",
  "2. ECs",
  
  # Skeletal lineage
  "15. Chondrocytes and fibroblasts",
  "3. Pre-osteoblasts",
  "12. Mature osteoblasts",
  
  # Hematopoietic stem/progenitors
  "13. HSCs",
  "10. GMPs DEFA1+",
  "14. GMPs CDC25C+",
  
  # Neutrophil lineage
  "4. Immature neutrophils",
  "8. Mature neutrophils IL1R2+",
  "17. Mature neutrophils LST1+",
  
  # Mononuclear myeloid
  "11. Monocytes and macrophages",
  
  # B-cell lineage
  "16. Pre B cells",
  "9. Immature B cells",
  "19. B cells",
  
  # T/NK lineage
  "5. T and NK cells",
  # Other
  "18. Multiplets"
)

group_labels <- c(  "Mesenchymal/stromal",
                    "Mesenchymal/stromal",
                    "Mesenchymal/stromal",
                    "Mesenchymal/stromal",
                    "Mesenchymal/stromal",
                    "Skeletal lineage",
                    "Skeletal lineage",
                    "Skeletal lineage",
                    "Hematopoietic progenitors",
                    "Hematopoietic progenitors",
                    "Hematopoietic progenitors",
                    "Neutrophils",
                    "Neutrophils",
                    "Neutrophils",
                    "Myeloid",
                    "B cells",
                    "B cells",
                    "B cells",
                    "T/NK cells",
                    "Other")

celltype_groups <- data.frame(
  Cluster_annotation = correct_order,
  group = group_labels)

celltype_groups$cell_type <- sub("^[0-9]+\\.\\s*", "", celltype_groups$Cluster_annotation)

tmp_meta.data = H.1@meta.data %>%
  dplyr::left_join(celltype_groups, by = "Cluster_annotation")
H.1@meta.data$cell_type = tmp_meta.data$cell_type
H.1@meta.data$cell_group = tmp_meta.data$group
H.1@meta.data$cell_type <- factor(x = H.1@meta.data$cell_type, levels = celltype_groups$cell_type)
H.1@meta.data$cell_group <- factor(x = H.1@meta.data$cell_group, levels = unique(celltype_groups$group))
celltype_groups$cell_type <- factor(x = celltype_groups$cell_type, levels = celltype_groups$cell_type)

human_genes <- c(
  "LEPR", "NES", # Bone marrow MSC
  "CD34", "CXCL12", "KITL", "ANGPT1", # Osteo-progenitor
  "RUNX2", "SP7", "DLX5", "SATB2", "rna_ALPL", "COL1A1", "COL1A2", # Early osteoblast
  "BGLAP", "BGLAP2", "IBSP", "SPARC", "SPP1", "PDPN", "DMP1", "PHEX", "IFITM5", # Mature osteoblast
  "CAPG", "MEPE", "SOST", "FGF23", "HYOU1" # Osteocyte
)

human_genes_found = human_genes[human_genes %in% rownames(H.1)]

H.1_plotdata = do.call(rbind, lapply(human_genes_found, function(x){
  tmp1 = VlnPlot(object = H.1, features = x, assay = "RNA", pt.size = 0.4, alpha = 0.1, group.by = "cell_type")
  tmp2 = tmp1@data
  colnames(tmp2) = c("gene", "cell_type")
  tmp2$gene_id = x
  return(tmp2)
})) %>% as.data.table() %>%
  dplyr::filter(cell_type != "Multiplets") %>%
  dplyr::left_join(celltype_groups, by = "cell_type") %>%
  dplyr::mutate(gene_noise = gene + rnorm(n(), sd = 1/100000))
H.1_plotdata$group <- factor(x = H.1_plotdata$group, levels = unique(celltype_groups$group))
H.1_plotdata$gene_id  = factor(x = H.1_plotdata$gene_id, levels = human_genes_found)

H.1_plotdata_out = ggplot(H.1_plotdata, aes(x = cell_type, y = gene, color = group)) + 
  theme_classic() +
  scale_color_manual(values = cbPalette) +
  rasterise(geom_jitter(size = 0.4, alpha = 0.05), dpi = 300) +
  geom_violin(trim = TRUE, scale = "width", linewidth = 0.1, color = "black", fill = NA) +
  facet_wrap(. ~ gene_id, scales = "free_y", ncol = 4, nrow = 6) +
  theme(axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        text = element_text(size = 6),
        strip.text = element_text(size = 6)) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Cell group") + xlab("Cell type") + ylab("Normalized expression")

rnaseq.0 = fread("revisions/gene_expression_data.csv") %>%
  dplyr::mutate(experiment = ifelse(experiment == "GSE293846 differentiation", "short-read RNA-seq", experiment),
                experiment = ifelse(experiment == "hFOB differentiation", "long-read RNA-seq", experiment))
rnaseq.0$gene <- factor(rnaseq.0$gene, levels = human_genes_found)
day_colors <- c(
  "0"  = "dodgerblue1",
  "2"  = "#b2df8a",
  "4"  = "#33a02c",
  "6"  = darken("#33a02c", amount = 0.15),
  "8"  = darken("#33a02c", amount = 0.30),
  "10" = darken("#33a02c", amount = 0.45)
)

gg1 = ggplot(subset(rnaseq.0, experiment == "short-read RNA-seq"), aes(x = factor(day), y = expression_count)) +
  theme_classic() +
  stat_summary(fun = mean, geom = "col",
               aes(fill = factor(day)), colour = "black", width = 0.7, alpha = 0.4) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.25, linewidth = 0.4,
               aes(ymin = after_stat(pmax(ymin, 0)))) +
  geom_jitter(width = 0.15, height = 0,
              size = 1, shape = 19, color = "black", alpha = 0.8) +
  facet_wrap(. ~ gene, scales = "free_y", nrow = 4) +
  scale_fill_manual(values = day_colors) +
  scale_color_manual(values = day_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  expand_limits(y = 0) +
  labs(x = "Short-read RNA-seq day", y = "Expression count") +
  theme(axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        text = element_text(size = 6),
        strip.text = element_text(size = 6),
        legend.position = "none")

day_colors <- c(
  "0"  = "dodgerblue1",
  "2"  = "#b2df8a",
  "4"  = "#33a02c",
  "10" = darken("#33a02c", amount = 0.45)
)

gg2 = ggplot(subset(rnaseq.0, experiment == "long-read RNA-seq"), aes(x = factor(day), y = expression_count)) +
  theme_classic() +
  stat_summary(fun = mean, geom = "col",
               aes(fill = factor(day)), colour = "black", width = 0.7, alpha = 0.4) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.25, linewidth = 0.4,
               aes(ymin = after_stat(pmax(ymin, 0)))) +
  geom_jitter(width = 0.15, height = 0,
              size = 1, shape = 19, color = "black", alpha = 0.8) +
  facet_wrap(. ~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = day_colors) +
  scale_color_manual(values = day_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  expand_limits(y = 0) +
  labs(x = "Long-read RNA-seq day", y = "Expression count") +
  theme(axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        text = element_text(size = 6),
        strip.text = element_text(size = 6),
        legend.position = "none")

plot_out = H.1_plotdata_out / gg1 / gg2 + plot_layout(heights = c(1.5, 1, 1))

ggsave(plot_out, filename = "revisions/260916_Human_GeneMarkers.pdf", device = "pdf", height = 11, width = 8.5, units = "in")

# ==========================================================================
# hMSC replication of CRE and TSS targeting (RT-qPCR, log2 fold-change)
# ==========================================================================
# Colour palette shared by both panels (Okabe-Ito)
locus_levels <- c("LINC02454","CXCL12","SMARCD3","TLN1","TRIP11","DAP3_YY1AP1")
locus_colors <- c(
  "LINC02454"    = "#0072B2",  # blue
  "CXCL12"       = "#E69F00",  # orange
  "SMARCD3"      = "#009E73",  # bluish green
  "TLN1"         = "#D55E00",  # vermillion
  "TRIP11"       = "#CC79A7",  # reddish purple
  "DAP3_YY1AP1"  = "#56B4E9"   # sky blue
)
locus_labels <- c(
  "LINC02454"    = "*LINC02454*",
  "CXCL12"       = "*CXCL12*",
  "SMARCD3"      = "*SMARCD3*",
  "TLN1"         = "*TLN1*",
  "TRIP11"       = "*TRIP11*",
  "DAP3_YY1AP1"  = "*DAP3*/*YY1AP1*"
)

# Axis label with a subscripted "2" in "log2"
log2_axis_lab <- function(prefix) {
  bquote(.(prefix) ~ "fold-change (" * log[2] * ")")
}

# ==========================================================================
# Panel A: hFOB vs hMSC effects
# ==========================================================================
cor_df <- tibble::tribble(
  ~Target,             ~Guide,      ~Locus,          ~hFOB_log2FC, ~hMSC_log2FC,
  "SMARCD3_E2-3",      "SNP151",    "SMARCD3",       -0.104697,     0.781952,
  "TRIP11_E7-9",       "TSS1",      "TRIP11",        -5.058894,    -2.362001,
  "LINC02454_E2-3",    "Upstream",  "LINC02454",     -6.058894,    -2.909396,
  "CXCL12_E2-3",       "SNP181",    "CXCL12",        -3.058894,    -0.969618,
  "TLN1_6479-6683",    "SNP165",    "TLN1",          -1.058894,    -0.915622,
  "DAP3_E9-11",        "TSS1",      "DAP3_YY1AP1",   -3.836501,    -1.140653,
  "SMARCD3_E1-3",      "SNP151",    "SMARCD3",       -4.836501,    -2.320304,
  "LINC02454_E2-3",    "SNP233",    "LINC02454",     -5.321928,    -1.831379,
  "HMGA2_E1_2-3",      "Upstream",  "LINC02454",     -0.152003,    -0.366269,
  "CXCL12_E2-E4",      "SNP181",    "CXCL12",        -3.473931,    -1.090091,
  "SMARCD3_994-1187",  "SNP151",    "SMARCD3",       -2.556393,    -1.103582,
  "YY1AP1_E9-11",      "SNP31",     "DAP3_YY1AP1",   -2.321928,    -0.523825,
  "DAP3_E9-11",        "SNP31",     "DAP3_YY1AP1",   -4.643856,    -0.204806,
  "YY1AP1_E9-11",      "TSS1",      "DAP3_YY1AP1",   -1.321928,    -0.351696,
  "HMGA2_E1_2-3",      "SNP233",    "LINC02454",     -0.074001,     0.070551,
  "CXCL12_E2-3",       "SNP179",    "CXCL12",        -1.251539,    -0.125624,
  "CXCL12_E2-E4",      "SNP179",    "CXCL12",        -1.643856,     0.011139
) %>%
  mutate(Locus = factor(Locus, levels = locus_levels))

gg_cor <- ggplot(cor_df, aes(x = hFOB_log2FC, y = hMSC_log2FC)) +
  theme_classic() +
  geom_hline(yintercept = 0, linewidth = 0.1, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.1, color = "black") +
  geom_point(aes(color = Locus), alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", formula = 'y ~ x', se = FALSE, color = "firebrick1", linewidth = 0.2) +
  stat_cor(method = "pearson", size = 2, color = "black") +
  scale_color_manual(values = locus_colors, drop = FALSE, labels = locus_labels) +
  labs(x = log2_axis_lab("hFOB"), y = log2_axis_lab("hMSC"), color = "Locus") +
  theme(
    axis.line       = element_line(linewidth = 0.1, color = "black"),
    axis.ticks      = element_line(linewidth = 0.1, color = "black"),
    text            = element_text(color = "black", size = 6),
    axis.text       = element_text(color = "black", size = 6),
    legend.position    = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background  = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.key.size    = unit(0.25, "cm"),
    legend.title       = element_text(size = 6),
    legend.text        = ggtext::element_markdown(size = 5.5),
    legend.margin      = margin(1, 2, 1, 2)
  )

# ==========================================================================
# Panel B: hMSC effects per gRNA and amplicon
# ==========================================================================
raw <- tribble(
  ~Target,             ~Guide,      ~log2FC,
  "CXCL12_E2-3",       "SNP179",     1.040158,
  "CXCL12_E2-3",       "SNP179",    -0.640301,
  "CXCL12_E2-3",       "SNP179",    -0.328024,
  "CXCL12_E2-3",       "SNP179",    -0.574329,
  "CXCL12_E2-3",       "SNP181",    -1.068612,
  "CXCL12_E2-3",       "SNP181",    -1.198834,
  "CXCL12_E2-3",       "SNP181",    -0.290255,
  "CXCL12_E2-3",       "SNP181",    -1.320769,
  "CXCL12_E2-E4",      "SNP179",     0.570646,
  "CXCL12_E2-E4",      "SNP179",    -0.074353,
  "CXCL12_E2-E4",      "SNP179",    -0.549533,
  "CXCL12_E2-E4",      "SNP179",     0.097798,
  "CXCL12_E2-E4",      "SNP181",    -2.279980,
  "CXCL12_E2-E4",      "SNP181",    -1.148863,
  "CXCL12_E2-E4",      "SNP181",    -0.234275,
  "CXCL12_E2-E4",      "SNP181",    -0.697245,
  "DAP3_E9-11",        "SNP31",     -0.525956,
  "DAP3_E9-11",        "SNP31",     -0.176106,
  "DAP3_E9-11",        "SNP31",      0.087644,
  "DAP3_E9-11",        "TSS1",      -1.617313,
  "DAP3_E9-11",        "TSS1",      -1.159326,
  "DAP3_E9-11",        "TSS1",      -0.645319,
  "HMGA2_E1_2-3",      "SNP233",    -0.016014,
  "HMGA2_E1_2-3",      "SNP233",     0.266283,
  "HMGA2_E1_2-3",      "SNP233",     0.269976,
  "HMGA2_E1_2-3",      "SNP233",    -0.238041,
  "HMGA2_E1_2-3",      "Upstream",  -0.149878,
  "HMGA2_E1_2-3",      "Upstream",  -0.619124,
  "HMGA2_E1_2-3",      "Upstream",  -0.108424,
  "HMGA2_E1_2-3",      "Upstream",  -0.587651,
  "LINC02454_E2-3",    "SNP233",    -0.695425,
  "LINC02454_E2-3",    "SNP233",    -1.721551,
  "LINC02454_E2-3",    "SNP233",    -1.216538,
  "LINC02454_E2-3",    "SNP233",    -3.692002,
  "LINC02454_E2-3",    "Upstream",  -3.209975,
  "LINC02454_E2-3",    "Upstream",  -2.954935,
  "LINC02454_E2-3",    "Upstream",  -1.878983,
  "LINC02454_E2-3",    "Upstream",  -3.593692,
  "SMARCD3_994-1187",  "SNP151",    -0.674301,
  "SMARCD3_994-1187",  "SNP151",    -1.883472,
  "SMARCD3_994-1187",  "SNP151",    -0.752971,
  "SMARCD3_E1-3",      "SNP151",    -1.886206,
  "SMARCD3_E1-3",      "SNP151",    -3.552103,
  "SMARCD3_E1-3",      "SNP151",    -1.522605,
  "SMARCD3_E2-3",      "SNP151",     0.777308,
  "SMARCD3_E2-3",      "SNP151",     0.778661,
  "SMARCD3_E2-3",      "SNP151",     0.789886,
  "TLN1_6479-6683",    "SNP165",    -1.031504,
  "TLN1_6479-6683",    "SNP165",    -1.102274,
  "TLN1_6479-6683",    "SNP165",    -0.613089,
  "TRIP11_E7-9",       "TSS1",      -2.505959,
  "TRIP11_E7-9",       "TSS1",      -2.511784,
  "TRIP11_E7-9",       "TSS1",      -2.068261,
  "YY1AP1_E9-11",      "SNP31",     -0.658182,
  "YY1AP1_E9-11",      "SNP31",     -1.032693,
  "YY1AP1_E9-11",      "SNP31",      0.119401,
  "YY1AP1_E9-11",      "TSS1",      -1.016165,
  "YY1AP1_E9-11",      "TSS1",      -0.524085,
  "YY1AP1_E9-11",      "TSS1",       0.485162
)

sig_lookup <- tribble(
  ~Target,             ~Guide,      ~sig_label,
  "SMARCD3_E2-3",      "SNP151",    "FDR significant",
  "TRIP11_E7-9",       "TSS1",      "FDR significant",
  "LINC02454_E2-3",    "Upstream",  "FDR significant",
  "CXCL12_E2-3",       "SNP181",    "Nominally significant",
  "TLN1_6479-6683",    "SNP165",    "Nominally significant",
  "DAP3_E9-11",        "TSS1",      "Trending",
  "SMARCD3_E1-3",      "SNP151",    "Trending",
  "LINC02454_E2-3",    "SNP233",    "Trending",
  "HMGA2_E1_2-3",      "Upstream",  "Trending",
  "CXCL12_E2-E4",      "SNP181",    "Trending",
  "SMARCD3_994-1187",  "SNP151",    "Not significant",
  "YY1AP1_E9-11",      "SNP31",     "Not significant",
  "DAP3_E9-11",        "SNP31",     "Not significant",
  "YY1AP1_E9-11",      "TSS1",      "Not significant",
  "HMGA2_E1_2-3",      "SNP233",    "Not significant",
  "CXCL12_E2-3",       "SNP179",    "Not significant",
  "CXCL12_E2-E4",      "SNP179",    "Not significant"
)

locus_map <- c(
  "SMARCD3_E2-3" = "SMARCD3", "SMARCD3_E1-3" = "SMARCD3", "SMARCD3_994-1187" = "SMARCD3",
  "TRIP11_E7-9" = "TRIP11",
  "LINC02454_E2-3" = "LINC02454", "HMGA2_E1_2-3" = "LINC02454",
  "CXCL12_E2-3" = "CXCL12", "CXCL12_E2-E4" = "CXCL12",
  "TLN1_6479-6683" = "TLN1",
  "DAP3_E9-11" = "DAP3_YY1AP1", "YY1AP1_E9-11" = "DAP3_YY1AP1"
)

guide_label <- function(target, guide) {
  case_when(
    guide == "SNP151" ~ "rs219233",
    guide == "SNP179" ~ "rs4948786",
    guide == "SNP181" ~ "rs10899923",
    guide == "SNP165" ~ "rs4142495",
    guide == "SNP233" ~ "rs74100689",
    guide == "SNP31"  ~ "rs678157",
    guide == "Upstream" ~ "rs74100689 (alt.)",
    guide == "TSS1" & target %in% c("DAP3_E9-11", "YY1AP1_E9-11") ~ "rs678157 (alt.)",
    TRUE ~ guide
  )
}

raw <- raw %>%
  left_join(sig_lookup, by = c("Target", "Guide")) %>%
  mutate(
    Locus       = factor(locus_map[Target], levels = locus_levels),
    Gene        = sub("_.*$", "", Target),
    Amplicon    = sub("^[^_]*_", "", Target),
    GuideLabel  = guide_label(Target, Guide),
    TargetLabel = paste0("*", Gene, "*", " (", Amplicon, ")"),
    Label       = paste0(GuideLabel, " \u2192 ", TargetLabel)
  )

label_order <- raw %>%
  distinct(Label, Locus, Gene, Guide) %>%
  mutate(is_alt = grepl("alt\\.", Label)) %>%
  arrange(Locus, Gene, is_alt, Guide) %>%
  pull(Label)

raw <- raw %>% mutate(Label = factor(Label, levels = rev(unique(label_order))))

summary_df <- raw %>%
  group_by(Label, Locus, sig_label) %>%
  summarise(
    n       = n(),
    mean_l2 = mean(log2FC),
    sd_l2   = sd(log2FC),
    sem_l2  = sd_l2 / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    sig_symbol = case_when(
      sig_label == "FDR significant"       ~ "*",
      sig_label == "Nominally significant" ~ "\u2020",
      sig_label == "Trending"              ~ "\u2021",
      TRUE                                  ~ ""
    ),
    star_x = if_else(mean_l2 >= 0, mean_l2 + sem_l2 + 0.1, mean_l2 - sem_l2 - 0.1)
  )

gg_forest <- ggplot(summary_df, aes(x = mean_l2, y = Label)) +
  theme_classic() +
  geom_vline(xintercept = 0, linewidth = 0.1, color = "darkgrey") +
  geom_col(aes(fill = Locus), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(xmin = mean_l2 - sem_l2, xmax = mean_l2 + sem_l2), width = 0.25, linewidth = 0.1) +
  scale_fill_manual(values = locus_colors, drop = FALSE) +
  geom_point(data = raw, aes(x = log2FC, y = Label), alpha = 0.8, size = 0.6, color = "black") +
  geom_text(aes(x = star_x, label = sig_symbol),
            position = position_nudge(y = 0.32), size = 2.2) +
  labs(
    x = log2_axis_lab("hMSC"), y = NULL,
    caption = "* FDR q < 0.05    \u2020 nominal p < 0.05    \u2021 trending p < 0.10"
  ) +
  theme(
    axis.line   = element_line(linewidth = 0.1, color = "black"),
    axis.ticks  = element_line(linewidth = 0.1, color = "black"),
    text        = element_text(color = "black", size = 6),
    axis.text.x = element_text(color = "black", size = 6),
    axis.text.y = ggtext::element_markdown(color = "black", size = 6),  # italicized gene names
    plot.caption = element_text(size = 5, hjust = 0),
    legend.position = "none"
  )

# ==========================================================================
# Combined figure
# ==========================================================================
combined <- gg_cor + gg_forest +
  plot_layout(widths = c(0.63, 0.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

ggsave("revisions/hMSC_correlation_and_forest.pdf", combined, width = 120, height = 70, units = "mm")
