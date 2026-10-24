# Violin plots of target gene expression in cells carrying gRNAs for a candidate CRE versus
# 1,500 randomly sampled cells without them (CXCL12, SMARCD3, DAP3/YY1AP1 and FHL3 CREs).

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
setwd(project_dir)
system("mkdir -p plots")

library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)

cDNA.pilot_full.umap <- readRDS("seurat/temp/umap_objects/cDNA.pilot_full.RDS")
cell_data_pilot_full <- fread("seurat/temp/cDNA.pilot_full.cell_data.txt")

GDO.pilot_full <- readRDS("seurat/temp/GDOminC.pilot_full.RDS")

# Per-cell metadata
cell_data_df <- as.data.frame(cell_data_pilot_full)
rownames(cell_data_df) <- cell_data_df$cell

matched_indices <- match(rownames(cDNA.pilot_full.umap@meta.data), cell_data_df$cell)
cell_data_matched <- cell_data_df[matched_indices, ]

cDNA.pilot_full.umap <- AddMetaData(cDNA.pilot_full.umap, metadata = cell_data_matched)

cDNA.pilot_full.umap <- NormalizeData(cDNA.pilot_full.umap, scale.factor = median(cell_data_pilot_full$nCount_RNA))

# Sum the (thresholded) counts of the gRNAs for each target, e.g. SNP-31-1/2/3 -> SNP-31
gdo_matrix <- GetAssayData(GDO.pilot_full, assay = "GDO", layer = "counts")

collapse_gRNAs_by_target <- function(gRNA_obj) {
  grna_names <- rownames(gRNA_obj)
  all_cells <- colnames(gRNA_obj)
  
  collapse_keys <- vapply(grna_names, function(name) {
    parts <- strsplit(name, "-")[[1]]
    if (length(parts) >= 2) {
      paste(parts[1:2], collapse = "-")
    } else {
      parts[1]
    }
  }, character(1))
  
  grna_groups <- split(seq_along(grna_names), collapse_keys)
  
  summed_list <- lapply(grna_groups, function(idxs) {
    colsum <- Matrix::colSums(gRNA_obj[idxs, , drop = FALSE])
    full_vector <- numeric(length(all_cells))
    names(full_vector) <- all_cells
    full_vector[names(colsum)] <- colsum
    Matrix::Matrix(full_vector, nrow = 1, dimnames = list(NULL, all_cells), sparse = TRUE)
  })
  
  collapsed_mat <- do.call(rbind, summed_list)
  rownames(collapsed_mat) <- names(grna_groups)
  return(collapsed_mat)
}

GDO.collapsed <- collapse_gRNAs_by_target(gdo_matrix)

# Expression of a gene in cells with >= min_count UMIs for the target's gRNAs vs 1,500 other cells
plot_locus_violin_by_perturbation <- function(seurat_obj = cDNA.pilot_full.umap, gRNA_obj = GDO.collapsed, gene, grna_prefix = NULL, min_count = 3) {
  
  if (!any(gene %in% rownames(seurat_obj))) {
    message(paste("Skipping:", gene, "not found in expression matrix"))
    return(NULL)
  }
  
  is_collapsed <- inherits(gRNA_obj, "dgCMatrix")
  
  shared_cells <- intersect(Cells(seurat_obj),
                            if (is_collapsed) colnames(gRNA_obj) else Cells(gRNA_obj))
  if (length(shared_cells) < 10) {
    message("Too few shared cells")
    return(NULL)
  }
  
  if (is_collapsed) {
    if (!(grna_prefix %in% rownames(gRNA_obj))) {
      message(paste("No matching row in collapsed gRNA object for:", grna_prefix))
      return(NULL)
    }
    gRNA_counts <- gRNA_obj[grna_prefix, shared_cells, drop = TRUE]
    perturbed_cells <- names(gRNA_counts)[gRNA_counts >= min_count]
  } else {
    matching_grnas <- grep(paste0("^", grna_prefix), rownames(gRNA_obj), value = TRUE)
    if (length(matching_grnas) == 0) {
      message(paste("No matching gRNAs found for prefix:", grna_prefix))
      return(NULL)
    }
    gRNA_counts <- FetchData(gRNA_obj, vars = matching_grnas)[shared_cells, , drop = FALSE]
    gRNA_sums <- rowSums(gRNA_counts)
    perturbed_cells <- names(gRNA_sums)[gRNA_sums >= min_count]
  }
  
  non_targeting_pool <- setdiff(shared_cells, perturbed_cells)
  set.seed(123); sampled_non_targeting <- if (length(non_targeting_pool) < 1500) {
    warning("Fewer than 1500 non-targeting cells available; using all available.")
    non_targeting_pool
  } else {
    set.seed(123); sample(non_targeting_pool, 1500)
  }
  
  final_cells <- c(perturbed_cells, sampled_non_targeting)
  temp_obj <- seurat_obj[, final_cells]
  
  perturb_status_vector <- ifelse(colnames(temp_obj) %in% perturbed_cells,
                                  "Perturbed", "Non-targeting")
  temp_obj@meta.data$perturb_status <- factor(perturb_status_vector,
                                              levels = c("Non-targeting", "Perturbed"))
  
  if (!"timepoint" %in% colnames(temp_obj@meta.data)) {
    message("Missing 'timepoint' metadata column")
    return(NULL)
  }
  
  gene_range = seq(1, length(gene))
  for(i in gene_range){
    gene_var = gene[i]
    temp_obj = AddMetaData(temp_obj, FetchData(temp_obj, layer = "counts", vars = gene[i]), col.name = gene_var)
  }
  temp_dt = temp_obj@meta.data
  temp_dt$perturb_label = as.factor(gsub("Non-targeting", "NT", gsub("Perturbed", grna_prefix, temp_dt$perturb_status)))
  
  scale_factor <- median(seurat_obj@meta.data$nCount_RNA)
  
  gene_index = which(colnames(temp_dt) == gene)
  colnames(temp_dt)[gene_index] = "gene_exp"
  
  set.seed(123); temp_dt_long = temp_dt %>%
    dplyr::mutate(norm_val = gene_exp/nCount_RNA,
                  scale_val = norm_val*scale_factor,
                  plot_val = log1p(scale_val)
                  ) %>%
    dplyr::select(-orig.ident, -nCount_RNA, -nFeature_RNA, -percent.mt, -nCount_GDO, -nFeature_GDO, -nCount_GDOminC, -nFeature_GDOminC, -CRISPRi, -study, -timepoint, -perturb_status, -gene_exp, -norm_val, -scale_val) %>%
    dplyr::relocate(cell, perturb_label) %>%
    tidyr::pivot_longer(-c(cell, perturb_label), values_to = "gex", names_to = "gene")  %>%
    dplyr::mutate(gex = gex + rnorm(n(), sd = 1/100000))
  temp_dt_long$perturb_label = factor(temp_dt_long$perturb_label, levels = c("NT", grna_prefix))
  temp_dt_long$gene = factor(temp_dt_long$gene, levels = gene)
  set.seed(123); gg1 = ggplot(temp_dt_long, aes(x = perturb_label, y = gex)) +
    theme_classic(base_size = 4) +
    geom_jitter(stroke = 0, alpha = 0.1, size = 0.4, color = "lightgrey", width = 0.4) +
    geom_violin(aes(color = perturb_label), linewidth = 0.2, fill = NA, width = 0.6, scale = "width") +
    scale_color_manual(values = c("darkgrey", "dodgerblue1")) +
    facet_wrap(~ gene, strip.position = "bottom") +
    stat_summary(fun = "mean", geom = "point", shape = 23, size = 0.8, color = "black", stroke = 0.2) +
    theme(legend.position = "none",
          axis.line = element_line(color = "black", linewidth = 0.1),
          axis.ticks = element_line(color = "black", linewidth = 0.1),
          axis.text = element_text(color = "black"),
          strip.background = element_blank(),
          strip.placement = "outside",
          plot.title = element_text(hjust = 0.5, size = 4)) +
    xlab("") + ylab("Gene expression") +
    ggtitle(grna_prefix)
  return(gg1)
}

# ---- CXCL12 CREs ----
plot_output_CXCL12 = lapply(c("SNP-178", "SNP-179", "SNP-180", "SNP-181"), function(g) {
  plot_locus_violin_by_perturbation(gene = "CXCL12", grna_prefix = g) + ylim(c(-0.2, 2.7))
}) %>% wrap_plots(nrow = 1)
ggsave(plot_output_CXCL12, filename = "plots/ViolinPlot_CXCL12.pdf", device = "pdf", height = 1.2, width = 2.4)

# ---- SMARCD3 CREs ----
plot_output_SMARCD3 = lapply(c("SNP-149", "SNP-150", "SNP-151", "SNP-152"), function(g) {
  plot_locus_violin_by_perturbation(gene = "SMARCD3", grna_prefix = g)
}) %>% wrap_plots(nrow = 1)
ggsave(plot_output_SMARCD3, filename = "plots/ViolinPlot_SMARCD3.pdf", device = "pdf", height = 1.2, width = 2.4)

# ---- DAP3/YY1AP1 bidirectional promoter CRE (SNP-31) ----
plot_output_DAP3 = plot_locus_violin_by_perturbation(gene = "DAP3", grna_prefix = "SNP-31") + ylim(c(-4.2, 3.5))
plot_output_YY1AP1 = plot_locus_violin_by_perturbation(gene = "YY1AP1", grna_prefix = "SNP-31") + ylim(c(-4.2, 3.5))
plot_output_SNP31 = plot_output_YY1AP1 + plot_output_DAP3 + plot_layout(nrow = 1)
ggsave(plot_output_SNP31, filename = "plots/ViolinPlot_SNP31.pdf", device = "pdf", height = 1.2, width = 1.6)

# ---- FHL3 CRE (SNP-17) ----
plot_output_FHL3 = plot_locus_violin_by_perturbation(gene = "FHL3", grna_prefix = "SNP-17") + ylim(c(-4.2, 3.4))
ggsave(plot_output_FHL3, filename = "plots/ViolinPlot_FHL3.pdf", device = "pdf", height = 1.2, width = 0.8)
