# Test each gRNA separately (rather than the three gRNAs per target together) for the
# significant cis effects (all cells) and for trans effects of selected targets, and compare
# effects between gRNAs targeting the same element (revision analyses).

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
setwd(project_dir)
system("mkdir -p revisions")

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  library(ondisc)
  library(sceptre)
  library(arrow)
  library(patchwork)
  library(ggpubr)
})

# QC'd UMI count matrices (all cells)
file_dir = paste0(project_dir, "/")
cDNA.pilot_full = readRDS(file.path(file_dir, "seurat/temp/cDNA.pilot_full.RDS"))
GDO.pilot_full = readRDS(file.path(file_dir, "seurat/temp/GDO.pilot_full.RDS"))
GDOminC.pilot_full = readRDS(file.path(file_dir, "seurat/temp/GDOminC.pilot_full.RDS"))

cell_data_pilot_full = fread(file.path(file_dir, "seurat/temp/cDNA.pilot_full.cell_data.txt"))
cDNA.pilot_full.umap = readRDS(file.path(file_dir, "seurat/temp/umap_objects/cDNA.pilot_full.RDS"))

# Cells carrying each gRNA (thresholded gRNA counts)
grna_cells = data.frame(
  grna_id = rownames(GDOminC.pilot_full),
  n_cells = rowSums(GetAssayData(GDOminC.pilot_full, layer = "counts")>0)) %>%
  as.data.table()
rownames(grna_cells) = NULL

# Gene TSS positions and gRNA targets
features <- fread(file.path(paste0(project_dir, "/cellranger/BMDSS-full-highnovax-1_submitjob/outs/filtered_feature_bc_matrix/features.tsv.gz")), header = F)
colnames(features) <- c("ensGene", "gene", "misc")
features_ = subset(features, misc == "Gene Expression" & gene != "CRISPRi")
biomart_genes <- fread(file.path(file_dir, "data/genes_hg38_ensembl110.txt"), header = T)
colnames(biomart_genes)[1] = "ensGene"
features_full <- features_ %>%
  dplyr::inner_join(biomart_genes, by = "ensGene") %>%
  dplyr::mutate(tss = ifelse(strand == 1, start_position, end_position))
gene_list = features_full %>%
  dplyr::select(gene, chromosome_name, tss) %>%
  dplyr::filter(chromosome_name %in% as.character(seq(1, 22))) %>%
  dplyr::mutate(chr_tss = as.numeric(chromosome_name)) %>%
  dplyr::select(-chromosome_name)
DT.lib <- fread(file.path(file_dir, "data/guides/230905_SSv3-BMD-hFOB-Library.txt"), fill = TRUE, na.strings = "")
colnames(DT.lib)[1] <- "gRNA"
tmp_gRNA_target = as.character(unlist(lapply(DT.lib$gRNA, function(x){strsplit(x, "__")[[1]][1]})))
tmp_gRNA_target = substr(tmp_gRNA_target, 0, nchar(tmp_gRNA_target)-2)
tmp_gRNA_target = ifelse(substr(tmp_gRNA_target, nchar(tmp_gRNA_target), nchar(tmp_gRNA_target)) == "-", substr(tmp_gRNA_target, 0, nchar(tmp_gRNA_target)-1), tmp_gRNA_target)
DT.lib$gRNA_target = tmp_gRNA_target
SNP_list = DT.lib %>%
  dplyr::filter(grepl("SNP-", gRNA)) %>%
  dplyr::select(gRNA, gRNA_target, snpid) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(chr = as.numeric(strsplit(snpid, split = ":")[[1]][1]),
                pos = as.numeric(strsplit(snpid, split = ":")[[1]][2])) %>%
  dplyr::ungroup()
POS_list = DT.lib %>%
  dplyr::filter(grepl("TSS", gRNA)) %>%
  dplyr::mutate(gene = gsub("-.*", "", gRNA_target)) %>%
  dplyr::left_join(gene_list, by = "gene") %>%
  dplyr::select(gRNA, gRNA_target, gene, chr_tss, tss) %>%
  dplyr::as_tibble()
NTC_list = DT.lib %>%
  dplyr::filter(grepl("NTC", gRNA)) %>%
  dplyr::select(gRNA, gRNA_target) %>%
  dplyr::mutate(gRNA_target = "non-targeting")
grna_target_data_frame = dplyr::bind_rows(
  SNP_list %>%
    dplyr::select(gRNA, gRNA_target, chr, pos) %>%
    dplyr::rename(grna_id = gRNA, grna_target = gRNA_target, start = pos) %>%
    dplyr::mutate(chr = paste0("chr", chr), end = start),
  POS_list %>%
    dplyr::select(gRNA, gRNA_target, chr_tss, tss) %>%
    dplyr::rename(grna_id = gRNA, grna_target = gRNA_target, chr = chr_tss, start = tss) %>%
    dplyr::mutate(chr = paste0("chr", chr), end = start),
  NTC_list %>%
    dplyr::select(gRNA, gRNA_target) %>%
    dplyr::rename(grna_id = gRNA, grna_target = gRNA_target))

# ---- Cis: each gRNA tested separately (singleton integration) ----
grna_target_data_frame_IndvGuides = grna_target_data_frame %>%
  dplyr::mutate(grna_target = ifelse(grna_target == "non-targeting", grna_target, grna_id))
counts <- GetAssayData(object = cDNA.pilot_full, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= ncol(counts)*0.05 # genes detected in >= 5% of cells
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = cDNA.pilot_full@meta.data)
cDNA.pilot_full_filtered = filtered_seurat

directory_to_write = paste0(project_dir, "/sceptre/260529_BMDSS_pilot_full_cis_pca_IndvGuides/ondisc")
input_cDNA = cDNA.pilot_full_filtered
input_PCs = cDNA.pilot_full.umap
input_GDO = GDO.pilot_full
input_celldata = cell_data_pilot_full

positive_control_pairs = dplyr::distinct(data.frame(grna_target = POS_list$gRNA, response_id = POS_list$gene)) %>%
  dplyr::filter(response_id %in% rownames(input_cDNA))
response_position_data_frame = gene_list %>%
  dplyr::filter(gene %in% rownames(input_cDNA)) %>%
  dplyr::rename(response_id = gene, chr = chr_tss, position = tss) %>%
  dplyr::mutate(chr = paste0("chr", chr)) %>%
  dplyr::relocate(response_id, chr)
response_matrix = LayerData(input_cDNA, layer = "counts")
grna_matrix = LayerData(input_GDO, layer = "counts")
extra_covariates = input_celldata %>%
  dplyr::left_join(as.data.table(input_PCs[["pca"]]@cell.embeddings) %>%
                     dplyr::mutate(cell = colnames(input_PCs)), by = "cell") %>%
  dplyr::rename(batch = orig.ident) %>%
  dplyr::select(batch, timepoint, CRISPRi, PC_1 , PC_2 , PC_3 , PC_4 , PC_5 , PC_6 , PC_7 , PC_8 , PC_9 , PC_10 ,
                PC_11 , PC_12 , PC_13 , PC_14 , PC_15 , PC_16 , PC_17 , PC_18 , PC_19 , PC_20 ,
                PC_21 , PC_22 , PC_23 , PC_24 , PC_25)
sceptre_object = import_data(response_matrix = response_matrix,
                             grna_matrix = grna_matrix,
                             grna_target_data_frame = grna_target_data_frame_IndvGuides,
                             moi = "high",
                             extra_covariates = extra_covariates,
                             use_ondisc = TRUE,
                             directory_to_write = directory_to_write)
discovery_pairs = dplyr::distinct(construct_cis_pairs(sceptre_object = sceptre_object,
                                                      distance_threshold = 1e6,
                                                      response_position_data_frame = response_position_data_frame))
formula_object = formula(~ log(response_n_nonzero) + log(response_n_umis) + response_p_mito + log(grna_n_umis) + log(grna_n_nonzero) + batch + timepoint + log1p(CRISPRi) +
                           PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10 +
                           PC_11 + PC_12 + PC_13 + PC_14 + PC_15 + PC_16 + PC_17 + PC_18 + PC_19 + PC_20 +
                           PC_21 + PC_22 + PC_23 + PC_24 + PC_25)
sceptre_object = set_analysis_parameters(sceptre_object = sceptre_object,
                                         discovery_pairs = discovery_pairs,
                                         positive_control_pairs = positive_control_pairs,
                                         side = "both",
                                         grna_integration_strategy = "singleton",
                                         formula_object = formula_object,
                                         resampling_mechanism = "permutations",
                                         multiple_testing_method = "BH",
                                         multiple_testing_alpha = 0.05)
write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object, directory_to_write = directory_to_write)

DC.0 = readRDS("sceptre/260529_BMDSS_pilot_full_cis_pca_IndvGuides/results_run_discovery_analysis.rds")

# Significant cis effects in the combined analysis (all time points)
SH.CRE = fread("summaries/250721_Sceptre-Results-CREs.txt")
SH.TSS = fread("summaries/250721_Sceptre-Results-TSSs.txt")

sig_hits = dplyr::bind_rows(
SH.CRE %>%
  dplyr::filter(timepoint == "dayA", sceptre_sig == 1) %>%
  dplyr::select(gene, grna_target),
SH.TSS %>%
  dplyr::filter(timepoint == "dayA", sceptre_sig == 1) %>%
  dplyr::select(gene, grna_target)
) %>%
  dplyr::distinct() %>%
  dplyr::mutate(gene_grna = paste0(gene, "_", grna_target))

DC.1 = DC.0 %>%
  dplyr::filter(pass_qc == TRUE) %>%
  dplyr::select(-grna_target, -pass_qc, -significant) %>%
  dplyr::left_join(grna_target_data_frame, by = "grna_id") %>%
  dplyr::rename(gene = response_id, sceptre_p = p_value, sceptre_log2fc = log_2_fold_change, chr_snp = chr, pos_snp = end) %>%
  dplyr::select(-start) %>%
  dplyr::left_join(gene_list, by = "gene") %>%
  dplyr::mutate(tss_dist = pos_snp - tss) %>%
  dplyr::relocate(gene, chr_tss, tss, tss_dist) %>%
  dplyr::filter(abs(tss_dist) < 1000000)

DC.2 = DC.1 %>%
  dplyr::mutate(gene_grna = paste0(gene, "_", grna_target)) %>%
  dplyr::filter(gene_grna %in% sig_hits$gene_grna) %>%
  dplyr::arrange(gene_grna, sceptre_p) %>%
  dplyr::left_join(grna_cells, by = "grna_id")

# Side-by-side results for gRNAs 1, 2 and 3 of each target
DC.3 = dplyr::full_join(DC.2 %>%
                    dplyr::filter(grepl("-1$", grna_id)), 
                  DC.2 %>%
                    dplyr::filter(grepl("-2$", grna_id)), 
                  by = "gene_grna", relationship = "many-to-many") %>%
  dplyr::select(gene_grna, gene.x, grna_id.x, n_nonzero_trt.x, n_cells.x, sceptre_p.x, sceptre_log2fc.x,
                grna_id.y, n_nonzero_trt.y, n_cells.y, sceptre_p.y, sceptre_log2fc.y) %>%
  dplyr::rename(gene = gene.x) %>%
  dplyr::full_join(DC.2 %>%
                      dplyr::filter(grepl("-3$", grna_id)) %>%
                      dplyr::select(gene_grna, grna_id, n_nonzero_trt, n_cells, sceptre_p, sceptre_log2fc), by = "gene_grna") %>%
  dplyr::rename(grna_id.z = grna_id, n_nonzero_trt.z = n_nonzero_trt, n_cells.z = n_cells, sceptre_p.z = sceptre_p, sceptre_log2fc.z = sceptre_log2fc)

# Candidate CREs: significance of the combined test and of each gRNA
DC.4 = DC.3 %>%
  dplyr::filter(!grepl("TSS", grna_id.x))

DC.5 = SH.CRE %>%
  dplyr::filter(timepoint == "dayA") %>%
  dplyr::mutate(gene_grna = paste0(gene, "_", grna_target)) %>%
  dplyr::select(gene_grna, grna_target, n_nonzero_trt, sceptre_p, sceptre_log2fc) %>%
  dplyr::inner_join(DC.4, by = "gene_grna") %>%
  dplyr::relocate(gene_grna, gene, grna_target)

DC.6 = DC.5 %>%
  pivot_longer(
    cols = starts_with("sceptre_p"),
    names_to = "name_var",
    values_to = "value_var"
  ) %>%
  dplyr::mutate(value_var = -log10(value_var)) %>%
  dplyr::mutate(value_var = ifelse(value_var > 16, 16, value_var))

fwrite(DC.5, file = "revisions/260709_IndvGuides_Comparison_CREs_LongTable.txt", col.names = T, row.names = F, quote = F, sep = "\t")
fwrite(DC.6, file = "revisions/260709_IndvGuides_Comparison_CREs_WideTable.txt", col.names = T, row.names = F, quote = F, sep = "\t")

gg1_cre = DC.6 %>%
  ggplot(aes(x = name_var, y = value_var, fill = name_var)) + 
  theme_classic(base_size = 6) +
  facet_wrap(. ~ gene_grna, nrow = 9) +
  geom_bar(stat = "identity", width = 0.4, color = NA, alpha = 0.8) + 
  geom_hline(yintercept = -log10(0.0009), linewidth = 0.1, color = "red") +
  geom_hline(yintercept = 0, linewidth = 0.1, color = "black") +
  scale_fill_manual(values = c("black", "skyblue3", "skyblue2", "skyblue1")) +
  scale_x_discrete(labels = c("All", "g1", "g2", "g3")) +
  xlab("gRNA tested") + ylab("Significance (-log10 p)") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.1),
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        legend.position = "none")

ggsave(filename = "revisions/IndvGuides_Comparison_CREs.pdf", plot = gg1_cre, height = 6, width = 8, device = "pdf")

# TSSs: the targeted gene only
DC.7 = DC.3 %>%
  dplyr::filter(grepl("TSS", grna_id.x))

DC.8 = SH.TSS %>%
  dplyr::filter(timepoint == "dayA") %>%
  dplyr::mutate(gene_grna = paste0(gene, "_", grna_target)) %>%
  dplyr::select(gene_grna, grna_target, n_nonzero_trt, sceptre_p, sceptre_log2fc) %>%
  dplyr::inner_join(DC.7, by = "gene_grna") %>%
  dplyr::relocate(gene_grna, gene, grna_target)

DC.9 = DC.8 %>%
  pivot_longer(
    cols = starts_with("sceptre_p"),
    names_to = "name_var",
    values_to = "value_var"
  ) %>%
  dplyr::mutate(value_var = -log10(value_var)) %>%
  dplyr::mutate(value_var = ifelse(value_var > 16, 16, value_var)) %>%
  dplyr::mutate(grna_genetar = sub("([^-]+).*", "\\1", grna_target)) %>%
  dplyr::filter(gene == grna_genetar)

fwrite(DC.8, file = "revisions/260709_IndvGuides_Comparison_TSSs_LongTable.txt", col.names = T, row.names = F, quote = F, sep = "\t")
fwrite(DC.9, file = "revisions/260709_IndvGuides_Comparison_TSSs_WideTable.txt", col.names = T, row.names = F, quote = F, sep = "\t")

gg2_tss = DC.9 %>%
  ggplot(aes(x = name_var, y = value_var, fill = name_var)) + 
  theme_classic(base_size = 6) +
  facet_wrap(. ~ gene_grna, nrow = 9) +
  geom_bar(stat = "identity", width = 0.4, color = NA, alpha = 0.8) + 
  geom_hline(yintercept = -log10(0.001), linewidth = 0.1, color = "red") +
  geom_hline(yintercept = 0, linewidth = 0.1, color = "black") +
  scale_fill_manual(values = c("black", "skyblue3", "skyblue2", "skyblue1")) +
  scale_x_discrete(labels = c("All", "g1", "g2", "g3")) +
  xlab("gRNA tested") + ylab("Significance (-log10 p)") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.1),
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        legend.position = "none")

ggsave(filename = "revisions/IndvGuides_Comparison_TSSs.pdf", plot = gg2_tss, height = 6, width = 4, device = "pdf")

# ---- Trans: each gRNA tested separately for targets with significant trans effects ----
counts <- GetAssayData(object = cDNA.pilot_full, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= ncol(counts)*0.2 # genes detected in >= 20% of cells
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = cDNA.pilot_full@meta.data)
cDNA.pilot_full_filtered = filtered_seurat

directory_to_write = "sceptre/260603_BMDSS_pilot_full_trans_pca_IndvGuides/ondisc"
input_cDNA = cDNA.pilot_full_filtered
input_PCs = cDNA.pilot_full.umap
input_GDO = GDO.pilot_full
input_celldata = cell_data_pilot_full

response_matrix = LayerData(input_cDNA, layer = "counts")
grna_matrix = LayerData(input_GDO, layer = "counts")
extra_covariates = input_celldata %>%
  dplyr::left_join(as.data.table(input_PCs[["pca"]]@cell.embeddings) %>%
                     dplyr::mutate(cell = colnames(input_PCs)), by = "cell") %>%
  dplyr::rename(batch = orig.ident) %>%
  dplyr::select(batch, timepoint, CRISPRi, PC_1 , PC_2 , PC_3 , PC_4 , PC_5 , PC_6 , PC_7 , PC_8 , PC_9 , PC_10 ,
                PC_11 , PC_12 , PC_13 , PC_14 , PC_15 , PC_16 , PC_17 , PC_18 , PC_19 , PC_20 ,
                PC_21 , PC_22 , PC_23 , PC_24 , PC_25)

# gRNAs for targets with significant trans effects, plus non-targeting gRNAs
DT.0 = fread("summaries/241002_BMDSS.dayA.trans-sig.txt")
all_hits = unique(DT.0$grna_target)
grna_matrix_hits = grna_matrix[rownames(grna_matrix) %in% c(paste0(all_hits, "-1"), paste0(all_hits, "-2"), paste0(all_hits, "-3")) | grepl("NTC", rownames(grna_matrix)),]

grna_target_data_frame_hits = subset(grna_target_data_frame, grna_id %in% rownames(grna_matrix_hits)) %>%
  dplyr::mutate(grna_target = ifelse(grna_target == "non-targeting", grna_target, grna_id))
response_matrix_noCRISPRi = response_matrix[rownames(response_matrix) != "CRISPRi",]
formula_object = formula(~ log(response_n_nonzero) + log(response_n_umis) + response_p_mito + log(grna_n_umis) + log(grna_n_nonzero) + batch + log1p(CRISPRi) + timepoint +
                           PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10 +
                           PC_11 + PC_12 + PC_13 + PC_14 + PC_15 + PC_16 + PC_17 + PC_18 + PC_19 + PC_20 +
                           PC_21 + PC_22 + PC_23 + PC_24 + PC_25)
sceptre_object = import_data(response_matrix = response_matrix_noCRISPRi, 
                             grna_matrix = grna_matrix_hits, 
                             grna_target_data_frame = grna_target_data_frame_hits, 
                             moi = "high", 
                             extra_covariates = extra_covariates, 
                             use_ondisc = TRUE, 
                             directory_to_write = directory_to_write)
sceptre_object = set_analysis_parameters(sceptre_object = sceptre_object,
                                         side = "both",
                                         grna_integration_strategy = "union",
                                         formula_object = formula_object,
                                         resampling_mechanism = "permutations")
write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object, directory_to_write = directory_to_write)

# Trans results for individual gRNAs (FDR 1% per gRNA)
DT.0 = collect(open_dataset(paste0(project_dir, "/sceptre/260618_BMDSS_pilot_full_trans_pca_massive_NTC_3_10_singleton/trans_results")))

GetSigNetworks = function(x){
  tmp0 = x %>%
    dplyr::distinct() %>%
    dplyr::filter(pass_qc == 1) %>%
    dplyr::group_by(grna_id) %>%
    dplyr::arrange(p_value) %>%
    dplyr::mutate(obs.p = -log10(p_value),
                  exp.p = -log10(1:length(p_value)/length(p_value)),
                  q_value = p.adjust(p_value, method = "BH"),
                  significant = ifelse(q_value < 0.01, 1, 0)) %>%
    dplyr::ungroup() %>%
    as.data.table()
}

DT.1 = GetSigNetworks(DT.0)

# Compare gRNA effects for the significant trans effects of each target
sig_hits = fread("summaries/241002_BMDSS.dayA.trans-sig.txt") %>%
dplyr::mutate(gene_grna = paste0(response_id, "_", grna_target))

DT.2 = DT.1 %>%
  dplyr::rename(gene = response_id, sceptre_p = p_value, sceptre_log2fc = log_2_fold_change) %>%
  dplyr::mutate(gene_grna = paste0(gene, "_", grna_target)) %>%
  dplyr::filter(gene_grna %in% sig_hits$gene_grna) %>%
  dplyr::arrange(gene_grna, sceptre_p) %>%
  dplyr::left_join(grna_cells, by = "grna_id")

DT.3 = dplyr::full_join(DT.2 %>%
                          dplyr::filter(grepl("-1$", grna_id)), 
                        DT.2 %>%
                          dplyr::filter(grepl("-2$", grna_id)), 
                        by = "gene_grna", relationship = "many-to-many") %>%
  dplyr::select(gene_grna, gene.x, grna_id.x, n_nonzero_trt.x, n_cells.x, sceptre_p.x, sceptre_log2fc.x,
                grna_id.y, n_nonzero_trt.y, n_cells.y, sceptre_p.y, sceptre_log2fc.y) %>%
  dplyr::rename(gene = gene.x) %>%
  dplyr::full_join(DT.2 %>%
                     dplyr::filter(grepl("-3$", grna_id)) %>%
                     dplyr::select(gene_grna, grna_id, n_nonzero_trt, n_cells, sceptre_p, sceptre_log2fc), by = "gene_grna") %>%
  dplyr::rename(grna_id.z = grna_id, n_nonzero_trt.z = n_nonzero_trt, n_cells.z = n_cells, sceptre_p.z = sceptre_p, sceptre_log2fc.z = sceptre_log2fc)

DT.4 = sig_hits %>%
  dplyr::rename(sceptre_p = p_value, sceptre_log2fc = log_2_fold_change) %>%
  dplyr::select(gene_grna, grna_target, n_nonzero_trt, sceptre_p, sceptre_log2fc) %>%
  dplyr::full_join(DT.3, by = "gene_grna") %>%
  dplyr::relocate(gene_grna, gene, grna_target)

my_style <- list(
  theme_classic(base_size = 6),
  geom_hline(yintercept = 0, linewidth = 0.1, color = "black"),
  geom_vline(xintercept = 0, linewidth = 0.1, color = "black"),
  geom_smooth(method = "lm", formula = 'y ~ x', se = FALSE, color = "red", linewidth = 0.2, fullrange = TRUE),
  theme(axis.line = element_line(colour = "black", linewidth = 0.1),
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        legend.position = "none"))

# Pearson correlation of gRNA fold-changes
cor_layer <- ggpubr::stat_cor(aes(label = after_stat(r.label)), method = "pearson", size = 2,
                              label.x.npc = "left", label.y.npc = "top", color = "black")

# gene name for the panel title, given the grna_target value
gene_label <- function(x) {
  dplyr::case_when(
    x == "SNP-31"  ~ "DAP3/YY1AP1",
    x == "SNP-219" ~ "CHD4",
    TRUE ~ sub("-TSS[0-9]+$", "", x)
  )
}

# axis label with a subscripted "2" in "log2"
log2_lab <- function(prefix) {
  bquote(.(prefix) ~ "fold-change (" * log[2] * ")")
}

# Fold-changes of the trans-regulated genes for each pair of gRNAs (cis target genes excluded)
plot_fxn = function(x, y){
  DT.4_plot = subset(DT.4, grna_target == x) %>%
    dplyr::filter(!(gene == "DAP3" & grna_target == "SNP-31"),
                  !(gene == "YY1AP1" & grna_target == "SNP-31"),
                  !(gene == "CHD4" & grna_target == "SNP-219"),
                  !(gene == "TRIP11" & grna_target == "TRIP11-TSS1"),
                  !(gene == "UBP1" & grna_target == "UBP1-TSS4"),
                  !(gene == "ZBTB4" & grna_target == "ZBTB4-TSS3"))
  lim_xmax = max(c(DT.4_plot$sceptre_log2fc.x, DT.4_plot$sceptre_log2fc.y, DT.4_plot$sceptre_log2fc.z), na.rm = T)*1.1
  lim_xmin = min(c(DT.4_plot$sceptre_log2fc.x, DT.4_plot$sceptre_log2fc.y, DT.4_plot$sceptre_log2fc.z), na.rm = T)*1.1
  lim_vars = c(lim_xmin, lim_xmax)
  color_var = y
  gene_lab = gene_label(x)
  
  gg_1v2 = ggplot(DT.4_plot, aes(x = sceptre_log2fc.x, y = sceptre_log2fc.y)) + geom_point(size = 1.2, alpha = 0.5, shape = 21, fill = color_var, stroke = NA) + my_style + cor_layer + xlab(log2_lab("g1")) + ylab(log2_lab("g2")) + lims(x = lim_vars, y = lim_vars) + ggtitle(paste0(gene_lab, ", g1 ~ g2"))
  gg_1v3 = ggplot(DT.4_plot, aes(x = sceptre_log2fc.x, y = sceptre_log2fc.z)) + geom_point(size = 1.2, alpha = 0.5, shape = 21, fill = color_var, stroke = NA) + my_style + cor_layer + xlab(log2_lab("g1")) + ylab(log2_lab("g3")) + lims(x = lim_vars, y = lim_vars) + ggtitle(paste0(gene_lab, ", g1 ~ g3"))
  gg_2v3 = ggplot(DT.4_plot, aes(x = sceptre_log2fc.y, y = sceptre_log2fc.z)) + geom_point(size = 1.2, alpha = 0.5, shape = 21, fill = color_var, stroke = NA) + my_style + cor_layer + xlab(log2_lab("g2")) + ylab(log2_lab("g3")) + lims(x = lim_vars, y = lim_vars) + ggtitle(paste0(gene_lab, ", g2 ~ g3"))
  
  gg_1v2 + gg_1v3 + gg_2v3 + plot_layout(nrow = 1)
}

gg_out_5 = plot_fxn("SNP-31", "#ffa000")
gg_out_3 = plot_fxn("SNP-219", "dodgerblue1")

gg_out_1 = plot_fxn("TRIP11-TSS1", "purple")
gg_out_2 = plot_fxn("UBP1-TSS4", "#0097a7")
gg_out_4 = plot_fxn("ZBTB4-TSS3", "firebrick1")

gg_panel = gg_out_1 / gg_out_2 / gg_out_3 / gg_out_4 / gg_out_5 +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold")))

ggsave(plot = gg_panel, filename = "revisions/IndvGuides_Comparison_Trans.pdf", height = 5*1.3, width = 3*1.35, device = "pdf")
