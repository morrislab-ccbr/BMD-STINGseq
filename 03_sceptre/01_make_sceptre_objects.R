# Build the ondisc-backed SCEPTRE objects analysed with sceptre-pipeline (02_run_sceptre.sh).
#   Rscript 01_make_sceptre_objects.R cis     # cis analyses (gRNA target vs genes within 1 Mb)
#   Rscript 01_make_sceptre_objects.R trans   # transcriptome-wide analyses; needs
#                                             # summaries/250114_gRNA-hits.txt (04_results/01_annotate_sceptre_results.R)
# All analyses adjust for lane (batch), CRISPRi transgene expression and 25 expression PCs.

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
setwd(project_dir)

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  library(ondisc)
  library(sceptre)
})

mode <- commandArgs(TRUE)[1]
pcs <- paste0("PC_", 1:25)

# ---- Genes (TSS positions) and gRNA targets ----
features <- fread("cellranger/BMDSS-full-highnovax-1_submitjob/outs/filtered_feature_bc_matrix/features.tsv.gz", header = F)
colnames(features) <- c("ensGene", "gene", "misc")
biomart_genes <- fread("data/genes_hg38_ensembl110.txt", header = T)
colnames(biomart_genes)[1] <- "ensGene"
gene_list <- subset(features, misc == "Gene Expression" & gene != "CRISPRi") %>%
  dplyr::inner_join(biomart_genes, by = "ensGene") %>%
  dplyr::mutate(tss = ifelse(strand == 1, start_position, end_position)) %>%
  dplyr::select(gene, chromosome_name, tss) %>%
  dplyr::filter(chromosome_name %in% as.character(seq(1, 22))) %>%
  dplyr::mutate(chr_tss = as.numeric(chromosome_name)) %>%
  dplyr::select(-chromosome_name)

# gRNA IDs are <target>-<1|2|3>; targets are SNP-<n> (candidate CREs), <gene>-TSS<n> or NTC
DT.lib <- fread("data/guides/230905_SSv3-BMD-hFOB-Library.txt", fill = TRUE, na.strings = "")
colnames(DT.lib)[1] <- "gRNA"
tmp_gRNA_target <- as.character(unlist(lapply(DT.lib$gRNA, function(x){strsplit(x, "__")[[1]][1]})))
tmp_gRNA_target <- substr(tmp_gRNA_target, 0, nchar(tmp_gRNA_target)-2)
tmp_gRNA_target <- ifelse(substr(tmp_gRNA_target, nchar(tmp_gRNA_target), nchar(tmp_gRNA_target)) == "-", substr(tmp_gRNA_target, 0, nchar(tmp_gRNA_target)-1), tmp_gRNA_target)
DT.lib$gRNA_target <- tmp_gRNA_target

SNP_list <- DT.lib %>%
  dplyr::filter(grepl("SNP-", gRNA)) %>%
  dplyr::select(gRNA, gRNA_target, snpid) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(chr = as.numeric(strsplit(snpid, split = ":")[[1]][1]),
                pos = as.numeric(strsplit(snpid, split = ":")[[1]][2])) %>%
  dplyr::ungroup()
POS_list <- DT.lib %>%
  dplyr::filter(grepl("TSS", gRNA)) %>%
  dplyr::mutate(gene = gsub("-.*", "", gRNA_target)) %>%
  dplyr::left_join(gene_list, by = "gene") %>%
  dplyr::select(gRNA, gRNA_target, gene, chr_tss, tss) %>%
  dplyr::as_tibble()
NTC_list <- DT.lib %>%
  dplyr::filter(grepl("NTC", gRNA)) %>%
  dplyr::select(gRNA, gRNA_target) %>%
  dplyr::mutate(gRNA_target = "non-targeting")

grna_target_data_frame <- dplyr::bind_rows(
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
    dplyr::rename(grna_id = gRNA, grna_target = gRNA_target)
)

# ---- Helpers ----
cell_data_pilot_full <- fread("seurat/temp/cDNA.pilot_full.cell_data.txt")

load_set <- function(label) {
  list(cDNA = readRDS(paste0("seurat/temp/cDNA.", label, ".RDS")),
       GDO = readRDS(paste0("seurat/temp/GDO.", label, ".RDS")),
       PCs = readRDS(paste0("seurat/temp/umap_objects/cDNA.", label, ".RDS")),
       cell_data = fread(paste0("seurat/temp/cDNA.", label, ".cell_data.txt")))
}

# Keep genes detected in at least the given fraction of cells
filter_genes <- function(obj, min_frac) {
  counts <- GetAssayData(object = obj, layer = "counts")
  keep_genes <- Matrix::rowSums(counts > 0) >= ncol(counts)*min_frac
  CreateSeuratObject(counts[keep_genes, ], meta.data = obj@meta.data)
}

# Covariates: lane, CRISPRi UMIs and 25 PCs (plus time point for sets spanning several days)
make_covariates <- function(set, timepoint) {
  cd <- set$cell_data
  if (timepoint && !"timepoint" %in% names(cd)) {
    cd <- dplyr::left_join(cd, cell_data_pilot_full %>% dplyr::select(cell, timepoint), by = "cell")
  }
  cd %>%
    dplyr::left_join(as.data.table(set$PCs[["pca"]]@cell.embeddings) %>%
                       dplyr::mutate(cell = colnames(set$PCs)), by = "cell") %>%
    dplyr::rename(batch = orig.ident) %>%
    dplyr::select(batch, CRISPRi, dplyr::any_of(if (timepoint) "timepoint" else character()), dplyr::all_of(pcs))
}

make_formula <- function(timepoint) {
  as.formula(paste("~ log(response_n_nonzero) + log(response_n_umis) + response_p_mito + log(grna_n_umis) + log(grna_n_nonzero) + batch + log1p(CRISPRi) +",
                   if (timepoint) "timepoint +" else "",
                   paste(pcs, collapse = " + ")))
}

# Cis: test each gRNA target against genes with a TSS within 1 Mb; TSS-targeting gRNAs vs their own gene are positive controls
make_cis_object <- function(directory_to_write, set, min_frac = 0.05, timepoint = FALSE) {
  system(paste0("mkdir -p ", directory_to_write))
  cDNA <- filter_genes(set$cDNA, min_frac)
  positive_control_pairs <- dplyr::distinct(data.frame(grna_target = POS_list$gRNA_target, response_id = POS_list$gene)) %>%
    dplyr::filter(response_id %in% rownames(cDNA))
  response_position_data_frame <- gene_list %>%
    dplyr::filter(gene %in% rownames(cDNA)) %>%
    dplyr::rename(response_id = gene, chr = chr_tss, position = tss) %>%
    dplyr::mutate(chr = paste0("chr", chr)) %>%
    dplyr::relocate(response_id, chr)
  sceptre_object <- import_data(response_matrix = LayerData(cDNA, layer = "counts"),
                                grna_matrix = LayerData(set$GDO, layer = "counts"),
                                grna_target_data_frame = grna_target_data_frame,
                                moi = "high",
                                extra_covariates = make_covariates(set, timepoint),
                                use_ondisc = TRUE,
                                directory_to_write = directory_to_write)
  discovery_pairs <- dplyr::distinct(construct_cis_pairs(sceptre_object = sceptre_object,
                                                         distance_threshold = 1e6,
                                                         response_position_data_frame = response_position_data_frame))
  sceptre_object <- set_analysis_parameters(sceptre_object = sceptre_object,
                                            discovery_pairs = discovery_pairs,
                                            positive_control_pairs = positive_control_pairs,
                                            side = "both",
                                            grna_integration_strategy = "union",
                                            formula_object = make_formula(timepoint),
                                            resampling_mechanism = "permutations",
                                            multiple_testing_method = "BH",
                                            multiple_testing_alpha = 0.05)
  write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object, directory_to_write = directory_to_write)
}

# Trans: test gRNA targets against all genes (CRISPRi transgene removed); optionally only gRNAs
# for targets with a significant cis effect (plus non-targeting gRNAs)
make_trans_object <- function(directory_to_write, set, min_frac, grna_targets = NULL) {
  system(paste0("mkdir -p ", directory_to_write))
  cDNA <- filter_genes(set$cDNA, min_frac)
  response_matrix <- LayerData(cDNA, layer = "counts")
  response_matrix <- response_matrix[rownames(response_matrix) != "CRISPRi",]
  grna_matrix <- LayerData(set$GDO, layer = "counts")
  grna_targets_df <- grna_target_data_frame
  if (!is.null(grna_targets)) {
    grna_matrix <- grna_matrix[rownames(grna_matrix) %in% c(paste0(grna_targets, "-1"), paste0(grna_targets, "-2"), paste0(grna_targets, "-3")) | grepl("NTC", rownames(grna_matrix)),]
    grna_targets_df <- subset(grna_target_data_frame, grna_target %in% grna_targets)
  }
  sceptre_object <- import_data(response_matrix = response_matrix, grna_matrix = grna_matrix, grna_target_data_frame = grna_targets_df,
                                moi = "high", extra_covariates = make_covariates(set, FALSE), use_ondisc = TRUE, directory_to_write = directory_to_write)
  sceptre_object <- set_analysis_parameters(sceptre_object = sceptre_object,
                                            side = "both",
                                            grna_integration_strategy = "union",
                                            formula_object = make_formula(FALSE),
                                            resampling_mechanism = "permutations")
  write_ondisc_backed_sceptre_object(sceptre_object = sceptre_object, directory_to_write = directory_to_write)
}

# ---- Cis analyses ----
if (mode == "cis") {
  make_cis_object("sceptre/240831_BMDSS_pilot_full_day0_cis_pca/ondisc", load_set("pilot_full_day0"))
  make_cis_object("sceptre/240831_BMDSS_full_day2_cis_pca/ondisc", load_set("day2"))
  make_cis_object("sceptre/240831_BMDSS_full_day4_cis_pca/ondisc", load_set("day4"))
  make_cis_object("sceptre/250210_BMDSS_full_day2_day4_cis_pca/ondisc", load_set("full_day2_day4"), timepoint = TRUE)
  make_cis_object("sceptre/240831_BMDSS_pilot_full_cis_pca25/ondisc", load_set("pilot_full"), timepoint = TRUE)

  # Downsampled to 10,000 cells per experiment and time point
  pilot_full <- list(cDNA = readRDS("seurat/temp/cDNA.pilot_full.RDS"), GDO = readRDS("seurat/temp/GDO.pilot_full.RDS"))
  downsample_sets <- list(pilot_day0 = quote(timepoint == "day0" & study == "pilot"),
                          full_day0 = quote(timepoint == "day0" & study == "full"),
                          full_day2 = quote(timepoint == "day2"),
                          full_day4 = quote(timepoint == "day4"))
  for (label in names(downsample_sets)) {
    set.seed(123); cell_data_ds <- cell_data_pilot_full %>%
      dplyr::filter(!!downsample_sets[[label]]) %>%
      dplyr::mutate(order_var = row_number()) %>%
      dplyr::sample_n(10000) %>%
      dplyr::arrange(order_var)
    cDNA_ds <- pilot_full$cDNA[, colnames(pilot_full$cDNA) %in% cell_data_ds$cell]
    ds_set <- list(cDNA = cDNA_ds,
                   GDO = pilot_full$GDO[, colnames(pilot_full$GDO) %in% cell_data_ds$cell],
                   PCs = NormalizeData(cDNA_ds, verbose = F) %>% FindVariableFeatures(verbose = F) %>% ScaleData(verbose = F) %>% RunPCA(verbose = F),
                   cell_data = cell_data_ds)
    make_cis_object(paste0("sceptre/250114_BMDSS_", label, "_ds/ondisc"), ds_set)
    gc()
  }
}

# ---- Trans analyses ----
if (mode == "trans") {
  all_hits <- fread("summaries/250114_gRNA-hits.txt", header = F)$V1
  make_trans_object("sceptre/240831_BMDSS_pilot_full_day0_trans_pca/ondisc", load_set("pilot_full_day0"), min_frac = 0.05, grna_targets = all_hits)
  make_trans_object("sceptre/250131_BMDSS_pilot_full_trans_pca_massive/ondisc", load_set("pilot_full"), min_frac = 0.1)
  make_trans_object("sceptre/250210_BMDSS_full_day2_day4_trans_pca_massive/ondisc", load_set("full_day2_day4"), min_frac = 0.1)
}
