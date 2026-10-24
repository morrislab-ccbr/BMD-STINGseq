# Combine QC'd lanes into the cell sets analysed with SCEPTRE and compute expression PCs.
#   pilot (day 0), full day 0 / day 2 / day 4 (split by hashtag), and the combined sets
#   pilot + full day 0, full day 2 + day 4, and pilot + full (all cells).

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
setwd(project_dir)
system("mkdir -p seurat/temp/umap_objects")

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
})

modalities <- c("cDNA", "GDO", "GDOminC", "HTO")

# Load one QC'd lane, restricting cDNA and HTO to cells that passed gRNA QC
load_lane <- function(name) {
  env <- new.env()
  for (m in modalities) load(paste0("seurat/temp/temp_", m, "-", name, ".RData"), envir = env)
  list(cDNA = env$cDNA.1[, intersect(colnames(env$cDNA.1), colnames(env$GDO.2))],
       GDO = env$GDO.2,
       GDOminC = env$GDO.minC,
       HTO = env$HTO.4[, intersect(colnames(env$GDO.2), colnames(env$HTO.4))])
}

# Merge Seurat objects (one list per modality) and join count layers
merge_objects <- function(objs, ids, project) {
  JoinLayers(merge(objs[[1]], y = objs[-1], add.cell.ids = ids, project = project))
}

# Keep genes detected in at least min_cells cells
filter_genes <- function(obj, min_cells) {
  counts <- GetAssayData(object = obj, layer = "counts")
  keep_genes <- Matrix::rowSums(counts > 0) >= min_cells
  CreateSeuratObject(counts[keep_genes, ], meta.data = obj@meta.data)
}

cell_table <- function(obj) {
  tbl <- as.data.table(obj@meta.data)
  tbl$cell <- rownames(obj@meta.data)
  tbl
}

# Per-cell covariates: cDNA metrics, gRNA counts (raw and thresholded), lane, CRISPRi transgene UMIs
make_cell_data <- function(merged) {
  cell_table(merged$cDNA) %>%
    dplyr::mutate(percent.mt = percent.mt/100) %>%
    dplyr::select(cell, nCount_RNA, nFeature_RNA, percent.mt) %>%
    dplyr::left_join(cell_table(merged$GDO), by = "cell") %>%
    dplyr::left_join(cell_table(merged$GDOminC) %>%
                       dplyr::select(-orig.ident) %>%
                       dplyr::rename(nCount_GDOminC = nCount_GDO, nFeature_GDOminC = nFeature_GDO),
                     by = "cell") %>%
    dplyr::mutate(orig.ident = gsub("_.*", "", cell))
}

save_set <- function(set, cell_data, label, prefix = "seurat/temp/") {
  for (m in modalities) saveRDS(object = set[[m]], file = paste0(prefix, m, ".", label, ".RDS"))
  fwrite(x = cell_data, file = paste0(prefix, "cDNA.", label, ".cell_data.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
}

read_set <- function(label) {
  set <- lapply(setNames(modalities, modalities), function(m) readRDS(paste0("seurat/temp/", m, ".", label, ".RDS")))
  set$cell_data <- fread(paste0("seurat/temp/cDNA.", label, ".cell_data.txt"))
  set
}

# ---- Pilot experiment (day 0) ----
lanes <- lapply(1:2, function(LN) load_lane(paste0(LN, "-pilot-midhigh2x")))
pilot <- lapply(setNames(modalities, modalities), function(m) merge_objects(lapply(lanes, `[[`, m), paste0("LN-", 1:2), "pilot-midhigh2x"))
pilot$cDNA <- filter_genes(pilot$cDNA, 200)
cell_data <- make_cell_data(pilot)
cell_data$CRISPRi <- FetchData(object = pilot$cDNA, layer = "counts", "CRISPRi")$CRISPRi
save_set(pilot, cell_data, "pilot")

# ---- Full experiment, split into days 0, 2 and 4 by hashtag ----
lanes <- lapply(1:6, function(LN) load_lane(paste0(LN, "-full-highnovax")))
full <- lapply(setNames(modalities, modalities), function(m) merge_objects(lapply(lanes, `[[`, m), paste0("LN-", 1:6), "full-highnovax"))
cell_data <- make_cell_data(full) %>%
  dplyr::left_join(cell_table(full$HTO) %>%
                     dplyr::select(cell, HTO_classification),
                   by = "cell") %>%
  dplyr::mutate(timepoint = ifelse(HTO_classification %in% c("HTO1", "HTO2"), "day0", ""),
                timepoint = ifelse(HTO_classification == "HTO3", "day2", timepoint),
                timepoint = ifelse(HTO_classification == "HTO4", "day4", timepoint))
cell_data$CRISPRi <- FetchData(object = full$cDNA, layer = "counts", "CRISPRi")$CRISPRi

for (day in c("day0", "day2", "day4")) {
  day_set <- lapply(full, function(obj) obj[, cell_data$timepoint == day])
  day_set$cDNA <- filter_genes(day_set$cDNA, 200)
  save_set(day_set, subset(cell_data, timepoint == day), day)
}
rm(lanes, pilot, full, day_set); gc()

# ---- Combined cell sets ----
pilot <- read_set("pilot")
days <- lapply(setNames(c("day0", "day2", "day4"), c("day0", "day2", "day4")), read_set)
prefix_cells <- function(cd, study) {
  cd %>% dplyr::mutate(cell = paste0(study, "_", cell), orig.ident = paste0(study, "_", orig.ident), study = study)
}
combine_sets <- function(sets, ids, project) {
  lapply(setNames(modalities, modalities), function(m) merge_objects(lapply(sets, `[[`, m), ids, project))
}

# Pilot + full day 0
cell_data <- prefix_cells(pilot$cell_data, "pilot") %>%
  dplyr::bind_rows(prefix_cells(days$day0$cell_data, "full")) %>%
  dplyr::select(-HTO_classification, -timepoint)
save_set(combine_sets(list(pilot, days$day0), c("pilot", "full"), "pilot-full-day0"), cell_data, "pilot_full_day0")

# Full day 2 + day 4
cell_data <- prefix_cells(days$day2$cell_data, "full") %>%
  dplyr::bind_rows(prefix_cells(days$day4$cell_data, "full")) %>%
  dplyr::select(-HTO_classification, -timepoint)
save_set(combine_sets(list(days$day2, days$day4), c("full", "full"), "full-day2-day4"), cell_data, "full_day2_day4")

# Pilot + full, all time points
cell_data <- prefix_cells(pilot$cell_data, "pilot") %>%
  dplyr::mutate(timepoint = "day0") %>%
  dplyr::bind_rows(prefix_cells(days$day0$cell_data, "full")) %>%
  dplyr::bind_rows(prefix_cells(days$day2$cell_data, "full")) %>%
  dplyr::bind_rows(prefix_cells(days$day4$cell_data, "full")) %>%
  dplyr::select(-HTO_classification)
save_set(combine_sets(list(pilot, days$day0, days$day2, days$day4), c("pilot", "full", "full", "full"), "pilot-full"), cell_data, "pilot_full")

# ---- Expression PCs (and UMAP) for each set analysed with SCEPTRE ----
for (label in c("day2", "day4", "pilot_full_day0", "full_day2_day4", "pilot_full")) {
  obj <- readRDS(paste0("seurat/temp/cDNA.", label, ".RDS"))
  obj <- NormalizeData(obj, verbose = F) %>% FindVariableFeatures(verbose = F) %>% ScaleData(verbose = F) %>% RunPCA(verbose = F) %>% RunUMAP(dims = 1:50, verbose = F)
  saveRDS(object = obj, file = paste0("seurat/temp/umap_objects/cDNA.", label, ".RDS"))
}
