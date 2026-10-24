# Per-lane QC of the Cell Ranger count matrices: filter cells on cDNA metrics, demultiplex
# hashtag oligos (HTOs) to keep singlets, and threshold gRNA (GDO) UMIs.
# Each 10x chip lane is processed separately and lanes are combined in 03_combine_lanes.R.

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
setwd(project_dir)
system("mkdir -p seurat/plots seurat/temp")

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
})

experiments <- list(
  pilot = list(lanes = 1:2, run = function(LN) paste0("BMDSS-pilot-midhigh2x-", LN), tag = "pilot-midhigh2x",
               mt_upper = 0.99, hto_labels = seq(0.90, 0.99, 0.01)),
  full  = list(lanes = 1:6, run = function(LN) paste0("BMDSS-full-highnovax-", LN, "_submitjob"), tag = "full-highnovax",
               mt_upper = 0.95, hto_labels = seq(0.98, 0.999, 0.001))
)

# Read one modality from a Cell Ranger matrix (1 = cDNA, 2 = HTO, 3 = GDO), removing the "-1" barcode suffix
read_modality <- function(run_id, i) {
  x <- Read10X(data.dir = paste0("cellranger/", run_id, "/outs/filtered_feature_bc_matrix/"))[[i]]
  colnames(x) <- substr(colnames(x), 1, nchar(colnames(x)) - 2)
  x
}

qc_lane <- function(cfg, LN) {
  run_id <- cfg$run(LN)
  name <- paste0(LN, "-", cfg$tag)
  project <- paste0(cfg$tag, "-", LN)

  # ---- cDNA: remove cells in the outer 1% tails of genes, UMIs and mitochondrial fraction ----
  cDNA.0 <- CreateSeuratObject(counts = read_modality(run_id, 1), project = project, min.cells = 3)
  cDNA.0[["percent.mt"]] <- PercentageFeatureSet(object = cDNA.0, pattern = "^MT-")
  cDNA.1 <- subset(cDNA.0, subset = nFeature_RNA > quantile(cDNA.0@meta.data$nFeature_RNA, 0.01) & nFeature_RNA < quantile(cDNA.0@meta.data$nFeature_RNA, 0.99) &
                     nCount_RNA > quantile(cDNA.0@meta.data$nCount_RNA, 0.01) & nCount_RNA < quantile(cDNA.0@meta.data$nCount_RNA, 0.99) &
                     percent.mt > quantile(cDNA.0@meta.data$percent.mt, 0.01) & percent.mt < quantile(cDNA.0@meta.data$percent.mt, cfg$mt_upper))

  mca <- NormalizeData(cDNA.1, verbose = F) %>% FindVariableFeatures(verbose = F) %>% ScaleData(verbose = F) %>% RunPCA(verbose = F) %>% RunUMAP(dims = 1:50, verbose = F)
  p1 <- VlnPlot(object = cDNA.0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  p2 <- VlnPlot(object = cDNA.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  u1 <- FeaturePlot(mca, features = c("nCount_RNA"), reduction = "umap", pt.size = 1, min.cutoff = "q01", max.cutoff = "q99") + ggtitle(paste0("cDNA-", name, " UMAP (", ncol(cDNA.1), " cells)")) + NoAxes()
  ggsave(filename = paste0("seurat/plots/cDNA-", name, "_preQC.png"), plot = p1, height = 5, width = 8, device = "png")
  ggsave(filename = paste0("seurat/plots/cDNA-", name, "_postQC.png"), plot = p2, height = 5, width = 8, device = "png")
  ggsave(filename = paste0("seurat/plots/cDNA-", name, "_postQC_UMAP.png"), plot = u1, device = "png", height = 6, width = 7)
  save(cDNA.1, file = paste0("seurat/temp/temp_cDNA-", name, ".RData"))

  # ---- HTO: demultiplex with the positive quantile that maximizes singlets ----
  HTO.0 <- CreateSeuratObject(counts = read_modality(run_id, 2), assay = "HTO", project = project, min.cells = 3)
  HTO.1 <- HTO.0[, intersect(colnames(cDNA.1), colnames(HTO.0))]
  HTO.2 <- subset(HTO.1, subset = nCount_HTO > quantile(HTO.1@meta.data$nCount_HTO, 0.01) & nCount_HTO < quantile(HTO.1@meta.data$nCount_HTO, 0.99))
  HTO.2 <- NormalizeData(HTO.2, normalization.method = "CLR", verbose = F)

  HTO.table <- data.table('Threshold' = cfg$hto_labels, 'Singlets' = sapply(seq(0.98,0.999,0.001), function(x){
    tmp <- HTODemux(HTO.2, assay = "HTO", positive.quantile = x, nsamples = 1000, verbose = F)
    table(tmp$HTO_classification.global)[3]
  }))
  custom_thresh <- (HTO.table %>% dplyr::arrange(desc(Singlets)))[1,]$Threshold
  HTO.2 <- HTODemux(HTO.2, positive.quantile = custom_thresh, nsamples = 1000, verbose = F)
  Idents(HTO.2) <- "hash.ID"

  h1 <- ggplot(HTO.table, aes(x = Threshold, y = Singlets)) + geom_point(size = 0.8) + ggtitle(paste0("HTO-", name, " demultiplex threshold"))
  p3 <- VlnPlot(object = HTO.1, features = c("nFeature_HTO", "nCount_HTO"), ncol = 2, pt.size = 0.1)
  p4 <- VlnPlot(object = HTO.2, group.by = "orig.ident", features = c("nFeature_HTO", "nCount_HTO"), ncol = 2, pt.size = 0.1)
  p5 <- RidgePlot(HTO.2, assay = "HTO", features = rownames(HTO.2[["HTO"]]), ncol = 4)
  ggsave(filename = paste0("seurat/plots/HTO-", name, "_preQC.png"), plot = p3, device = "png", height = 5, width = 5)
  ggsave(filename = paste0("seurat/plots/HTO-", name, "_postQC.png"), plot = p4, device = "png", height = 5, width = 5)
  ggsave(filename = paste0("seurat/plots/HTO-", name, "_postQC_RidgePlot.png"), plot = p5, device = "png", height = 10, width = 16)
  ggsave(filename = paste0("seurat/plots/HTO-", name, "_postQC_DemuxThresh.png"), plot = h1, device = "png", height = 3, width = 6)

  HTO.3 <- subset(HTO.2, idents = "Negative", invert = TRUE)
  HTO.3 <- ScaleData(HTO.3, verbose = F)
  HTO.3 <- RunPCA(HTO.3, features = rownames(HTO.3), reduction.name = "HPCA", reduction.key = "HPC_", approx = F, verbose = F)
  HTO.3 <- RunUMAP(HTO.3, dims = 1:4, reduction = "HPCA", reduction.name = "HUMAP", reduction.key = "HUMAP_", verbose = F)
  u2 <- DimPlot(HTO.3, reduction = "HUMAP", group.by = "hash.ID") + ggtitle(paste0("HTO-", name, " UMAP (", ncol(HTO.3), " cells)")) + NoAxes()
  ggsave(filename = paste0("seurat/plots/HTO-", name, "_postQC_UMAP.png"), plot = u2, device = "png", height = 6, width = 7)

  Idents(HTO.3) <- "HTO_classification.global"
  HTO.4 <- subset(HTO.3, idents = "Singlet")
  save(HTO.4, file = paste0("seurat/temp/temp_HTO-", name, ".RData"))

  # ---- GDO: singlets only; zero counts below Cell Ranger's per-gRNA UMI threshold (minimum 3) ----
  GDO.0 <- CreateSeuratObject(counts = read_modality(run_id, 3), assay = "GDO", project, min.cells = 3)
  GDO.1 <- GDO.0[, intersect(colnames(HTO.4), colnames(GDO.0))]
  GDO.2 <- subset(GDO.1, subset = nCount_GDO < quantile(GDO.1@meta.data$nCount_GDO, 0.99) & nCount_GDO > quantile(GDO.1@meta.data$nCount_GDO, 0.01))
  GDO.mtx <- as.matrix(GetAssayData(GDO.2))

  GDO.thresholds <- fread(paste0("cellranger/", run_id, "/outs/crispr_analysis/protospacer_umi_thresholds.csv")) %>%
    dplyr::mutate(Protospacer = gsub("_", "-", Protospacer))
  GDO.mtx.minC <- do.call(rbind, lapply(rownames(GDO.mtx), function(x){
    if(x %in% GDO.thresholds$Protospacer) UMI_threshold = subset(GDO.thresholds, Protospacer == x)$`UMI threshold` else UMI_threshold = 3
    if(UMI_threshold < 3) UMI_threshold = 3
    GDO.mtx.sub = GDO.mtx[x, ]
    GDO.mtx.sub[GDO.mtx.sub < UMI_threshold] = 0
    return(GDO.mtx.sub)
  }))
  rownames(GDO.mtx.minC) <- rownames(GDO.mtx)
  GDO.minC <- CreateSeuratObject(counts = GDO.mtx.minC, assay = "GDO", project = paste0(project, "-minC"))

  p6 <- VlnPlot(object = GDO.1, features = c("nFeature_GDO", "nCount_GDO"), ncol = 2, pt.size = 0.1)
  p7 <- VlnPlot(object = GDO.2, features = c("nFeature_GDO", "nCount_GDO"), ncol = 2, pt.size = 0.1)
  p8 <- VlnPlot(object = GDO.minC, features = c("nFeature_GDO", "nCount_GDO"), ncol = 2, pt.size = 0.1)
  ggsave(filename = paste0("seurat/plots/GDO-", name, "_preQC.png"), plot = p6, device = "png", height = 5, width = 5)
  ggsave(filename = paste0("seurat/plots/GDO-", name, "_postQC.png"), plot = p7, device = "png", height = 5, width = 5)
  ggsave(filename = paste0("seurat/plots/GDO-", name, "_postQC_minC.png"), plot = p8, device = "png", height = 5, width = 5)
  save(GDO.2, file = paste0("seurat/temp/temp_GDO-", name, ".RData"))
  save(GDO.minC, file = paste0("seurat/temp/temp_GDOminC-", name, ".RData"))
}

for (cfg in experiments) {
  for (LN in cfg$lanes) qc_lane(cfg, LN)
}
