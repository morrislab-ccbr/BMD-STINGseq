# Summarise hFOB short-read (Kallisto) and long-read (SQANTI3) RNA-seq across osteoblast
# differentiation (days 0, 2, 4, 10) and select the isoforms used to define target TSSs.

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
setwd(project_dir)

library(tidyverse)
library(data.table)
library(biomaRt)
library(patchwork)

timepoints <- list("0" = c("a", "c"), "2" = c("a", "b", "c"), "4" = c("a", "b", "c"), "10" = c("a", "b", "c"))

# ---- Short-read RNA-seq: read Kallisto isoform estimates per replicate ----
kallisto <- list()
for (tp in names(timepoints)) {
  for (rep in timepoints[[tp]]) {
    id <- paste0("t", tp, rep)
    kallisto[[id]] <- fread(paste0("data/hFOB_RNA/Kallisto_results/", id, "_iso"))
    colnames(kallisto[[id]]) <- c("ensTranscript", id)
  }
}
join_all <- function(tables) Reduce(function(a, b) dplyr::full_join(a, b, by = "ensTranscript"), tables)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 100)
biomart_transcripts <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "ensembl_gene_id_version",
                                            "ensembl_transcript_id", "ensembl_transcript_id_version",
                                            "chromosome_name", "start_position", "end_position", "strand",
                                            "transcript_start", "transcript_end", "transcription_start_site", "transcript_length", "transcript_version",
                                            "gene_biotype", "transcript_biotype", "description"),
                             filters = "ensembl_transcript_id_version",
                             values = unique(join_all(kallisto)$ensTranscript),
                             mart = ensembl)
transcript_anno <- biomart_transcripts %>% dplyr::rename(ensTranscript = ensembl_transcript_id_version)

# ---- Isoform table with replicate means per time point ----
RA.0 <- join_all(lapply(names(timepoints), function(tp) {
  reps <- paste0("t", tp, timepoints[[tp]])
  out <- join_all(kallisto[reps])
  out[[paste0("t", tp)]] <- rowMeans(as.data.frame(out)[, reps, drop = FALSE], na.rm = TRUE)
  out
})) %>% dplyr::ungroup()

RA.2 <- RA.0 %>%
  dplyr::filter(t0 != 0 | t2 != 0 | t4 != 0 | t10 != 0) %>%
  dplyr::filter(!grepl("_", ensTranscript)) %>%
  dplyr::left_join(transcript_anno, by = "ensTranscript")
fwrite(RA.2, file = "data/hFOB_RNA/hFOB_RNA_Isoforms.txt", col.names = T, row.names = F, quote = T, sep = "\t")

# ---- Gene- and isoform-level expression per time point ----
# Gene expression = mean across replicates of the summed isoform counts
collapse_timepoint <- function(tp) {
  tcols <- paste0("t", tp, timepoints[[tp]])
  gcols <- paste0("g", tp, timepoints[[tp]])
  out <- join_all(kallisto[tcols]) %>%
    dplyr::inner_join(transcript_anno, by = "ensTranscript") %>%
    dplyr::rename(ensGene = ensembl_gene_id) %>%
    dplyr::group_by(ensGene)
  for (i in seq_along(tcols)) out <- dplyr::mutate(out, !!gcols[i] := sum(.data[[tcols[i]]]))
  out <- as.data.table(dplyr::ungroup(out))
  out[[paste0("g", tp)]] <- Reduce(`+`, out[, ..gcols]) / length(gcols)
  out[[paste0("t", tp)]] <- Reduce(`+`, out[, ..tcols]) / length(tcols)
  out %>% dplyr::select(-dplyr::all_of(c(tcols, gcols)), -ensembl_gene_id_version, -ensembl_transcript_id, -transcript_version)
}
per_tp <- lapply(setNames(names(timepoints), names(timepoints)), collapse_timepoint)

RA.0 <- dplyr::bind_rows(lapply(names(per_tp), function(tp) per_tp[[tp]] %>% dplyr::select(-paste0("g", tp), -paste0("t", tp)))) %>%
  dplyr::distinct()
for (tp in names(per_tp)) {
  RA.0 <- dplyr::left_join(RA.0, per_tp[[tp]] %>% dplyr::select(ensTranscript, paste0("g", tp), paste0("t", tp)), by = "ensTranscript")
}

# Direction of change between consecutive time points
trend <- function(a, b) ifelse(b == a, "o", ifelse(b > a, "+", "-"))
RA.1 <- RA.0
RA.1$t_pattern <- with(RA.0, paste0(trend(t0, t2), trend(t2, t4), trend(t4, t10)))
RA.1$g_pattern <- with(RA.0, paste0(trend(g0, g2), trend(g2, g4), trend(g4, g10)))

RA.2 <- RA.1 %>%
  dplyr::filter(!(t0 == 0 & t2 == 0 & t4 == 0 & t10 == 0))
fwrite(RA.2, file = "data/hFOB_RNA/hFOB_RNA.isoforms.txt", col.names = T, row.names = F, quote = T, sep = "\t")

# ---- Top isoforms per gene (those covering 80% of expression) at each time point ----
to_tss_bed <- function(x) {
  x %>%
    dplyr::mutate(chr = paste0("chr", chromosome_name), start = transcription_start_site - 1, end = transcription_start_site) %>%
    dplyr::select(chr, start, end, ensTranscript)
}
qc_plots <- list()
fills <- c("0" = "dodgerblue1", "2" = "dodgerblue2", "4" = "dodgerblue3", "10" = "dodgerblue4")
for (tp in names(timepoints)) {
  g <- paste0("g", tp); t <- paste0("t", tp)
  top <- RA.2 %>%
    dplyr::filter(.data[[g]] > 0) %>%
    dplyr::group_by(ensGene) %>%
    dplyr::arrange(desc(.data[[t]])) %>%
    dplyr::mutate(n_transcripts = length(ensTranscript)) %>%
    dplyr::mutate(sum_transcripts = sum(.data[[t]])) %>%
    dplyr::filter(sum_transcripts > 0) %>%
    dplyr::mutate(frac_transcripts = .data[[t]] / sum_transcripts) %>%
    dplyr::mutate(cumsum_transcripts = cumsum(frac_transcripts),
                  cumsum_lag = lag(cumsum_transcripts, default = 0)) %>%
    ungroup() %>% as.data.table() %>%
    dplyr::filter(cumsum_lag <= 0.8)
  fwrite(top, file = paste0("data/hFOB_RNA/hFOB_RNA.isoforms.top_t", tp, ".hg38.txt"))
  fwrite(to_tss_bed(subset(top, transcript_biotype == "protein_coding")), file = paste0("data/hFOB_RNA/hFOB_RNA.isoforms.top_t", tp, ".protein.hg38.bed"))
  fwrite(to_tss_bed(subset(top, transcript_biotype %in% c("lncRNA", "miRNA"))), file = paste0("data/hFOB_RNA/hFOB_RNA.isoforms.top_t", tp, ".ncRNA.hg38.bed"))

  qc_plots[[tp]] <- subset(top, transcript_biotype %in% c("protein_coding", "lncRNA", "miRNA")) %>%
    ggplot(aes(x = log2(.data[[t]]))) +
    theme_classic() +
    geom_histogram(bins = 30, color = "black", fill = fills[[tp]]) +
    ylab("Isoforms") + xlab("Expression (log2)") +
    facet_grid(. ~ transcript_biotype) +
    ggtitle(paste("Day", tp))
}
ggsave(wrap_plots(qc_plots, ncol = 1), filename = "data/plots/hFOB_RNA_postQC.png", device = "png", height = 8, width = 6)

# ---- Long-read RNA-seq: isoforms covering 80% of each gene's expression ----
DT.0 <- fread("data/hFOB_RNA/Long_read_RNAseq_full/SQANTI3_results_full_classification.txt")
DT.1 <- fread("data/hFOB_RNA/counts_cupcake.mapped_fl_count.txt")
GP.0 <- fread("data/hFOB_RNA/SQANTI_output/SQANTI3_results_full_corrected.genePred")

GP.1 <- GP.0 %>%
  dplyr::mutate(isoform = V1, chr = V2, start = V4, end = V5, strand = V3) %>%
  dplyr::mutate(tss = ifelse(V3 == "-", V5, V4)) %>%
  dplyr::select(isoform, chr, start, end, strand, tss)

colnames(DT.1)[1] <- "isoform"
colnames(DT.1) <- sub("_3p--NEB_5p$", "", colnames(DT.1))

# Keep isoforms detected in every replicate of at least one time point
DT.2 <- DT.1 %>%
  dplyr::filter((t0a > 0 & t0c > 0) |
                (t2a > 0 & t2b > 0 & t2c > 0) |
                (t4a > 0 & t4b > 0 & t4c > 0) |
                (t10a > 0 & t10b > 0 & t10c > 0)) %>%
  dplyr::mutate(t0 = (t0a + t0c) / 2,
                t2 = (t2a + t2b + t2c) / 3,
                t4 = (t4a + t4b + t4c) / 3,
                t10 = (t10a + t10b + t10c) / 3)

DT.3 <- DT.2 %>%
  dplyr::select(isoform, t0, t2, t4, t10) %>%
  dplyr::left_join(DT.0, by = "isoform") %>%
  dplyr::filter(chrom %in% paste0("chr", c(1:22, "Y", "X")))

DT.4 <- DT.3 %>% dplyr::group_by(associated_gene)
for (t in paste0("t", names(timepoints))) DT.4 <- dplyr::mutate(DT.4, !!paste0(t, "_sum") := sum(.data[[t]]))
for (t in paste0("t", names(timepoints))) DT.4 <- dplyr::mutate(DT.4, !!paste0(t, "_prop") := .data[[t]] / .data[[paste0(t, "_sum")]])
for (tp in names(timepoints)) {
  t <- paste0("t", tp)
  DT.4 <- DT.4 %>%
    dplyr::arrange(desc(.data[[paste0(t, "_prop")]])) %>%
    dplyr::mutate(!!paste0(t, "_runsum") := cumsum(.data[[paste0(t, "_prop")]]),
                  !!paste0(t, "_lag") := lag(.data[[paste0(t, "_runsum")]], default = 0))
}
DT.4 <- DT.4 %>% dplyr::ungroup() %>% as.data.table()

DT.5 <- DT.4 %>%
  dplyr::select(-strand) %>%
  dplyr::group_by(associated_gene) %>%
  dplyr::filter(t0_lag <= 0.8 | t2_lag <= 0.8 | t4_lag <= 0.8 | t10_lag <= 0.8) %>%
  dplyr::ungroup() %>% as.data.table() %>%
  dplyr::left_join(GP.1, by = "isoform")

system("mkdir -p data/hFOB_RNA/longreads_process")
fwrite(DT.5, file = "data/hFOB_RNA/longreads_process/hFOB_longReads_Top80p.txt", col.names = T, row.names = F, quote = F, sep = "\t")
for (tp in names(timepoints)) {
  t <- paste0("t", tp)
  lr_tp <- DT.5 %>%
    dplyr::select(associated_gene, associated_transcript, isoform, chr, start, end, strand, tss, dplyr::all_of(c(t, paste0(t, "_prop")))) %>%
    dplyr::filter(.data[[t]] >= 1) %>%
    dplyr::arrange(desc(.data[[t]]))
  lr_bed <- lr_tp %>%
    dplyr::mutate(tss_start = tss - 1) %>%
    dplyr::select(chr, tss_start, tss, isoform) %>%
    dplyr::distinct() %>%
    dplyr::arrange(chr, tss_start)
  fwrite(lr_tp, file = paste0("data/hFOB_RNA/longreads_process/hFOB_longReads_Top80p_", t, ".txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  fwrite(lr_bed, file = paste0("data/hFOB_RNA/longreads_process/hFOB_longReads_Top80p_", t, ".bed"), col.names = F, row.names = F, quote = F, sep = "\t")
}
