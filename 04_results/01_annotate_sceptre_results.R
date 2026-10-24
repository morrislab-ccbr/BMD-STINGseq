# Annotate SCEPTRE results:
#   - cis results for candidate CREs and TSSs at each time point (summaries/250721_*)
#   - positive-control summary figure
#   - per-time-point cis summaries used for Figure 2 and the gRNA list for trans testing
#   - transcriptome-wide (trans) results across all time points (summaries/250219_*)

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
setwd(project_dir)
system("mkdir -p summaries plots")

library(tidyverse)
library(data.table)
library(patchwork)
library(arrow)
library(ggbeeswarm)

# ---- Library design and target annotation ----
# gRNA library; gRNA target = gRNA ID without the -1/-2/-3 suffix
DT.lib <- fread("data/guides/230905_SSv3-BMD-hFOB-Library.txt", fill = TRUE, na.strings = "")
colnames(DT.lib)[1] <- "gRNA"
tmp_gRNA_target = as.character(unlist(lapply(DT.lib$gRNA, function(x){strsplit(x, "__")[[1]][1]})))
tmp_gRNA_target = substr(tmp_gRNA_target, 0, nchar(tmp_gRNA_target)-2)
tmp_gRNA_target = ifelse(substr(tmp_gRNA_target, nchar(tmp_gRNA_target), nchar(tmp_gRNA_target)) == "-", substr(tmp_gRNA_target, 0, nchar(tmp_gRNA_target)-1), tmp_gRNA_target)
DT.lib$gRNA_target = tmp_gRNA_target

# Candidate CRE variants, with the summed PIP of each variant and the redundant variants
# (within 1 kb) it represents, per credible set
DG.0 = fread("data/guides/230830_hFOB-SSv1-TargetCRE.txt")
DG.1 = fread("data/guides/230830_SNPs-preLibrary.txt")
DP.0 = fread("data/guides/230825_hFOB-SSv1-TargetTSS_preFilter.txt")
DP.1 = fread("data/guides/230830_TSSs-preLibrary.txt")
redundant_snps = unique(DG.1$redundant_variants)[-1]
redundant_snps_pipsums = as.data.table(do.call(rbind, lapply(redundant_snps, function(x){
  var1 = strsplit(x, ",")[[1]]
  var2 = subset(DG.1, redundant_variants == x)$snpid
  tmp1 = dplyr::bind_rows(subset(DG.0, snpid %in% var1), subset(DG.0, snpid %in% var2)) %>%
    dplyr::select(snpid, susie_locus, susie_cs, susie_pip)
  tmp2 = tmp1 %>%
    dplyr::group_by(susie_locus, susie_cs) %>%
    dplyr::mutate(pipsum = sum(susie_pip),
                  pipsum_n = length(snpid)) %>%
    dplyr::ungroup() %>%
    dplyr::select(susie_locus, susie_cs, pipsum, pipsum_n) %>%
    dplyr::arrange(susie_locus, susie_cs) %>%
    dplyr::distinct()
  tmp2$redundant_variants = x
  out = tmp2
  return(out)})))
DG.2 = DG.1 %>%
  dplyr::left_join(redundant_snps_pipsums, by = "redundant_variants") %>%
  dplyr::mutate(susie_locus.y = ifelse(redundant_variants == "", susie_locus.x, susie_locus.y),
                susie_cs.y = ifelse(redundant_variants == "", susie_cs.x, susie_cs.y),
                pipsum = ifelse(redundant_variants == "", susie_pip, pipsum),
                pipsum_n = ifelse(redundant_variants == "", 1, pipsum_n)) %>%
  dplyr::distinct()
DG.3 = DG.2 %>%
  dplyr::select(snpid, rsid, chr, pos_hg38, p, susie_locus.x, susie_cs.x, n_cs, susie_pip, pipsum, susie_cs.y, pipsum_n, ATACp_d0, ATACs_d0_dist, ATACp_d2, ATACs_d2_dist, ATACp_d4, ATACs_d4_dist) %>%
  dplyr::group_by(snpid) %>%
  dplyr::arrange(susie_cs.y) %>%
  dplyr::mutate(cs_list = paste(susie_cs.y, collapse = ","),
                pipsum_list = paste(signif(pipsum, 2), collapse = ","),
                pipsum_n_list = paste(pipsum_n, collapse = ",")) %>%
  dplyr::mutate(pipmax = max(pipsum)) %>%
  dplyr::rename(locus_main = susie_locus.x, cs_main = susie_cs.x, cs_n = n_cs) %>%
  dplyr::relocate(snpid, rsid, chr, pos_hg38, p, locus_main, cs_main, cs_n, susie_pip, pipmax, cs_list, pipsum_n_list, pipsum_list) %>%
  dplyr::select(-susie_cs.y, -pipsum, -pipsum_n) %>%
  dplyr::filter(row_number() == 1) %>%
  as.data.table()
DG.4 = DG.3 %>%
  dplyr::left_join(DT.lib %>%
                     dplyr::select(snpid, gRNA_target) %>%
                     dplyr::distinct(), by = "snpid") %>%
  dplyr::filter(!is.na(gRNA_target))

# TSS targets and the protein-altering variants that nominated them
DP.2 = DP.1 %>%
  dplyr::select(gene_id, gene) %>%
  dplyr::rename(grna_target = gene_id)
DP.3 = DP.0 %>%
  dplyr::select(snpid, rsid, chr, pos_hg38, p, susie_locus, susie_cs, susie_pip, susie_n, vep_genes, ATACp_d0, ATACp_d2, ATACp_d4, ATACs_d0_dist, ATACs_d2_dist, ATACs_d4_dist) %>%
  dplyr::rename(snp_chr = chr, snp_pos = pos_hg38, gwas_p = p, gwas_locus = susie_locus, cs = susie_cs, cs_n = susie_n, gene = vep_genes)
DP.4 = DP.2 %>%
  dplyr::left_join(DP.3, by = "gene")
snp_list = DG.4 %>%
  dplyr::rename(snp_chr = chr, snp_pos = pos_hg38, gwas_p = p, gwas_locus = locus_main, cs = cs_main, grna_target = gRNA_target) %>%
  dplyr::select(-pipmax) %>%
  dplyr::bind_rows(DP.4) %>%
  dplyr::rename(pc_target = gene) %>%
  dplyr::distinct()

# Gene TSS positions (Ensembl 110)
features <- fread(paste0("cellranger/BMDSS-full-highnovax-1_submitjob/outs/filtered_feature_bc_matrix/features.tsv.gz"), header = F)
colnames(features) <- c("ensGene", "gene", "misc")
features_ = subset(features, misc == "Gene Expression" & ensGene != "CRISPRi")
biomart_genes <- fread("data/genes_hg38_ensembl110.txt", header = T)
colnames(biomart_genes)[1] = "ensGene"
features_full <- features_ %>%
  dplyr::inner_join(biomart_genes, by = "ensGene") %>%
  dplyr::mutate(tss = ifelse(strand == 1, start_position, end_position))
gene_list = features_full %>%
  dplyr::select(gene, chromosome_name, tss) %>%
  dplyr::filter(chromosome_name %in% as.character(seq(1, 22))) %>%
  dplyr::mutate(chr_tss = as.numeric(chromosome_name)) %>%
  dplyr::select(-chromosome_name)

# Short- and long-read isoform expression
lr.1 = fread("data/hFOB_RNA/hFOB_RNA_Isoforms.txt")
lr.3 = fread("data/hFOB_RNA/counts_cupcake.mapped_fl_count.txt")

redundant_tss = DP.1$redundant_tss[DP.1$redundant_tss != ""]

# ---- Cis results (SCEPTRE discovery analyses) ----
cis_runs <- c(day0 = "240831_BMDSS_pilot_full_day0_cis_pca",
              day2 = "240831_BMDSS_full_day2_cis_pca",
              day4 = "240831_BMDSS_full_day4_cis_pca",
              dayD = "250210_BMDSS_full_day2_day4_cis_pca",
              dayA = "240831_BMDSS_pilot_full_cis_pca25")

# Add target annotations and keep gene-target pairs within 1 Mb
SimpleAnnotation = function(x){
  out = x %>%
    dplyr::filter(pass_qc == TRUE) %>%
    dplyr::select(-pass_qc, -significant) %>%
    dplyr::rename(gene = response_id, sceptre_p = p_value, sceptre_log2fc = log_2_fold_change) %>%
    dplyr::left_join(snp_list, by = "grna_target") %>%
    dplyr::filter(!pc_target %in% c("FAM3C", "CPED1", "ING3", "TSPAN12")) %>% # WNT16-CPED1 locus TSSs
    dplyr::left_join(gene_list, by = "gene") %>%
    dplyr::mutate(tss_dist = snp_pos - tss) %>%
    dplyr::relocate(gene, chr_tss, tss, tss_dist, pc_target) %>%
    dplyr::filter(abs(tss_dist) < 1000000)
  return(out)}

# Expected p-values and BH-adjusted significance (FDR 5%), separately for CRE and TSS targets
add_significance = function(x){
  x %>%
    dplyr::arrange(sceptre_p) %>%
    dplyr::mutate(obs.p = -log10(sceptre_p),
                  exp.p = -log10(1:length(sceptre_p)/length(sceptre_p)),
                  sceptre_q = p.adjust(sceptre_p, method = "BH"),
                  sceptre_sig = ifelse(sceptre_q < 0.05, 1, 0))
}
cis <- lapply(cis_runs, function(run) SimpleAnnotation(readRDS(file.path("sceptre", run, "results_run_discovery_analysis.rds"))))
CRE <- lapply(cis, function(x) x %>% dplyr::filter(is.na(pc_target)) %>% dplyr::select(-pc_target) %>% add_significance())
TSS <- lapply(cis, function(x) x %>% dplyr::filter(!is.na(pc_target)) %>% add_significance())

# Non-targeting gRNA results, downsampled to the number of TSS tests for QQ plots
AddExpectedPvalue_NT = function(x, y){
  set.seed(1234); out = x %>%
    dplyr::filter(pass_qc == TRUE) %>%
    dplyr::sample_n(nrow(y)) %>%
    dplyr::arrange(p_value) %>%
    dplyr::mutate(obs.p = -log10(p_value),
                  exp.p = -log10(1:length(p_value)/length(p_value)),
                  q_value = p.adjust(p_value, method = "BH"),
                  significant = ifelse(q_value < 0.05, 1, 0))
  return(out)}
NT <- lapply(names(cis_runs), function(tp) AddExpectedPvalue_NT(readRDS(file.path("sceptre", cis_runs[[tp]], "results_run_calibration_check.rds")), TSS[[tp]]))
names(NT) <- names(cis_runs)

stack_timepoints <- function(x) dplyr::bind_rows(lapply(names(x), function(tp) x[[tp]] %>% dplyr::mutate(timepoint = tp)))
PC_export <- stack_timepoints(TSS)
DC_export <- stack_timepoints(CRE)
fwrite(PC_export, file = "summaries/250721_Sceptre-Results-TSSs.txt", col.names = T, row.names = F, quote = T, sep = "\t")
fwrite(DC_export, file = "summaries/250721_Sceptre-Results-CREs.txt", col.names = T, row.names = F, quote = T, sep = "\t")

# ---- Positive controls (TSS-targeting gRNAs vs their own gene): summary figure ----
full_qqplot_pc = function(x, y, z, color_value){
  y$gene_label = ""
  y$gene_label[y$gene == y$pc_target & y$sceptre_sig == 1] = y$gene[y$gene == y$pc_target & y$sceptre_sig == 1]
  out = ggplot() +
    theme_classic(base_size = 6) +
    geom_point(data = x, mapping = aes(x = exp.p, y = obs.p), color = "darkgrey", stroke = 0, alpha = 0.8, size = 0.8) +
    geom_point(data = subset(y, gene_label == ""), mapping = aes(x = exp.p, y = obs.p), color = "lightgrey", stroke = 0, alpha = 0.8, size = 0.8) +
    geom_point(data = subset(y, gene_label != ""), mapping = aes(x = exp.p, y = obs.p), color = color_value, stroke = 0, alpha = 0.8, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.1) +
    xlab("Expected p (-log10)") + ylab("Observed p (-log10)") +
    theme(axis.line = element_line(color = "black", linewidth = 0.1),
          axis.ticks = element_line(color = "black", linewidth = 0.1),
          axis.text = element_text(color = "black")) +
    ggtitle(z)
  return(out)
}

qq_titles <- c(day0 = "Day 0", day2 = "Day 2", day4 = "Day 4", dayD = "Day 2 + 4", dayA = "All")
qq_colors <- c(day0 = "dodgerblue1", day2 = "#b2df8a", day4 = "#33a02c", dayD = "darkgreen", dayA = "#008080")
qq <- lapply(names(cis_runs), function(tp) full_qqplot_pc(NT[[tp]], TSS[[tp]], qq_titles[[tp]], qq_colors[[tp]]))
names(qq) <- names(cis_runs)

PC_summary = PC_export %>%
  dplyr::select(gene, pc_target, grna_target, sceptre_sig, timepoint, sceptre_log2fc, sceptre_p) %>%
  as.data.table() %>%
  dplyr::mutate(tss_target = ifelse(pc_target == gene, 1, 0))

# Relative expression of each targeted TSS isoform (including the redundant isoforms it represents)
lr.5 = lr.1 %>%
  dplyr::select(ensTranscript, t0, t2, t4) %>%
  dplyr::rename(isoform = ensTranscript)
lr.4 = lr.3 %>%
  dplyr::mutate(t0 = (`t0a_3p--NEB_5p` + `t0c_3p--NEB_5p`) / 2,
                t2 = (`t2a_3p--NEB_5p` + `t2b_3p--NEB_5p` + `t2c_3p--NEB_5p`) / 2,
                t4 = (`t4a_3p--NEB_5p` + `t4b_3p--NEB_5p` + `t4c_3p--NEB_5p`) / 2) %>%
  dplyr::select(id, t0, t2, t4) %>%
  dplyr::rename(isoform = id)
lr.6 = dplyr::bind_rows(lr.4, lr.5)

redundant_tss_df = as.data.frame(do.call(rbind, lapply(redundant_tss, function(x){
  out = as.data.table(unlist(strsplit(x, split = ",")))
  out$redundant_tss = x
  out2 = out %>%
    dplyr::rename(isoform = V1)
  return(out2)
}))) %>% 
  dplyr::left_join(lr.6, by = "isoform") %>%
  dplyr::group_by(redundant_tss) %>%
  dplyr::mutate(t0_sum = sum(t0),
                t2_sum = sum(t2),
                t4_sum = sum(t4)) %>%
  dplyr::select(-isoform, -t0, -t2, -t4) %>%
  as.data.table()

PC_summary_dayA_sig = PC_summary %>%
  dplyr::filter(timepoint == "dayA", gene == pc_target) %>%
  dplyr::arrange(pc_target, grna_target)

PC_anno = PC_summary_dayA_sig %>%
  dplyr::left_join(DP.1 %>% 
                     dplyr::select(gene_id, isoform, chr, tss_hg38, redundant_tss) %>%
                     dplyr::rename(grna_target = gene_id), by = "grna_target") %>%
  dplyr::left_join(lr.6, by = "isoform") %>%
  dplyr::left_join(redundant_tss_df, by = "redundant_tss") %>%
  dplyr::mutate(t0_sum = ifelse(is.na(t0_sum), 0, t0_sum),
                t2_sum = ifelse(is.na(t2_sum), 0, t2_sum),
                t4_sum = ifelse(is.na(t4_sum), 0, t4_sum)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(t0_all = t0 + t0_sum, 
                t2_all = t2 + t2_sum, 
                t4_all = t4 + t4_sum) %>%
  dplyr::mutate(t_avg = (t0_all + t2_all + t4_all) / 3) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(g0 = sum(t0_all),
                g2 = sum(t2_all),
                g4 = sum(t4_all)) %>%
  dplyr::mutate(t0_ratio = t0_all / g0,
                t2_ratio = t2_all / g2,
                t4_ratio = t4_all / g4) %>%
  dplyr::mutate(t_ratio_avg = ((t0_ratio + t2_ratio + t4_ratio) / 3) * 100) %>%
  dplyr::ungroup() %>%
  dplyr::select(-t0, -t2, -t4, -t0_sum, -t2_sum, -t4_sum, -redundant_tss) %>%
  as.data.table()

gg_pos_v2 = ggplot(PC_anno, aes(x = as.factor(sceptre_sig), y = pc_target)) +
  theme_classic(base_size = 4) +
  geom_rect(data = PC_summary_dayA_sig %>% 
              dplyr::distinct(pc_target, .keep_all = T), aes(fill = pc_target), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
  scale_fill_manual(values = rep(c("grey", "white"), 9)) +
  geom_beeswarm(data = PC_anno, aes(color = t_ratio_avg), priority = 'random', size = 1.2, stroke = 0, pch = 19, cex = 4) +
  scale_color_gradient(low = "darkgrey", high = "#d7191c", limits = c(0, 100)) +
  facet_wrap(. ~ pc_target, scales = "free_y", ncol = 1, strip.position = "left") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1),
        axis.text = element_text(color = "black"),
        legend.spacing.y = unit(0.2, 'cm')) +
  guides(fill = "none", color = guide_legend(override.aes = list(size = 1.5), reverse = TRUE))

gg_panel2 = (free(gg_pos_v2) + qq$day0 + qq$day2) / (qq$day4 + qq$dayD + qq$dayA) + plot_layout(widths = c(1, 1, 1))

ggsave(gg_panel2, filename = "plots/251217_ControlSummary.PanelA.pdf", device = cairo_pdf, height = 4, width = 6)

# ---- Per-time-point cis summaries for candidate CREs (Figure 2) and gRNAs with cis effects ----
# Raw SCEPTRE output, filtered to candidate CRE (SNP) targets passing QC
cis_raw <- lapply(cis_runs, function(run) {
  readRDS(file.path("sceptre", run, "results_run_discovery_analysis.rds")) %>%
    dplyr::filter(grepl("SNP", grna_target)) %>%
    dplyr::arrange(p_value) %>%
    dplyr::filter(pass_qc == TRUE) %>%
    dplyr::mutate(obs.p = -log10(p_value),
                  exp.p = -log10(1:length(p_value)/length(p_value)),
                  q_value = p.adjust(p_value, method = "BH"),
                  significant = ifelse(q_value < 0.05, 1, 0))
})

# Positive controls from the SCEPTRE power check, significant at FDR 5%
power_check <- lapply(cis_runs, function(run) {
  readRDS(file.path("sceptre", run, "results_run_power_check.rds")) %>%
    dplyr::filter(pass_qc == TRUE) %>%
    dplyr::arrange(p_value) %>%
    dplyr::mutate(q_value = p.adjust(p_value, method = "BH"),
                  significant = ifelse(q_value < 0.05, 1, 0))
})

features_full$tss_position = with(features_full, ifelse(strand == 1, start_position, end_position))

# Annotate variants and genes; recompute significance for pairs within 1 Mb
annotate_targets = function(x){
  x %>%
    dplyr::filter(!grepl("TSS", grna_target), pass_qc == TRUE) %>%
    dplyr::left_join(DG.4 %>% dplyr::rename(grna_target = gRNA_target), by = "grna_target") %>%
    dplyr::left_join(features_full %>% dplyr::rename(response_id = gene), by = "response_id") %>%
    dplyr::select(-description, -misc) %>%
    dplyr::mutate(dist_gene = pos_hg38 - tss_position) %>%
    dplyr::filter(abs(dist_gene) < 1000000) %>%
    dplyr::arrange(p_value) %>%
    dplyr::mutate(obs.p = -log10(p_value),
                  exp.p = -log10(1:length(p_value)/length(p_value)),
                  q_value = p.adjust(p_value, method = "BH"),
                  significant = ifelse(q_value < 0.05, 1, 0))
}
atac_cols <- list(day0 = c("ATACp_d2", "ATACs_d2_dist", "ATACp_d4", "ATACs_d4_dist"),
                  day2 = c("ATACp_d0", "ATACs_d0_dist", "ATACp_d4", "ATACs_d4_dist"),
                  day4 = c("ATACp_d2", "ATACs_d2_dist", "ATACp_d0", "ATACs_d0_dist"),
                  dayD = c("ATACp_d0", "ATACs_d0_dist"),
                  dayA = character())
cis_anno <- lapply(names(cis_raw), function(tp) annotate_targets(cis_raw[[tp]]) %>% dplyr::select(-dplyr::all_of(atac_cols[[tp]])))
names(cis_anno) <- names(cis_raw)

# gRNA targets with a significant cis effect at any time point (used to build the day 0 trans object)
grnas_hits = c(unlist(lapply(cis_anno, function(x) subset(x, significant == 1)$grna_target)),
               unlist(lapply(power_check, function(x) as.character(subset(x, significant == 1)$grna_target)))) %>%
  sort() %>% unique() %>% as.data.frame()
fwrite(grnas_hits, file = "summaries/250114_gRNA-hits.txt", col.names = F, row.names = F, quote = F)

# Nearest gene TSS by long-read RNA-seq: the most highly expressed isoform and the closest isoform
transcripts <- lapply(c(t0 = 0, t2 = 2, t4 = 4), function(tp) {
  fread(paste0("data/hFOB_RNA/longreads_process/hFOB_longReads_Top80p_t", tp, ".txt")) %>%
    dplyr::rename(ensGene = associated_gene, ensTranscript = associated_transcript, isoform_chr = chr, isoform_tss = tss) %>%
    dplyr::select(ensGene, ensTranscript, isoform, isoform_chr, isoform_tss, paste0("t", tp)) %>%
    dplyr::mutate(ensGene = gsub("\\..*","", ensGene))
})
add_isoform_tss = function(x, tp){
  t <- paste0("t", tp)
  for (which in c("top", "closest")) {
    sfx <- paste0(".", tp, "_", which)
    tr <- transcripts[[t]] %>%
      dplyr::rename(!!paste0("ensT", sfx) := ensTranscript, !!paste0("isoform", sfx) := isoform,
                    !!paste0("isoform_chr", sfx) := isoform_chr, !!paste0("isoform_tss", sfx) := isoform_tss,
                    !!paste0(t, "_", which) := !!t)
    dist <- paste0("dist_", t, "_", which)
    x <- x %>%
      dplyr::left_join(tr, by = "ensGene", relationship = "many-to-many") %>%
      dplyr::mutate(!!dist := pos_hg38 - .data[[paste0("isoform_tss", sfx)]]) %>%
      dplyr::group_by(gene_grna)
    x <- if (which == "top") dplyr::arrange(x, desc(.data[[paste0(t, "_top")]])) else dplyr::arrange(x, .data[[dist]])
    x <- x %>% dplyr::filter(row_number() == 1) %>% dplyr::ungroup()
  }
  x
}
prep_isoforms = function(x){
  x %>%
    dplyr::select(-gene_biotype, -hgnc_symbol, -tss, -pass_qc, -obs.p, -exp.p) %>%
    dplyr::mutate(gene_grna = paste0(response_id, "_", grna_target),
                  dist_gene = pos_hg38 - tss_position)
}

# CREs within 1 kb of the target gene TSS are promoter CREs (pCRE), otherwise enhancer CREs (eCRE)
annotate_cres = function(x){
  x %>%
    dplyr::group_by(grna_target) %>%
    dplyr::mutate(cre_class = "no effect",
                  cre_class = ifelse(abs(dist_gene) > 1000 & significant == 1, "eCRE", cre_class),
                  cre_class = ifelse(abs(dist_gene) <= 1000 & significant == 1, "pCRE", cre_class)) %>%
    dplyr::ungroup() %>%
    as.data.table()
}

DC.d0.4 = cis_anno$day0 %>% prep_isoforms() %>% add_isoform_tss(0) %>% dplyr::arrange(grna_target, p_value) %>% as.data.table() %>% annotate_cres()
DC.d2.4 = cis_anno$day2 %>% prep_isoforms() %>% add_isoform_tss(2) %>% dplyr::arrange(grna_target, p_value) %>% as.data.table() %>% annotate_cres()
DC.d4.4 = cis_anno$day4 %>% prep_isoforms() %>% add_isoform_tss(4) %>% dplyr::arrange(grna_target, p_value) %>% as.data.table() %>% annotate_cres()
DC.dA.4 = cis_anno$dayA %>% prep_isoforms() %>% add_isoform_tss(0) %>% add_isoform_tss(2) %>% add_isoform_tss(4) %>% dplyr::arrange(grna_target, p_value) %>% as.data.table() %>% annotate_cres()

system("mkdir -p summaries/241211")
fwrite(DC.dA.4, file = "summaries/241211/241211_BMDSS-dA-cis.txt", col.names = T, row.names = F, quote = T, sep = "\t")
fwrite(DC.d0.4, file = "summaries/241211/241211_BMDSS-d0-cis.txt", col.names = T, row.names = F, quote = T, sep = "\t")
fwrite(DC.d2.4, file = "summaries/241211/241211_BMDSS-d2-cis.txt", col.names = T, row.names = F, quote = T, sep = "\t")
fwrite(DC.d4.4, file = "summaries/241211/241211_BMDSS-d4-cis.txt", col.names = T, row.names = F, quote = T, sep = "\t")

# ---- Trans results (all time points) ----
DT.dA.0 = collect(open_dataset(paste0(project_dir, "/sceptre/250526_BMDSS_pilot_full_trans_pca_massive/trans_results/")))

GetSigNetworks = function(x){
  tmp0 = x %>%
    dplyr::distinct() %>%
    dplyr::group_by(grna_target) %>%
    dplyr::arrange(p_value) %>%
    dplyr::mutate(obs.p = -log10(p_value),
                  exp.p = -log10(1:length(p_value)/length(p_value)),
                  q_value = p.adjust(p_value, method = "BH"),
                  significant = ifelse(q_value < 0.01, 1, 0)) %>%
    dplyr::ungroup() %>%
    as.data.table()
}

DT.dA.1 = GetSigNetworks(DT.dA.0)
fwrite(DT.dA.1, file = "summaries/250219_BMDSS_dayAll_trans.txt", col.names = T, row.names = F, quote = F, sep = "\t")
