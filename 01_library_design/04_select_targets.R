# Select STING-seq targets from fine-mapped eBMD credible sets:
#   - TSSs of transcription factors carrying a protein-altering credible-set variant
#   - candidate CREs: credible-set variants in hFOB ATAC-seq peaks, away from expressed TSSs,
#     at signals near a transcription factor or miRNA gene (plus all candidate CREs at the
#     WNT16-CPED1 locus)

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
external_dir <- Sys.getenv("EXTERNAL_DIR", unset = "path/to/external_inputs")
setwd(project_dir)

library(tidyverse)
library(data.table)
library(biomaRt)

DT.0 = fread("data/UKB-eBMD-SuSiE-SumStats.230707.txt")

# Protein-altering variants (VEP MODERATE or HIGH impact: missense, stop gained/lost, etc.)
PC.0 = fread("data/vep/UKB-eBMD-SuSiE.230707.EA-NEA.vep.txt")
PC.1 = PC.0 %>%
  dplyr::filter(BIOTYPE == "protein_coding", EXON != "-", IMPACT %in% c("HIGH", "MODERATE"))

PC.2 = PC.1 %>%
  dplyr::rename(snpid = `#Uploaded_variation`,
                vep_gene = SYMBOL,
                vep_transcript = Feature) %>%
  dplyr::select(snpid, vep_gene, vep_transcript)

PC.3 = PC.2 %>%
  dplyr::group_by(snpid) %>%
  dplyr::filter(vep_gene != "-") %>%
  dplyr::mutate(vep_genes = paste(unique(vep_gene), collapse = "__")) %>%
  dplyr::distinct(snpid, .keep_all = T) %>%
  dplyr::select(snpid, vep_genes) %>%
  dplyr::ungroup()

# Transcription factors (Lambert et al. 2018, plus an additional curated TF list)
TF.0 = fread("data/human_TFs/DatabaseExtract_v_1.01.csv")
TF.1 = TF.0 %>%
  dplyr::rename(gene_id = `HGNC symbol`, tf_assessment = `TF assessment`) %>%
  dplyr::select(gene_id, tf_assessment)
TF_alt = fread(paste0(external_dir, "/STINGseq_v2/merged_interactions.txt"))

big_TFlist = c(TF.1$gene_id, TF_alt$TF) %>% unique() %>% sort()

A.0 = fread("data/hFOB_ATAC/hFOB_dayAll_peaksAll.hg38.txt")

LR.0 = fread("data/hFOB_RNA/longreads_process/hFOB_longReads_Top80p.txt")
LR.0_d0 = fread("data/gwas_anno/UKB-eBMD.lrRNA-D0.closest3.hg38.txt")
LR.0_d2 = fread("data/gwas_anno/UKB-eBMD.lrRNA-D2.closest3.hg38.txt")
LR.0_d4 = fread("data/gwas_anno/UKB-eBMD.lrRNA-D4.closest3.hg38.txt")
LR.0_d10 = fread("data/gwas_anno/UKB-eBMD.lrRNA-D10.closest3.hg38.txt")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 104) # Ensembl release used to annotate the long-read RNA-seq
transcript_ids = unique(LR.0$associated_transcript)
transcript_ids = transcript_ids[!grepl("_", transcript_ids)]
transcript_ids = transcript_ids[!grepl("novel", transcript_ids)]

biomart_transcripts <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "ensembl_gene_id_version",
                                            "ensembl_transcript_id", "ensembl_transcript_id_version",
                                            "chromosome_name", "start_position", "end_position", "strand",
                                            "transcript_start", "transcript_end", "transcription_start_site", "transcript_length", "transcript_version",
                                            "gene_biotype", "transcript_biotype", "description"), 
                             filters = "ensembl_transcript_id_version", 
                             values = transcript_ids, 
                             mart = ensembl)

lr_biomart = biomart_transcripts %>%
  dplyr::rename(lr_ensTranscript = ensembl_transcript_id_version, lr_gene = hgnc_symbol) %>%
  dplyr::select(lr_ensTranscript, lr_gene) %>%
  dplyr::distinct()

LR.1 = LR.0 %>%
  dplyr::rename(lr_ensTranscript = associated_transcript) %>%
  dplyr::left_join(lr_biomart)

# Expression tiers for long-read isoforms at each time point
LR.0_d0_key = LR.0 %>%
  dplyr::select(isoform, associated_transcript, associated_gene, t0) %>%
  dplyr::rename(lr_isoform = isoform, lr_ensTranscript = associated_transcript, lr_ensGene = associated_gene, lr_d0 = t0) %>%
  dplyr::mutate(lr_tier_d0 = 0,
                lr_tier_d0 = ifelse(lr_d0 >= 0.5, 1, lr_tier_d0),
                lr_tier_d0 = ifelse(lr_d0 >= 10, 2, lr_tier_d0),
                lr_tier_d0 = ifelse(lr_d0 >= 1000, 3, lr_tier_d0))
LR.0_d2_key = LR.0 %>%
  dplyr::select(isoform, associated_transcript, associated_gene, t2) %>%
  dplyr::rename(lr_isoform = isoform, lr_ensTranscript = associated_transcript, lr_ensGene = associated_gene, lr_d2 = t2) %>%
  dplyr::mutate(lr_tier_d2 = 0,
                lr_tier_d2 = ifelse(lr_d2 >= 0.5, 1, lr_tier_d2),
                lr_tier_d2 = ifelse(lr_d2 >= 10, 2, lr_tier_d2),
                lr_tier_d2 = ifelse(lr_d2 >= 1000, 3, lr_tier_d2))
LR.0_d4_key = LR.0 %>%
  dplyr::select(isoform, associated_transcript, associated_gene, t4) %>%
  dplyr::rename(lr_isoform = isoform, lr_ensTranscript = associated_transcript, lr_ensGene = associated_gene, lr_d4 = t4) %>%
  dplyr::mutate(lr_tier_d4 = 0,
                lr_tier_d4 = ifelse(lr_d4 >= 0.5, 1, lr_tier_d4),
                lr_tier_d4 = ifelse(lr_d4 >= 10, 2, lr_tier_d4),
                lr_tier_d4 = ifelse(lr_d4 >= 1000, 3, lr_tier_d4))
LR.0_d10_key = LR.0 %>%
  dplyr::select(isoform, associated_transcript, associated_gene, t10) %>%
  dplyr::rename(lr_isoform = isoform, lr_ensTranscript = associated_transcript, lr_ensGene = associated_gene, lr_d10 = t10) %>%
  dplyr::mutate(lr_tier_d10 = 0,
                lr_tier_d10 = ifelse(lr_d10 >= 0.5, 1, lr_tier_d10),
                lr_tier_d10 = ifelse(lr_d10 >= 10, 2, lr_tier_d10),
                lr_tier_d10 = ifelse(lr_d10 >= 1000, 3, lr_tier_d10))

# Consensus ATAC-seq peaks across time points
A.1 = A.0 %>%
  dplyr::arrange(V6, V7) %>%
  dplyr::group_by(V4) %>%
  dplyr::mutate(n_peaks = length(unique(V9)),
                n_times = sum(unique(grepl("day0", V9))) + sum(unique(grepl("day2", V9))) + sum(unique(grepl("day4", V9))) + sum(unique(grepl("day10", V9))),
                peak_dist = V7 - lag(V8),
                peak_size = V8 - V7,
                consensus_size = V3 - V2) %>%
  dplyr::ungroup() %>% as.data.table()

A.2 = A.1 %>%
  dplyr::select(-V5, -V6, -V7, -V8, -V10)
colnames(A.2)[1:5] = c("merge_chr", "merge_start", "merge_end", "ATACp_all", "ATACp_time")

# Nearest ATAC-seq peak/summit and long-read TSSs per time point
Ap.d0.0 = fread("data/gwas_anno/UKB-eBMD.ATACp-D0.closest.hg38.txt")
As.d0.0 = fread("data/gwas_anno/UKB-eBMD.ATACs-D0.closest.hg38.txt") 

Ap.d0.1 = Ap.d0.0 %>%
  dplyr::mutate(ATACp_d0_coords = paste0(V5, ":", V6, "-", V7)) %>%
  dplyr::rename(snpid = V4, ATACp_d0 = V8) %>%
  dplyr::select(snpid, ATACp_d0, ATACp_d0_coords) %>%
  dplyr::distinct(snpid, .keep_all = T)
As.d0.1 = As.d0.0 %>%
  dplyr::rename(snpid = V4, ATACs_d0_pos = V7, ATACs_d0 = V8, ATACs_d0_dist = V9) %>%
  dplyr::select(snpid, ATACs_d0, ATACs_d0_pos, ATACs_d0_dist) %>%
  dplyr::distinct(snpid, .keep_all = T)

Rl.d0.1 = LR.0_d0 %>% # three closest TSSs per variant
  dplyr::rename(snpid = V4, lr_d0_chr = V5, lr_d0_tss = V7, lr_isoform = V8, lr_d0_dist = V9) %>%
  dplyr::left_join(LR.0_d0_key, by = "lr_isoform") %>%
  dplyr::left_join(lr_biomart, by = "lr_ensTranscript") %>%
  dplyr::rename(lr_d0_symbol = lr_gene, lr_d0_isoform = lr_isoform, lr_d0_tier = lr_tier_d0, lr_d0_ensTranscript = lr_ensTranscript) %>%
  dplyr::select(snpid, lr_d0_isoform, lr_d0_ensTranscript, lr_d0_symbol, lr_d0_tss, lr_d0_dist, lr_d0_tier, lr_d0)

Ap.d2.0 = fread("data/gwas_anno/UKB-eBMD.ATACp-D2.closest.hg38.txt")
As.d2.0 = fread("data/gwas_anno/UKB-eBMD.ATACs-D2.closest.hg38.txt") 

Ap.d2.1 = Ap.d2.0 %>%
  dplyr::mutate(ATACp_d2_coords = paste0(V5, ":", V6, "-", V7)) %>%
  dplyr::rename(snpid = V4, ATACp_d2 = V8) %>%
  dplyr::select(snpid, ATACp_d2, ATACp_d2_coords) %>%
  dplyr::distinct(snpid, .keep_all = T)
As.d2.1 = As.d2.0 %>%
  dplyr::rename(snpid = V4, ATACs_d2_pos = V7, ATACs_d2 = V8, ATACs_d2_dist = V9) %>%
  dplyr::select(snpid, ATACs_d2, ATACs_d2_pos, ATACs_d2_dist) %>%
  dplyr::distinct(snpid, .keep_all = T)

Rl.d2.1 = LR.0_d2 %>%
  dplyr::rename(snpid = V4, lr_d2_chr = V5, lr_d2_tss = V7, lr_isoform = V8, lr_d2_dist = V9) %>%
  dplyr::left_join(LR.0_d2_key, by = "lr_isoform") %>%
  dplyr::left_join(lr_biomart, by = "lr_ensTranscript") %>%
  dplyr::rename(lr_d2_symbol = lr_gene, lr_d2_isoform = lr_isoform, lr_d2_tier = lr_tier_d2, lr_d2_ensTranscript = lr_ensTranscript) %>%
  dplyr::select(snpid, lr_d2_isoform, lr_d2_ensTranscript, lr_d2_symbol, lr_d2_tss, lr_d2_dist, lr_d2_tier, lr_d2)

Ap.d4.0 = fread("data/gwas_anno/UKB-eBMD.ATACp-D4.closest.hg38.txt")
As.d4.0 = fread("data/gwas_anno/UKB-eBMD.ATACs-D4.closest.hg38.txt") 

Ap.d4.1 = Ap.d4.0 %>%
  dplyr::mutate(ATACp_d4_coords = paste0(V5, ":", V6, "-", V7)) %>%
  dplyr::rename(snpid = V4, ATACp_d4 = V8) %>%
  dplyr::select(snpid, ATACp_d4, ATACp_d4_coords) %>%
  dplyr::distinct(snpid, .keep_all = T)
As.d4.1 = As.d4.0 %>%
  dplyr::rename(snpid = V4, ATACs_d4_pos = V7, ATACs_d4 = V8, ATACs_d4_dist = V9) %>%
  dplyr::select(snpid, ATACs_d4, ATACs_d4_pos, ATACs_d4_dist) %>%
  dplyr::distinct(snpid, .keep_all = T)

Rl.d4.1 = LR.0_d4 %>%
  dplyr::rename(snpid = V4, lr_d4_chr = V5, lr_d4_tss = V7, lr_isoform = V8, lr_d4_dist = V9) %>%
  dplyr::left_join(LR.0_d4_key, by = "lr_isoform") %>%
  dplyr::left_join(lr_biomart, by = "lr_ensTranscript") %>%
  dplyr::rename(lr_d4_symbol = lr_gene, lr_d4_isoform = lr_isoform, lr_d4_tier = lr_tier_d4, lr_d4_ensTranscript = lr_ensTranscript) %>%
  dplyr::select(snpid, lr_d4_isoform, lr_d4_ensTranscript, lr_d4_symbol, lr_d4_tss, lr_d4_dist, lr_d4_tier, lr_d4)

Ap.d10.0 = fread("data/gwas_anno/UKB-eBMD.ATACp-D10.closest.hg38.txt")
As.d10.0 = fread("data/gwas_anno/UKB-eBMD.ATACs-D10.closest.hg38.txt") 

Ap.d10.1 = Ap.d10.0 %>%
  dplyr::mutate(ATACp_d10_coords = paste0(V5, ":", V6, "-", V7)) %>%
  dplyr::rename(snpid = V4, ATACp_d10 = V8) %>%
  dplyr::select(snpid, ATACp_d10, ATACp_d10_coords) %>%
  dplyr::distinct(snpid, .keep_all = T)
As.d10.1 = As.d10.0 %>%
  dplyr::rename(snpid = V4, ATACs_d10_pos = V7, ATACs_d10 = V8, ATACs_d10_dist = V9) %>%
  dplyr::select(snpid, ATACs_d10, ATACs_d10_pos, ATACs_d10_dist) %>%
  dplyr::distinct(snpid, .keep_all = T)

Rl.d10.1 = LR.0_d10 %>%
  dplyr::rename(snpid = V4, lr_d10_chr = V5, lr_d10_tss = V7, lr_isoform = V8, lr_d10_dist = V9) %>%
  dplyr::left_join(LR.0_d10_key, by = "lr_isoform") %>%
  dplyr::left_join(lr_biomart, by = "lr_ensTranscript") %>%
  dplyr::rename(lr_d10_symbol = lr_gene, lr_d10_isoform = lr_isoform, lr_d10_tier = lr_tier_d10, lr_d10_ensTranscript = lr_ensTranscript) %>%
  dplyr::select(snpid, lr_d10_isoform, lr_d10_ensTranscript, lr_d10_symbol, lr_d10_tss, lr_d10_dist, lr_d10_tier, lr_d10)

# Per-variant annotation across time points
DT.1 = DT.0 %>%
  dplyr::group_by(susie_locus, susie_cs) %>%
  dplyr::mutate(any_sig = ifelse(any(p <= 6.6e-9), 1, 0)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(susie_locus) %>%
  dplyr::mutate(n_cs = length(unique(susie_cs))) %>%
  dplyr::ungroup() %>% as.data.table() %>%
  dplyr::left_join(PC.3, by = "snpid") 

DT.2 = DT.1 %>%
  dplyr::left_join(Ap.d0.1, by = "snpid") %>%
  dplyr::left_join(As.d0.1, by = "snpid") %>%
  dplyr::left_join(Rl.d0.1, by = "snpid") %>% 
  dplyr::group_by(snpid) %>%
  dplyr::mutate(tss_d0_1kb_any3 = 0, tss_d0_1kb_any3 = ifelse(any(lr_d0_tier >= 2 & lr_d0_dist <= 1000), 1, 0),
                tf_d0_any3 = 0, tf_d0_any3 = ifelse(any(lr_d0_symbol %in% big_TFlist), 1, 0),
                mir_d0_any3 = 0, mir_d0_any3 = ifelse(any(grepl("MIR", lr_d0_symbol)), 1, 0)) %>%
  dplyr::arrange(lr_d0_dist) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>% 
  dplyr::left_join(Ap.d2.1, by = "snpid") %>%
  dplyr::left_join(As.d2.1, by = "snpid") %>%
  dplyr::left_join(Rl.d2.1, by = "snpid") %>% 
  dplyr::group_by(snpid) %>%
  dplyr::mutate(tss_d2_1kb_any3 = 0, tss_d2_1kb_any3 = ifelse(any(lr_d2_tier >= 2 & lr_d2_dist <= 1000), 1, 0),
                tf_d2_any3 = 0, tf_d2_any3 = ifelse(any(lr_d2_symbol %in% big_TFlist), 1, 0),
                mir_d2_any3 = 0, mir_d2_any3 = ifelse(any(grepl("MIR", lr_d2_symbol)), 1, 0)) %>%
  dplyr::arrange(lr_d2_dist) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>% 
  dplyr::left_join(Ap.d4.1, by = "snpid") %>%
  dplyr::left_join(As.d4.1, by = "snpid") %>%
  dplyr::left_join(Rl.d4.1, by = "snpid") %>% 
  dplyr::group_by(snpid) %>%
  dplyr::mutate(tss_d4_1kb_any3 = 0, tss_d4_1kb_any3 = ifelse(any(lr_d4_tier >= 2 & lr_d4_dist <= 1000), 1, 0),
                tf_d4_any3 = 0, tf_d4_any3 = ifelse(any(lr_d4_symbol %in% big_TFlist), 1, 0),
                mir_d4_any3 = 0, mir_d4_any3 = ifelse(any(grepl("MIR", lr_d4_symbol)), 1, 0)) %>%
  dplyr::arrange(lr_d4_dist) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>% 
  dplyr::left_join(Ap.d10.1, by = "snpid") %>%
  dplyr::left_join(As.d10.1, by = "snpid") %>%
  dplyr::left_join(Rl.d10.1, by = "snpid") %>% 
  dplyr::group_by(snpid) %>%
  dplyr::mutate(tss_d10_1kb_any3 = 0, tss_d10_1kb_any3 = ifelse(any(lr_d10_tier >= 2 & lr_d10_dist <= 1000), 1, 0),
                tf_d10_any3 = 0, tf_d10_any3 = ifelse(any(lr_d10_symbol %in% TF.1$gene_id), 1, 0),
                mir_d10_any3 = 0, mir_d10_any3 = ifelse(any(grepl("MIR", lr_d10_symbol)), 1, 0)) %>%
  dplyr::arrange(lr_d10_dist) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>% as.data.table() %>%
  dplyr::distinct()

# Use the isoform ID as the symbol for novel transcripts
DT.2_ = DT.2 %>%
  dplyr::mutate(lr_d0_symbol = ifelse(lr_d0_ensTranscript == "novel", lr_d0_isoform, lr_d0_symbol),
                lr_d2_symbol = ifelse(lr_d2_ensTranscript == "novel", lr_d2_isoform, lr_d2_symbol),
                lr_d4_symbol = ifelse(lr_d4_ensTranscript == "novel", lr_d4_isoform, lr_d4_symbol),
                lr_d10_symbol = ifelse(lr_d10_ensTranscript == "novel", lr_d10_isoform, lr_d10_symbol))

# Variants with SuSiE PIP >= 0.01
DT.3 = DT.2_ %>%
  dplyr::filter(susie_pip >= 0.01) %>%
  dplyr::select(snpid, rsid, chr, pos_hg38, p, susie_locus, susie_cs, n_cs, susie_pip, susie_rank, susie_n, vep_genes, 
                lr_d0_symbol, lr_d2_symbol, lr_d4_symbol, lr_d10_symbol,
                lr_d0_dist, lr_d2_dist, lr_d4_dist, lr_d10_dist,
                tss_d0_1kb_any3, tss_d2_1kb_any3, tss_d4_1kb_any3, tss_d10_1kb_any3,
                tf_d0_any3, tf_d2_any3, tf_d4_any3, tf_d10_any3,
                mir_d0_any3, mir_d2_any3, mir_d4_any3, mir_d10_any3,
                ATACp_d0, ATACp_d2, ATACp_d4, ATACp_d10,
                ATACs_d0_dist, ATACs_d2_dist, ATACs_d4_dist, ATACs_d10_dist) %>%
  dplyr::mutate(susie_uniq = paste0(susie_locus, "__", susie_cs))

# Signals with a protein-altering variant
protein_coding = DT.3 %>%
  dplyr::filter(!is.na(vep_genes), vep_genes != "", susie_pip >= 0.01) %>%
  dplyr::select(susie_locus, susie_cs, susie_uniq, vep_genes) %>%
  dplyr::group_by(susie_uniq) %>%
  dplyr::mutate(vep_genes_signal = paste(unique(vep_genes), collapse = "__")) %>%
  dplyr::ungroup() %>% as.data.table() %>%
  dplyr::select(-vep_genes) %>%
  dplyr::distinct()

# ... where the affected gene is a transcription factor (TSS targets)
DT.3_pc_tf = DT.3 %>%
  dplyr::filter(susie_uniq %in% protein_coding$susie_uniq, susie_pip >= 0.01) %>%
  dplyr::rowwise() %>%
  dplyr::filter(c(strsplit(vep_genes, "__")[[1]])[1] %in% big_TFlist) %>%
  dplyr::ungroup() %>% as.data.table()

# Signals without protein-altering variants that have a TF or miRNA among the three nearest
# TSSs, no expressed TSS within 1 kb, and accessible chromatin at every time point
DT.4_tf_mir = DT.3 %>%
  dplyr::filter(!susie_uniq %in% protein_coding$susie_uniq, susie_pip >= 0.01) %>%
  dplyr::filter(tss_d0_1kb_any3 == 0, tss_d2_1kb_any3 == 0, tss_d4_1kb_any3 == 0, tss_d10_1kb_any3 == 0,
                tf_d0_any3 == 1 | tf_d2_any3 == 1 | tf_d4_any3 == 1 | tf_d10_any3 == 1 |
                  mir_d0_any3 == 1 | mir_d2_any3 == 1 | mir_d4_any3 == 1 | mir_d10_any3 == 1) %>%
  dplyr::filter(
    ATACs_d0_dist < 1000 | !is.na(ATACp_d0), ATACs_d2_dist < 1000 | !is.na(ATACp_d2), ATACs_d4_dist < 1000 | !is.na(ATACp_d4), ATACs_d10_dist < 1000 | !is.na(ATACp_d10)) %>%
  as.data.table() %>%
  dplyr::mutate(ATACp_time = ATACp_d0,
                ATACp_time = ifelse(is.na(ATACp_time), ATACp_d2, ATACp_time),
                ATACp_time = ifelse(is.na(ATACp_time), ATACp_d4, ATACp_time),
                ATACp_time = ifelse(is.na(ATACp_time), ATACp_d10, ATACp_time)) %>%
  dplyr::left_join(A.2, by = "ATACp_time")

# Candidate CREs: all accessible variants in those credible sets that are > 1 kb from a TSS
DT.6_tf_mir = DT.3 %>%
  dplyr::filter(susie_uniq %in% DT.4_tf_mir$susie_uniq, susie_pip >= 0.01) %>%
  dplyr::filter(
    ATACs_d0_dist < 1000 | !is.na(ATACp_d0), ATACs_d2_dist < 1000 | !is.na(ATACp_d2), ATACs_d4_dist < 1000 | !is.na(ATACp_d4), ATACs_d10_dist < 1000 | !is.na(ATACp_d10),
    tss_d0_1kb_any3 == 0, tss_d2_1kb_any3 == 0, tss_d4_1kb_any3 == 0, tss_d10_1kb_any3 == 0
    ) %>%
  dplyr::mutate(ATACp_time = ATACp_d0,
                ATACp_time = ifelse(is.na(ATACp_time), ATACp_d2, ATACp_time),
                ATACp_time = ifelse(is.na(ATACp_time), ATACp_d4, ATACp_time),
                ATACp_time = ifelse(is.na(ATACp_time), ATACp_d10, ATACp_time)) %>%
  dplyr::left_join(A.2, by = "ATACp_time")

DT.3_pc_tf_genes = LR.1 %>%
  dplyr::select(isoform, lr_ensTranscript, t0, t2, t4, t10, chr, tss, strand, lr_gene) %>%
  dplyr::rename(lr_chr = chr, lr_tss = tss, lr_strand = strand) %>%
  dplyr::filter(lr_gene %in% DT.3_pc_tf$vep_genes) %>%
  dplyr::arrange(lr_gene, desc(t0)) %>%
  dplyr::distinct() %>%
  dplyr::left_join(DT.3_pc_tf %>%
                     dplyr::rename(lr_gene = vep_genes) %>%
                     dplyr::select(lr_gene, snpid, rsid, chr, pos_hg38, p, susie_locus, susie_cs, n_cs, susie_pip, susie_rank, susie_n), by = "lr_gene") %>%
  dplyr::relocate(snpid, rsid, chr, pos_hg38, p, susie_locus, susie_cs, n_cs, susie_pip, susie_rank, susie_n) %>%
  dplyr::filter(lr_ensTranscript %in% PC.1$Feature) %>%
  dplyr::group_by(lr_gene) %>%
  dplyr::arrange(lr_tss) %>%
  dplyr::mutate(start_lag = lr_tss - lag(lr_tss)) %>%
  dplyr::ungroup() %>% as.data.table()

# All candidate CREs at the WNT16-CPED1 locus (locus 235)
DT.4_wnt16 = DT.3 %>%
  dplyr::filter(susie_locus == 235, susie_pip >= 0.01, tss_d0_1kb_any3 == 0, tss_d2_1kb_any3 == 0, tss_d4_1kb_any3 == 0, tss_d10_1kb_any3 == 0) %>%
  dplyr::filter(ATACs_d0_dist < 1000 | !is.na(ATACp_d0), ATACs_d2_dist < 1000 | !is.na(ATACp_d2), ATACs_d4_dist < 1000 | !is.na(ATACp_d4), ATACs_d10_dist < 1000 | !is.na(ATACp_d10)) %>%
  dplyr::mutate(ATACp_time = ATACp_d0,
                ATACp_time = ifelse(is.na(ATACp_time), ATACp_d2, ATACp_time),
                ATACp_time = ifelse(is.na(ATACp_time), ATACp_d4, ATACp_time),
                ATACp_time = ifelse(is.na(ATACp_time), ATACp_d10, ATACp_time)) %>%
  dplyr::left_join(A.2, by = "ATACp_time") %>%
  as.data.table()

DT.7 = dplyr::bind_rows(
  DT.6_tf_mir, DT.4_wnt16
)

# Target lists for gRNA design
system("mkdir -p data/guides")
fwrite(DT.7, file = "data/guides/230830_hFOB-SSv1-TargetCRE.txt", col.names = T, row.names = F, quote = F, sep = "\t")
fwrite(DT.3_pc_tf_genes, file = "data/guides/230825_hFOB-SSv1-TargetTSS.txt", col.names = T, row.names = F, quote = F, sep = "\t")
fwrite(DT.3_pc_tf, file = "data/guides/230825_hFOB-SSv1-TargetTSS_preFilter.txt", col.names = T, row.names = F, quote = F, sep = "\t")
