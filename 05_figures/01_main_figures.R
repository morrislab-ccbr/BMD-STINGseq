# Main figure panels:
#   Figure 1: eBMD GWAS Manhattan plot highlighting variants within 500 kb of STING-seq targets,
#             and an ideogram of the GWAS loci and targets
#   Figure 2: QQ and volcano plots of candidate CRE cis effects (all time points combined)

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
setwd(project_dir)
system("mkdir -p figures")

library(data.table)
library(tidyverse)
library(ggrastr)
library(GenomicRanges)
library(RIdeogram)
library(parallel)
ncores = detectCores()

file_dir = paste0(project_dir, "/data/")

# ==========================================================================
# Figure 1
# ==========================================================================
MH.0 = fread(file.path(file_dir, "Morris2019/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz"), nThread = ncores)

# p-values below R's numerical limit are stored as 0; replace with the smallest non-zero value
temp_pvals = as.numeric(MH.0$P.NI) 
temp_pvals_ = temp_pvals
repl_pval = min(temp_pvals_[temp_pvals > 0])
temp_pvals_[temp_pvals == 0] = repl_pval
MH.0$P.NI_exact = MH.0$P.NI
MH.0$P.NI = temp_pvals_

# hg38 positions from 02_finemap_gwas.sh
snp_hg38 = dplyr::bind_rows(lapply(seq(1,22), function(x){
  pos_hg38 = fread(paste0(file_dir, "Morris2019/liftover/UKB-eBMD.chr", x, ".hg38.bed")) %>%
    dplyr::mutate(RSID = V4, BP_hg38 = V3) %>%
    dplyr::select(RSID, BP_hg38)
  return(pos_hg38)
}))

MH.1 = MH.0 %>%
  dplyr::left_join(snp_hg38, by = "RSID", relationship = "many-to-many")

MH.2 = MH.1 %>%
  dplyr::select(RSID, SNPID, CHR, BP_hg38, P.NI) %>%
  dplyr::filter(!is.na(BP_hg38)) %>%
  dplyr::mutate(P.NI_neglog10 = -log10(P.NI)) %>%
  dplyr::rename(rsid = RSID, snpid = SNPID, chr_hg38 = CHR, pos_hg38 = BP_hg38, gwas_p = P.NI, gwas_p_neglog10 = P.NI_neglog10)

# STING-seq targets
DG.0 = fread(file.path(file_dir, "guides/230830_hFOB-SSv1-TargetCRE.txt"))
DP.0 = fread(file.path(file_dir, "guides/230825_hFOB-SSv1-TargetTSS_preFilter.txt"))
DP.1 = fread(file.path(file_dir, "guides/230830_TSSs-preLibrary.txt"))

DT.lib <- fread(file.path(file_dir,"guides/230905_SSv3-BMD-hFOB-Library.txt"), fill = TRUE, na.strings = "")
colnames(DT.lib)[1] <- "gRNA"
tmp_gRNA_target = as.character(unlist(lapply(DT.lib$gRNA, function(x){strsplit(x, "__")[[1]][1]})))
tmp_gRNA_target = substr(tmp_gRNA_target, 0, nchar(tmp_gRNA_target)-2)
tmp_gRNA_target = ifelse(substr(tmp_gRNA_target, nchar(tmp_gRNA_target), nchar(tmp_gRNA_target)) == "-", substr(tmp_gRNA_target, 0, nchar(tmp_gRNA_target)-1), tmp_gRNA_target)
DT.lib$gRNA_target = tmp_gRNA_target

# Highlight variants (p < 0.1) within 500 kb of a target
MH.3 = MH.2 %>%
  dplyr::filter(gwas_p_neglog10 > 1) %>%
  dplyr::distinct(rsid, .keep_all = T)

DP_targets = subset(DP.0, vep_genes %in% DP.1$gene) %>%
  dplyr::select(rsid, chr, pos_hg38)

DT_targets = DG.0 %>%
  dplyr::filter(snpid %in% DT.lib$snpid) %>%
  dplyr::select(rsid, chr, pos_hg38)

target_list_window = dplyr::bind_rows(DT_targets, DP_targets) %>%
  dplyr::mutate(wstart = pos_hg38-500000, wend = pos_hg38+500000) %>%
  dplyr::arrange(chr, pos_hg38)

mh_ranges = GRanges(seqnames = paste("chr", MH.3$chr_hg38),
                    ranges = IRanges(start = MH.3$pos_hg38-1,
                                     end = MH.3$pos_hg38),
                    names = MH.3$rsid)

target_ranges = GRanges(seqnames = paste("chr", target_list_window$chr),
                    ranges = IRanges(start = target_list_window$wstart-1,
                                     end = target_list_window$wend))

intersect_window = GenomicRanges::intersect(mh_ranges, target_ranges)

intersect_df = data.frame('chr' = intersect_window@seqnames, intersect_window@ranges) %>%
  dplyr::mutate(chr = as.numeric(gsub("chr", "", as.character(chr)))) %>%
  dplyr::mutate(snp_pos = paste0(chr, "_", end))

color_snps = (MH.3 %>%
                dplyr::mutate(snp_pos = paste0(chr_hg38, "_", pos_hg38)) %>%
                dplyr::filter(snp_pos %in% intersect_df$snp_pos))$rsid

MH.4 = MH.3 %>%
  dplyr::group_by(chr_hg38) %>%
  dplyr::arrange(pos_hg38) %>%
  dplyr::summarise(chr_len = max(pos_hg38)) %>%
  dplyr::mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  dplyr::left_join(MH.3, ., by=c("chr_hg38" = "chr_hg38")) %>%
  dplyr::arrange(chr_hg38, pos_hg38) %>%
  dplyr::mutate(BPtotal = pos_hg38 + tot) %>%
  dplyr::mutate(is_highlight = ifelse(rsid %in% color_snps, "yes", "no"))

# Manhattan plot
axisdf = MH.4 %>%
  dplyr::group_by(chr_hg38) %>%
  dplyr::summarize(center=( max(BPtotal) + min(BPtotal) ) / 2 )

m1 = ggplot(MH.4, aes(x = BPtotal, y = gwas_p_neglog10)) +
  rasterise(geom_point(aes(color = as.factor(chr_hg38)), alpha = 0.8, size = 0.3, stroke = NA, color = "lightgrey"), dpi = 600) +
  rasterize(geom_point(data = subset(MH.4, is_highlight == "yes"), size = 0.6, stroke = NA, color = "skyblue1"), dpi = 600) +
  scale_x_continuous(label = axisdf$chr_hg38, breaks = axisdf$center ) +
  scale_y_continuous(limits = c(0, max(MH.4$gwas_p_neglog10)+10), expand = c(0, 0)) +
  theme_classic(base_size = 6) +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("Chromosome") + ylab("Significance (-log10 p-value)")

ggsave(m1, filename = file.path(project_dir, "figures/draft_MH_v2.pdf"), device = cairo_pdf, height = 1.5, width = 3.5)

# ---- Ideogram: significant GWAS loci (+/- 500 kb) and STING-seq targets ----
data(human_karyotype, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")
human_karyotype = human_karyotype[1:22,] %>% # autosomes only (sex chromosomes were not fine-mapped)
  dplyr::mutate(End = End + 500000)

# Lead variant per GWAS locus
GW.0 = fread(file.path(file_dir, "Morris2019/Morris2019.ST2.txt"))
GW.1 = GW.0 %>%
  dplyr::group_by(L.BIN) %>%
  dplyr::mutate(Z = β/SE) %>%
  dplyr::arrange(desc(abs(Z))) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
GW.2 = subset(MH.2, rsid %in% GW.1$RSID) %>%
  dplyr::group_by(rsid) %>%
  dplyr::arrange(desc(gwas_p_neglog10)) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(chr_hg38, pos_hg38) %>%
  as.data.table()

gwas_hits = GW.2 %>%
  dplyr::mutate(Chr = as.character(chr_hg38), Start = pos_hg38 - 1 - 500000, End = pos_hg38 + 500000, Value = 1) %>%
  dplyr::select(Chr, Start, End, Value)

# RIdeogram's example markers (Random_RNAs_500) are appended so that the marker labels render,
# and the list is truncated to the first 385 rows
screen_targets = dplyr::bind_rows(DT_targets %>%
                                    dplyr::mutate(Type = "Candidate CRE targets", Shape = "circle", color = "ff3030"),
                                  DP_targets %>%
                                    dplyr::mutate(Type = "Direct TSS targets", Shape = "triangle", color = "5D3FD3")) %>%
  dplyr::mutate(Chr = as.character(chr), Start = pos_hg38 - 1, End = pos_hg38) %>%
  dplyr::select(Type, Shape, Chr, Start, End, color) %>%
  dplyr::bind_rows(Random_RNAs_500) %>%
  dplyr::filter(Type != "miRNA", Type != "tRNA") %>%
  head(385)

setwd(paste0(project_dir, "/figures/"))
ideogram(karyotype = human_karyotype, overlaid = gwas_hits, label = screen_targets, label_type = "marker", colorset1 = c("white", "dodgerblue1"))
convertSVG("chromosome.svg", device = "png")

# ==========================================================================
# Figure 2
# ==========================================================================
DC.dA.0 = fread(file.path(project_dir, "summaries/241211/241211_BMDSS-dA-cis.txt"))

nt.0 = readRDS(file.path(project_dir, "sceptre/240831_BMDSS_pilot_full_cis_pca25/results_run_calibration_check.rds"))

# Non-targeting gRNA results downsampled to the number of CRE tests
set.seed(1234); nt.1 = nt.0 %>%
  dplyr::sample_n(nrow(DC.dA.0)) %>%
  dplyr::arrange(p_value) %>%
  dplyr::mutate(q_value = p.adjust(p_value, method = "BH"),
                significant = ifelse(q_value < 0.05, "TRUE", "FALSE"),
                obs.p = -log10(p_value),
                exp.p = -log10(1:length(p_value)/length(p_value)))
qq.0 = DC.dA.0 %>%
  dplyr::arrange(p_value) %>%
  dplyr::mutate(obs.p = -log10(p_value),
                exp.p = -log10(1:length(p_value)/length(p_value))) %>%
  dplyr::mutate(obs.p = ifelse(obs.p > 16, 16, obs.p)) # cap at p = 1e-16 for display

gg.qq = ggplot() +
  theme_classic(base_size = 6) +
  geom_point(data = nt.1, mapping = aes(x = exp.p, y = obs.p), color = "lightgrey", stroke = 0, alpha = 0.4, size = 0.8) +
  geom_point(data = qq.0, mapping = aes(x = exp.p, y = obs.p), color = "dodgerblue1", stroke = 0, alpha = 0.8, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.1) +
  xlab("Expected p (-log10)") + ylab("Observed p (-log10)") +
  theme(axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", size = 0.1),
        axis.text = element_text(color = "black"))

gg.volcano = ggplot() +
  theme_classic(base_size = 6) + 
  geom_point(data = subset(qq.0, cre_class == "no effect"), aes(x = log_2_fold_change, y = obs.p), color = "lightgrey", stroke = 0, alpha = 0.4, size = 0.8) +
  geom_point(data = subset(qq.0, cre_class == "eCRE"), aes(x = log_2_fold_change, y = obs.p), color = "dodgerblue1", stroke = 0, alpha = 0.8, size = 0.8) +
  geom_point(data = subset(qq.0, cre_class == "pCRE"), aes(x = log_2_fold_change, y = obs.p), color = "dodgerblue4", stroke = 0, alpha = 0.8, size = 0.8) +
  xlim(c(-1.7, 1.7)) +
  xlab("Effect size (log2 fold change)") + ylab("Significance (-log10 p)") +
  theme(axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", size = 0.1),
        axis.text = element_text(color = "black"))

ggsave(gg.qq, filename = file.path(project_dir, "figures/draft_QQ.pdf"), height = 1.5, width = 1.5, device = cairo_pdf)
ggsave(gg.volcano, filename = file.path(project_dir, "figures/draft_Volcano.pdf"), height = 1.5, width = 1.5, device = cairo_pdf)
