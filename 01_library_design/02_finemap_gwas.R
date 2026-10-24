# SuSiE fine-mapping of the Morris et al. 2019 UK Biobank eBMD GWAS.
# Called by 02_finemap_gwas.sh in three steps:
#   Rscript 02_finemap_gwas.R define_loci     # define +/- 500 kb regions around each locus
#   Rscript 02_finemap_gwas.R susie <locus>   # fine-map one locus (SLURM array task)
#   Rscript 02_finemap_gwas.R collect         # combine credible sets across loci

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
setwd(project_dir)

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

# BOLT-LMM p-values below R's numerical limit are stored as 0; replace with the smallest non-zero value
read_gwas <- function() {
  gwas <- fread("data/Morris2019/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz")
  pvals <- as.numeric(gwas$P.NI)
  pvals_fixed <- pvals
  pvals_fixed[pvals == 0] <- min(pvals[pvals > 0])
  gwas$P.NI_exact <- gwas$P.NI
  gwas$P.NI <- pvals_fixed
  gwas
}

# ---- Step 1: locus definitions ----
define_loci <- function() {
  GW.0 <- read_gwas()
  GW.0$Z <- with(GW.0, BETA/SE)

  # Independent autosomal signals (Supplementary Table 2 of Morris et al. 2019)
  DT.3 <- fread("data/Morris2019/Morris2019.ST2.txt") %>%
    dplyr::filter(CHR != 23) %>%
    dplyr::select(L.BIN, RSID) %>%
    dplyr::left_join(GW.0, by = "RSID") %>%
    dplyr::arrange(CHR, BP) %>%
    dplyr::group_by(L.BIN) %>%
    dplyr::mutate(LOWER = min(BP)-500000,
                  LOWER = ifelse(LOWER < 0, 1, LOWER),
                  UPPER = max(BP)+500000,
                  RANGE = paste0(LOWER, ".", UPPER)) %>%
    dplyr::arrange(P.NI, desc(abs(Z))) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(L.BIN, .keep_all = T) %>%
    dplyr::select(L.BIN, RSID, CHR, LOWER, UPPER, RANGE) %>%
    dplyr::arrange(L.BIN, CHR, LOWER)

  fwrite(DT.3, file = "data/Morris2019/UKB-eBMD-SigRegions.txt", col.names = F, row.names = F, quote = F, sep = "\t")
}

# ---- Step 2: SuSiE for one locus, using the LD matrix computed by plink ----
run_susie <- function(locus) {
  suppressPackageStartupMessages({
    library(coloc)
    library(susieR)
    library(patchwork)
    library(see)
    library(ggrepel)
  })

  maha <- read_gwas()
  sig_regions <- fread("data/Morris2019/UKB-eBMD-SigRegions.txt")
  colnames(sig_regions) <- c("locus", "RSID", "CHR", "lower", "upper", "region_rounded")

  CHR_var <- sig_regions$CHR[locus]
  LOWER_var <- sig_regions$lower[locus]
  UPPER_var <- sig_regions$upper[locus]
  REGION_var <- paste(LOWER_var, UPPER_var, sep = ".")

  snp_region <- maha %>%
    dplyr::filter(CHR == CHR_var, BP >= LOWER_var, BP <= UPPER_var)

  ld_prefix <- paste0("data/Morris2019/plink/", CHR_var, "/", REGION_var, "/", REGION_var)
  LD <- fread(paste0(ld_prefix, ".ld"))
  BIM <- fread(paste0(ld_prefix, ".bim"))
  colnames(BIM)[1:6] <- c("chr", "rsid", "dk", "pos", "alt", "ref")
  BIM$SNP <- paste(BIM$chr, BIM$pos, BIM$alt, BIM$ref, sep = ":")

  # Remove SNPs missing from the summary statistics, duplicates and NAs
  colnames(LD) <- BIM$SNP
  LD$SNP <- BIM$SNP
  LD <- LD[apply(LD, 1, function(x) sum(is.na(x))) < nrow(LD)]
  LD <- LD[, colnames(LD) %in% c(LD$SNP, "SNP"), with = F]
  LD <- LD[complete.cases(LD)]
  LD <- LD[, colnames(LD) %in% c(LD$SNP, "SNP"), with = F]
  LD <- LD[LD$SNP %in% snp_region$SNPID,]
  LD_gwas <- LD[!duplicated(LD$SNP),]
  BIM_gwas <- dplyr::distinct(BIM, SNP, .keep_all = T)

  int1 <- intersect(snp_region$SNPID, LD_gwas$SNP)
  LD_gwas2 <- LD_gwas[LD_gwas$SNP %in% int1,]
  LD_gwas2 <- LD_gwas2[, colnames(LD_gwas2) %in% int1, with = F]
  LD_gwas_coloc <- as.matrix(LD_gwas2)
  dimnames(LD_gwas_coloc)[[1]] <- dimnames(LD_gwas_coloc)[[2]]

  BIM_gwas_ <- BIM_gwas %>%
    dplyr::left_join(snp_region %>%
                       dplyr::rename(SNP = SNPID), by = "SNP") %>%
    dplyr::filter(SNP %in% colnames(LD_gwas2))

  gwas_region <- snp_region[SNPID %in% colnames(LD_gwas2)] %>%
    mutate(varbeta = SE^2, z = BETA/SE) %>%
    arrange(SNPID)

  D2 <- list(snp = gwas_region$SNPID,
             position = gwas_region$BP,
             beta = setNames(gwas_region$BETA, gwas_region$SNPID),
             varbeta = setNames(gwas_region$varbeta, gwas_region$SNPID),
             MAF = setNames(gwas_region$EAF, gwas_region$SNPID),
             N = gwas_region$N[1],
             type = "cc",
             LD = LD_gwas_coloc)

  # Run with the default L = 10; if no credible set is found, retry with L = 25
  S2 <- try(runsusie(D2))
  sets <- S2$sets
  fine_res <- data.frame()

  if (is.null(sets$cs)) {
    warning("No credible sets found")
    cat("Increasing maximum number of causal variants in region from 10 to 25/n")
    S2 <- try(runsusie(D2, L=25))
    sets <- S2$sets

    if (is.null(sets$cs)) {
      warning("Still no credible sets found after raising L=25")
      pr <- as.data.frame(cbind(S2$pip, t(S2$lbf_variable))) %>%
        mutate(SNP = names(S2$pip), Chr = CHR_var, region = REGION_var, cs = 0)
      colnames(pr) <- c("pip", paste("lbf_cs",seq(1,25), sep=""), "SNP", "Chr", "region", "cs")
      fine_res <- bind_rows(fine_res, pr)
    } else {
      cat("Converged when L raised to 25/n")
      for (i in 1:length(sets$cs_index)) {
        pr <- as.data.frame(cbind(S2$pip[sets$cs[[i]]], t(S2$lbf_variable[,sets$cs[[i]]]))) %>%
          mutate(SNP = names(S2$pip[sets$cs[[i]]]), Chr = CHR_var, region = REGION_var, cs = i)
        colnames(pr) <- c("pip", paste("lbf_cs",seq(1,25), sep=""), "SNP", "Chr", "region", "cs")
        fine_res <- bind_rows(fine_res, pr)
      }
    }
  } else {
    for (i in 1:length(sets$cs_index)) {
      pr <- as.data.frame(cbind(S2$pip[sets$cs[[i]]],
                                t(S2$lbf_variable[,sets$cs[[i]]]))) %>%
        dplyr::mutate(SNP = names(S2$pip[sets$cs[[i]]]),
                      Chr = CHR_var,
                      region = REGION_var,
                      cs = sets$cs_index[i])
      empt <- data.frame(matrix(ncol = 15, nrow = 1))
      pr <- dplyr::bind_cols(pr, empt)
      colnames(pr) <- c("pip", paste0("lbf_cs", seq(1, 10)), "SNP", "Chr", "region", "cs", paste0("lbf_cs", seq(11, 25)))
      fine_res <- dplyr::bind_rows(fine_res, pr)
    }
  }

  # Annotate with the gene reported in the original GWAS and hg38 positions
  GW.0 <- fread("data/Morris2019/Morris2019.ST2.txt") %>%
    dplyr::select(RSID, A.GENE) %>%
    dplyr::distinct()
  pos_hg38 <- fread(paste0("data/Morris2019/liftover/UKB-eBMD.chr", CHR_var, ".hg38.bed")) %>%
    dplyr::mutate(RSID = V4, BP_hg38 = V3) %>%
    dplyr::select(RSID, BP_hg38)
  fine_res_ <- fine_res %>%
    dplyr::mutate(locus = locus) %>%
    dplyr::mutate(region = gsub("\\.", "-", region)) %>%
    dplyr::relocate(SNP, Chr, region, cs, locus) %>%
    dplyr::rename(chr = Chr)
  rownames(fine_res_) <- NULL

  DT.1 <- BIM_gwas_ %>%
    dplyr::select(-chr, -rsid, -dk, -pos, -alt, -ref) %>%
    dplyr::rename(BP_hg19 = BP) %>%
    dplyr::left_join(pos_hg38, by = "RSID") %>%
    dplyr::relocate(SNP, RSID, CHR, BP_hg19, BP_hg38) %>%
    dplyr::left_join(fine_res_, by = "SNP") %>%
    dplyr::mutate(pip = ifelse(is.na(pip), 0, pip)) %>%
    dplyr::group_by(cs) %>%
    dplyr::arrange(desc(pip)) %>%
    dplyr::mutate(in_cs = ifelse(cs > 0, 1, 0),
                  top_variant = ifelse(cs > 0 & row_number() == 1, 1, 0),
                  top_label = ifelse(top_variant == 1, RSID, "")) %>%
    dplyr::left_join(GW.0, by = "RSID") %>%
    as.data.table() %>%
    dplyr::arrange(BP_hg38)

  # Two non-significant SNPs at locus 30 have incorrect hg38 positions after liftover
  if (locus == 30) DT.1 <- subset(DT.1, BP_hg38 > 121039511)

  fwrite(x = DT.1, file = paste0("data/Morris2019/susie/eBMD_SuSiE_locus", locus, ".txt"),
         sep = "\t", row.names = F, quote = F, col.names = T)

  # Locus plot: GWAS significance and SuSiE PIPs
  lead_snp <- sig_regions$RSID[locus]
  lead_snp_gene <- subset(GW.0, RSID == lead_snp)$A.GENE
  n_cs <- length(unique(subset(DT.1, cs > 0)$cs))
  total_pip <- sum(DT.1$pip)

  set.seed(1234); g1 <- ggplot(DT.1, aes(x = BP_hg38/1000000, y = -log10(P.NI), label = top_label)) +
    theme_classic() +
    geom_point_borderless(alpha = 0.2) +
    geom_text_repel(min.segment.length = 0, box.padding = 1, force = 100, segment.alpha = 0.8, segment.size = 0.25, segment.color = "darkgrey") +
    geom_point_borderless(data = subset(DT.1, in_cs == 1), mapping = aes(x = BP_hg38/1000000, y= -log10(P.NI), color = as.factor(cs)), size = 2, alpha = 0.5) +
    ylab("Significance (-log10 p)") +
    ggtitle(paste0(lead_snp, " (locus ", locus, ", ", lead_snp_gene, ")")) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(color = "Credible set") + geom_hline(yintercept = -log10(6.6e-9), color = "red", linetype = "dashed", size = 0.25)

  set.seed(1234); g2 <- ggplot(DT.1, aes(x = BP_hg38/1000000, y = pip, label = top_label)) +
    theme_classic() +
    geom_point_borderless(alpha = 0.2) +
    geom_text_repel(min.segment.length = 0, box.padding = 1, force = 10, segment.alpha = 0.8, segment.size = 0.25, segment.color = "darkgrey") +
    geom_point_borderless(data = subset(DT.1, in_cs == 1), mapping = aes(x = BP_hg38/1000000, y= pip, color = as.factor(cs)), size = 2, alpha = 0.5) +
    xlab(paste0("Chromosome ", CHR_var, ", hg38 position (Mb)")) + ylab("SuSiE PIP") +
    ggtitle(paste0(n_cs, " credible sets (PIP sum = ", signif(total_pip, digits = 2), ")")) +
    labs(color = "Credible set") +
    geom_hline(yintercept = 0.01, color = "red", linetype = "dashed", size = 0.25)

  ggsave(plot = g1 / g2, filename = paste0("data/Morris2019/susie_plots/eBMD_SuSiE_locus", locus, ".png"), device = png, height = 6, width = 7)
}

# ---- Step 3: variants in 95% credible sets across all loci ----
collect <- function() {
  SU.1 <- fread("data/Morris2019/UKB-eBMD-SuSiE.230707.txt") %>%
    dplyr::distinct(RSID, .keep_all = T) %>% # one multi-allelic SNP (rs10755214)
    dplyr::filter(in_cs == 1) %>%
    dplyr::select(-chr) %>%
    dplyr::rename(snpid = SNP, rsid = RSID, chr = CHR, pos_hg19 = BP_hg19, pos_hg38 = BP_hg38, ea = EA, nea = NEA, eaf = EAF, beta = BETA, beta_se = SE, p = P.NI, p_exact = P.NI_exact,
                  susie_cs = cs, susie_locus = locus, susie_pip = pip) %>%
    dplyr::select(snpid, rsid, chr, pos_hg19, pos_hg38, ea, nea, eaf, beta, beta_se, p, p_exact, susie_locus, susie_cs, susie_pip) %>%
    dplyr::group_by(susie_locus, susie_cs) %>%
    dplyr::arrange(desc(susie_pip)) %>%
    dplyr::mutate(susie_rank = row_number(), susie_n = length(unique(snpid)), max_pip = max(susie_pip)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(susie_locus, susie_cs, susie_rank)

  fwrite(SU.1, file = "data/UKB-eBMD-SuSiE-SumStats.230707.txt", col.names = T, row.names = F, quote = F, sep = "\t")
}

args <- commandArgs(TRUE)
switch(args[1],
       define_loci = define_loci(),
       susie = run_susie(as.numeric(args[2])),
       collect = collect(),
       stop("Usage: Rscript 02_finemap_gwas.R define_loci | susie <locus> | collect"))
