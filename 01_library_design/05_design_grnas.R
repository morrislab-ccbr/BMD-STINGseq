# Design candidate CRISPRi gRNAs for the selected target variants (candidate CREs) and TSSs.
# Targets closer than 1 kb (the CRISPRi window) are collapsed, keeping one representative and
# recording the others as redundant. Output gRNAs are scored with FlashFry (06_score_grnas.sh).

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
setwd(project_dir)

library(data.table)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)

# ---- Helper functions ----

# All 20-nt spacers with an NGG PAM on either strand of a sequence
GetAllGuidesInRegion <- function(regionToSearch, PAMpattern) {

  GetGuidesFromMatch_TopStrand <- function(matchPattern_Output){
    PAMends = end(matchPattern_Output); indsToKeep = which(PAMends >= 23)
    if (length(indsToKeep) == 0) {return(NULL)}
    PAMends_Keep = PAMends[indsToKeep]
    return(data.frame("GUIDE_START" = (PAMends_Keep - 22),
                      "GUIDE_END" = (PAMends_Keep - 3),
                      "GUIDE_STRAND" = "+",
                      "GUIDE_SEQ" = sapply(PAMends_Keep, function(x){return(matchPattern_Output %>% subject %>% as.character %>% substr(start=(x-22),stop=(x)))}),
                      "PAM" = (matchPattern_Output[indsToKeep,] %>% as.data.frame %>% as.matrix)) %>% mutate_if(is.factor, as.character) %>% mutate(PAM=x) %>% dplyr::select(-x))}

  GetGuidesFromMatch_BottomStrand <- function(matchPattern_Output){
    PAMends = end(matchPattern_Output)
    indsToKeep = which(PAMends <= ((matchPattern_Output %>% subject %>% length)-20))
    if (length(indsToKeep) == 0) {return(NULL)}; PAMends_Keep = PAMends[indsToKeep]
    return(data.frame("GUIDE_START" = (PAMends_Keep + 20),
                      "GUIDE_END" = (PAMends_Keep + 1),
                      "GUIDE_STRAND" = "-",
                      "GUIDE_SEQ" = sapply(PAMends_Keep, function(x){return(matchPattern_Output %>% subject %>% as.character %>% substr(start=(x-2),stop=(x+20)) %>% DNAString %>% reverseComplement %>% as.character)}),
                      "PAM" = (matchPattern_Output[indsToKeep,] %>% reverseComplement %>% as.data.frame %>% as.matrix)) %>% mutate_if(is.factor, as.character) %>% mutate(PAM=x) %>% dplyr::select(-x))}

  temp <- bind_rows(GetGuidesFromMatch_TopStrand(matchPattern(PAMpattern, regionToSearch, fixed=F)),
                    GetGuidesFromMatch_BottomStrand(matchPattern(reverseComplement(PAMpattern), regionToSearch, fixed=F))) %>%
    mutate(Design_PAM = as.character(PAMpattern))

  return(temp)}

# Number each run of consecutive targets that lie within 1 kb of the previous target
number_blast_radius <- function(rel_pos) {
  cumsum(c(TRUE, rel_pos[-1] > 1000))
}

# Flag targets within 1 kb of a neighbouring target (blastrad_pass = 1)
flag_blast_radius <- function(x, group_col) {
  near_prev <- x$lag_pos < 1000
  near_prev[is.na(near_prev)] <- FALSE
  flag <- near_prev + c(near_prev[-1], FALSE)
  flag[flag > 1] <- 1
  x$blastrad_pass <- flag
  x$blastrad_group <- paste0(x[[group_col]], "__", flag)
  x
}

# Iteratively resolve flagged targets: within each group keep the top-ranked target (rank_expr),
# record the next one as redundant and drop it, then re-check until no targets are within 1 kb.
collapse_blast_radius <- function(x, group_col, id_col, pos_col, rank_expr) {
  x <- flag_blast_radius(x, group_col)
  keep <- list(subset(x, blastrad_pass == 0))
  pass <- x; k <- 0
  while (any(pass$blastrad_pass == 1)) {
    redundant_col <- paste0("redundant_variant", strrep("_", k))
    pass <- subset(pass, blastrad_pass == 1) %>%
      dplyr::arrange(blastrad_group, !!rank_expr) %>%
      dplyr::group_by(blastrad_group) %>%
      dplyr::mutate(keep_ = ifelse(row_number() == 1, 1, 0)) %>%
      dplyr::mutate(!!redundant_col := ifelse(keep_ == 1, .data[[id_col]][2], "")) %>%
      dplyr::filter(row_number() != 2) %>%
      dplyr::arrange(chr, .data[[pos_col]]) %>%
      dplyr::group_by(chr) %>%
      dplyr::mutate(lag_pos = .data[[pos_col]] - lag(.data[[pos_col]])) %>%
      dplyr::ungroup() %>%
      as.data.table() %>%
      flag_blast_radius(group_col)
    keep[[length(keep) + 1]] <- subset(pass, blastrad_pass == 0)
    k <- k + 1
  }
  dplyr::bind_rows(keep) %>%
    dplyr::mutate_if(is.character, ~replace_na(.,"")) %>%
    dplyr::arrange(chr, .data[[pos_col]]) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(lag_pos = .data[[pos_col]] - lag(.data[[pos_col]])) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(redundant = paste(sort(unique(c_across(matches("^redundant_variant_*$")))), collapse = ",")) %>%
    dplyr::mutate(redundant = gsub('^\\,|\\,$', '', redundant)) %>%
    dplyr::ungroup() %>% as.data.table() %>%
    dplyr::select(-matches("^redundant_variant_*$")) %>%
    flag_blast_radius(group_col)
}

# Count NGG PAMs (either strand) in each target window
count_ngg <- function(chr, start, end) {
  region <- getSeq(Hsapiens, names = chr, start = start, end = end)
  countPattern(DNAString("NGG"), region, fixed = FALSE) + countPattern(reverseComplement(DNAString("NGG")), region, fixed = FALSE)
}

# Remove spacers with homopolymers (AAAAA, CCCCC, GGGGG), U6 terminators (TTTT) or duplicates;
# return gRNA coordinates and FlashFry input (spacer with random 10-nt flanks for the Doench 2014 score)
guide_tables <- function(guides, targets, id_col) {
  spacers <- guides$GUIDE_SEQ
  remove <- sort(unique(c(grep("AAAAA", spacers), grep("CCCCC", spacers), grep("GGGGG", spacers),
                          grep("TTTT", spacers), which(duplicated(spacers)))))
  kept <- guides[-remove,]
  positions <- kept %>%
    left_join(targets, by = id_col) %>%
    dplyr::mutate(FlashFryID = gsub(":", "-", .data[[id_col]]),
                  gRNA_chr = TARGET_CHR,
                  gRNA_start = ifelse(GUIDE_STRAND == "+", as.character(TARGET_START + GUIDE_START - 1), as.character(TARGET_START + GUIDE_END - 1)),
                  gRNA_end = ifelse(GUIDE_STRAND == "+", as.character(TARGET_START + GUIDE_END - 1), as.character(TARGET_START + GUIDE_START - 1)),
                  gRNA_strand = GUIDE_STRAND,
                  gRNA_seq_wPAM = substr(GUIDE_SEQ, 0, 20)) %>%
    dplyr::select(all_of(id_col), FlashFryID, gRNA_chr, gRNA_start, gRNA_end, gRNA_strand, gRNA_seq_wPAM)
  to_score <- kept %>%
    left_join(targets %>% dplyr::select(TARGET_CHR, all_of(id_col)), by = id_col) %>%
    dplyr::mutate(FlashFryID = gsub(":", "-", .data[[id_col]]),
                  FlashFrySEQ = paste0("TTCGTACAAA", GUIDE_SEQ, "TGTCCGCACT"))
  list(positions = positions, to_score = to_score)
}

# ---- Candidate CRE variants ----
DT.0 <- fread("data/guides/230830_hFOB-SSv1-TargetCRE.txt")

DT.1 <- DT.0 %>%
  dplyr::group_by(susie_locus) %>% # credible sets can overlap within a locus
  dplyr::arrange(pos_hg38) %>%
  dplyr::mutate(lag_pos = pos_hg38 - lag(pos_hg38)) %>%
  dplyr::mutate(rel_pos = pos_hg38 - lag(pos_hg38, default = pos_hg38[1])) %>%
  dplyr::ungroup() %>% as.data.table()
DT.1$blastrad <- paste0(DT.1$chr, "__", number_blast_radius(DT.1$rel_pos))

# Among nearby variants keep the one closest (on average) to an ATAC-seq summit
DT.2 <- DT.1 %>%
  dplyr::mutate(ATACs_avg_dist = (ATACs_d0_dist + ATACs_d2_dist + ATACs_d4_dist + ATACs_d10_dist) / 4) %>%
  collapse_blast_radius(group_col = "susie_locus", id_col = "snpid", pos_col = "pos_hg38", rank_expr = quo(ATACs_avg_dist)) %>%
  dplyr::rename(redundant_variants = redundant)

# gRNAs within 150 bp of each variant
DT.3 <- DT.2 %>%
  dplyr::select(snpid, chr, pos_hg38) %>%
  dplyr::arrange(chr, pos_hg38) %>%
  dplyr::mutate(TARGET_CHR = paste0("chr", chr),
                TARGET_START = pos_hg38-(150),
                TARGET_END = pos_hg38+(150)) %>%
  dplyr::group_by(snpid) %>%
  dplyr::mutate(NGG = count_ngg(TARGET_CHR, TARGET_START, TARGET_END)) %>%
  dplyr::filter(NGG > 0) %>%
  dplyr::ungroup() %>% as.data.table()

DT.4 <- DT.3 %>%
  dplyr::group_by(snpid) %>%
  do((GetAllGuidesInRegion(regionToSearch = getSeq(Hsapiens, names=.$TARGET_CHR, start=.$TARGET_START, end=.$TARGET_END), PAMpattern = DNAString("NGG")))) %>%
  dplyr::ungroup() %>% as.data.table()
snp_guides <- guide_tables(DT.4, DT.3, "snpid")

# ---- Target TSSs ----
PC.0 <- fread("data/guides/230825_hFOB-SSv1-TargetTSS.txt")
PC.0_tmp <- fread("data/guides/230825_hFOB-SSv1-TargetTSS_preFilter.txt")

# Short-read isoforms that make up 80% of gene expression at any time point
SR.0 <- fread("data/hFOB_RNA/hFOB_RNA.isoforms.txt")
keep_transcripts <- unlist(lapply(c(0, 2, 4, 10), function(tp) unique(fread(paste0("data/hFOB_RNA/hFOB_RNA.isoforms.top_t", tp, ".hg38.txt"))$ensTranscript)))
SR.1 <- subset(SR.0, ensTranscript %in% keep_transcripts)
SR.2 <- subset(SR.1, hgnc_symbol %in% unique(c(PC.0$lr_gene, PC.0_tmp$vep_genes)))

# Long-read isoforms for the WNT16-CPED1 locus genes and TSPAN12
LR.0 <- fread("data/hFOB_RNA/longreads_process/hFOB_longReads_Top80p.txt")
LR.1 <- dplyr::bind_rows(lapply(c("WNT16", "FAM3C", "ING3", "CPED1", "TSPAN12"), function(g) {
  subset(LR.0, grepl(unique(subset(SR.1, hgnc_symbol == g)$ensGene), associated_gene)) %>%
    dplyr::mutate(lr_gene = g)
}))

PC.1 <- PC.0 %>%
  dplyr::mutate(lr_avg_exp = (t0 + t2 + t4 + t10) / 4) %>%
  dplyr::select(lr_gene, isoform, lr_chr, lr_tss, lr_strand, lr_avg_exp) %>%
  dplyr::bind_rows(LR.1 %>%
                     dplyr::mutate(lr_avg_exp = (t0 + t2 + t4 + t10) / 4) %>%
                     dplyr::select(lr_gene, isoform, chr, tss, strand, lr_avg_exp) %>%
                     dplyr::rename(lr_chr = chr, lr_tss = tss, lr_strand = strand)) %>%
  dplyr::distinct() %>%
  dplyr::rename(gene = lr_gene, chr = lr_chr, tss_hg38 = lr_tss, strand = lr_strand, avg_exp = lr_avg_exp) %>%
  dplyr::bind_rows(SR.2 %>%
                     dplyr::mutate(avg_exp = (t0 + t2 + t4 + t10) / 4, strand = ifelse(strand == 1, "+", "-")) %>%
                     dplyr::select(hgnc_symbol, ensTranscript, chromosome_name, transcription_start_site, strand, avg_exp) %>%
                     dplyr::rename(gene = hgnc_symbol, isoform = ensTranscript, chr = chromosome_name, tss_hg38 = transcription_start_site)) %>%
  dplyr::mutate(chr = gsub("chr", "", chr),
                chr = paste0("chr", chr)) %>%
  dplyr::distinct() %>%
  dplyr::group_by(chr) %>%
  dplyr::arrange(tss_hg38) %>%
  dplyr::mutate(lag_pos = tss_hg38 - lag(tss_hg38)) %>%
  dplyr::mutate(rel_pos = tss_hg38 - lag(tss_hg38, default = tss_hg38[1])) %>%
  dplyr::ungroup() %>% as.data.table()
PC.1$blastrad <- paste0(PC.1$gene, "__", number_blast_radius(PC.1$rel_pos))

# Among nearby TSSs keep the most highly expressed isoform
PC.2 <- PC.1 %>%
  collapse_blast_radius(group_col = "gene", id_col = "isoform", pos_col = "tss_hg38", rank_expr = quo(desc(avg_exp))) %>%
  dplyr::rename(redundant_tss = redundant)

# gRNAs within 100 bp downstream of each TSS
PC.3_tmp <- PC.2 %>%
  dplyr::group_by(gene) %>%
  dplyr::arrange(tss_hg38) %>%
  dplyr::mutate(gene_id = paste0(gene, "-TSS", row_number())) %>%
  dplyr::ungroup()
PC.3 <- PC.3_tmp %>%
  dplyr::select(gene_id, chr, tss_hg38, strand) %>%
  dplyr::arrange(chr, tss_hg38) %>%
  dplyr::mutate(TARGET_CHR = chr,
                TARGET_START = ifelse(strand == "+", tss_hg38, tss_hg38 - 100),
                TARGET_END = ifelse(strand == "+", tss_hg38 + 100, tss_hg38)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(NGG = count_ngg(TARGET_CHR, TARGET_START, TARGET_END)) %>%
  dplyr::filter(NGG > 0) %>%
  dplyr::ungroup() %>% as.data.table()

PC.4 <- PC.3 %>%
  dplyr::group_by(gene_id) %>%
  do((GetAllGuidesInRegion(regionToSearch = getSeq(Hsapiens, names=.$TARGET_CHR, start=.$TARGET_START, end=.$TARGET_END), PAMpattern = DNAString("NGG")))) %>%
  dplyr::ungroup() %>% as.data.table()
tss_guides <- guide_tables(PC.4, PC.3, "gene_id")

# ---- Export ----
gRNA_positions <- dplyr::bind_rows(snp_guides$positions, tss_guides$positions %>% dplyr::rename(snpid = gene_id))
gRNA_to_score <- dplyr::bind_rows(snp_guides$to_score, tss_guides$to_score %>% dplyr::rename(snpid = gene_id))

fwrite(gRNA_positions, file = "data/guides/230830_gRNA-positions.txt", col.names = T, row.names = F, quote = T, sep = "\t")
fwrite(gRNA_to_score, file = "data/guides/230830_gRNA-toScore.txt", col.names = T, row.names = F, quote = T, sep = "\t")
fwrite(DT.2, file = "data/guides/230830_SNPs-preLibrary.txt", col.names = T, row.names = F, quote = F, sep = "\t")
fwrite(PC.3_tmp, file = "data/guides/230830_TSSs-preLibrary.txt", col.names = T, row.names = F, quote = F, sep = "\t")
