# Select the final STING-seq v3 gRNA library: three gRNAs per target ranked by FlashFry
# Hsu 2013 specificity, plus CD55-targeting (MOI) and non-targeting control gRNAs.

project_dir <- Sys.getenv("PROJECT_DIR", unset = "path/to/230126_BMD-STINGseq")
external_dir <- Sys.getenv("EXTERNAL_DIR", unset = "path/to/external_inputs")
setwd(project_dir)

library(data.table)
library(tidyverse)

# Top three gRNAs per target by Hsu 2013 score; drop targets with fewer than three
top3 <- function(x) {
  x %>%
    dplyr::group_by(snpid) %>%
    dplyr::arrange(desc(Hsu2013)) %>%
    dplyr::filter(row_number() <= 3) %>%
    dplyr::mutate(n_guides = length(snpid)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n_guides == 3) %>%
    dplyr::select(-n_guides) %>%
    as.data.table()
}

DT.0 = fread("data/guides/230830_gRNA-positions.txt")
DT.1 = fread("data/guides/230830_gRNA-toScore.txt")
DT.2 = fread("data/guides/230830_gRNA-Scored.txt")
DT.3 = fread("data/guides/230830_SNPs-preLibrary.txt")

# Remove gRNAs with poly-T or multiple perfect genomic matches; keep Hsu 2013 >= 50 and <= 500 off-targets
DT.4 <- DT.2 %>%
  dplyr::filter(dangerous_polyT == "NONE", dangerous_in_genome == "IN_GENOME=1") %>%
  dplyr::select(ID, guide, Doench2014OnTarget, Hsu2013, DoenchCFD_specificityscore, otCount) %>%
  dplyr::mutate(GUIDE_ID = paste0(ID, "_", guide))

DT.5 <- DT.4 %>%
  dplyr::group_by(ID) %>%
  dplyr::filter(Hsu2013 >= 50, otCount <= 500) %>%
  dplyr::arrange(desc(Hsu2013)) %>%
  dplyr::mutate(n_guides = length(ID)) %>%
  dplyr::ungroup()

# Add PAM sequences and genomic positions to the scored gRNAs
DT.1_ = DT.1 %>%
  dplyr::select(FlashFryID, GUIDE_SEQ, PAM) %>%
  dplyr::rename(ID = FlashFryID, gRNA_seq_wPAM = GUIDE_SEQ) %>%
  dplyr::mutate(gRNA_seq_woPAM = substr(gRNA_seq_wPAM, 1, 20)) %>%
  dplyr::mutate(ID__gRNA_seq_woPAM = paste0(ID, "__", gRNA_seq_woPAM))
DT.0_ = DT.0 %>%
  dplyr::rename(ID = FlashFryID, gRNA_seq_woPAM = gRNA_seq_wPAM) %>%
  dplyr::mutate(ID__gRNA_seq_woPAM = paste0(ID, "__", gRNA_seq_woPAM)) %>%
  dplyr::select(-ID, -gRNA_seq_woPAM) %>%
  dplyr::left_join(DT.1_, by = "ID__gRNA_seq_woPAM") %>%
  dplyr::select(-ID__gRNA_seq_woPAM, -ID) %>%
  dplyr::mutate(snpid__gRNA_seq_wPAM = paste0(snpid, "__", gRNA_seq_wPAM)) %>%
  dplyr::select(-snpid, -gRNA_seq_wPAM)

# Candidate CRE (variant) gRNAs
DT.6 = DT.5 %>%
  dplyr::filter(!grepl("TSS", ID)) %>%
  dplyr::mutate(snpid = gsub("-", ":", ID)) %>%
  dplyr::rename(gRNA_seq_wPAM = guide) %>%
  dplyr::select(-ID, -GUIDE_ID) %>%
  dplyr::mutate(snpid__gRNA_seq_wPAM = paste0(snpid, "__", gRNA_seq_wPAM)) %>%
  dplyr::inner_join(DT.0_, by = "snpid__gRNA_seq_wPAM") %>%
  dplyr::select(-snpid__gRNA_seq_wPAM, -Doench2014OnTarget, -DoenchCFD_specificityscore, -n_guides, -gRNA_seq_wPAM) %>%
  as.data.table()

DT.7 = top3(DT.6)

# TSS gRNAs
PC.2 = DT.5 %>%
  dplyr::filter(grepl("TSS", ID)) %>%
  dplyr::rename(snpid = ID, gRNA_seq_wPAM = guide) %>%
  dplyr::select(-GUIDE_ID) %>%
  dplyr::mutate(snpid__gRNA_seq_wPAM = paste0(snpid, "__", gRNA_seq_wPAM)) %>%
  dplyr::inner_join(DT.0_, by = "snpid__gRNA_seq_wPAM") %>%
  dplyr::select(-snpid__gRNA_seq_wPAM, -Doench2014OnTarget, -DoenchCFD_specificityscore, -n_guides, -gRNA_seq_wPAM) %>%
  as.data.table()
PC.3 = top3(PC.2)

# Name gRNAs: variants are numbered by genomic position (SNP-1 ... SNP-n); TSSs by gene and TSS
snp_numbers = DT.7 %>%
  dplyr::select(snpid) %>%
  dplyr::left_join(DT.3 %>%
                     dplyr::select(snpid, chr, pos_hg38), by = "snpid") %>%
  dplyr::distinct() %>%
  dplyr::arrange(chr, pos_hg38) %>%
  dplyr::mutate(snp_index = paste0("SNP-", row_number())) %>%
  dplyr::select(-chr, -pos_hg38)

snp_numbers$snp_index = factor(snp_numbers$snp_index, levels = unique(snp_numbers$snp_index))

LB.0 = DT.7 %>%
  dplyr::left_join(snp_numbers, by = "snpid") %>%
  dplyr::arrange(snp_index) %>%
  dplyr::group_by(snp_index) %>%
  dplyr::arrange(gRNA_start) %>%
  dplyr::mutate(gRNA_index = paste0(snp_index, "-", row_number())) %>%
  dplyr::ungroup() %>% as.data.table() %>%
  dplyr::relocate(gRNA_index, snp_index, snpid, gRNA_chr, gRNA_start, gRNA_end, gRNA_strand, gRNA_seq_woPAM, PAM, Hsu2013) %>%
  dplyr::arrange(snp_index)

LB.1 = PC.3 %>%
  dplyr::arrange(gRNA_start) %>%
  dplyr::group_by(snpid) %>%
  dplyr::mutate(gRNA_index = paste0(snpid, "-", row_number())) %>%
  dplyr::arrange(snpid) %>%
  dplyr::ungroup() %>% as.data.table() %>%
  dplyr::relocate(gRNA_index, snpid, gRNA_chr, gRNA_start, gRNA_end, gRNA_strand, gRNA_seq_woPAM, PAM, Hsu2013)

LB.2 = dplyr::bind_rows(
  LB.0, LB.1
)

# Controls: 22 CD55-targeting gRNAs from the STING-seq v2 library (MOI controls) and
# 30 GeCKO v2 non-targeting gRNAs
SS.0 = fread(paste0(external_dir, "/STINGseq_v2/210413.STINGseq_v2_library.txt"))
SS.1 = subset(SS.0, grepl("CD55", ConstID))
SS.2 = subset(SS.1, ConstID %in% c("CD55-MOI-27", "CD55-MOI-41", "CD55-MOI-42", "CD55-MOI-20", "CD55-MOI-17", "CD55-MOI-38", "CD55-MOI-37", "CD55-MOI-12", "CD55-MOI-4", "CD55-MOI-2", "CD55-MOI-3",
                                   "CD55-MOI-9", "CD55-MOI-16", "CD55-MOI-25", "CD55-MOI-10", "CD55-MOI-22", "CD55-MOI-43", "CD55-MOI-15", "CD55-MOI-19", "CD55-MOI-36", "CD55-MOI-34", "CD55-MOI-39"))

GV.0 = fread(paste0(external_dir, "/sanjanalab/GeCKO_human_ntcguides.txt"), header = F) %>%
  dplyr::filter(!grepl("TTTT", V1), !grepl("CCCCC", V1), !grepl("AAAAA", V1), !grepl("GGGGG", V1)) %>%
  dplyr::rename(gRNA_seq_woPAM = V1)
set.seed(1234); GV.1 = GV.0 %>%
  dplyr::sample_n(30) %>%
  dplyr::mutate(gRNA_index = paste0("NTC-GV2-", row_number()),
                snpid = "NTC-GV2")

LB.3 = dplyr::bind_rows(LB.2,
                        SS.2 %>%
                          dplyr::rename(gRNA_index = ConstID, gRNA_seq_woPAM = guide_noPAM) %>%
                          dplyr::select(-ID) %>%
                          dplyr::mutate(snpid = "CD55-MOI"),
                        GV.1) 

# Oligo with flanking sequences for cloning into the gRNA backbone
cas9_EF_beforeGuideForOligo = 'GGAAAGGACGAAACACCG'
cas9_EF_afterGuideForOligo = 'GTTTAAGAGCTATGCTGGAAACAGC'

LB.3$fullOligo = paste0(cas9_EF_beforeGuideForOligo, LB.3$gRNA_seq_woPAM, cas9_EF_afterGuideForOligo)

fwrite(LB.3, file = "data/guides/230905_SSv3-BMD-hFOB-Library.txt", col.names = T, row.names = F, quote = F, sep = "\t")
