# BMD STING-seq

Analysis code for:

> Kahlani M, Mesner L, *et al.* Identification of target genes and regulatory networks for bone
> mineral density GWAS loci through systematic targeting and inhibition. *Nature Genetics* (2026).
> DOI: TBD

STING-seq links fine-mapped noncoding GWAS variants to their target genes and downstream
regulatory networks with pooled CRISPRi and single-cell RNA-seq. Here it was applied to bone
mineral density (BMD) loci in human fetal osteoblast (hFOB) cells across osteogenic
differentiation.

## Pipeline

Run the scripts in order within each folder. Shell scripts that launch SLURM array jobs
(`02_finemap_gwas.sh`, `06_score_grnas.sh`) submit themselves for the per-locus or per-target steps.

| Step | Script | What it does |
|---|---|---|
| **1. Library design** | `01_library_design/01_process_rna.R` | hFOB short-read and long-read RNA-seq: expressed isoforms and TSSs per time point |
| | `01_library_design/02_finemap_gwas.sh` (+ `.R`) | SuSiE fine-mapping of eBMD GWAS loci with UK Biobank LD; VEP input |
| | `01_library_design/03_annotate_variants.sh` | hFOB ATAC-seq peaks and bigWigs; nearest peaks and TSSs for fine-mapped variants |
| | `01_library_design/04_select_targets.R` | Selects candidate CREs and TSS targets |
| | `01_library_design/05_design_grnas.R` | Designs candidate gRNAs |
| | `01_library_design/06_score_grnas.sh` | Scores gRNAs with FlashFry |
| | `01_library_design/07_finalize_library.R` | Final gRNA library with controls (Tables S2B–D) |
| | `01_library_design/08_plasmid_library_qc.sh` | Plasmid library sequencing QC |
| **2. Single-cell processing** | `02_singlecell/01_cellranger.sh` | FASTQs, CRISPRi-augmented reference, Cell Ranger count |
| | `02_singlecell/02_qc_lanes.R` | Per-lane cell QC, hashtag demultiplexing and gRNA UMI thresholds |
| | `02_singlecell/03_combine_lanes.R` | Merges lanes into the analysed cell sets; expression PCs |
| **3. SCEPTRE** | `03_sceptre/01_make_sceptre_objects.R` | Builds cis (`cis`) and trans (`trans`) SCEPTRE objects |
| | `03_sceptre/02_run_sceptre.sh` | Runs sceptre-pipeline (Nextflow) for one analysis |
| **4. Results** | `04_results/01_annotate_sceptre_results.R` | Annotated cis and trans results, positive-control summary |
| | `04_results/02_individual_grnas.R` | Individual-gRNA analyses (Figs S6, S8, S12) |
| **5. Figures** | `05_figures/01_main_figures.R` | Figures 1 and 2 |
| | `05_figures/02_violin_plots.R` | Single-cell expression violin plots for selected CREs |
| | `05_figures/03_revision_figures.R` | ATAC-seq PCA, marker gene expression (Fig. S1) and hMSC replication (Fig. 3B–C) |

SCEPTRE analyses (`02_run_sceptre.sh <analysis>`):

| Analysis | Cells | Run directory |
|---|---|---|
| `pilot_full_day0` | Day 0 (pilot + full) | `sceptre/240831_BMDSS_pilot_full_day0_cis_pca` |
| `full_day2`, `full_day4` | Day 2, day 4 | `sceptre/240831_BMDSS_full_day{2,4}_cis_pca` |
| `full_day2_day4` | Days 2 + 4 | `sceptre/250210_BMDSS_full_day2_day4_cis_pca` |
| `pilot_full` | All cells | `sceptre/240831_BMDSS_pilot_full_cis_pca25` |
| `*_ds` | 10,000 cells per experiment and day | `sceptre/250114_BMDSS_*_ds` |
| `pilot_full_day0_trans`, `pilot_full_trans`, `full_day2_day4_trans` | Trans analyses | see script |
| `indv_cis`, `indv_trans` | Individual gRNAs | `sceptre/2605*/2606*_..._IndvGuides` |

The trans results used in the paper were run by a co-author with these settings and are read
from `sceptre/250526_BMDSS_pilot_full_trans_pca_massive`, `sceptre/250526_BMDSS_full_day2_day4_trans_pca_massive`,
`sceptre/250622_240831_BMDSS_pilot_full_day0_trans_pca` and
`sceptre/260618_BMDSS_pilot_full_trans_pca_massive_NTC_3_10_singleton`.

Script names were shortened during clean-up; the dated originals (e.g. `250626_AnnotateSceptreResults.R`)
are available on request.

## Running the code

Paths are relative to one project directory with subfolders `data/`, `cellranger/`, `seurat/`,
`sceptre/`, `summaries/`, `plots/`, `figures/` and `revisions/`. Set two environment variables:

| Variable | Meaning |
|---|---|
| `PROJECT_DIR` | Project directory; every script works from here |
| `EXTERNAL_DIR` | Inputs produced outside this project: hFOB ATAC-seq peak calls and nf-core output, UK Biobank imputed genotypes, UCSC liftOver chains, and reference lists from earlier STING-seq work |

Other inputs expected under `data/`: the Morris *et al.* 2019 eBMD summary statistics
(`data/Morris2019/`), hFOB RNA-seq quantifications (`data/hFOB_RNA/`), Ensembl VEP output
(`data/vep/`), the Lambert *et al.* 2018 transcription factor list (`data/human_TFs/`) and Ensembl
110 gene coordinates (`data/genes_hg38_ensembl110.txt`).

UK Biobank genotypes (used only to build LD matrices in `02_finemap_gwas.sh`) require an approved
UK Biobank application and are not distributed here.

## Software

R 4.0–4.5 with data.table, tidyverse, Seurat v5.3.1, sceptre, ondisc and susieR v0.14.2;
10x Cell Ranger v8.0.1; sceptre-pipeline run with Nextflow; FlashFry v1.15; PLINK; bedtools;
UCSC liftOver; Ensembl VEP.

## Code available on request

The following analyses in the paper are not included here and are available from the
corresponding authors on request:

- ATAC-seq read processing and peak calling (Genrich, nf-core/atacseq)
- Consensus non-negative matrix factorization (cNMF) and gene-program association tests
- LD score regression (LDSC)
- Allele-specific chromatin accessibility (WASP, GATK ASEReadCounter)
- GWAS–eQTL colocalization (fastENLOC)
- hMSC bulk RNA-seq processing and differential expression (nf-core/rnaseq, DESeq2)
- Arrayed CRISPRi/CRISPRa, siRNA, qPCR and mineralization analyses

## Data availability

- hFOB ATAC-seq: GEO GSE318537
- BMD STING-seq single-cell data: GEO GSE319266
- hFOB whole-genome sequencing: dbGaP phs004654.v1.p1
- hFOB short-read and long-read RNA-seq: GEO GSE186922 and GSE224588
- hMSC RNA-seq: TBD

## License

MIT; see `LICENSE`.
