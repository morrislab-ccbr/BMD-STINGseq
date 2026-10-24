#!/bin/bash
# Generate UMI count matrices for the BMD STING-seq experiments with 10x Cell Ranger.
#   pilot: 2 10x chip lanes, three NextSeq runs (one Mid-, two High-Output)
#   full:  6 10x chip lanes, one NextSeq High-Output run and one NovaSeq X run
# Each library has gene expression (cDNA), hashtag oligo (HTO) and gRNA (GDO) reads.
#
# Usage: bash 01_cellranger.sh            prepare FASTQs, references and library files
#        bash 01_cellranger.sh count <experiment> <lane>   run Cell Ranger count for one lane
#
# Cell Ranger mkfastq/mkref v7.1.0; count v8.0.1.

PROJECT_DIR="${PROJECT_DIR:-path/to/230126_BMD-STINGseq}"

# ---- Cell Ranger count for one chip lane (submit as a job; 12 cores, 124 GB) ----
if [[ "$1" == "count" ]]; then
  EXP=$2; S=$3
  if [[ "${EXP}" == "pilot" ]]; then ID="BMDSS-pilot-midhigh2x-${S}"; else ID="BMDSS-full-highnovax-${S}_submitjob"; fi
  cd ${PROJECT_DIR}/cellranger
  cellranger count --feature-ref=../data/${EXP}_seq/all_features.csv \
   --transcriptome=refdata-gex-GRCh38-2024-A_CRISPRi \
   --id=${ID} \
   --libraries=../data/${EXP}_seq/cDNA-${S}_library_${EXP}.csv \
   --create-bam=true --localcores=12 --localmem=120
  exit 0
fi

# Feature reference (HTO + gRNA) and per-lane library CSVs for one experiment
make_library_files() {
  local EXP=$1; shift
  cd ${PROJECT_DIR}/data/${EXP}_seq

  awk 'BEGIN{FS="\t"; OFS=","; print "id","name","read","pattern","sequence","feature_type"}{
   print $1,$1,"R2","5PNNNNNNNNNN(BC)",$2,"Antibody Capture"}' HTO_barcodes.txt > HTO_features.csv

  # gRNA capture sequences are reverse complemented; target gene fields are placeholders (CD55)
  # because Cell Ranger's gRNA-target analysis is not used
  paste -d ',' <(cut -f 1 GDO_sequences.txt) <(cut -f 2 GDO_sequences.txt | tr ACGTacgt TGCAtgca | rev) > GDO_barcodes.txt
  awk 'BEGIN{FS=","; OFS=","; print "id","name","read","pattern","sequence","feature_type","target_gene_id","target_gene_name"}{
   print $1,$1,"R2","GCATAGCTCTTAAAC(BC)",$2,"CRISPR Guide Capture","ENSG00000196352","CD55"}' GDO_barcodes.txt > GDO_features.csv

  cat HTO_features.csv <(cut -d ',' -f 1-6 GDO_features.csv | tail -n +2) > all_features.csv

  for S in "$@"; do
   {
    echo "fastqs,sample,library_type"
    echo "${PROJECT_DIR}/data/${EXP}_seq/fastq_cat,cDNA-${S},Gene Expression"
    echo "${PROJECT_DIR}/data/${EXP}_seq/fastq_cat,HTO-${S},Antibody Capture"
    echo "${PROJECT_DIR}/data/${EXP}_seq/fastq_cat,GDO-${S},CRISPR Guide Capture"
   } > cDNA-${S}_library_${EXP}.csv
  done
}

# Concatenate mkfastq output per library; each sequencing run is labelled as a separate "lane"
concat_fastq() {
  local FLOWCELL=$1 LANE=$2; shift 2
  for G in I1 I2 R1 R2; do
   for M in cDNA GDO HTO; do
    for S in "$@"; do
     cat fastq/outs/fastq_path/${FLOWCELL}/${M}-${S}_S*_L00*_${G}_001.fastq.gz > fastq_cat/${M}-${S}_S1_${LANE}_${G}_001.fastq.gz
  done;done;done
}

# ---- Pilot experiment ----
cd ${PROJECT_DIR}/data/pilot_seq
for run in 240325_NS500146_0040_AHMMCJAFX5 240409_NS500146_0041_AHWJKKBGXT 240411_NS500146_0042_AH5YK2BGXW; do
 cellranger mkfastq --run=${PROJECT_DIR}/data/pilot_seq/${run} --id=fastq --csv=cellranger-pilot-info.csv --jobmode=local --localcores=10 --localmem=96
done
mkdir -p fastq_cat
concat_fastq HMMCJAFX5 L001 1 2
concat_fastq HWJKKBGXT L002 1 2
concat_fastq H5YK2BGXW L003 1 2
make_library_files pilot 1 2

# ---- Full experiment ----
cd ${PROJECT_DIR}/data/full_seq
cellranger mkfastq --run=${PROJECT_DIR}/data/full_seq/240624_NS500146_0047_AH5YL2BGXW --id=fastq --csv=cellranger-full-info.csv --jobmode=local --localcores=8 --localmem=64
mkdir -p fastq_cat
concat_fastq H5YL2BGXW L001 1 2 3 4 5 6

# NovaSeq X FASTQs (provider samples 1-6 = cDNA, 7-12 = GDO, 13-18 = HTO for lanes 1-6)
raw=${PROJECT_DIR}/data/full_seq/AN00020436_10X_RawData_Outs
for G in I1 I2 R1 R2; do
 for M in $(seq 1 18); do
  S=$(( (M - 1) % 6 + 1 ))
  if (( M <= 6 )); then LIB=cDNA; elif (( M <= 12 )); then LIB=GDO; else LIB=HTO; fi
  ln -s ${raw}/Farber_ANQ00005202_${M}/22F3KYLT4/Farber_ANQ00005202_${M}_S${M}_L005_${G}_001.fastq.gz fastq_cat/${LIB}-${S}_S1_L002_${G}_001.fastq.gz
 done
done
make_library_files full 1 2 3 4 5 6

# ---- Reference: GRCh38 (10x 2024-A) plus the KRAB-dCas9-MeCP2 (CRISPRi) transgene ----
mkdir -p ${PROJECT_DIR}/cellranger/refs ${PROJECT_DIR}/logs/cellranger
cd ${PROJECT_DIR}/cellranger/refs
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz

cd ${PROJECT_DIR}/cellranger
cp refs/refdata-gex-GRCh38-2024-A/fasta/genome.fa genome_CRISPRi.fa
cat CRISPRi.fa >> genome_CRISPRi.fa
zcat refs/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz > genes_CRISPRi.gtf
cat CRISPRi.gtf >> genes_CRISPRi.gtf
cellranger mkref --genome=refdata-gex-GRCh38-2024-A_CRISPRi \
 --fasta=genome_CRISPRi.fa \
 --genes=genes_CRISPRi.gtf

# ---- Count each lane ----
echo "Submit: bash 01_cellranger.sh count pilot <1-2>  and  bash 01_cellranger.sh count full <1-6>"
