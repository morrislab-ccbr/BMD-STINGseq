#!/bin/bash
#SBATCH --job-name=FlashFry
#SBATCH --cpus-per-task=2
#SBATCH --mem=4g
#SBATCH --output=logs/flashfry/%A_%a.out
#SBATCH --error=logs/flashfry/%A_%a.out
#
# Score candidate gRNAs (from 05_design_grnas.R) for on-target efficiency and off-target
# specificity with FlashFry v1.15 against hg38, one SLURM array task per target.
#
# Usage: bash 06_score_grnas.sh

PROJECT_DIR="${PROJECT_DIR:-path/to/230126_BMD-STINGseq}"
export SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
FF="${PROJECT_DIR}/download"
GUIDES="230830_gRNA-toScore.txt"

module purge
module load java/1.8
cd ${PROJECT_DIR}/data/guides

# ---- Array task: score all gRNAs for one target ----
if [[ "$1" == "score" ]]; then
  n=$(sed -n "${SLURM_ARRAY_TASK_ID}p" FlashFry/temp/var_list_tars | cut -f 2)

  sed -e 's/"//g' ${GUIDES} | awk -v SNP=${n} '{if ($9 == SNP) print ">"$9"_"$5"\n"$10}' > FlashFry/fasta/${n}.fa

  java -Xmx4g -jar ${FF}/FlashFry-assembly-1.15.jar \
   discover \
   --database ${FF}/hg38_cas9ngg_database \
   --fasta FlashFry/fasta/${n}.fa \
   --output FlashFry/output/${n}.output

  java -Xmx4g -jar ${FF}/FlashFry-assembly-1.15.jar \
   score \
   --input FlashFry/output/${n}.output \
   --output FlashFry/scored/${n}.output.scored \
   --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013 \
   --database ${FF}/hg38_cas9ngg_database
  exit 0
fi

# ---- Build the FlashFry hg38 SpCas9 (NGG) database ----
mkdir -p ${FF}/tmp
cd ${FF}
wget https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
java -Xmx4g -jar FlashFry-assembly-1.15.jar \
 index \
 --tmpLocation ./tmp \
 --database hg38_cas9ngg_database \
 --reference hg38.fa.gz \
 --enzyme spcas9ngg

# ---- Score gRNAs, one array task per target (12,276 gRNAs for 389 targets) ----
cd ${PROJECT_DIR}/data/guides
mkdir -p FlashFry/temp FlashFry/output FlashFry/scored FlashFry/fasta ${PROJECT_DIR}/logs/flashfry

tail -n +2 ${GUIDES} | awk '{print $8"\t"$9}' | sed -e 's/"//g' | sed -e 's/chr//g' | sort | uniq > FlashFry/temp/var_list_tars
n_targets=$(wc -l < FlashFry/temp/var_list_tars)

cd ${PROJECT_DIR}
sbatch --wait --array=1-${n_targets} ${SCRIPT_DIR}/06_score_grnas.sh score

cd ${PROJECT_DIR}/data/guides
cat <(head -n 1 -q FlashFry/scored/* | uniq | sed -e 's/contig/ID\tguide/') <(tail -n +2 -q FlashFry/scored/*.output.scored | sed -e 's/_/\t/1' | awk '{if ($2 == $5) print $0}') > 230830_gRNA-Scored.txt
