#!/bin/bash
#SBATCH --job-name=finemap
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --output=logs/finemap/%A_%a.out
#SBATCH --error=logs/finemap/%A_%a.out
#
# Fine-map eBMD GWAS loci (Morris et al. 2019) with SuSiE, using LD from 50,000 randomly
# selected UK Biobank participants (data/Morris2019/plink/Random_50k.seed42.txt).
#
# Usage: bash 02_finemap_gwas.sh
# The same script is submitted as the per-locus SLURM array tasks (plink | susie).

PROJECT_DIR="${PROJECT_DIR:-path/to/230126_BMD-STINGseq}"
EXTERNAL_DIR="${EXTERNAL_DIR:-path/to/external_inputs}"
export SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd ${PROJECT_DIR}

# ---- Array task: LD matrix for one locus ----
if [[ "$1" == "plink" ]]; then
  read LOCUS RSID CHR LOWER UPPER RANGE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" data/Morris2019/UKB-eBMD-SigRegions.txt)
  cd data/Morris2019/plink
  mkdir -p ${CHR}/${RANGE}

  module purge
  module load plink/2.0
  plink --bgen ${EXTERNAL_DIR}/UKBiobank_imputed_bgen/UKWB-ImpQc-Chr${CHR}.bgen \
   --sample ${EXTERNAL_DIR}/UKBiobank_imputed_bgen/UKWB-ImpQc-Chr${CHR}.sample \
   --keep Random_50k.seed42.txt \
   --memory 60000 \
   --threads 8 \
   --chr ${CHR} \
   --from-bp ${LOWER} \
   --to-bp ${UPPER} \
   --make-bed \
   --out ./${CHR}/${RANGE}/${RANGE}

  module purge
  module load plink/1.90b6.24
  plink --bfile ./${CHR}/${RANGE}/${RANGE} \
   --memory 60000 \
   --threads 8 \
   --keep-allele-order \
   --r square \
   --out ./${CHR}/${RANGE}/${RANGE}
  exit 0
fi

# ---- Array task: SuSiE for one locus ----
if [[ "$1" == "susie" ]]; then
  module purge
  module load R/4.0.2
  Rscript ${SCRIPT_DIR}/02_finemap_gwas.R susie ${SLURM_ARRAY_TASK_ID}
  exit 0
fi

# ---- Pipeline ----
mkdir -p logs/finemap data/Morris2019/susie data/Morris2019/susie_plots data/Morris2019/plink data/Morris2019/liftover

# Lift GWAS variants from hg19 to hg38
module purge
module load liftover
cd data/Morris2019
for i in {1..22}; do
 zcat Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz | tail -n +2 | awk -v CHR=${i} 'OFS="\t"{ if ($3 == CHR) print "chr"$3, $4-1, $4, $2 }' | sort -k 1,1 -k 2,2n > liftover/UKB-eBMD.chr${i}.hg19.bed
 liftOver liftover/UKB-eBMD.chr${i}.hg19.bed ${EXTERNAL_DIR}/ucsc/hg19ToHg38.over.chain.gz liftover/UKB-eBMD.chr${i}.hg38.bed liftover/unlifted.chr${i}.bed
done
cd ${PROJECT_DIR}

# Define loci, compute LD (515 loci), fine-map (510 autosomal loci with LD)
Rscript ${SCRIPT_DIR}/02_finemap_gwas.R define_loci
sbatch --wait --array=1-515 ${SCRIPT_DIR}/02_finemap_gwas.sh plink
sbatch --wait --array=1-510 ${SCRIPT_DIR}/02_finemap_gwas.sh susie

# Combine loci and keep variants in 95% credible sets
cd data/Morris2019
cat <(head -n 1 susie/eBMD_SuSiE_locus9.txt) <(tail -n +2 -q susie/eBMD_SuSiE_locus*.txt) > UKB-eBMD-SuSiE.230707.txt
cd ${PROJECT_DIR}
Rscript ${SCRIPT_DIR}/02_finemap_gwas.R collect

# Export credible-set variants as VCF for Ensembl VEP (output: data/vep/UKB-eBMD-SuSiE.230707.EA-NEA.vep.txt)
mkdir -p data/vep
tail -n +2 data/UKB-eBMD-SuSiE-SumStats.230707.txt | awk 'OFS=" "{print $3, $5, $1, $6, $7, ".", ".", "."}' | sort -k 1,1n -k 2,2n | uniq > data/vep/UKB-eBMD-SuSiE.230707.EA-NEA.vcf
