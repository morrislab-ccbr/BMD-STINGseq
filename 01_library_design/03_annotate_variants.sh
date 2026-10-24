#!/bin/bash
# Process hFOB ATAC-seq peaks and annotate fine-mapped eBMD variants with the nearest
# ATAC-seq peaks/summits and long-read RNA-seq TSSs at each time point (days 0, 2, 4, 10).
#
# Inputs: data/hFOB_ATAC/hFOB_day{0,2,4,10}_peaks.hg38.txt (MACS2 narrowPeak calls, hg38) and
#         ${EXTERNAL_DIR}/hFOB_ATAC_peaks/peakcall/*.hg38.bedgraph;
#         outputs of 01_process_rna.R and 02_finemap_gwas.sh.
# Requires bedtools and UCSC bedGraphToBigWig.

PROJECT_DIR="${PROJECT_DIR:-path/to/230126_BMD-STINGseq}"
EXTERNAL_DIR="${EXTERNAL_DIR:-path/to/external_inputs}"

module purge
module load bedtools

# ---- ATAC-seq peaks, consensus peaks and summits ----
cd ${PROJECT_DIR}/data/hFOB_ATAC
for i in 0 2 4 10; do
 awk '$1 !~ /_/' hFOB_day${i}_peaks.hg38.txt | awk 'OFS="\t"{print $1,$2,$3,$4}' | sed -e "s/peak_/day${i}_peak_/g" | sort -k 1,1 -k 2,2n > hFOB_day${i}_peaks.hg38.bed
 awk '$1 !~ /_/' hFOB_day${i}_peaks.hg38.txt | awk 'OFS="\t"{print $1,$2+$10-1,$2+$10,$4}' | sed -e "s/peak_/day${i}_summit_/g" | sort -k 1,1 -k 2,2n > hFOB_day${i}_summits.hg38.bed
done

# Merge peaks across time points into consensus peaks, then link back to the per-day peaks
bedtools merge -i <(cat hFOB_day*_peaks.hg38.bed | sort -k 1,1 -k 2,2n) | awk 'OFS="\t"{print $0, "dayAll_peak_"NR}' > hFOB_dayAll_peaks.hg38.bed
bedtools intersect -a hFOB_dayAll_peaks.hg38.bed -b hFOB_day0_peaks.hg38.bed hFOB_day2_peaks.hg38.bed hFOB_day4_peaks.hg38.bed hFOB_day10_peaks.hg38.bed -wao -names day0 day2 day4 day10 > hFOB_dayAll_peaksAll.hg38.txt

# ---- ATAC-seq bigWig tracks ----
mkdir -p bedgraphs bigwigs
cp ${EXTERNAL_DIR}/hFOB_ATAC_peaks/peakcall/*.hg38.bedgraph bedgraphs/.
n=1
for i in 0 2 4 10; do
 tail -n +2 bedgraphs/${n}_day${i}.hg38.bedgraph | awk 'OFS="\t"{print $1, $2, $3, $8}' | sort -k 1,1 -k 2,2n > bedgraphs/hFOB_day${i}_ATAC.hg38.sorted.bedgraph
 bedGraphToBigWig bedgraphs/hFOB_day${i}_ATAC.hg38.sorted.bedgraph http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes bigwigs/hFOB_day${i}_ATAC.hg38.bw
 n=$((n + 1))
done

# ---- Annotate fine-mapped variants ----
cd ${PROJECT_DIR}/data
tail -n +2 UKB-eBMD-SuSiE-SumStats.230707.txt | awk 'OFS="\t"{print "chr"$3, $5-1, $5, $1}' | sort -k 1,1 -k 2,2n | uniq > UKB-eBMD.hg38.bed

mkdir -p gwas_anno
for i in 0 2 4 10; do
 bedtools closest -a UKB-eBMD.hg38.bed -b hFOB_ATAC/hFOB_day${i}_summits.hg38.bed -d > gwas_anno/UKB-eBMD.ATACs-D${i}.closest.hg38.txt
 bedtools intersect -a UKB-eBMD.hg38.bed -b hFOB_ATAC/hFOB_day${i}_peaks.hg38.bed -wb > gwas_anno/UKB-eBMD.ATACp-D${i}.closest.hg38.txt
 bedtools closest -a UKB-eBMD.hg38.bed -b <(sort -k 1,1 -k 2,2n hFOB_RNA/longreads_process/hFOB_longReads_Top80p_t${i}.bed) -d -k 3 > gwas_anno/UKB-eBMD.lrRNA-D${i}.closest3.hg38.txt
done
