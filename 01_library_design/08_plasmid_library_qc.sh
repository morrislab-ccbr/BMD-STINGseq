#!/bin/bash
# Sequence the cloned STING-seq v3 plasmid library and check gRNA representation
# (recovery and skew, the ratio of the 90th to 10th percentile of gRNA counts).

PROJECT_DIR="${PROJECT_DIR:-path/to/230126_BMD-STINGseq}"
cd ${PROJECT_DIR}/data
mkdir -p basespace basespace_fastq basespace_fastq_qc

# Download the MiSeq run and demultiplex
module purge
module load basespace/0.10.7
bs-cp -v https://basespace.illumina.com/Run/271775511 basespace

module purge
module load bcl2fastq/2.20.0.422
module load fastqc/0.11.4
bcl2fastq -R basespace/ \
 --sample-sheet basespace/SampleSheet.csv \
 --output-dir basespace_fastq/
fastqc basespace_fastq/*.fastq.gz -o basespace_fastq_qc/

# Trim the plasmid sequences flanking the spacer
cd basespace_fastq
module purge
module load cutadapt
cutadapt -g CGAAACACCG -o Undetermined_S0_L001_R1_001.cut1.fastq Undetermined_S0_L001_R1_001.fastq.gz --discard-untrimmed -m 20 -e 0.2 > cutadapt1.log.txt
cutadapt -a GTTTAAGAGC -o Undetermined_S0_L001_R1_001.cut2.fastq Undetermined_S0_L001_R1_001.cut1.fastq --discard-untrimmed -m 20 -e 0.2 > cutadapt2.log.txt

# Align spacers to the library (unique hits, <= 1 mismatch, forward strand only)
module purge
module load python/2.7.8
module load bowtie/1.1.2
tail -n +2 ../guides/230905_SSv3-BMD-hFOB-Library.txt | awk 'FS="\t"{print ">"$1"\n"$8}' > STINGseq_v3_library.fa
bowtie-build STINGseq_v3_library.fa STINGseq_v3_library.lib > bowtie_build.log.txt
bowtie -m 1 -v 1 --norc STINGseq_v3_library.lib -q Undetermined_S0_L001_R1_001.cut2.fastq Undetermined_S0_L001_R1_001.cut2.fastq.aligned --un Undetermined_S0_L001_R1_001.cut2.fastq.unaligned 2> bowtie_log.txt

# Library recovery and skew
Rscript - <<'EOF'
library(data.table)
library_size <- nrow(fread("../guides/230905_SSv3-BMD-hFOB-Library.txt"))
counts <- as.data.table(table(fread("Undetermined_S0_L001_R1_001.cut2.fastq.aligned")$V3))
pct <- quantile(counts$N, c(0.1, 0.9))
cat(sprintf("gRNAs detected: %d of %d (%.1f%%)\nSkew (90th / 10th percentile): %.2f\n",
            nrow(counts), library_size, 100 * nrow(counts) / library_size, pct[2] / pct[1]))
EOF
