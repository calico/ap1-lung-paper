#!/bin/bash
#SBATCH --job-name=atac-process
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --error=%x.%A_%a.error
#SBATCH --output=%x.%A_%a.out
#SBATCH --array=0-5

# ATAC-seq processing: raw fastq -> sorted, indexed BAM
#
# Usage:  sbatch atac_process.sh <fastq_dir> <output_dir>
#
# Runs as a 6-element job array, one sample per array task.
# Requires: cutadapt 2.8, bowtie2, samtools

set -euo pipefail

FQ_DIR=$1
OUT_DIR=$2

# mm10 bowtie2 index prefix
INDEX=/home/alyssa/references/mm10/mm10

# Nextera transposase adapter (Tn5 mosaic end)
ADAPTER=CTGTCTCTTATACACATCT

# The 6 samples in this experiment: 3 EGFP, 3 AFOS
SAMPLES=(1-0_S129 2-0_S130 3-0_S131 4-0_S132 5-0_S133 6-0_S134)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

module load bowtie2
module load samtools

mkdir -p "$OUT_DIR"/{trimmed,sam,bam}
cd "$OUT_DIR"

echo "=== $SAMPLE ==="

# 1. Trim Nextera adapters, quality trim both ends at q10, drop pairs < 30 bp
cutadapt \
    -q 10,10 \
    -a "$ADAPTER" \
    -A "$ADAPTER" \
    --minimum-length=30 \
    -o trimmed/"$SAMPLE".out.1.fastq \
    -p trimmed/"$SAMPLE".out.2.fastq \
    "$FQ_DIR/${SAMPLE}_R1_001.fastq" \
    "$FQ_DIR/${SAMPLE}_R2_001.fastq"

# 2. Align to mm10, allowing fragments up to 1 kb
bowtie2 \
    -x "$INDEX" \
    -1 trimmed/"$SAMPLE".out.1.fastq \
    -2 trimmed/"$SAMPLE".out.2.fastq \
    -S sam/"$SAMPLE".sam \
    -X 1000

# 3. SAM -> BAM, sort, index
samtools view -S -b sam/"$SAMPLE".sam > bam/"$SAMPLE".bam
samtools sort bam/"$SAMPLE".bam -o bam/"$SAMPLE".sorted.bam
samtools index -b bam/"$SAMPLE".sorted.bam

# Drop the intermediates once the sorted BAM exists
rm sam/"$SAMPLE".sam bam/"$SAMPLE".bam

echo "=== $SAMPLE done ==="
