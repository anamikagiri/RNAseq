#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --output=logs/%x-%j.log
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --partition=medium
#SBATCH --time=07:00:00

STAR=/data/cephfs-1/home/users/agiri_m/work/miniforge3/bin/STAR
FASTA=/data/cephfs-1/home/users/agiri_m/work/RNAseq/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa
GTF=/data/cephfs-1/home/users/agiri_m/work/RNAseq/references/Homo_sapiens.GRCh38.110.gtf 
OUT=/data/cephfs-1/home/users/agiri_m/work/RNAseq/references/star_index

srun "$STAR" --runThreadN 2 \
  --runMode genomeGenerate \
  --genomeDir "$OUT" \
  --genomeFastaFiles "$FASTA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 100

wait
echo "`date`: STAR index is generated"
