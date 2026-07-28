#!/usr/bin/env bash
# Step 0 - Raw data download and Cell Ranger processing (bash)
# Auto-derived from the methods paper (Part 1). Run OUTSIDE the container, on a
# host / server / HPC where Cell Ranger 9.0.1 is installed -- the
# matigara/scrnaseq image does NOT bundle Cell Ranger. Run every command from
# the cloned repository directory and write outputs into it, so they are visible
# inside the container (which mounts that directory as /workspace).
#
# This is a reference script: the Cell Ranger download URL is signed and
# time-limited, so it must be pasted by hand (see the first block). Review each
# block before running.

# ==============================================================================
# Install Cell Ranger 9.0.1
# ==============================================================================
# Proprietary 10x Genomics software. Get the signed curl command for 9.0.1 from
# https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions
# (log in, accept the license, copy its curl). The URL carries an Expires=
# timestamp tied to your account, so it cannot be hard-coded here.
#
# curl -o cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=...&Key-Pair-Id=...&Signature=..."

tar -xzvf cellranger-9.0.1.tar.gz
export PATH="$PWD/cellranger-9.0.1/bin:$PATH"
cellranger --version

# ==============================================================================
# Download FASTQ (CRA010863) and rename to the Cell Ranger convention
# ==============================================================================
# Public data from Han et al., 2023 (Nature Plants), BioProject PRJCA016521,
# run accession CRA010863: ScWT = CRR775298, Scpifq = CRR775297.
mkdir -p /workspace/fastq/CRA010863

wget -c -P /workspace/fastq/CRA010863 \
  ftp://download.big.ac.cn/gsa2/CRA010863/CRR775298/CRR775298_f1.fastq.gz
wget -c -P /workspace/fastq/CRA010863 \
  ftp://download.big.ac.cn/gsa2/CRA010863/CRR775298/CRR775298_r2.fastq.gz
wget -c -P /workspace/fastq/CRA010863 \
  ftp://download.big.ac.cn/gsa2/CRA010863/CRR775297/CRR775297_f1.fastq.gz
wget -c -P /workspace/fastq/CRA010863 \
  ftp://download.big.ac.cn/gsa2/CRA010863/CRR775297/CRR775297_r2.fastq.gz

cd /workspace/fastq/CRA010863
ln -s CRR775298_f1.fastq.gz ScWT_S1_L001_R1_001.fastq.gz
ln -s CRR775298_r2.fastq.gz ScWT_S1_L001_R2_001.fastq.gz
ln -s CRR775297_f1.fastq.gz Scpifq_S1_L001_R1_001.fastq.gz
ln -s CRR775297_r2.fastq.gz Scpifq_S1_L001_R2_001.fastq.gz

# ==============================================================================
# Build the Cell Ranger reference (TAIR10.1 + Araport11)
# ==============================================================================
# The genome and annotation ship compressed under data/, so no download is
# needed. Run these from the cloned repository root.
cd /workspace
gunzip -k data/Arabidopsis_thaliana_TAIR10.1_Araport11_genome.fasta.gz
gunzip -k data/Arabidopsis_thaliana_TAIR10.1_Araport11_annotation.gtf.gz

cellranger mkref \
  --genome=Ara \
  --fasta=data/Arabidopsis_thaliana_TAIR10.1_Araport11_genome.fasta \
  --genes=data/Arabidopsis_thaliana_TAIR10.1_Araport11_annotation.gtf

# ==============================================================================
# Cell Ranger count for each sample
# ==============================================================================
# Use the resulting filtered_feature_bc_matrix/ path in the sample manifest in
# Step 1, Section 1.
cellranger count --localcores=80 --id=ScWT --fastqs=/workspace/fastq/CRA010863 \
  --sample=ScWT --transcriptome=Ara --create-bam false

cellranger count --localcores=80 --id=Scpifq --fastqs=/workspace/fastq/CRA010863 \
  --sample=Scpifq --transcriptome=Ara --create-bam false
