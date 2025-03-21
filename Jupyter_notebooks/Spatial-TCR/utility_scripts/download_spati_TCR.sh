#!/bin/bash
# Paciente 15
fastq-dump --split-files --gzip --outdir /data/scratch/LAB/enric/zz_Spatial_Seq_Hudosn/Patient_15/Spatial_TCR/Raw_data/ SRR15052442 &

# Paciente 16
fastq-dump --split-files --gzip --outdir /data/scratch/LAB/enric/zz_Spatial_Seq_Hudosn/Patient_16/Spatial_TCR/Raw_data/ SRR15052441 &

# Paciente 19
fastq-dump --split-files --gzip --outdir /data/scratch/LAB/enric/zz_Spatial_Seq_Hudosn/Patient_19/Spatial_TCR/Raw_data/ SRR15052440 &

# Paciente 24
fastq-dump --split-files --gzip --outdir /data/scratch/LAB/enric/zz_Spatial_Seq_Hudosn/Patient_24/Spatial_TCR/Raw_data/ SRR15052439 &

# Paciente 26
fastq-dump --split-files --gzip --outdir /data/scratch/LAB/enric/zz_Spatial_Seq_Hudosn/Patient_26/Spatial_TCR/Raw_data/ SRR15052438 &

# Paciente 27
fastq-dump --split-files --gzip --outdir /data/scratch/LAB/enric/zz_Spatial_Seq_Hudosn/Patient_27/Spatial_TCR/Raw_data/ SRR15052437 &
