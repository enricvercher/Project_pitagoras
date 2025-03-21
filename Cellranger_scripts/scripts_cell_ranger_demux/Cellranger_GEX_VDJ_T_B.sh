#!/bin/bash
#SBATCH --job-name=CellRanger_single_P14
#SBATCH --output=CellRanger_single_P14_%j.out
#SBATCH --error=CellRanger_single_P14_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --partition=intel_std

# Cargar el módulo de singularity
module load singularity/3.4.1

# Ejecutar cellranger para el paciente P14
singularity exec --bind /data:/data /data/scratch/LAB/enric/TFM_enric/Contenedores/cellranger.sif \
cellranger multi --id=P14-pancreas \
--csv=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Scripts/scripts_cellranger_8/scripts_cell_ranger_demux/P14.csv


