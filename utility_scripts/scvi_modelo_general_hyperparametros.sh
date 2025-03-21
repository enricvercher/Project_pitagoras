#!/bin/bash
#SBATCH --job-name=integracion_scanpy
#SBATCH --cpus-per-task=9
#SBATCH --time=24:00:00
#SBATCH --mem=45GB
#SBATCH --partition=intel_std

# Cargar el entorno Conda
source ~/miniforge3/etc/profile.d/conda.sh

# Activar el entorno Conda
conda activate scvi-env

# Ejecutar el script Python
python /data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/scripts_complementarios/scvi_tuning.py

conda deactivate
