#!/bin/bash
#SBATCH --job-name=integracion_scanpy
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --mem=32GB
#SBATCH --partition=intel_std

# Cargar el entorno Conda
source ~/miniforge3/etc/profile.d/conda.sh

# Activar el entorno Conda
conda activate single_cell_3.10

# Ejecutar el script Python
python /data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/scripts_complementarios/run_scanpy_integration.py
