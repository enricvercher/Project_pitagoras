#!/bin/bash
#SBATCH --job-name=CellRanger_multi_VDJ
#SBATCH --output=CellRanger_multi_VDJ_%A_%a.out
#SBATCH --error=CellRanger_multi_VDJ_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --array=0-1
#SBATCH --partition=intel_std

# Cargar el módulo de singularity
module load singularity/3.4.1

# Definir los comandos en un array de Bash, uno para cada paciente
commands=(
    "singularity exec --bind /data:/data /data/scratch/LAB/enric/TFM_enric/Contenedores/cellranger.sif cellranger multi --id=PT_49-final --csv=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Scripts/scripts_cellranger_8/scripts_cell_ranger_demux/PIT_49.csv"
    "singularity exec --bind /data:/data /data/scratch/LAB/enric/TFM_enric/Contenedores/cellranger.sif cellranger multi --id=PT_50-final --csv=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Scripts/scripts_cellranger_8/scripts_cell_ranger_demux/PIT_50.csv"
)

# Ejecutar el comando correspondiente al índice de este trabajo en el array
eval "${commands[$SLURM_ARRAY_TASK_ID]}"

