#!/bin/bash
#SBATCH --job-name=CellRanger_multi_VDJ
#SBATCH --output=CellRanger_multi_VDJ_%A_%a.out
#SBATCH --error=CellRanger_multi_VDJ_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --array=0-5
#SBATCH --partition=intel_std

# Cargar el módulo de singularity
module load singularity/3.4.1

# Definir los comandos en un array de Bash
commands=(
    "singularity exec --bind /data:/data /data/scratch/LAB/enric/TFM_enric/Contenedores/cellranger.sif cellranger multi --id=PIT_18-final --csv=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Scripts/scripts_cellranger_8/scripts_hcc_hashing/scripts_cell_ranger_hashing/PIT_18_final_Multiplex_config.csv"
    "singularity exec --bind /data:/data /data/scratch/LAB/enric/TFM_enric/Contenedores/cellranger.sif cellranger multi --id=PIT_41-final --csv=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Scripts/scripts_cellranger_8/scripts_hcc_hashing/scripts_cell_ranger_hashing/PIT_41_final_Multiplex_config.csv"
)
# Ejecutar el comando correspondiente al índice de este trabajo en el array
eval "${commands[$SLURM_ARRAY_TASK_ID]}"

