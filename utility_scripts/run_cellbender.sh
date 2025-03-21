#!/bin/bash
#SBATCH --job-name=cellbender
#SBATCH --partition=intel_std
#SBATCH --cpus-per-task=5
#SBATCH --mem=30G             # Aumentamos a 30 GB por tarea
#SBATCH --time=3-00:00:00
#SBATCH --array=0-4%3         # Ejecutar un máximo de 3 tareas en paralelo

# Inicializar Conda
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate cb

# Crear el directorio de salida si no existe
output_dir="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/matrices_corregidas_cellbender/"
mkdir -p "$output_dir"

# Rutas de los archivos sin filtrar proporcionados por CellRanger
raw_files=(
    "/data/scratch/LAB/enric/Proyecto_pitagoras/PT_14-final_2/outs/multi/count/raw_feature_bc_matrix.h5"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/PT_17-final_2/outs/multi/count/raw_feature_bc_matrix.h5"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/PT_20-final_2/outs/multi/count/raw_feature_bc_matrix.h5"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/PT_22-final_2/outs/multi/count/raw_feature_bc_matrix.h5"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/PT_28-final_2/outs/multi/count/raw_feature_bc_matrix.h5"
)

# Nombres base para los pacientes
pt_names=(
    "PT_14"
    "PT_17"
    "PT_20"
    "PT_22"
    "PT_28"
)

# Obtener el índice del array actual
i=$SLURM_ARRAY_TASK_ID

# Procesar el archivo correspondiente al índice
file="${raw_files[$i]}"
pt_name="${pt_names[$i]}"

cellbender remove-background \
    --input "$file" \
    --output "$output_dir/${pt_name}_denoised.h5" \

# Desactivar el entorno Conda
conda deactivate
