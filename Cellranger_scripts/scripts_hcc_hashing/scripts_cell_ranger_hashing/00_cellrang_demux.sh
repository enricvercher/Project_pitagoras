#!/bin/bash

#############################################################################################################################
## Single Cell Multi (CellRanger)
## evercheh@nasertic.es
#############################################################################################################################

## Initial SBATCH commands (In these lines we define the parameters to run the job in the cluster)
#SBATCH --job-name=CellRanger_multi
#SBATCH --mail-type=END
#SBATCH --mail-user=evercheh@nasertic.es
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH -p gpu_a100

## Lines to demultiplex the samples
module load singularity/3.4.1

# Ejecutar cellranger multi usando el contenedor de Singularity
singularity exec --nv --bind /data:/data /data/scratch/LAB/enric/TFM_enric/Contenedores/cellranger.sif \
cellranger multi \
--id=hashing_PIT_18_41 \
--csv=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Scripts/scripts_cellranger_8/scripts_hcc_hashing/multi_config_CSV/0_hashing_demux_config.csv \
--output-dir=/data/scratch/LAB/enric/Proyecto_pitagoras/hashing_PIT_18_41