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

singularity exec --nv --bind /data:/data /data/scratch/LAB/enric/TFM_enric/Contenedores/cellranger.sif \
cellranger mkfastq \
--id=Pitagoras_HCC_demux \
--run=/data/scratch/LAB/enric/Proyecto_pitagoras/241127_A01215_0166_BHLN5NDRX5 \
--samplesheet=/data/scratch/LAB/enric/Proyecto_pitagoras/241127_A01215_0166_BHLN5NDRX5/SampleSheet_final_modified.csv \
--output-dir=/data/scratch/LAB/temp_demultiplex/241127_sc_PITAGORAS_HCC_demux/ 