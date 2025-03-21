#!/bin/bash
#SBATCH --job-name=crear_loom
#SBATCH --partition=gpu_a100
##SBATCH --cores-per-socket=3
#SBATCH --cpus-per-task=6
#SBATCH --gpus=a100:2
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --array=0-9
#SBATCH --output=/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/creation_loom.log


##############################
#NO FUNCIONA CON EL OUTPUT DE CELLRANGER 8.
##############################
# Cargar el entorno Conda
source ~/miniforge3/etc/profile.d/conda.sh

# Activar el entorno Conda
conda activate velocyto


# Definir la ruta del archivo de máscara y del archivo de anotación del genoma
MASK_FILE="/data/scratch/LAB/enric/Referencias_para_cellranger/repeat_mask_velocyto/repeat_mask.gtf"
GTF_FILE="/data/scratch/LAB/enric/Referencias_para_cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"

# Lista de rutas de las muestras
SAMPLES=(
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/P14-pancreas/outs/per_sample_outs/P14-pancreas/"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PT_17-final_2/"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PT_20-final_2/"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PT_22-final_2/"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PT_28-final_2/"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PT_49-final/"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PT_50-final/"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PIT_18-final/"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PIT_41-final/"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/P14-pancreas/"
)

# Ejecutar el comando correspondiente al índice de este trabajo en el array
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Procesando $SAMPLE"
velocyto run10x -m "$MASK_FILE" "$SAMPLE" "$GTF_FILE"

velocyto run10x -m /data/scratch/LAB/enric/Referencias_para_cellranger/repeat_mask_velocyto/repeat_mask.gtf /data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PT_14-final_2_para_loom /data/scratch/LAB/enric/Referencias_para_cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz

import loompy as lp
# Ruta al directorio de salida de Cell Ranger (el que contiene 'outs')
indir = "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PT_14-final_2/"

# Genoma humano
genome = "GRCh38"

# Crear el archivo .loom
loom_file = lp.create_from_cellranger(indir, outdir=None, genome=genome)

# Ruta al archivo .loom
loom_file_path = "/data/scratch/LAB/enric/Proyecto_pitagoras/Resultados_CellRanger/PT_14-final_2_para_loom/PT_14-final_2_para_loom.loom"

# Abrir el archivo .loom
with lp.connect(loom_file_path) as ds:
    # Verificar las dimensiones de los datos
    print("Número de células:", ds.shape[1])
    print("Número de genes:", ds.shape[0])

    # Verificar los metadatos de las células
    print("Metadatos de las células:", ds.ca)

    # Verificar los metadatos de los genes
    print("Metadatos de los genes:", ds.ra)

    # Verificar los datos de expresión génica
    print("Datos de expresión génica (primera célula, primeros 5 genes):", ds[:, 0][:5])