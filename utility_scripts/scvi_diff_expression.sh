#!/bin/bash
#SBATCH --job-name=diff_expresion
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --partition=intel_std
#SBATCH --output=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/analisis_diferencia_expresion/scvi_DE_leiden.log

# Cargar el módulo de Singularity
module load singularity/3.4.1

# Definir las rutas
INPUT_FILE="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/firmas_transcriptomicas/Firmas_AUC_ATLAS.h5ad"
MODEL_PATH="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/mod_scvi_all_genes_recommended/vae_model"
OUTPUT_FILE="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/analisis_diferencia_expresion/markers_leiden_res1.pkl"
CONTAINER_PATH="/data/scratch/LAB/enric/TFM_enric/Contenedores/scvi-tools.sif"

# Crear el directorio de salida si no existe
mkdir -p $(dirname "${OUTPUT_FILE}")

# Ejecutar el código dentro del contenedor de Singularity
singularity exec --nv --bind /data:/data $CONTAINER_PATH \
python -c "
import anndata as ad
import scvi
import torch
import pandas as pd

# Verificar el uso de la GPU
print('CUDA disponible:', torch.cuda.is_available())
if torch.cuda.is_available():
    print('Usando GPU:', torch.cuda.get_device_name(0))
else:
    print('Usando CPU')

# Cargar el objeto AnnData
adata = ad.read_h5ad('${INPUT_FILE}')

# Cambiar 'X' para usar 'counts_soupx_crude'
adata.X = adata.layers['counts_soupx_crude'].copy()

# Verificar que 'leiden' existe en los datos cargados
print('Columnas disponibles en adata.obs:', adata.obs.columns)
if 'leiden' not in adata.obs.columns:
    raise KeyError('leiden no está presente en adata.obs')

# Cargar el modelo entrenado desde el contenedor
model = scvi.model.SCVI.load('${MODEL_PATH}')

# Asociar el objeto AnnData al modelo cargado
model.adata = adata

# Configurar SCVI con la clave de agrupamiento 
scvi.model.SCVI.setup_anndata(adata, labels_key='leiden')

# Realizar la expresión diferencial
markers_scvi = model.differential_expression(groupby='leiden')

# Guardar los resultados de la expresión diferencial en un archivo
markers_scvi.to_pickle('${OUTPUT_FILE}')
"

# Comprobar si el archivo de salida se generó correctamente
if [ ! -f "${OUTPUT_FILE}" ]; then
    echo "Error: Los resultados de la expresión diferencial no se encontraron en ${OUTPUT_FILE}"
    exit 1
fi