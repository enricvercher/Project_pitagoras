#!/bin/bash
#SBATCH --job-name=integracion_scvi
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=32GB
#SBATCH --partition=intel_std
#SBATCH --output=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/modelo_integrar_scvi/modelo_int_scvi.log

###################si esta disponible la GPU################  
##SBATCH --job-name=scvi_integrar
##SBATCH --partition=gpu_a100
##SBATCH --cores-per-socket=5
##SBATCH --cpus-per-task=5
##SBATCH --gpus=a100:1
##SBATCH --mem=20G
##SBATCH --time=3-00:00:00
##SBATCH --output=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/modelo_integracion_scvi/modelo_int_scvi.log

# Cargar el módulo de Singularity
module load singularity/3.4.1

# Definir las rutas
INPUT_PATH="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/processed_data_for_scvi/adata_scvi_3000_genes.h5ad"
OUTPUT_DIR="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/modelo_integrar_scvi"
CONTAINER_PATH="/data/scratch/LAB/enric/TFM_enric/Contenedores/scvi-tools.sif"

# Crear el directorio de salida
mkdir -p $OUTPUT_DIR

# Ejecutar el código dentro del contenedor de Singularity
singularity exec --nv --bind /data:/data $CONTAINER_PATH \
python -c '
import anndata as ad
import scvi

# Cargar el objeto AnnData preprocesado
adata = ad.read_h5ad("'"${INPUT_PATH}"'")

# Configurar SCVI con la capa "counts_soupx_crude" y las covariables
scvi.model.SCVI.setup_anndata(adata, layer="counts_soupx_crude",
                             categorical_covariate_keys=["Sample"],
                             continuous_covariate_keys=["pct_counts_mt", "total_counts", "pct_counts_ribo"])

# Inicializar y entrenar el modelo SCVI por defecto
model_default = scvi.model.SCVI(adata)
model_default.train()  

# Inicializar y entrenar el modelo SCVI con parámetros recomendados
model_recommended = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model_recommended.train()  

# Guardar los modelos entrenados
model_default.save("'"${OUTPUT_DIR}"'/vae_model_default/", save_anndata=True, overwrite=True)
model_recommended.save("'"${OUTPUT_DIR}"'/vae_model_recommended/", save_anndata=True, overwrite=True)
'

# Comprobar si el entrenamiento se completó correctamente
if [ ! -f "${OUTPUT_DIR}/vae_model_default/model.pt" ]; then
    echo "Error: El archivo del modelo SCVI por defecto no se encontró en ${OUTPUT_DIR}/vae_model_default/model.pt"
    exit 1
fi

if [ ! -f "${OUTPUT_DIR}/vae_model_recommended/model.pt" ]; then
    echo "Error: El archivo del modelo SCVI recomendado no se encontró en ${OUTPUT_DIR}/vae_model_recommended/model.pt"
    exit 1
fi