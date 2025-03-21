#!/bin/bash
#SBATCH --job-name=scvi_integrar
#SBATCH --partition=gpu_a100
#SBATCH --cores-per-socket=12
#SBATCH --cpus-per-task=12
#SBATCH --gpus=a100:1
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00
#SBATCH --output=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/mod_scvi_all_genes_recommended_10pt/modelo_all_recommended_10pt.log

##################si esta disponible la GPU################ 

##SBATCH --job-name=integracion_scanpy
##SBATCH --cpus-per-task=8
##SBATCH --time=6:00:00
##SBATCH --mem=32GB
##SBATCH --partition=intel_std
##SBATCH --output=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/mod_scvi_all_genes_recommended/modelo_all_recommended.log


# Cargar el módulo de Singularity
module load singularity/3.4.1

# Definir las rutas
INPUT_PATH="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/04_datos_concatenados/adata_conTCR_GEX_firmas.h5ad"
OUTPUT_DIR="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/mod_scvi_all_genes_recommended_10pt_conTCR"
CONTAINER_PATH="/data/scratch/LAB/enric/TFM_enric/Contenedores/scvi-tools.sif"  

# Crear el directorio de salida
mkdir -p $OUTPUT_DIR

# Ejecutar el código dentro del contenedor de Singularity
singularity exec --nv --bind /data:/data $CONTAINER_PATH \
python -c "
import anndata as ad
import torch
import scvi

# Ajustar la precisión de la multiplicación de matrices para aprovechar Tensor Cores
torch.set_float32_matmul_precision('medium')  # Puedes cambiar a 'high' si prefieres más rendimiento

# Comprobar si CUDA está disponible y se está utilizando
print('CUDA disponible:', torch.cuda.is_available())
print('Dispositivo:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU')

# Cargar el objeto AnnData preprocesado
adata = ad.read_h5ad('${INPUT_PATH}')

epochs = 50

# Configurar SCVI con la capa 'counts_soupx_crude' y las covariables
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts_soupx_crude',
    categorical_covariate_keys=['Sample'],
    continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'],
    labels_key='ATLAS_TIL'  # Importante para scANVI
)

# Inicializar el modelo SCVI con los parametros recomendados.
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood='nb')

# Entrenar el modelo
#model.train()  
model.train(check_val_every_n_epoch=1, max_epochs=epochs, early_stopping=True)

# Guardar el modelo entrenado en la carpeta especificada
model.save('${OUTPUT_DIR}/vae_model/', save_anndata=True, overwrite=True)
"

# Comprobar si el entrenamiento se completó correctamente
if [ ! -f "${OUTPUT_DIR}/vae_model/model.pt" ]; then
    echo "Error: El archivo del modelo SCVI no se encontró en ${OUTPUT_DIR}/vae_model/model.pt"
    exit 1
fi
