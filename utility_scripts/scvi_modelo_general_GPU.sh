#!/bin/bash
#SBATCH --job-name=training_solo_combined_array
#SBATCH --partition=gpu_a100
#SBATCH --cores-per-socket=5
#SBATCH --cpus-per-task=5
#SBATCH --gpus=a100:1
#SBATCH --mem=20G
#SBATCH --time=3-00:00:00
#SBATCH --output=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/modelo_general_scvi/modelo_general.log        

# Cargar el módulo de Singularity
module load singularity/3.4.1

# Definir las rutas
INPUT_PATH="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/processed_data_for_scvi/adata_mixto_processed.h5ad"
OUTPUT_DIR="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/modelo_general_scvi"
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

# Configurar SCVI con la capa 'counts_soupx_crude' y las covariables
scvi.model.SCVI.setup_anndata(adata, layer='counts_soupx_crude',
                             categorical_covariate_keys=['Sample'],
                             continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'])

# Inicializar el modelo SCVI
model = scvi.model.SCVI(adata)

# Entrenar el modelo
model.train()  

# Guardar el modelo entrenado en la carpeta especificada
model.save('${OUTPUT_DIR}/vae_model/', save_anndata=True, overwrite=True)
"

# Comprobar si el entrenamiento se completó correctamente
if [ ! -f "${OUTPUT_DIR}/vae_model/model.pt" ]; then
    echo "Error: El archivo del modelo SCVI no se encontró en ${OUTPUT_DIR}/vae_model/model.pt"
    exit 1
fi
