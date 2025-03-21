#!/bin/bash
#SBATCH --job-name=scvi_integrar
#SBATCH --partition=gpu_a100
#SBATCH --cores-per-socket=5
#SBATCH --cpus-per-task=5
#SBATCH --gpus=a100:1
#SBATCH --mem=20G
#SBATCH --time=3-00:00:00
#SBATCH --output=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/mod_sANVI_all_genes_recommended_10pt_conTCR/modelo_SCANVI_10pt_conTCR.log

##################si esta disponible la GPU################ 

##SBATCH --job-name=integracion_scanpy
##SBATCH --cpus-per-task=8
##SBATCH --time=6:00:00
##SBATCH --mem=32GB
##SBATCH --partition=intel_std
##SBATCH --output=/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/mod_scvi_all_genes_recommended/modelo_all_recommended.log


# -------------------------------------------------------------------------
# 1) Cargar el módulo de Singularity (opcional, según tu entorno)
# -------------------------------------------------------------------------
module load singularity/3.4.1

# -------------------------------------------------------------------------
# 2) Definir rutas de entrada/salida y contenedor
# -------------------------------------------------------------------------
INPUT_ADATA="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/04_datos_concatenados/adata_conTCR_GEX_firmas.h5ad"
SCVI_MODEL_PATH="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/mod_scvi_all_genes_recommended_10pt_conTCR/vae_model"
OUTPUT_DIR="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/modelos/mod_sANVI_all_genes_recommended_10pt_conTCR"
CONTAINER_PATH="/data/scratch/LAB/enric/TFM_enric/Contenedores/scvi-tools.sif"

# Crear directorio de salida
mkdir -p "$OUTPUT_DIR"

# -------------------------------------------------------------------------
# 3) Ejecutar el código Python dentro del contenedor con GPU
# -------------------------------------------------------------------------
singularity exec --nv --bind /data:/data "$CONTAINER_PATH" \
python -c "
import anndata as ad
import torch
import scvi
import pandas as pd

# A) Comprobar CUDA
print('CUDA disponible:', torch.cuda.is_available())
if torch.cuda.is_available():
    print('Dispositivo:', torch.cuda.get_device_name(0))
else:
    print('Usando CPU')

# B) Cargar el AnnData
adata = ad.read_h5ad('${INPUT_ADATA}')
print('Dimensiones adata:', adata.shape)

# C) Reemplazar los NaN de la columna ATLAS_TIL por 'unlabeled'
col_etiquetas = 'ATLAS_TIL'
unlabeled_cat = 'unlabeled'

if col_etiquetas not in adata.obs.columns:
    raise ValueError(f'La columna {col_etiquetas} no existe en adata.obs')

# Corrección del warning de dtype
if isinstance(adata.obs[col_etiquetas].dtype, pd.CategoricalDtype):
    # Verificar si 'unlabeled' ya existe antes de añadirla
    if unlabeled_cat not in adata.obs[col_etiquetas].cat.categories:
        adata.obs[col_etiquetas] = adata.obs[col_etiquetas].cat.add_categories([unlabeled_cat])
    
    adata.obs[col_etiquetas] = adata.obs[col_etiquetas].fillna(unlabeled_cat)  # Sin inplace
else:
    adata.obs[col_etiquetas] = adata.obs[col_etiquetas].astype(str).fillna(unlabeled_cat)

n_unlabeled = (adata.obs[col_etiquetas] == unlabeled_cat).sum()
print(f'Celdas con etiqueta \"{unlabeled_cat}\":', n_unlabeled)

# D) Configurar scVI incluyendo labels_key
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts_soupx_crude',
    categorical_covariate_keys=['Sample'],
    continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'],
    labels_key='ATLAS_TIL'  # Añadido para asociar las etiquetas
)

# E) Cargar el modelo scVI pre-entrenado
scvi_model = scvi.model.SCVI.load('${SCVI_MODEL_PATH}', adata=adata)

# F) Crear un modelo scANVI desde scVI (añadiendo dropout_rate aquí)
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    unlabeled_category='unlabeled',
    dropout_rate=0.2  # Aplicar dropout aquí
)

# G) Entrenar scANVI (sin dropout_rate en plan_kwargs)
scanvi_model.train(
    max_epochs=50,
    check_val_every_n_epoch=1,
    early_stopping=True,
    plan_kwargs={
        'lr': 1e-3,
        'weight_decay': 1e-4
    }
)

# H) Guardar el modelo SCANVI y el AnnData
scanvi_model.save('${OUTPUT_DIR}/', overwrite=True, save_anndata=True)
print(f'SCANVI guardado en ${OUTPUT_DIR}/')
"

# -------------------------------------------------------------------------
# 4) Verificar que el modelo se guardó
# -------------------------------------------------------------------------
if [ ! -f "${OUTPUT_DIR}/model.pt" ]; then
    echo 'Error: No se encontró ${OUTPUT_DIR}/model.pt'
    exit 1
fi

echo 'Entrenamiento SCANVI completado con éxito. Ubicación del modelo: ${OUTPUT_DIR}/'
