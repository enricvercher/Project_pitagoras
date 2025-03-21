#!/bin/bash
#SBATCH --job-name=training_solo_
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --mem=32GB
#SBATCH --partition=intel_std
#SBATCH --array=0-4  # Actualizado para 5 muestras

# Cargar el módulo de Singularity
module load singularity/3.4.1

# Definir las rutas a los archivos preprocesados
INPUT_FILES=(
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/pred_dob_Soupx_clean/process_for_doublets/PT_14_preprocessed_SoupX_for_doublets.h5ad"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/pred_dob_Soupx_clean/process_for_doublets/PT_17_preprocessed_SoupX_for_doublets.h5ad"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/pred_dob_Soupx_clean/process_for_doublets/PT_20_preprocessed_SoupX_for_doublets.h5ad"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/pred_dob_Soupx_clean/process_for_doublets/PT_22_preprocessed_SoupX_for_doublets.h5ad"
    "/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/pred_dob_Soupx_clean/process_for_doublets/PT_28_preprocessed_SoupX_for_doublets.h5ad"
)

PATIENTS_NAMES=("PT_14"
"PT_17"
"PT_20"
"PT_22"
"PT_28"
)

INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}
PATIENT_NAME=${PATIENTS_NAMES[$SLURM_ARRAY_TASK_ID]}

# Cambiar la ruta de salida para los modelos
OUTPUT_MODEL_DIR="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/pred_dob_Soupx_clean/modelos/${PATIENT_NAME}"
OUTPUT_PREDICTION_DIR="/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/pred_dob_Soupx_clean/predicciones_dobletes/${PATIENT_NAME}"

# Crear los directorios de salida si no existen
mkdir -p $OUTPUT_MODEL_DIR
mkdir -p $OUTPUT_PREDICTION_DIR

# Ejecutar el código dentro del contenedor de Singularity
singularity exec --nv --bind /data:/data /data/scratch/LAB/enric/TFM_enric/Contenedores/scvi-tools.sif \
python -c "
import scvi
import pandas as pd

# Cargar los datos preprocesados del paciente
adata = scvi.data.read_h5ad('${INPUT_FILE}')

# Configurar SCVI y entrenar el modelo
scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)
vae.train()

# Guardar el modelo entrenado
vae.save('${OUTPUT_MODEL_DIR}/vae_model/', save_anndata=True, overwrite=True)

# Cargar el modelo entrenado para usar SOLO
vae = scvi.model.SCVI.load('${OUTPUT_MODEL_DIR}/vae_model/')

# Crear y entrenar el modelo SOLO
solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train()

# Realizar predicciones
df = solo.predict()
df['prediction'] = solo.predict(soft=False)

# Guardar las predicciones en un archivo CSV
df.to_csv('${OUTPUT_PREDICTION_DIR}/solo_predictions_soupX.csv', index=True)
"

# Comprobar si SOLO se ejecutó correctamente
if [ ! -f "${OUTPUT_PREDICTION_DIR}/solo_predictions_soupX.csv" ]; then
    echo "Error: El archivo de predicciones SOLO no se encontró en ${OUTPUT_PREDICTION_DIR}/solo_predictions_soupX.csv"
    exit 1
fi
