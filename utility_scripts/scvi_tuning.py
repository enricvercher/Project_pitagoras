import os
import scanpy as sc
import numpy as np
import scvi
from ray import tune
from scvi import autotune

# Leer los datos
adata = sc.read_h5ad("/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/datos_integrados_scvi/datos_busq_hiperparams.h5ad")

# Configuración del modelo SCVI
model_cls = scvi.model.SCVI
model_cls.setup_anndata(adata, layer='counts_soupx_crude',
                        categorical_covariate_keys=['Sample'],
                        continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'])

# Definir el espacio de búsqueda
search_space = {
    "n_hidden": tune.choice([92, 128, 192, 256]),
    "n_latent": tune.choice([10, 20, 30, 40, 50, 60]),
    "n_layers": tune.choice([1, 2, 3]),
    "lr": tune.loguniform(1e-4, 1e-2),
    "gene_likelihood": tune.choice(["nb", "zinb"])
}

# Ejecutar el autotune
results = autotune.run_autotune(
    model_cls,
    data=adata,
    mode="min",  # Minimizar la métrica de pérdida de validación
    metrics="validation_loss",
    search_space=search_space,
    num_samples=25,  # Número de muestras de hiperparámetros
    resources={'cpu': 9, 'gpu': 0},  # Solo usar CPU
)

# Inicializar el valor más bajo de validación
best_vl = 10000
best_i = 0

# Iterar a través de los resultados
for i, res in enumerate(results.result_grid):
    vl = res.metrics['validation_loss']  # Accede a la métrica de pérdida de validación
    if vl < best_vl:
        best_vl = vl  # Actualiza si es menor que el mejor valor actual
        best_i = i  # Guarda el índice del mejor resultado

# Obtener los mejores hiperparámetros
best_result = results.result_grid[best_i]
print('Mejores hiperparámetros encontrados:', best_result)

# Guardar los resultados en un archivo
output_file = "/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/scvi_tuning_results.txt"
with open(output_file, "w") as f:
    f.write(f"Mejores hiperparámetros: {best_result}\n")
