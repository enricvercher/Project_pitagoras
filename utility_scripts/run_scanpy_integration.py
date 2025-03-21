import scanpy as sc
import os
import scipy

# Definir la ruta de salida
output_dir = "/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/integracion_scanpy"
os.makedirs(output_dir, exist_ok=True)

# Función para procesar los datos
def process_data(adata, sample_name):
    # Asegurarte de que adata.X contiene los datos crudos desde counts_soupx_crude
    if 'counts_soupx_crude' in adata.layers:
        adata.X = adata.layers['counts_soupx_crude'].copy()  # Asignar capa cruda a X
        print(f"Asignando counts_soupx_crude a X para el sample {sample_name}")
    else:
        raise ValueError(f"La capa 'counts_soupx_crude' no está disponible en adata para {sample_name}.")
    
    # Normalización
    sc.pp.normalize_total(adata, target_sum=1e4)  # Normalizar cada célula a 10,000 UMI
    
    # Aplicar la transformación logarítmica
    sc.pp.log1p(adata)  # Aplicar log-transformación (log1p)
    
    # Identificación de genes altamente variables
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key = 'Sample')
    
    adata.raw = adata  # Guardar los datos sin procesar antes de hacer más cambios
    adata = adata[:, adata.var.highly_variable]  # Filtrar genes altamente variables

    # Regresión y escalado
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])  # Eliminar efectos de total counts y genes mitocondriales
    sc.pp.scale(adata, max_value=10)  # Escalar los datos
    
    # PCA, vecinos y UMAP
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(adata, resolution=0.25)
    sc.tl.umap(adata)

    # Guardar el objeto procesado
    adata.write(os.path.join(output_dir, f"{sample_name}_processed.h5ad"))

# Cargar y procesar los datos de cada paciente
pt_14 = sc.read_h5ad("/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/filtered_matrices_qc/pt_14_filtered_QC_mixto.h5ad")
pt_17 = sc.read_h5ad("/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/filtered_matrices_qc/pt_17_filtered_QC_mixto.h5ad")
pt_20 = sc.read_h5ad("/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/filtered_matrices_qc/pt_20_filtered_QC_mixto.h5ad")
pt_22 = sc.read_h5ad("/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/filtered_matrices_qc/pt_22_filtered_QC_mixto.h5ad")
pt_28 = sc.read_h5ad("/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Results/filtered_matrices_qc/pt_28_filtered_QC_mixto.h5ad")

# Procesar los datos asegurando la capa cruda se use en X
process_data(pt_14, "pt_14")
process_data(pt_17, "pt_17")
process_data(pt_20, "pt_20")
process_data(pt_22, "pt_22")
process_data(pt_28, "pt_28")

print("Procesamiento completado.")
