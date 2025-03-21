# Identification of Tumor-Specific TCRs Based on the Transcriptomic Signature and Tumor Localization of TILs

This repository contains an ongoing single-cell analysis focused on identifying tumor-specific TCRs by integrating transcriptomic signatures and tumor localization of tumor-infiltrating lymphocytes (TILs). The study is based on **10 patients** with different tumor types:

- **5 ovarian cancer patients**
- **4 hepatocellular carcinoma (HCC) patients**
- **1 pancreatic ductal adenocarcinoma (PDAC) patient**

**⚠️ This analysis is still in progress and subject to continuous updates.**

---

## 📂 Repository Structure

### 🔹 **Jupyter Notebooks**
The main analysis is structured in multiple Jupyter notebooks, covering various aspects of **single-cell TCR and transcriptomic data analysis**, including:

- **Quality control and filtering**
  - `00_Ambient_RNA_Correction.ipynb`
  - `00_Analysis_Doublet_SOLO.ipynb`
  - `00_Analysis_sc_QC_filtering.ipynb`

- **Evaluation of Single-Cell Automatic Annotation and Integration Methods**
  - `01_Analysis_sc_integration_cluster.ipynb`
  - `02_Analysis_sc_integration_scANVI.ipynb`
  - `02_Celltypist_Annotation.ipynb`
  - `03_Annotation_signatures_SH-TCR_10pts.ipynb`
  - `03_Atlas_ProjecTIL.ipynb`

- **TCR repertoire and clonality analysis, avidity assessment and trajectory inference**
  - `05_Analysis_Repertoire_clonality.ipynb`
  - `05_Study_Avidity_TCR.ipynb`
  - `05_Study_TCRdist3.ipynb`
  - `05_Study_TCR_Trajectory.ipynb`

- **GEX & TCR integration**
  - `06_Analysis_GEX_VDJ_signatures_integration.ipynb`
  - `07_Analysis_GEX_VDJ_KNN.ipynb`
  - `08_Analysis_Leiden_DE.ipynb`
  - `09_Final_Study_GEX_VDJ.ipynb`

- **Figures & Data export**
  - `10_Figures_ppt_final_presentation.ipynb`
    
- **Interplay between scanpy and Seurat objects**
  - `Z_Transform_from_scanpy_to_Seurat_R.ipynb`

---

### 🔹 **Spatial-TCR Analysis**
The `Spatial-TCR/` subdirectory contains ongoing analyses integrating **spatial transcriptomics and TCR sequencing**, aiming to **map tumor-reactive TCRs to their spatial location within the tumor microenvironment**.

This analysis is based on **Hudson et al.** *Distinct phenotypic states and spatial distribution of CD8+ T cell clonotypes in human brain metastases* (PMID: 35584630) and serves as a **technical validation and pipeline setup** for spatial TCR analysis.

- **Notebooks:**
  - `00_Analysis_spatial_hudson.ipynb` → Processing spatial transcriptomics datasets.
  - `01_Analysis_TCR_Repertoire_Seurat.ipynb` → Integrating spatial data with Seurat.
    
- **Raw data processing scripts** (inside `utility_scripts/`):
  - `download_scRNAseq.sh`
  - `download_spati_scRNAseq.sh`
  - `download_spati_TCR.sh`

---

### 🔹 **Cell Ranger Pipeline**
The `Cellranger_scripts/` directory contains scripts used to preprocess **single-cell RNA-seq and TCR-seq data** from 10 patients using **Cell Ranger** for demultiplexing, alignment, and feature counting.

---

### 🔹 **TCR Repertoire Analysis (Immcantation)**
The `Immcantation_repertoire_scripts/` directory includes **sample sheets and commands** for running the **Immcantation framework**, a powerful pipeline for TCR/BCR repertoire analysis.

---

## 🚧 **Work in Progress**
This project is still under active development. Some analyses are being optimized, and additional validation is required to confirm the tumor specificity of the identified TCRs. Future updates will include:

- **Refinement of single-cell annotation based on differentially expressed genes (DEGs) and literature-derived signatures**
- **Additional spatial transcriptomic datasets**
- **Functional avidity study using different methods and TCRdist-based physicochemical classification of TCRs**

---

## 📢 **Contributors & Contact**
For any questions or collaborations, please reach out via GitHub Issues or contact **evercheh@nasertic.es**.
