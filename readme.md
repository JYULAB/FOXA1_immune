# Code for FOXA1 Paper
Citation:

This repository includes the main pieces of code used to support the conclusions in the paper. The following code were used in a SLURM High Performance Computing environment, using a combination of bash scripting, R (4.0.3) and python >=3.6. Detailed versions of softwares used are listed in Methods section of the paper.

The scripts are organized by technology/analysis type. If any scripts are to be run in order, prefix numbers are provided in the script name. Otherwise, each script is mostly standalone. Some files are provided here, such as the final cells passing filtering that were used in Multiome analysis. Expected output of each script are the indicated figures from the paper.
 
### ChIP-seq

- [**command.sh**](ChIP-seq/command.sh): creating heatmaps, Fig 4D, 4E, 4H, 5E, 7D, S5E, S7E, S9A, S9B (same tool is used for ATAC-seq heatmaps).

### scRNA-seq

- [**01_d0_single.R**](Multiome/01_d0_single.R),  [**02_d14_single.R**](Multiome/02_d14_single.R),  [**03_d21_single.R**](Multiome/03_d21_single.R): Individual timepoint processing and filtering
- [**04_GEX_integration.R**](Multiome/04_GEX_integration.R): Integration and analysis of scRNA-Seq, Fig 3C, D, E, F, Fig S3F 
- [**05_ATAC_integration.R**](Multiome/05_ATAC_integration.R): Integration and analysis of scATAC-Seq, Fig 3H, I, J
- [**clinical_patients_intergration.R**](Multiome/clinical_patients_intergration.R): Integration with clinical patients and pseudotime analysis, Fig 3D, E
- [**RNA_velocity.py**](Multiome/RNA_velocity.py), [**RNA_velocity_metadata.py**](Multiome/RNA_velocity_metadata.py), Fig 3G
- [**Multiome_cells.csv**](Multiome/Multiome_cells.csv): cell barcodes that passed filtering

### Spatial_Transcriptomics

#### CosMx
- [**FXOA1_stain_score_quantification.py**](/Spatial_Transcriptomics/CosMx_analyses/FXOA1_stain_score_quantification.py): Quantify FOXA1 stain score
- [**Tumor_infiltration_score.py**](/Spatial_Transcriptomics/CosMx_analyses/Tumor_infiltration_score.py): Calculate tumor infiltration score

#### TESLA
- [**TESLA_analyses,py**](/Spatial_Transcriptomics/TESLA_analyses/TESLA_analyses.py): Tumor Edge Structure and Lymphocyte multi-level Annotation Analysis

### miscellaneous_code

- [**gene_correlation_for_all_tables_4ppt_ggplot.R**](miscellaneous_code/gene_correlation_for_all_tables_4ppt_ggplot.R): R code to draw scatter plots with prostate cancer patient data, Fig S6A, S7B  
- [**command.sh**](miscellaneous_code/command.sh): generating super enhancer plots, Fig S9C, S9E
- [**TF_volcano.R**](miscellaneous_code/TF_volcano.R): creating TF volcano plots, Fig 5A
