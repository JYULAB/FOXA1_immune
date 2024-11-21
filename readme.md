# Code for FOXA1 Paper
Citation:

This repository includes the main pieces of code used to support the conclusions in the paper. The following code were used in a SLURM High Performance Computing environment, using a combination of bash scripting, R (4.0.3) and python >=3.6. Detailed versions of softwares used are listed in Methods section of the paper.

The scripts are organized by technology/analysis type. If any scripts are to be run in order, prefix numbers are provided in the script name. Otherwise, each script is mostly standalone. Some files are provided here, such as the final cells passing filtering that were used in Multiome analysis. Expected output of each script are the indicated figures from the paper.
 
### ChIP-seq

- [**command.sh**](ChIP-seq/command.sh): creating heatmaps, Fig 4D, 4E, 4H, 5E, 7D, S5E, S7E, S9A, S9B (same tool is used for ATAC-seq heatmaps).

### scRNA-seq

- [**clinical_patients_intergration.R**](Multiome/clinical_patients_intergration.R): Integration with clinical patients and pseudotime analysis, Fig 3D, E


### Spatial_Transcriptomics

#### CosMx
- [**FXOA1_stain_score_quantification.py**](/Spatial_Transcriptomics/CosMx_analyses/FXOA1_stain_score_quantification.py): Quantify FOXA1 stain score, Fig 7
- [**Tumor_infiltration_score.py**](/Spatial_Transcriptomics/CosMx_analyses/Tumor_infiltration_score.py): Calculate tumor infiltration score, Fig 7

#### TESLA
- [**TESLA_analyses,py**](/Spatial_Transcriptomics/TESLA_analyses/TESLA_analyses,py): Tumor Edge Structure and Lymphocyte multi-level Annotation Analysis, Fig 5

### miscellaneous_code

- [**TF_volcano.R**](miscellaneous_code/TF_volcano.R): creating TF volcano plots, Fig 5A
