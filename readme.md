# Code for FOXA1 Paper
Citation:

The analyses were conducted in a High-Performance Computing (HPC) environment using a combination of Bash scripting, R (v4.0.3), and Python (v3.6 or later). Detailed software versions are described in the Methods section of the paper.

The scripts are organized by technology or analysis type. For workflows requiring sequential execution, scripts are prefixed with numbers indicating the order. Otherwise, most scripts are standalone. Each script's expected output corresponds to specific figures in the paper.
 
### ChIP-seq

- [**heatmap.sh**](ChIP-seq/heatmap.sh): Creating heatmaps, Fig 3
- [**mapping.sh**](ChIP-seq/mapping.sh): Map raw reads
- [**peakcalling.sh**](ChIP-seq/peakcalling.sh): Call peaks
- [**spike_in_downsampling.sh**](ChIP-seq/spike_in_downsampling.sh): Down sample with Spike-in
- [**venn.R**](ChIP-seq/venn.R): Draw venn diagram, Fig 3
    
### Spatial_Transcriptomics

#### CosMx
- [**FXOA1_stain_score_quantification.py**](Spatial_Transcriptomics/CosMx_analyses/FXOA1_stain_score_quantification.py): Quantify FOXA1 stain score, Fig 7
- [**Tumor_infiltration_score.py**](Spatial_Transcriptomics/CosMx_analyses/Tumor_infiltration_score.py): Calculate tumor infiltration score, Fig 7

#### TESLA
- [**TESLA_analyses,py**](/Spatial_Transcriptomics/TESLA_analyses/TESLA_analyses,py): Tumor Edge Structure and Lymphocyte multi-level Annotation Analysis, Fig 5

### scRNA-seq

- [**scRNA-seq_processing.R**](scRNA-seq/scRNA-seq_processing.R), [**scRNA-seq_integration.R**](scRNA-seq/scRNA-seq_integration.R), processing, integrating and annotation of scRNA-Seq samples, Fig 2, 3, 6
- [**10x_scRNAseq_analysis_whole_prostate.R**](scRNA-seq/10x_scRNAseq_analysis_whole_prostate.R): Analyses of the 10x scRNA-seq data with seurat, Fig 2, 3, 6
- [**trajectory_analysis.R**](scRNA-seq/trajectory_analysis.R): Trajectory analysis of epithelial populations Fig 2
- [**CellChat_analysis.R**](scRNA-seq/CellChat_analysis.R): Analysis of 10x whole tumor scRNA-seq data with CellChat, Fig 6
- [**Pip-seq_processing_processing.R**](scRNA-seq/Pip-seq_processing_processing.R): Processing, integrating and annotation of PIP-Seq data, Fig 4
- [**pipseq_analysis_CD45.R**](scRNA-seq/pipseq_analysis_CD45.R): Analyses of the CD45+ sorted PIP-Seq data, Fig 4 
