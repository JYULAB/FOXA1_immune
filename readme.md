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
- [**FXOA1_stain_score_quantification.py**](/Spatial_Transcriptomics/CosMx_analyses/FXOA1_stain_score_quantification.py): Quantify FOXA1 stain score, Fig 7
- [**Tumor_infiltration_score.py**](/Spatial_Transcriptomics/CosMx_analyses/Tumor_infiltration_score.py): Calculate tumor infiltration score, Fig 7

#### TESLA
- [**TESLA_analyses,py**](/Spatial_Transcriptomics/TESLA_analyses/TESLA_analyses,py): Tumor Edge Structure and Lymphocyte multi-level Annotation Analysis, Fig 5

### scRNA-seq

- [**seq_analysis_whole_prostate.R**](/scRNA-seq/seq_analysis_whole_prostate.R): Analyses of the 10x scRNA-seq data with seurat, Fig 2, 3, 6
- [**Cellchat_analysis.R**](/scRNA-seq/Cellchat_analysis.R): Cellchat analysis of 10x whole tumor scRNA-seq data with CellChat, Fig 6
- [**pipseq_analysis_CD45.R**](/scRNA-seq/pipseq_analysis_CD45.R): Analyses of the CD45+ sorted PIP-Seq data, Fig 4 
