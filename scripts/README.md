# Analysis Scripts and Notebooks

This folder contains a collection of R scripts and Python notebooks used for analyzing data from the PVAT cell atlas. Each file is dedicated to a specific type of analysis ranging from data export to visual representation through heatmaps and dotplots.

## File Descriptions

### Python Notebooks

- **adata_to_expr_matrix.ipynb**
  - **Description:** This notebook handles the export of original data from the PVAT cell atlas, specifically excluding BAT (Brown Adipose Tissue) data.

- **pvat_col_to_itg_lliana.ipynb**
  - **Description:** Contains analysis using CellChat and CellPhoneDB to explore cell-cell communication.

- **dotplots.ipynb**
  - **Description:** Generates heatmap visualizations and dotplots of the interaction data derived from the previous analyses.

### R Scripts

- **pvat_collagen_integrin_nichenet.r**
  - **Description:** Performs NicheNet analysis to explore collagen-integrin interactions in the PVAT, focusing on cell-cell communication inferences (CCI).
