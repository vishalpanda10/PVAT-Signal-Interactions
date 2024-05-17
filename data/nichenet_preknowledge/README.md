# NicheNet Preknowledge Databases

This folder contain details about the NicheNet preknowledge databases used in our analyses. These files contain essential data for performing ligand-target interaction predictions and other related computational tasks in the study of cell-cell communication.

## Files

1. **ligand_target_matrix_nsga2r_final_mouse.rds**
   - Description: Contains the ligand-target interaction matrix for mouse, crucial for predicting potential ligand-receptor interactions.

2. **signaling_network_mouse_21122021.rds**
   - Description: Represents the signaling network data for mouse, detailing the pathways and interactions relevant to signal transduction.

3. **weighted_networks_nsga2r_final_mouse.rds**
   - Description: Includes weighted network data for mouse, which is used to weigh the significance of different signaling pathways based on the context.

4. **gr_network_mouse_21122021.rds**
   - Description: Contains gene regulation networks for mouse, which are essential for understanding how ligands influence gene expression.

5. **lr_network_mouse_21122021.rds**
   - Description: Involves ligand-receptor pairs and their interactions for mouse, foundational for the identification of communication links between cells.

6. **ligand_tf_matrix_nsga2r_final_mouse.rds**
   - Description: Features a matrix of ligand-transcription factor interactions for mouse, providing insights into the regulatory mechanisms affected by ligands.

## Usage

These files are utilized in various scripts and analyses throughout the project to facilitate a deeper understanding of cellular communication mechanisms, supporting both hypothesis generation and testing. Ensure that the paths to these files are correctly configured in your scripts as per the repository structure.
