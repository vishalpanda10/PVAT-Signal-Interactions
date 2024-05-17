# Cell-Cell Communication Analysis in PVAT

## Introduction
This repository contains analysis of cell-cell communication (CCC) in perivascular adipose tissue (PVAT), focusing on the role of integrins in establishing tone after stress relaxation. Our research utilizes a comprehensive CCC analysis using several analytical tools to investigate the interactions between various cell types within PVAT.

## Methods Summary
We employed the Ligand-receptor ANAnalysis Framework (LIANA) to conduct CellPhoneDB and CellChat analyses, alongside an independent run of NicheNet. These tools helped us identify key ligand-receptor interactions based on the expression levels of the ligands and receptors, backed by extensive knowledgebases:
- **CellPhoneDB** and **CellChat** infer CCC based on ligand and receptor expressions, using permutation-based statistical methods to highlight significant interactions.
- **NicheNet** offers a more detailed approach by considering the downstream responses resulting from known ligand-receptor interactions, using the top 700 highly expressed genes to infer signaling activity measured by AUROC and AUPRC scores.

## Key Findings
1. **Integrin's Role in PVAT:** Integrins significantly contribute to the mechanical tone of PVAT following stress relaxation, particularly through interactions involving collagen.
2. **Collagen-Integrin Interactions:** Both CellPhoneDB and CellChat confirm interactions predominantly involving collagen and integrin, with significant contributions across several cell types, including adipocytes, ASPCs, and SMCs/Pericytes.
3. **Consensus and Unique Findings:**
   - NicheNet identified unique integrin interactions (Itgav) and corroborated findings with CellPhoneDB for Itga5 specific to ASPCs.
   - CellPhoneDB alone identified potential activation of integrin Itga9 and Itga11 in ASPCs.
   - A consensus among the tools highlighted ASPCs as crucial in mediating collagen-integrin interactions, supported by the significant number of interactions identified.

## Visualization of Findings
The repository includes detailed visualizations of the CCC analysis:
- **Heatmaps and Dot Plots:** Illustrate the intensity and significance of the collagen-integrin interactions across different cell types.
- **Directional Interaction Diagrams:** Provide insights into the directional nature of these interactions, highlighting the role of collagen as ligands emanating from various cell sources impacting integrin receptors.

## Supplementary Materials
Detailed tables and supplementary figures (Supplemental tables 1-3 and Figures 9A, 9B, and 9C) provide exhaustive lists of significant interactions and a comprehensive overview of the inferred CCC among PVAT cell types.

## Conclusion
Our analysis underscores the complex nature of cell-cell communication in PVAT, with a specific focus on the mechanistic roles of collagen and integrin. This detailed examination aids in understanding the cellular interactions that contribute to the structural and functional integrity of PVAT.

