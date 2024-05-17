# Load necessary libraries
library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(dplyr)
library(circlize)
library(networkD3)

# Load exported PVAT metadata
pvat_metadata = read.csv('../data/exported_data/pvat_metadata.csv', row.names=1)

# Load exported PVAT cell by gene expression matrix
pvat_data = read.csv('../data/exported_data/pvat_cell_by_gene_matrix.csv', row.names=1)

# Load preknowledge databases from nichenet
lr_network = readRDS("../data/nichenet_preknowledge/lr_network_mouse_21122021.rds")
ligand_target_matrix = readRDS("../data/nichenet_preknowledge/ligand_target_matrix_nsga2r_final_mouse.rds")
weighted_networks = readRDS("../data/nichenet_preknowledge/weighted_networks_nsga2r_final_mouse.rds")
ligand_tf_matrix = readRDS("../data/nichenet_preknowledge/ligand_tf_matrix_nsga2r_final_mouse.rds")
sig_network = readRDS("../data/nichenet_preknowledge/signaling_network_mouse_21122021.rds")
gr_network = readRDS("../data/nichenet_preknowledge/gr_network_mouse_21122021.rds")

# Join LR signaling weight with LR network data
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
weighted_networks_lr = weighted_networks_lr[, -c(4,5)]  # Removing unused columns

# Create a Seurat object from PVAT cell by gene matrix and metadata
seuratObj = CreateSeuratObject(counts = t(pvat_data), meta.data = pvat_metadata)

# Construct a list containing the Seurat object and both the original and transposed expression matrices
cell_data <- list(seuratObj = seuratObj,
                  expression_matrix = t(pvat_data), 
                  transposed_expression_matrix = pvat_data)

# List nicheNet analysis components using mouse preknowledge cell-cell communication data
nichenet_networks <- list(
  lr_network = lr_network,
  ligand_target_matrix = ligand_target_matrix,
  weighted_networks = weighted_networks,
  ligand_tf_matrix = ligand_tf_matrix,
  sig_network = sig_network,
  gr_network = gr_network,
  weighted_networks_lr = weighted_networks_lr
)

# Function to identify genes with high expression levels across cells
highly_expressed_genes <- function(expression_matrix, cell_ids) {
    transposed_expression_matrix = t(expression_matrix)
    genes <- transposed_expression_matrix[cell_ids, ] %>%
             apply(2, function(x) {10 * (2^x - 1)}) %>%
             apply(2, function(x) {log2(mean(x) + 1)}) %>%
             .[. >= 2] %>%
             names()
    return(genes)
}

# Function to generate a heatmap for active ligand-target interactions
ligand_target_heatmap <- function(active_ligand_target_links_df, ligand_target_matrix, best_upstream_ligands) {
    active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
    order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
    order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # Prepare gene names for heatmap visualization
    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # Prepare gene names for heatmap visualization
    vis_ligand_target = active_ligand_target_links[order_targets, order_ligands] %>% t()
    options(repr.plot.width=20, repr.plot.height=10)
    return(vis_ligand_target)
}

# Construct and order a heatmap for top ligand-receptor interactions
ligand_receptor_heatmap <- function(lr_network, best_upstream_ligands, expressed_receptors, weighted_networks_lr) {
    lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from, to)
    best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
    lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
    lr_network_top_df = lr_network_top_df_large %>% spread("from", "weight", fill = 0)
    lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
    dist_receptors = dist(lr_network_top_matrix, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
    order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
    order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
    vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
    rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
    colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
    options(repr.plot.width=20, repr.plot.height=10)
    return(vis_ligand_receptor_network)
}

# Analyze ligand-receptor activities between specified source and target cell types using nicheNet networks and cell data
ligand_receptor_activites <- function(source_celltype, target_celltype, cell_data, nichenet_networks) {
    # Assign local variables for data
    lr_network = nichenet_networks$lr_network
    ligand_target_matrix = nichenet_networks$ligand_target_matrix
    weighted_networks = nichenet_networks$weighted_networks
    ligand_tf_matrix = nichenet_networks$ligand_tf_matrix
    sig_network = nichenet_networks$sig_network
    gr_network = nichenet_networks$gr_network
    weighted_networks_lr = nichenet_networks$weighted_networks_lr
    seuratObj = cell_data$seuratObj
    expression_matrix = cell_data$expression_matrix
    transposed_expression_matrix = cell_data$transposed_expression_matrix

    n = 700  # Top number of target genes

    source_cell_ids = WhichCells(seuratObj, expression = celltype == source_celltype)
    target_cell_ids = WhichCells(seuratObj, expression = celltype == target_celltype)

    source_genes = highly_expressed_genes(expression_matrix, source_cell_ids)
    target_genes = highly_expressed_genes(expression_matrix, target_cell_ids)

    Idents(seuratObj) = seuratObj$celltype

    background_expressed_genes = target_genes %>% .[. %in% rownames(ligand_target_matrix)]
    ligands = lr_network %>% pull(from) %>% unique()
    expressed_ligands = intersect(ligands, source_genes)
    receptors = lr_network %>% pull(to) %>% unique()
    expressed_receptors = intersect(receptors, target_genes)
    potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

    ligand_activities = predict_ligand_activities(geneset = (target_genes[0:n]), background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
    ligand_activities = ligand_activities %>% arrange(-aupr_corrected)
    best_upstream_ligands = ligand_activities %>% arrange(-aupr_corrected) %>% pull(test_ligand)
    best_upstream_ligands = best_upstream_ligands %>% intersect(expressed_ligands)

    ligand_expression_tbl = tibble(
      ligand = best_upstream_ligands, 
      source = transposed_expression_matrix[source_cell_ids, best_upstream_ligands] %>% apply(2, function(x) {10 * (2 ** x - 1)}) %>% apply(2, function(x) {log2(mean(x) + 1)})
    )

    general_ligands = ligand_expression_tbl %>% pull(ligand)

    ligand_type_indication_df = tibble(
      ligand_type = c(rep(source_celltype, times = general_ligands %>% length())),
      ligand = c(general_ligands)
    )

    active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = target_genes[0:n], ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
    active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = target_celltype) %>% inner_join(ligand_type_indication_df)
    active_ligand_target_links_df = na.omit(active_ligand_target_links_df)

    p_ligand_target_network = ligand_target_heatmap(active_ligand_target_links_df, ligand_target_matrix, best_upstream_ligands)
    p_ligand_receptor_network = ligand_receptor_heatmap(lr_network, best_upstream_ligands, expressed_receptors, weighted_networks_lr)

    return(list(ligand_activities = ligand_activities, 
                p_ligand_target_network = p_ligand_target_network, 
                p_ligand_receptor_network = p_ligand_receptor_network))
}

# Clean up environment by removing large variables
rm(pvat_data, pvat_metadata, seuratObj, lr_network, ligand_target_matrix, weighted_networks, ligand_tf_matrix, sig_network, gr_network)

# Define cell types
cell_types <- c("Adipocytes", "Adipocytes_2", "ASPC", "Endothelial cells",
                "Immune cells", "Mesothelium",
                "Neuronal-like cells", "SMCs & Pericytes")

# Initialize progress bar for batch analysis
pb <- txtProgressBar(min = 0, max = length(cell_types)^2, style = 3, width = 50, char = "=")

i = 1
cell_atlas_cci = list()

# Loop through each cell type combination
for (source_cell in cell_types) {
  for (target_cell in cell_types) {
    key = paste(gsub(" ", "", source_cell), "2", gsub(" ", "", target_cell), sep="")
    cell_atlas_cci[[key]] <- tryCatch({
      ligand_receptor_activites(source_cell, target_cell, cell_data, nichenet_networks)
    }, error = function(e) {
      cat("Error with", source_cell, "to", target_cell, ":", e$message, "\n")
      NULL
    })

    i = i + 1
    setTxtProgressBar(pb, i)
  }
}
close(pb)

# Convert interactions list to data frames if applicable
names_interactions <- names(cell_atlas_cci)
for (interaction_name in names_interactions) {
  interaction <- cell_atlas_cci[[interaction_name]]

  matrix_names <- c("ligand_activities", "p_ligand_receptor_network", "p_ligand_target_network")
  for (matrix_name in matrix_names) {
    tryCatch({
      if (!is.null(interaction[[matrix_name]])) {
        interaction[[matrix_name]] <- as.data.frame(interaction[[matrix_name]])
      }
    }, error = function(e) {
      cat("Issue converting", matrix_name, "in", interaction_name, "to a data.frame:", e$message, "\n")
    })
  }
  cell_atlas_cci[[interaction_name]] <- interaction
}

# Save the results
saveRDS(cell_atlas_cci, 'pvat_cci_autocrine_included.rds')
