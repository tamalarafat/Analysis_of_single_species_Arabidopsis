# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Load data tables
###

###
# Protoplasting induced genes
###

# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Col0.leaf.vs.protoplast_table_2pseudorep_final_August_2022.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID


###
# Orthologues table
###

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")


###
# WT A. thaliana
###

# Load data - WT OX 1st Experiment (leaf 5 and 6)
COL_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_1ST_2/filtered_feature_bc_matrix/")

# Load data - WT OX 2nd Experiment (leaf 6 and 7)
COL_data_2E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col_RNA_2nd_ALL_2/filtered_feature_bc_matrix/")

# Load data - WT OX 3rd Experiment (leaf 5 and 6)
COL_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_3rd_ALL/filtered_feature_bc_matrix/")

# Load data - WT OX 7th Experiment (leaf 6 and 7)
COL_data_5E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_5th_ALL_2/filtered_feature_bc_matrix/")

# All gene IDs - Arabidopsis Thaliana
thaliana_genes = rownames(COL_data_1E)

# extracting the Cardamine IDs that are present in orthologues table 
thaliana_ortho_genes = as.character(ortho_table$A.thaliana.TAIR10)

# not all the thaliana ids are present in the ortho data - 51 thaliana genes are missing in the thaliana data
thaliana_ortho_genes = intersect(thaliana_genes, thaliana_ortho_genes)

# Let's subset the data with the ortho genes
COL_data_1E <- COL_data_1E[thaliana_ortho_genes, ]
COL_data_2E <- COL_data_2E[thaliana_ortho_genes, ]
COL_data_3E <- COL_data_3E[thaliana_ortho_genes, ]
COL_data_5E <- COL_data_5E[thaliana_ortho_genes, ]

# All gene IDs - all datasets
ortho_genes = rownames(COL_data_1E)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(ortho_genes, PP_genes)

##### Remove the protoplasting induced genes
COL_data_1E <- COL_data_1E[genes_to_keep, ]
COL_data_2E <- COL_data_2E[genes_to_keep, ]
COL_data_3E <- COL_data_3E[genes_to_keep, ]
COL_data_5E <- COL_data_5E[genes_to_keep, ]


###
# COL - 1 E
###

# First replicate - COL 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
COL_1E <- CreateSeuratObject(counts = COL_data_1E, project = "COL_1E", min.cells = 10, min.features = 200)

# Add metadata information to the seurat object
COL_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-1", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_1E <- subset(COL_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_1E[["percent.mt"]] <- PercentageFeatureSet(COL_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_1E[["percent.pt"]] <- PercentageFeatureSet(COL_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_1E <- subset(COL_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_1E <- NormalizeData(COL_1E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
COL_1E <- FindVariableFeatures(COL_1E, selection.method = "vst", nfeatures = 2000)


###
# COL - 2 E
###

# First replicate - COL 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
COL_2E <- CreateSeuratObject(counts = COL_data_2E, project = "COL_2E", min.cells = 27, min.features = 200)

# Add metadata information to the seurat object
COL_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-2", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_2E <- subset(COL_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_2E[["percent.mt"]] <- PercentageFeatureSet(COL_2E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_2E[["percent.pt"]] <- PercentageFeatureSet(COL_2E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_2E <- subset(COL_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_2E <- NormalizeData(COL_2E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
COL_2E <- FindVariableFeatures(COL_2E, selection.method = "vst", nfeatures = 2000)


###
# COL - 3 E
###

# First replicate - COL 3E - total cells 4100; filter out genes that are not detected in at least 8 cells
COL_3E <- CreateSeuratObject(counts = COL_data_3E, project = "COL_3E", min.cells = 6, min.features = 200)

# Add metadata information to the seurat object
COL_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-3", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_3E <- subset(COL_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_3E[["percent.mt"]] <- PercentageFeatureSet(COL_3E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_3E[["percent.pt"]] <- PercentageFeatureSet(COL_3E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_3E <- subset(COL_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_3E <- NormalizeData(COL_3E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
COL_3E <- FindVariableFeatures(COL_3E, selection.method = "vst", nfeatures = 2000)


###
# COL - 5 E
###

# First replicate - COL 5E - total cells 8420; filter out genes that are not detected in at least 18 cells
COL_5E <- CreateSeuratObject(counts = COL_data_5E, project = "COL_5E", min.cells = 17, min.features = 200)

# Add metadata information to the seurat object
COL_5E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-5", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_5E <- subset(COL_5E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_5E[["percent.mt"]] <- PercentageFeatureSet(COL_5E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_5E[["percent.pt"]] <- PercentageFeatureSet(COL_5E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_5E <- subset(COL_5E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_5E <- NormalizeData(COL_5E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
COL_5E <- FindVariableFeatures(COL_5E, selection.method = "vst", nfeatures = 2000)

# Integration of the replicates - find features from the datasets to anchor cells from different sources
anchFeatures <- SelectIntegrationFeatures(object.list = list(COL_1E, COL_2E, COL_3E, COL_5E))

fileGenerator(anchFeatures, "Seurat_HVG_standard.txt")

ingAnchors <- FindIntegrationAnchors(object.list = list(COL_1E, COL_2E, COL_3E, COL_5E), dims = 1:50, anchor.features = anchFeatures)

# To keep records of all the genes in the integrated assay, create a feature set with all of the genes from different replicates.
features_integrated <- unique(rownames(COL_1E), rownames(COL_2E), rownames(COL_3E), rownames(COL_5E))

# Integrate the replicates
integrated.data <- IntegrateData(anchorset = ingAnchors, dims = 1:50, verbose = T, features.to.integrate = features_integrated)

# Setting the default assay to "integrated"
DefaultAssay(integrated.data) <- "integrated"

# Gene level scaling - standardization
integrated.data <- ScaleData(integrated.data, verbose = FALSE)

# Run PCA
integrated.data <- RunPCA(integrated.data, npcs = 50, verbose = FALSE)

# Run UMAP and tSNE
integrated.data <- RunUMAP(integrated.data, reduction = "pca", dims = 1:50, n.components = 2)

integrated.data <- RunTSNE(integrated.data, reduction = "pca", dims = 1:50, dim.embed = 2)

# Find neighbours and clusters
integrated.data <- FindNeighbors(integrated.data, reduction = "pca", dims = 1:50)

for (i in seq(0.1, 1.2, 0.1)) {
  integrated.data <- FindClusters(integrated.data, resolution = i, n.start = 50, n.iter = 50)
}

save(integrated.data, file = "integrated_wt_arabidopsis_seurat.RData")

writeLines(capture.output(sessionInfo()), "Session_info_integrated_wt_arabidopsis_seurat.txt")
