#Integration of scRNAseq and scATACseq using WNN.
# Load necessary libraries
library(Seurat)
library(Signac)
library(ggplot2)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Define parameters for the analysis
min_cells <- 10      # Minimum number of cells expressing a gene or region
min_features <- 200   # Minimum number of features per cell

# Load the multi-modal data for multiple samples (TEA-seq: RNA, ATAC, Protein)
# Example assumes three samples with both RNA and ATAC data.
sample1_rna <- Read10X_h5("sample1/filtered_feature_bc_matrix.h5")
sample1_atac <- Read10X_h5("sample1/atac_matrix.h5")
sample2_rna <- Read10X_h5("sample2/filtered_feature_bc_matrix.h5")
sample2_atac <- Read10X_h5("sample2/atac_matrix.h5")
sample3_rna <- Read10X_h5("sample3/filtered_feature_bc_matrix.h5")
sample3_atac <- Read10X_h5("sample3/atac_matrix.h5")

# Create Seurat objects for each sample
seurat_obj1 <- CreateSeuratObject(counts = sample1_rna, assay = "RNA")
seurat_obj1[["percent.mt"]] <- PercentageFeatureSet(seurat_obj1, pattern = "^MT-") 
seurat_obj2 <- CreateSeuratObject(counts = sample2_rna, assay = "RNA")
seurat_obj2[["percent.mt"]] <- PercentageFeatureSet(seurat_obj2, pattern = "^MT-")
seurat_obj3 <- CreateSeuratObject(counts = sample3_rna, assay = "RNA")
seurat_obj3[["percent.mt"]] <- PercentageFeatureSet(seurat_obj3, pattern = "^MT-")

# Plot the QC metrics
pdf("Vlnplots.pdf", width = 15, height = 8)
VlnPlot(seurat_obj1, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3) 
VlnPlot(seurat_obj2, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3) 
VlnPlot(seurat_obj3, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3) 
dev.off()

#Filter low quality cells
seurat_obj1 <- subset(seurat_obj1, subset = nFeature_RNA < 700 & nCount_RNA < 1500 & percent.mt < 25)
seurat_obj2 <- subset(seurat_obj2, subset = nFeature_RNA < 2000 & nCount_RNA < 4000 & percent.mt < 15)
seurat_obj3 <- subset(seurat_obj3, subset = nFeature_RNA < 1800 & nCount_RNA < 4000 & percent.mt < 25)


# Add ATAC data to each Seurat object
seurat_obj1[["ATAC"]] <- CreateChromatinAssay(counts = sample1_atac, min.cells = min_cells)
seurat_obj2[["ATAC"]] <- CreateChromatinAssay(counts = sample2_atac, min.cells = min_cells)
seurat_obj3[["ATAC"]] <- CreateChromatinAssay(counts = sample3_atac, min.cells = min_cells)

# Merge the datasets into one Seurat object
seurat_obj <- merge(seurat_obj1, y = list(seurat_obj2, seurat_obj3), add.cell.ids = c("Sample1", "Sample2", "Sample3"))

# Preprocessing RNA Data (Normalization, Feature Selection, PCA)
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 50)

# Preprocessing ATAC Data (TF-IDF, LSI)
DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = "q0")
seurat_obj <- RunTFIDF(seurat_obj)
seurat_obj <- RunSVD(seurat_obj, reduction.key = "LSI_", n = 50)

# Perform Weighted Nearest Neighbor (WNN) Analysis
# This step combines RNA and ATAC (and optionally protein) modalities.
seurat_obj <- FindMultiModalNeighbors(
  seurat_obj, 
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 1:30)
)

# Run UMAP based on WNN graph
seurat_obj <- RunUMAP(seurat_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap")

# Find clusters using the WNN graph
seurat_obj <- FindClusters(seurat_obj, graph.name = "wsnn", algorithm = 1, resolution = 0.8)

# Annotate the clusters using public data
reference <- LoadH5Seurat("../pbmc_multimodal.h5seurat")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE,raster=FALSE) + NoLegend()

anchors <- FindTransferAnchors(
  reference = reference,
  query = seurat_obj,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
seurat_obj <- MapQuery(
  anchorset = anchors,
  query = seurat_obj,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca"
)

query.nn <- Seurat:::NNHelper(query = seurat_obj[['ref.spca']]@cell.embeddings , 
                              data =  reference[['spca']]@cell.embeddings , 
                              method = "annoy",
                              k = 20, 
                              metric = "cosine"
)
seurat_obj[['ref.umap']]  <- RunUMAP(object = query.nn,  reduction.model = reference[[ "wnn.umap" ]] )


# Visualization of UMAP with WNN
pdf("WNN_UMAP_plot.pdf", width = 13, height = 10)
DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "seurat_clusters", pt.size = 4,label = FALSE, label.size = 3, repel = TRUE) + ggtitle("WNN Clusters")
DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "orig.ident", pt.size = 4,label = FALSE, label.size = 3, repel = TRUE) + ggtitle("Samples")
DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "TEMRA_level",pt.size = 4, label = FALSE, label.size = 3, repel = TRUE) + ggtitle("TEMRA_level")
DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "Condition",pt.size = 4, label = FALSE, label.size = 3, repel = TRUE) + ggtitle("Condition")
DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "celltype.l2", pt.size = 4,label = FALSE, label.size = 3, repel = TRUE) + ggtitle("Annotation") 
dev.off()

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(seurat_obj) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(seurat_obj))) %in% main.chroms)
seurat_obj[["ATAC"]] <- subset(seurat_obj[["ATAC"]], features = rownames(seurat_obj[["ATAC"]])[keep.peaks])
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(seurat_obj), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
seurat_obj <- SetAssayData(seurat_obj, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

seurat_obj <- RunChromVAR(
  object = seurat_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Determine differentially expressed Markers and Motifs
markers_rna <- presto:::wilcoxauc.Seurat(X = seurat_obj, group_by = 'seurat_clusters', assay = 'data', seurat_assay = 'RNA')
saveRDS(markers_rna,"Differentially_expressed_markers_WNN.rds")
markers_motifs <- presto:::wilcoxauc.Seurat(X = seurat_obj, group_by = 'seurat_clusters', assay = 'data', seurat_assay = 'chromvar')
saveRDS(markers_motifs,"Differentially_expressed_motifs_WNN.rds")

# Save Seurat object 
saveRDS(seurat_obj, file = "seurat_tea_wnn_analysis.rds")

