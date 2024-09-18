# scRNA_scATAC
Integration of scRNAseq and scATACseq data using Weighted Nearest Neighbour (WNN) approach.

--> Loading Multi-Modal Data:

Load RNA and ATAC data using Read10X_h5. For each sample, CreateSeuratObject is used to create a Seurat object from the RNA data.
The ATAC data is added to the Seurat object using CreateChromatinAssay from the Signac package.

--> Merging Seurat Objects:

Seurat objects from multiple samples are merged into a single object using the merge function, which allows you to analyze all samples together.

--> Preprocessing RNA Data:

RNA data is normalized, variable features are identified, and PCA is run using Seurat’s standard RNA workflow.

--> Preprocessing ATAC Data:

ATAC data is preprocessed by running the TF-IDF transformation and performing dimensionality reduction using Latent Semantic Indexing (LSI).

--> WNN Integration:

The FindMultiModalNeighbors function is used to compute the Weighted Nearest Neighbors (WNN) graph, which integrates the RNA and ATAC data into a single structure.
Clustering and UMAP:

UMAP is run using the WNN graph, and clusters are identified based on the WNN graph (wsnn). The clusters are visualized using DimPlot. Differentially expressed genes and motifs are determined.
