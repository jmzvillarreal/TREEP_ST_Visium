# TREEP_ST_Visium
# Hands-on session for Visium Spatial Transcriptomics

1️⃣ Install conda environment, activate it and check what is installed:
```bash
conda env create -f ST_omics-env.yml
conda activate ST_omics-env
conda list
```
2️⃣ Download scRNAseq reference data:
```bash
wget https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds
```
3️⃣ Open R console:
```bash
R
```
4️⃣ Seurat analysis (inspired in: https://satijalab.org/seurat/articles/spatial_vignette.html)


```r

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)
library(SingleCellExperiment)

### Load 10X Visium data:

InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

## Data preprocessing

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

## Data normalization

# run normalization to store sctransform residuals for all genes
brain <- SCTransform(brain, assay = "Spatial", return.only.var.genes = FALSE, verbose = TRUE)
# also run standard log normalization for comparison
brain <- NormalizeData(brain, verbose = FALSE, assay = "Spatial")

## Dimensionality reduction, clustering, and visualization 

brain <- RunPCA(brain, assay = "SCT", verbose = TRUE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = TRUE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

# Highlight clusters

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3, 5, 8)), facet.highlight = TRUE, ncol = 3)
    
## Identification of Spatially Variable Features

# Differentially expressed genes

de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:4], alpha = c(0.1, 1), ncol = 4)

# Moran's Index spatial correlation

brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
    selection.method = "moransi")

# Visualize top 6

top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))

## Subset out anatomical regions 

# Subset to clusters of interest
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))

# Extract and attach spatial coordinates from anterior1 image
centroids <- cortex[["anterior1"]]@boundaries$centroids
coords <- setNames(as.data.frame(centroids@coords), c("x", "y"))
rownames(coords) <- centroids@cells
cortex$x <- coords[colnames(cortex), "x"]
cortex$y <- coords[colnames(cortex), "y"]

# Perform spatial subsetting using image position Subset to upper right quadrant
cortex <- subset(cortex, y < 2719 | x > 7835, invert = TRUE)

# Further remove cells in lower right quadrant
m <- (9395 - 5747)/(4960 - 7715)
b <- 5747 - m * 7715
cortex <- subset(cortex, y > m * x + b, invert = TRUE)


# Visualize the subsetted data on full image vs cropped image
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2



### Integration with single-cell data

allen_reference <- readRDS("allen_cortex.rds")

# Subset to 25 cells per subclass BEFORE processing

cells_per_type <- 25

cell_ids <- unlist(lapply(
  split(colnames(allen_reference), allen_reference$subclass),
  function(x) {
    if(length(x) >= cells_per_type){
      sample(x, cells_per_type)
    } else {
      NULL
    }
  }
))

allen_reference <- subset(allen_reference, cells = cell_ids)

# Preprocess reference 
allen_reference <- SCTransform(allen_reference, verbose = TRUE)
allen_reference <- RunPCA(allen_reference, verbose = TRUE)
allen_reference <- RunUMAP(allen_reference, dims = 1:30)

DimPlot(allen_reference, group.by = "subclass", label = TRUE, reduction = "umap")


# Build RCTD reference

counts <- allen_reference@assays$RNA$counts

cell_types <- allen_reference$subclass
cell_types <- gsub("/", "_", cell_types)
cell_types <- factor(cell_types)

names(cell_types) <- colnames(allen_reference)

nUMI <- allen_reference$nCount_RNA
names(nUMI) <- colnames(allen_reference)

reference <- Reference(counts, cell_types, nUMI)


# Prepare spatial data

coords <- GetTissueCoordinates(brain)

counts <- as.matrix(brain@assays$Spatial$counts)

stopifnot(all(rownames(coords) == colnames(counts)))

spatial <- SpatialRNA(
  coords = coords[, c("x", "y")],
  counts = counts,
  require_int = FALSE
)


# Run RCTD

rctd <- create.RCTD(
  spatialRNA = spatial,
  reference = reference,
  max_cores = 4
)

rctd <- run.RCTD(rctd, doublet_mode = "full")

# Add results to Seurat object


brain_RCTD <- subset(brain, cells = rownames(rctd@results$weights))

results <- rctd@results$weights
allen_subclass <- as.data.frame(results)

brain_RCTD <- AddMetaData(brain_RCTD, metadata = allen_subclass)

allen_subclass_dominant <- colnames(results)[apply(results, 1, which.max)]

dominant_df <- data.frame(allen_subclass_dominant)
rownames(dominant_df) <- rownames(results)

brain_RCTD <- AddMetaData(brain_RCTD, metadata = dominant_df)

SpatialDimPlot(
  brain_RCTD,
  group.by = "allen_subclass_dominant",
  image.alpha = 0,
  pt.size.factor = 2.5
)

```
