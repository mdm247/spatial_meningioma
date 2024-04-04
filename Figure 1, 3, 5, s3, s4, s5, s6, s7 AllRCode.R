# Script for Figure 1 (c, d), Figure 3 (b, c), Figure 5 (c, d), Supplementary Figure 3 (c, d), Supplementary Figure 4 (a, b, f, g, h), Supplementary Figure 5b, Supplementary Figure 6a, Supplementary Figure 7 (d, e, i)
#########################################################################################################################################################################
### Individual tumor spatial assay and feature plot
Spatial_Meningioma_number <- Load10X_Spatial(
  data.dir="./Spatial_Meningioma_01_WithProbeSet/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial")

feature_plot <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "gene", pt.size.factor = 1.6, alpha = c(0.5, 1), image.alpha = 0, min.cutoff = 0, max.cutoff='q90')
#########################################################################################################################################################################




#########################################################################################################################################################################
### Normalization, batch correction, and clustering
library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
library(harmony)
library(future)
library(clustree)

counts = Read10X(data.dir="./Spatial_Meningioma_01_WithProbeSet/filtered_feature_bc_matrix/")
a = CreateSeuratObject(counts=counts, project="A", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_02_WithProbeSet/filtered_feature_bc_matrix/")
b = CreateSeuratObject(counts=counts, project="B", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_03_WithProbeSet/filtered_feature_bc_matrix/")
c = CreateSeuratObject(counts=counts, project="C", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_04_WithProbeSet/filtered_feature_bc_matrix/")
d = CreateSeuratObject(counts=counts, project="D", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_05_WithProbeSet/filtered_feature_bc_matrix/")
e = CreateSeuratObject(counts=counts, project="E", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_06_WithProbeSet/filtered_feature_bc_matrix/")
f = CreateSeuratObject(counts=counts, project="F", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_07_WithProbeSet/filtered_feature_bc_matrix/")
g = CreateSeuratObject(counts=counts, project="G", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_08_WithProbeSet/filtered_feature_bc_matrix/")
h = CreateSeuratObject(counts=counts, project="H", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_09_WithProbeSet/filtered_feature_bc_matrix/")
i = CreateSeuratObject(counts=counts, project="I", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_10_WithProbeSet/filtered_feature_bc_matrix/")
j = CreateSeuratObject(counts=counts, project="J", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_11_WithProbeSet/filtered_feature_bc_matrix/")
k = CreateSeuratObject(counts=counts, project="K", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_12_WithProbeSet/filtered_feature_bc_matrix/")
l = CreateSeuratObject(counts=counts, project="L", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_13_WithProbeSet/filtered_feature_bc_matrix/")
m = CreateSeuratObject(counts=counts, project="M", min.cells=1, min.features=10)

counts = Read10X(data.dir="./Spatial_Meningioma_14_WithProbeSet/filtered_feature_bc_matrix/")
n = CreateSeuratObject(counts=counts, project="N", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_15_WithProbeSet/filtered_feature_bc_matrix/")
o = CreateSeuratObject(counts=counts, project="O", min.cells=3, min.features=100)

counts = Read10X(data.dir="./Spatial_Meningioma_16_WithProbeSet/filtered_feature_bc_matrix/")
p = CreateSeuratObject(counts=counts, project="P", min.cells=3, min.features=100)

Spatial_Meningioma_aggr <- merge(x=a, y = c(b, c, d, e, f, g, h, i, j, k, l, m, n, o, p), add.cell.ids=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"))

VlnPlot(Spatial_Meningioma_aggr, features="nFeature_RNA", pt.size=0, group.by="orig.ident")
Spatial_Meningioma_aggr = subset(Spatial_Meningioma_aggr, subset = nFeature_RNA>1500 & nFeature_RNA<9500)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = Spatial_Meningioma_aggr), value = TRUE, ignore.case = TRUE)
Spatial_Meningioma_aggr[["percent.mt"]] <- PercentageFeatureSet(Spatial_Meningioma_aggr, pattern = "^MT-")
VlnPlot(Spatial_Meningioma_aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Spatial_Meningioma_aggr = subset(Spatial_Meningioma_aggr, subset = percent.mt<10)
#Spatial_Meningioma_aggr <- NormalizeData(Spatial_Meningioma_aggr, normalization.method = "LogNormalize", scale.factor = 10000)
VlnPlot(Spatial_Meningioma_aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plan("multisession", workers=16)
options(future.globals.maxSize = 12000 * 1024^2)
options(error = function() {traceback(2, max.lines=100); if(!interactive()) quit(save="no", status=1, runLast=T)})

Spatial_Meningioma_aggr = SCTransform(Spatial_Meningioma_aggr, vars.to.regress="nCount_RNA")
Spatial_Meningioma_aggr = RunPCA(Spatial_Meningioma_aggr, npcs=50)
Spatial_Meningioma_aggr = RunHarmony(Spatial_Meningioma_aggr, group.by.vars="orig.ident", assay.use="SCT")
Spatial_Meningioma_aggr = RunUMAP(Spatial_Meningioma_aggr, dims=1:50, reduction="harmony", min.dist=0.01)
DimPlot(Spatial_Meningioma_aggr, reduction="umap", group.by="orig.ident")
Spatial_Meningioma_aggr = FindNeighbors(Spatial_Meningioma_aggr, dims=1:30, reduction="harmony")
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=0.7)
DimPlot(Spatial_Meningioma_aggr, reduction="umap", label=T)

#Cluster analysis can be performed using clustree - run through iterations on FindClusters and append those to the Spatial_Meningioma_aggr Seuratobject
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=5)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=2)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=1.2)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=1.0)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=0.8)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=0.7)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=0.6)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=0.5)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=0.4)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=0.2)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=0.1)
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=0.0)

#view clustree
clustree(Spatial_Meningioma_aggr, prefix = "SCT_snn_res.")
clustree(Spatial_Meningioma_aggr, prefix = "SCT_snn_res.", node_colour = "sc3_stability")
clustree(Spatial_Meningioma_aggr, prefix = "SCT_snn_res.", layout = "sugiyama")

Spatial_Meningioma_aggr[['SCT_snn_res.5']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.2']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.1.2']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.1']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.0.8']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.0.7']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.0.6']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.0.5']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.0.4']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.0.2']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.0.1']] <- NULL
Spatial_Meningioma_aggr[['SCT_snn_res.0']] <- NULL

#run with optimal clustering resolution
Spatial_Meningioma_aggr = FindClusters(Spatial_Meningioma_aggr, resolution=n)
DimPlot(Spatial_Meningioma_aggr, reduction="umap", label=T)
#########################################################################################################################################################################



#########################################################################################################################################################################
# cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Spatial_Meningioma_aggr <- CellCycleScoring(Spatial_Meningioma_aggr, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(Spatial_Meningioma_aggr[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(Spatial_Meningioma_aggr, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
Spatial_Meningioma_aggr <- RunPCA(Spatial_Meningioma_aggr, features = c(s.genes, g2m.genes))
DimPlot(Spatial_Meningioma_aggr, order = c('G2M', 'S', 'G1'), cols = c('black', 'red', 'blue'), pt.size = 1)
#########################################################################################################################################################################



#########################################################################################################################################################################
###DGE
Spatial_Meningioma_aggr <- SetIdent(Spatial_Meningioma_aggr, value = "seurat_clusters")

# run PrepSCTFindMarkers to ensure fixed value is set properly
Spatial_Meningioma_aggr <- PrepSCTFindMarkers(Spatial_Meningioma_aggr, verbose = F)

# Getting DEGs for each cluster vs everything else
markers_clust <- FindAllMarkers(Spatial_Meningioma_aggr, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
markers_clust %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(x = markers_clust, file = "Spatial_Meningioma_aggr_clusters_DEGs.csv")
#########################################################################################################################################################################



#########################################################################################################################################################################
###feature plot
library("SCpubr")

featuresCustom = c("gene")
scpubr_plot <- do_FeaturePlot(
  Spatial_Meningioma_aggr,
  featuresCustom,
  viridis_color_map = "A",
  split.by = 'orig.ident'
)
ggsave(filename = "Spatial_Meningioma_aggr_feature_plot_custom_byOrig.png", plot = scpubr_plot, width = 80, height = 50, limitsize = FALSE)
#########################################################################################################################################################################



#########################################################################################################################################################################
### inferCNV
library(infercnv)
library(SpatialInferCNV)
library(tidyverse)

#Setting source dataset
read_file("./filtered_feature_bc_matrix.h5")

# Reading the 10X Visium count data, and appending the string "spatial_menin_10X" before each barcode with ImportCountData()
spatial_menin_ENSBMLID_Counts <- ImportCountData("spatial_menin_10X", "./filtered_feature_bc_matrix.h5")

# Showing the first part of the spatial_menin_ENSBMLID_Counts dataframe 
head(spatial_menin_ENSBMLID_Counts)

#Importing Histological Annotation Data
# Reading the histological annotation file created above, and appending the string "spatial_menin_10X" before each barcode with ImportHistologicalAnnotations()
spatial_menin_10x_Histology <- ImportHistologicalAnnotations("spatial_menin_10X", "./spatial_menin_10x_HistologyAnnotations.csv")

# Showing the first part of the spatial_menin_10x_Histology dataframe 
head(spatial_menin_10x_Histology)

#Merging histological annotation data with count data
# Merging annotations, and counts, and applying QC with MergingCountAndAnnotationData()
spatial_menin_10x_Joined_Counts <- MergingCountAndAnnotationData("spatial_menin_10X",spatial_menin_10x_Histology, spatial_menin_ENSBMLID_Counts)

# Setting the "Genes" column as the rownames
spatial_menin_10x_Joined_Counts <- spatial_menin_10x_Joined_Counts %>% column_to_rownames(var = "Genes")

# Removing the pre-merged count file as it is no longer needed
rm(spatial_menin_ENSBMLID_Counts)

# Showing the first part of the spatial_menin_10x_Histology dataframe 
head(spatial_menin_10x_Joined_Counts)

#Selecting finalized annotations for export for use in infercnv::run
dim(spatial_menin_10x_Histology)

#Running FinalAnnotations() to return only the barcodes in the count file
FinalAnnotationsForExport <- FinalAnnotations(spatial_menin_10x_Histology, spatial_menin_10x_Joined_Counts)

#Diplaying the dimensions of the final annotation:
dim(FinalAnnotationsForExport)

# Write the spatial_menin_10x_Joined_Counts.tsv to the local directory
write.table(spatial_menin_10x_Joined_Counts, "spatial_menin_10x_Joined_Counts.tsv", sep = "\t")

# Write the spatial_menin_10xFinalAnnotationsForExport.tsv to the local directory
write.table(FinalAnnotationsForExport, "spatial_menin_10xFinalAnnotationsForExport.tsv", 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

# Create the inferCNV object from the 3 files, with/without the reference group and mitochondria
spatial_menin_10xCancer_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix="./spatial_menin_10x_Joined_Counts.tsv", 
                                                                 gene_order_file="./siCNV_GeneOrderFile.tsv",
                                                                 annotations_file="./spatial_menin_10xFinalAnnotationsForExport.tsv",
                                                                 delim="\t",
                                                                 ref_group_names=NULL, 
                                                                 chr_exclude = c("chrM")) 
# Run inferCNV
spatial_menin_10xCancer_infCNV = infercnv::run(spatial_menin_10xCancer_infCNV, 
                                               cutoff=0.1, #(see infercnv::run documentation)
                                               out_dir="./InferCNVrun_outputs", 
                                               cluster_by_groups=FALSE, #unsupervised analysis
                                               HMM = FALSE, 
                                               denoise=TRUE)
#########################################################################################################################################################################




#########################################################################################################################################################################
### deconvolution

library(SeuratData)
library(SeuratDisk)
library(SCDC)
library(scater)
library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
library(harmony)
library(future)
library(clustree)
library(SpatialInferCNV)
library(tidyverse)
library(infercnv)
library(NGCHM)
library(infercnvNGCHM)
library(ggpubr)

msc_NG <- readRDS('single_cell_10x_meningioma.rds')
msc_NG[[]]
msc_NG$original <- msc_NG$seurat_clusters

list1 <- msc_NG@meta.data[["original"]]
list1
list2 <- gsub("10", "Active monocytes", list1)
list3 <- gsub("11", "Glia", list2)
list4 <- gsub("12", "SSTR2A meningioma cells", list3)
list5 <- gsub("13", "Monocytes", list4)
list6 <- gsub("14", "CD4 T cells", list5)
list7 <- gsub("0", "G1 phase meningioma cells", list6)
list8 <- gsub("1", "ECM remodeling meningioma cells", list7)
list9 <- gsub("2", "Cycling S phase meningioma cells", list8)
list10 <- gsub("3", "CD163 monocytes", list9)
list11 <- gsub("4", "Cycling G2M phase meningioma cells", list10)
list12 <- gsub("5", "CD8 T cells", list11)
list13 <- gsub("6", "Pericytes", list12)
list14 <- gsub("7", "Vascular cells", list13)
list15 <- gsub("8", "Cycling G2M phase meningioma cells", list14)
list16 <- gsub("9", "Meningioma cells", list15)
list16

# check number of cells per subclass
table(msc_NG$customclassif)

# select 200 cells per subclass
Idents(msc_NG) <- msc_NG$customclassif
msc_NG <- subset(msc_NG, cells = WhichCells(msc_NG, downsample = 200))

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(msc_NG, label = TRUE)

msc_NG@active.assay = "RNA"

markers_sc <- FindAllMarkers(msc_NG, only.pos = TRUE, max.cells.per.ident = 1000, assay = "RNA")

# Filter for genes that are also present in the ST data
markers_sc <- markers_sc[markers_sc$gene %in% rownames(Spatial_Meningioma_number), ]

# Select top n genes per cluster, select top by first p-value, then absolute diff in pct, then quota of pct.
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +
                                                               1))
top50 <- markers_sc %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
m_feats <- unique(as.character(top50$gene))

eset_SC <- ExpressionSet(assayData = as.matrix(msc_NG@assays$RNA@counts[m_feats,
]), phenoData = AnnotatedDataFrame(msc_NG@meta.data))
eset_ST <- ExpressionSet(assayData = as.matrix(Spatial_Meningioma_number@assays$Spatial@counts[m_feats,
]), phenoData = AnnotatedDataFrame(Spatial_Meningioma_number@meta.data))

deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "customclassif",
                                     ct.sub = as.character(unique(eset_SC$customclassif)))
table_deconv <- deconvolution_crc$prop.est.mvw
Spatial_Meningioma_number@assays[["SCDC"]] <- CreateAssayObject(data = t(deconvolution_crc$prop.est.mvw))

if (length(Spatial_Meningioma_number@assays$SCDC@key) == 0) {
  Spatial_Meningioma_number@assays$SCDC@key = "scdc_"
}

DefaultAssay(Spatial_Meningioma_number) <- "SCDC"
p1 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "G1 phase meningioma cells", pt.size.factor = 1.6,
                         crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p2 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "Cycling G2M phase meningioma cells", pt.size.factor = 1.6,
                         crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p3 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "Cycling G2M phase meningioma cells", pt.size.factor = 1.6,
                         crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p4 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "Cycling S phase meningioma cells", pt.size.factor = 1.6,
                         crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p5 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "Glia", pt.size.factor = 1.6,
                         crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p6 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "Pericytes", pt.size.factor = 1.6,
                         crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p7 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "Vascular cells", pt.size.factor = 1.6,
                         crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p8 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "Active monocytes", pt.size.factor = 1.6,
                         crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p9 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "Monocytes", pt.size.factor = 1.6,
                         crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p10 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "CD163 monocytes", pt.size.factor = 1.6,
                          crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p11 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "CD4 T cells", pt.size.factor = 1.6,
                          crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p12 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "CD8 T cells", pt.size.factor = 1.6,
                          crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p13 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "ECM remodeling meningioma cells", pt.size.factor = 1.6,
                          crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p14 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "SSTR2A meningioma cells", pt.size.factor = 1.6,
                          crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
p15 <- SpatialFeaturePlot(Spatial_Meningioma_number, features = "Meningioma cells", pt.size.factor = 1.6,
                          crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
pcombined = (p1 | p2 | p3| p4) / (p5 | p6 | p7 | p11) / (p8 | p9 | p10 | p12) / (p13 | p14 | p15)
#########################################################################################################################################################################