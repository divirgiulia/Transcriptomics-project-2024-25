setwd("/Users/giuliadivirgilio/Desktop/BCG/I anno/II semestre/GENOMICS AND TRANSCRIPTOMICS/Transcriptomics/Pavesi_proj")
getwd
library(dplyr)
#install.packages("Seurat")
library(Seurat)
library(patchwork)
library(ggplot2)
load("SRA701877_SRS3279690.sparse.RData")

sm@Dimnames[[1]] <- sub("[_].*","", sm@Dimnames[[1]])
head(sm)

anyDuplicated(rownames(sm))
rownames(sm) <- make.unique(rownames(sm))

PanIslets <- CreateSeuratObject(counts = sm,
                              project = "PanIslets", min.cells = 3,
                              min.features = 200)
PanIslets
head(PanIslets@meta.data)
meta_data <- PanIslets@meta.data
head(colnames(PanIslets))

################### QUALITY CONTROL ################### 
# reads that map to the mitochondrial genome. Dying cells are affected by mitochondrial contamination
grep("^MT-",rownames(PanIslets),value = TRUE)
# The [[ operator can add columns to object metadata. This is a great place to store additional info/data
PanIslets[["percent.mt"]] <- PercentageFeatureSet(PanIslets, pattern = "^MT-")
grep("^RP[LS]",rownames(PanIslets),value = TRUE) # ribosomal protein genes
PanIslets[["percent.rbp"]] <- PercentageFeatureSet(PanIslets, pattern = "^RP[LS]")
# first 5 rows of the table
head(PanIslets@meta.data, 5)
# Visualize QC metrics as violin plots - also adding the RPL genes
# here we refer to genes as features
VlnPlot(PanIslets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4)
#without dots:

VlnPlot(PanIslets, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), 
        ncol = 4, 
        pt.size = 0, 
        cols = c("#78c679")) 

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
meta_data <- PanIslets@meta.data
#plot1 <- FeatureScatter(PanIslets, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(PanIslets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# with ggplot
cor1 <- cor(meta_data$nCount_RNA, meta_data$percent.mt, method = "pearson")
cor2 <- cor(meta_data$nCount_RNA, meta_data$nFeature_RNA, method = "pearson")

# Format the correlation text
cor_label1 <- paste0("Pearson's r = ", round(cor1, 2))
cor_label2 <- paste0("Pearson's r = ", round(cor2, 2))
plot1 <- ggplot(meta_data, aes(x = nCount_RNA, y = percent.mt)) +
  geom_point(color = "#addd8e", alpha = 0.6, size = 1) +
  theme_minimal() +
  labs(title = "nCount_RNA vs percent.mt", subtitle = cor_label1,
       x = "nCount_RNA", y = "percent.mt") +
  theme(
    plot.title = element_text(face = "bold"),           # Bold title
  )

plot2 <- ggplot(meta_data, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(color = "#78c679", alpha = 0.6, size = 1) +
  theme_minimal() +
  labs(title = "nCount_RNA vs nFeature_RNA", subtitle = cor_label2,
       x = "nCount_RNA", y = "nFeature_RNA") +
  theme(
    plot.title = element_text(face = "bold"),           # Bold title
  )

cor3 <- cor(meta_data$nCount_RNA, meta_data$percent.rbp, method = "pearson")
cor_label3 <- paste0("Pearson's r = ", round(cor3, 2))
plot3 <-  ggplot(meta_data, aes(x = nCount_RNA, y = percent.rbp)) +
  geom_point(color = "#31a354", alpha = 0.6, size = 1) +
  theme_minimal() +
  labs(title = "nCount_RNA vs percent.rbp", subtitle = cor_label3,
       x = "nCount_RNA", y = "percent.rbp") +
  theme(
    plot.title = element_text(face = "bold"),           # Bold title
  )
plot1 + plot2 + plot3

PanIslets <- subset(PanIslets, subset = nFeature_RNA > 100 & nFeature_RNA < 5600 & percent.mt < 20)
PanIslets



########### NORMALIZATION ########### 
# 10x data are usually just transformed into counts per 10,000 reads
# (more human readable than counts per million). 
# But, the final “expression estimate” used for downstream analyses 
# is given by the log of the normalized counts
PanIslets <- NormalizeData(PanIslets, normalization.method = "LogNormalize", scale.factor = 10000)
PanIslets@assays #23480 genes and 2492 cells

PanIslets@assays$RNA
#raw counts are here
#LayerData(pbmc, assay = "RNA", layer = "counts")
#pbmc[["RNA"]]$counts
#normalized counts are here - two possible notations
#LayerData(pbmc, assay = "RNA", layer = "data")
#pbmc[["RNA"]]$data
#normalized counts after scaling are here
#LayerData(pbmc, assay = "RNA", layer = "scale.data")
#pbmc[["RNA"]]$scale.data

apply(PanIslets[["RNA"]]$data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
# top 50 genes that have the highest mean expression across our cells
head(gene.expression, n=50)
VlnPlot(PanIslets, features = c("RPS2","GAPDH"))
cc.genes.updated.2019

CellCycleScoring(PanIslets, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> PanIslets     
PanIslets[[]]

# Each cell is a point in a n-dimensional space, where n is the number of genes considered. 
# The closer two points, the more similar are the transcriptomes of the corresponding cells
#the default method -vst- computes (or better, estimates) the mean-variance relationship of each gene, and chooses the 2000 genes with hte highest variance. 
PanIslets <- FindVariableFeatures(PanIslets, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(PanIslets), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(PanIslets)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Add title and display
MostVarGenes <- plot1 + plot2
MostVarGenes + ggtitle("Top 10 Variable Features in PanIslets")

all.genes <- rownames(PanIslets)
# scaling data
PanIslets <- ScaleData(PanIslets, features = all.genes)

PanIslets@assays$RNA
# top 10 PPY, PRSS1, CPA2, CTRB1, REG3A, PNLIP, REG1B, SPINK1, PLA2G1B, CCL2 


################ DIMENSIONALITY REDUCTUION ################
PanIslets <- RunPCA(PanIslets, features = VariableFeatures(object = PanIslets))
# Examine and visualize PCA results a few different ways
print(PanIslets[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(PanIslets, dims = 1:2, reduction = "pca")
DimPlot(PanIslets, reduction = "pca")

#with ndims we can choose how many PC to plot
ElbowPlot(PanIslets, ndims = 25)

pc.touse <- (PanIslets$pca@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse>=0.70))
pc.touse




################# CLUSTERING ################# 
PanIslets <- FindNeighbors(PanIslets, dims = 1:20)
PanIslets <- FindClusters(PanIslets, resolution = 1.0)

# Look at cluster IDs of the first 5 cells
head(Idents(PanIslets), 5)
head(PanIslets[[]],5)
DimPlot(PanIslets, reduction = "pca")
DimPlot(PanIslets,reduction="pca", dims=c(1,2))


PanIslets <- RunTSNE(PanIslets, dims=1:15)
DimPlot(PanIslets, reduction = "tsne")

PanIslets <- RunUMAP(PanIslets, dims = 1:10)
DimPlot(PanIslets, reduction = "umap")

# PCA has already been done before this
# Use 1:20 PCs consistently
PanIslets <- FindNeighbors(PanIslets, dims = 1:20)
PanIslets <- FindClusters(PanIslets, resolution = 0.8)

# UMAP and t-SNE use same PC space
PanIslets <- RunTSNE(PanIslets, dims = 1:15)
PanIslets <- RunUMAP(PanIslets, dims = 1:20)

# Plot results
DimPlot(PanIslets, reduction = "umap", label = FALSE)  # label clusters
DimPlot(PanIslets, reduction = "tsne", label = FALSE)
DimPlot(PanIslets, reduction = "pca", dims = c(1, 2), label = FALSE)

#######

VlnPlot(PanIslets,features="nCount_RNA")
VlnPlot(PanIslets,features="nFeature_RNA")
VlnPlot(PanIslets,features="percent.mt")
VlnPlot(PanIslets,features="percent.rbp")


library(wesanderson)
PanIslets@meta.data %>%
  group_by(seurat_clusters, Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = seurat_clusters, y = percent, fill = Phase)) +
  geom_col() +
  scale_fill_manual(values = wes_palette("Darjeeling1", n = length(unique(PanIslets$Phase)), type = "discrete")) +
  ggtitle("Percentage of Cell Cycle Phases per Cluster") +
  theme_minimal()



############ FINDING MARKER GENES and assigning cell types to clusters #########
markers <- FindAllMarkers(PanIslets, 
                          only.pos = TRUE,     # only return upregulated genes
                          min.pct = 0.25,      # gene must be expressed in at least 25% of cells
                          logfc.threshold = 0.25)  # log fold change threshold
markers %>% filter(cluster == 0) %>% top_n(10, avg_log2FC) # 
markers %>% filter(cluster == 1) %>% top_n(10, avg_log2FC) # beta cells
markers %>% filter(cluster == 2) %>% top_n(10, avg_log2FC) # alpha
markers %>% filter(cluster == 3) %>% top_n(10, avg_log2FC) # ductal cells
markers %>% filter(cluster == 4) %>% top_n(10, avg_log2FC)
markers %>% filter(cluster == 5) %>% top_n(10, avg_log2FC) # acinar cells CELA3B
markers %>% filter(cluster == 6) %>% top_n(10, avg_log2FC) # delta SST, PPY(gammaPP)
markers %>% filter(cluster == 7) %>% top_n(10, avg_log2FC) # pancreatic stellar cells VIM
markers %>% filter(cluster == 8) %>% top_n(10, avg_log2FC) # 
markers %>% filter(cluster == 9) %>% top_n(10, avg_log2FC) 
markers %>% filter(cluster == 10) %>% top_n(10, avg_log2FC)
markers %>% filter(cluster == 11) %>% top_n(10, avg_log2FC) # endothelial cells CDH5

# DotPlot of top markers
DotPlot(PanIslets, features = unique(top_markers$gene)) + RotatedAxis()
# Violin plot of a single gene
VlnPlot(PanIslets, features = "INS")
# UMAP with marker overlay
FeaturePlot(PanIslets, features = c("INS", "GCG", "SST", "KRT19"))
FeaturePlot(PanIslets, features = c("INS", "GCG", "SST", "PDGFRB", "KRT19"))
FeaturePlot(PanIslets, features = c("INS", "GCG", "SST", "PPY", "KRT19", "PECAM1", "COL1A1"))

table(Idents(PanIslets))


# find all markers of cluster 2 versus all the others
cluster2.markers <- FindMarkers(PanIslets, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
head(cluster2.markers, n = 5)

cluster2_01.markers <- FindMarkers(PanIslets, ident.1 = 2, ident.2 = c(0, 1), min.pct = 0.25)
head(cluster2_01.markers, n = 5)

#we return only genes "over expressed", found in at least 25% of the cells, and with a logFC threshold of at least 0.25
PanIslets.markers <- FindAllMarkers(PanIslets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

PanIslets.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
PanIslets.markers %>%
  group_by(cluster) %>%
  slice_min(n = 2, order_by = p_val_adj )


PanIslets.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(PanIslets, features = top10$gene) + NoLegend()


cluster0.markers <- FindMarkers(PanIslets, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster0.markers <- cluster0.markers[order(-cluster0.markers$avg_log2FC),]
head(cluster0.markers, n = 10)

cluster1.markers <- FindMarkers(PanIslets, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster1.markers <- cluster1.markers[order(-cluster1.markers$avg_log2FC),]
head(cluster1.markers, n = 10)

cluster2.markers <- FindMarkers(PanIslets, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster2.markers <- cluster2.markers[order(-cluster2.markers$avg_log2FC),]
head(cluster2.markers, n = 10)

cluster3.markers <- FindMarkers(PanIslets, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster3.markers <- cluster3.markers[order(-cluster3.markers$avg_log2FC),]
head(cluster3.markers, n = 10)


cluster4.markers <- FindMarkers(PanIslets, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster4.markers <- cluster4.markers[order(-cluster4.markers$avg_log2FC),]
head(cluster3.markers, n = 10)

cluster5.markers <- FindMarkers(PanIslets, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster5.markers <- cluster5.markers[order(-cluster5.markers$avg_log2FC),]
head(cluster5.markers, n = 10)

cluster6.markers <- FindMarkers(PanIslets, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster6.markers <- cluster6.markers[order(-cluster6.markers$avg_log2FC),]
head(cluster6.markers, n = 10)

cluster7.markers <- FindMarkers(PanIslets, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster7.markers <- cluster7.markers[order(-cluster7.markers$avg_log2FC),]
head(cluster7.markers, n = 10)






