setwd("/Users/giuliadivirgilio/Desktop/BCG/I anno/II semestre/GENOMICS AND TRANSCRIPTOMICS/Transcriptomics/Pavesi_proj")
getwd
library(dplyr)
#install.packages("Seurat")
library(Seurat)
library(patchwork)
library(ggplot2)
load("SRA850958_SRS4386125.sparse.RData")

sm@Dimnames[[1]] <- sub("[_].*","", sm@Dimnames[[1]])
head(sm)

anyDuplicated(rownames(sm))
rownames(sm) <- make.unique(rownames(sm))


bednucleus <- CreateSeuratObject(counts = sm,
                                project = "bednucleus", min.cells = 3,
                                min.features = 200)
bednucleus # 9045 cells
head(bednucleus@meta.data)
meta_data <- bednucleus@meta.data
head(colnames(bednucleus))

################### QUALITY CONTROL ################### 
# reads that map to the mitochondrial genome. Dying cells are affected by mitochondrial contamination
grep("^mt-",rownames(bednucleus),value = TRUE)
# The [[ operator can add columns to object metadata. This is a great place to store additional info/data
bednucleus[["percent.mt"]] <- PercentageFeatureSet(bednucleus, pattern = "^mt-")
grep("^Rp[ls]",rownames(bednucleus),value = TRUE) # ribosomal protein genes
bednucleus[["percent.rbp"]] <- PercentageFeatureSet(bednucleus, pattern = "^Rp[ls]")
# first 5 rows of the table
head(bednucleus@meta.data, 5)
# Visualize QC metrics as violin plots - also adding the RPL genes
# here we refer to genes as features
VlnPlot(bednucleus, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4)
#without dots:

VlnPlot(bednucleus, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), 
        ncol = 4, 
        pt.size = 0, 
        cols = c("#78c679")) 

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
meta_data <- bednucleus@meta.data
plot1 <- FeatureScatter(bednucleus, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bednucleus, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2


bednucleus <- subset(bednucleus, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 2.0)
bednucleus

# 8870 cells

########### NORMALIZATION ########### 
# 10x data are usually just transformed into counts per 10,000 reads
# (more human readable than counts per million). 
# But, the final “expression estimate” used for downstream analyses 
# is given by the log of the normalized counts
bednucleus <- NormalizeData(bednucleus, normalization.method = "LogNormalize", scale.factor = 10000)
bednucleus@assays #27823 genes and 7954 cells

bednucleus@assays$RNA
#raw counts are here
#LayerData(pbmc, assay = "RNA", layer = "counts")
#pbmc[["RNA"]]$counts
#normalized counts are here - two possible notations
#LayerData(pbmc, assay = "RNA", layer = "data")
#pbmc[["RNA"]]$data
#normalized counts after scaling are here
#LayerData(pbmc, assay = "RNA", layer = "scale.data")
#pbmc[["RNA"]]$scale.data

apply(bednucleus[["RNA"]]$data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
# top 50 genes that have the highest mean expression across our cells
head(gene.expression, n=50)
VlnPlot(bednucleus, features = c("Malat1","MCM5"))
cc.genes

library(biomaRt)
cc.genes.mouse <- list(
  s.genes = c(
    "Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1",
    "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1",
    "Mlf1ip", "Hells", "Rfc2", "Rpa2", "Nap1l1", "Rad51ap1",
    "Gmnn", "Wdr76", "Slbp", "Cdc45", "Cdc6", "Exo1", "Tipin",
    "Dscc1", "Blm", "Casp8ap2", "Ubr7", "Pola1", "Chaf1b",
    "Brip1", "E2f8"
  ),
  g2m.genes = c(
    "Hmgb2", "Cdc20", "Cdca3", "Bub1", "Kif11", "Cdc25c",
    "Nusap1", "Ckap2l", "Ckap2", "Aurkb", "Birc5", "Tpx2",
    "Top2a", "Nek2", "G2e3", "Gas2l3", "Cenpf", "Cenpe",
    "Ccnf", "Cdc6", "Mki67", "Tmpo", "Cks2", "Nuf2", "Cks1b",
    "Ccnb2", "Top2a", "Ccnb1", "Cdca2", "Cdca3", "Cdca8",
    "Ect2", "Kif20b", "Kif2c", "Kif4", "Anp32e", "Tubb4b"
  )
)

CellCycleScoring(
  bednucleus,
  s.features = cc.genes.mouse$s.genes,
  g2m.features = cc.genes.mouse$g2m.genes,
  set.ident = TRUE
) -> bednucleus


bednucleus[[]]  # View metadata (with S.Score, G2M.Score, Phase)


# Each cell is a point in a n-dimensional space, where n is the number of genes considered. 
# The closer two points, the more similar are the transcriptomes of the corresponding cells
#the default method -vst- computes (or better, estimates) the mean-variance relationship of each gene, and chooses the 2000 genes with hte highest variance. 
bednucleus <- FindVariableFeatures(bednucleus, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(bednucleus), 10)
top10 

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(bednucleus)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(bednucleus)
# scaling data
bednucleus <- ScaleData(bednucleus, features = all.genes)

bednucleus@assays$RNA
# top 10  Sst, Klf2, Gal, Trh, Vtn, Npy, Nts, Ly6c1, C1qa, Ly6a


################ DIMENSIONALITY REDUCTUION ################
bednucleus <- RunPCA(bednucleus, features = VariableFeatures(object = bednucleus))
# Examine and visualize PCA results a few different ways
print(bednucleus[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(bednucleus, dims = 1:2, reduction = "pca")
DimPlot(bednucleus, reduction = "pca")

#with ndims we can choose how many PC to plot
ElbowPlot(bednucleus, ndims = 25)

pc.touse <- (bednucleus$pca@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse>=0.75))
pc.touse
# 13 PCs for 75% var
# 18 PCs for 80% var


################# CLUSTERING ################# 
bednucleus <- FindNeighbors(bednucleus, dims = 1:13)
bednucleus <- FindClusters(bednucleus, resolution = 0.9)



# Look at cluster IDs of the first 5 cells
head(Idents(bednucleus), 5)
head(bednucleus[[]],5)
DimPlot(bednucleus, reduction = "pca")
DimPlot(bednucleus,reduction="pca", dims=c(1,2))


bednucleus <- RunTSNE(bednucleus, dims=1:13)
DimPlot(bednucleus, reduction = "tsne")

bednucleus <- RunUMAP(bednucleus, dims = 1:13)
DimPlot(bednucleus, reduction = "umap")

# PCA has already been done before this
# Use 1:20 PCs consistently
bednucleus <- FindNeighbors(bednucleus, dims = 1:13)
bednucleus <- FindClusters(bednucleus, resolution = 0.9)

# UMAP and t-SNE use same PC space
bednucleus <- RunTSNE(bednucleus, dims = 1:18)
bednucleus <- RunUMAP(bednucleus, dims = 1:13)

# Plot results
DimPlot(bednucleus, reduction = "umap", label = FALSE)  # label clusters
DimPlot(bednucleus, reduction = "tsne", label = FALSE)
DimPlot(bednucleus, reduction = "pca", dims = c(1, 2), label = FALSE)

#######
# After clustering (e.g., using Seurat's FindClusters)
table(Idents(bednucleus))


VlnPlot(bednucleus,features="nCount_RNA")
VlnPlot(bednucleus,features="nFeature_RNA")
VlnPlot(bednucleus,features="percent.mt")
VlnPlot(bednucleus,features="percent.rbp")


library(wesanderson)
bednucleus@meta.data %>%
  group_by(seurat_clusters, Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = seurat_clusters, y = percent, fill = Phase)) +
  geom_col() +
  scale_fill_manual(values = wes_palette("Darjeeling1", n = length(unique(bednucleus$Phase)), type = "discrete")) +
  ggtitle("Percentage of Cell Cycle Phases per Cluster") +
  theme_minimal()



############ FINDING MARKER GENES and assigning cell types to clusters #########
markers <- FindAllMarkers(bednucleus, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)

markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top.markers <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
subset(markers, cluster %in% c(0, 1,9, 11)) # neurons
subset(markers, cluster %in% c(3, 5)) # oligodendrocytes
FeaturePlot(bednucleus, features = c("Slc1a3", "Gja1", "F3", "Plpp3"))
subset(markers, cluster %in% c(15, 2)) # astrocytes
#subset(markers, cluster %in% c(9, 11)) # neurons
subset(markers, cluster %in% c(12, 13)) # microglia


#markers <- FindAllMarkers(bednucleus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Check markers for cluster 0 and 3
top0 <- markers %>% filter(cluster == 0) %>% arrange(desc(avg_log2FC)) # neurons 
top1 <- markers %>% filter(cluster == 1) %>% arrange(desc(avg_log2FC)) # astrocytes
top2 <- markers %>% filter(cluster == 2) %>% arrange(desc(avg_log2FC)) # oligo
top3 <- markers %>% filter(cluster == 3) %>% arrange(desc(avg_log2FC)) # oligo
top4 <- markers %>% filter(cluster == 4) %>% arrange(desc(avg_log2FC)) # neurons?
top5 <- markers %>% filter(cluster == 5) %>% arrange(desc(avg_log2FC)) # inhibitory neurons
top6 <- markers %>% filter(cluster == 6) %>% arrange(desc(avg_log2FC)) # exitatory neurons
top7 <- markers %>% filter(cluster == 7) %>% arrange(desc(avg_log2FC)) # OPC
top8 <- markers %>% filter(cluster == 8) %>% arrange(desc(avg_log2FC)) # neurons (inhibitory)
top9 <- markers %>% filter(cluster == 9) %>% arrange(desc(avg_log2FC)) # neurons inhibitory (pretty sure)
top10 <- markers %>% filter(cluster == 10) %>% arrange(desc(avg_log2FC)) # exitatory neurons
top11 <- markers %>% filter(cluster == 11) %>% arrange(desc(avg_log2FC)) # inhibitory neurons
top12 <- markers %>% filter(cluster == 12) %>% arrange(desc(avg_log2FC)) # microglia
top13 <- markers %>% filter(cluster == 13) %>% arrange(desc(avg_log2FC)) # inhibitory neuronal
top14 <- markers %>% filter(cluster == 14) %>% arrange(desc(avg_log2FC)) # olygodendrocytes
top15 <- markers %>% filter(cluster == 15) %>% arrange(desc(avg_log2FC)) # excitatory neurons
top16 <- markers %>% filter(cluster == 16) %>% arrange(desc(avg_log2FC)) # neurons
top17 <- markers %>% filter(cluster == 17) %>% arrange(desc(avg_log2FC)) # endothelial
top18 <- markers %>% filter(cluster == 18) %>% arrange(desc(avg_log2FC)) # astrocytes
top19 <- markers %>% filter(cluster == 19) %>% arrange(desc(avg_log2FC)) # endothelial OK

intersect(top9$gene[1:20], top11$gene[1:20]) #RP23-245A10.3

top2

FeaturePlot(bednucleus, features = c(
  "Snap25",       # neurons
  "Camk2a",      # excitatory neurons
  "Gad1",         # inhibitory neurons
  "Slc1a2",       # astrocytes
  "Plp1",         # mature oligodendrocytes
  "Pdgfra",       # OPCs
  "Flt1",       # endothelial
  "Cx3cr1"      # microglia
))

FeaturePlot(bednucleus, features = c(
  "Slc7a10",
  'Bdnf',
  "Tbr1",
  "Slc17a6"# exitatory neurons
))
FeaturePlot(bednucleus, features = c(
  "Gad1",
  'Gad2',
  "Lamp5",
  "Npy"# inhibitory neurons
))


markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(bednucleus, features = top10$gene) + NoLegend()






 #######
markers <- FindAllMarkers(bednucleus18_08, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
markers %>%  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(bednucleus, features = top10$gene) + NoLegend()

VlnPlot(bednucleus, features = c( "Snap25",       # neurons
                                  "Camk2a",      # excitatory neurons
                                  "Gad1",         # inhibitory neurons
                                  "Slc1a2",       # astrocytes
                                  "Plp1",         # mature oligodendrocytes
                                  "Pdgfra",       # OPCs
                                  "Flt1",       # endothelial
                                  "Cx3cr1"  ))
VlnPlot(bednucleus, 
        features = c("Snap25",       # neurons
                     "Camk2a",      # excitatory neurons
                     "Gad1",         # inhibitory neurons
                     "Slc1a2",       # astrocytes
                     "Plp1",         # mature oligodendrocytes
                     "Pdgfra",       # OPCs
                     "Flt1",       # endothelial
                     "Cx3cr1"), 
        ncol = 4, 
        pt.size = 0) 
FeaturePlot(bednucleus, features = c("Snap25",       # neurons
                                     "Camk2a",      # excitatory neurons
                                     "Gad1",         # inhibitory neurons
                                     "Slc1a2",       # astrocytes
                                     "Plp1",         # mature oligodendrocytes
                                     "Pdgfra",       # OPCs
                                     "Flt1",       # endothelial
                                     "Cx3cr1"))
########


# DotPlot of top markers
DotPlot(bednucleus, features = unique(top10$gene)) + RotatedAxis()


table(Idents(bednucleus))
# find all markers of cluster 2 versus all the others
cluster2.markers <- FindMarkers(bednucleus, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
head(cluster2.markers, n = 10) # oligodendrocytes
cluster1.markers <- FindMarkers(bednucleus, ident.1 = 1, min.pct = 0.25, test.use = "wilcox")
head(cluster1.markers, n = 10) # astrocytes
cluster0.markers <- FindMarkers(bednucleus, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
head(cluster0.markers, n = 10)

#######################

# Current cluster identities
Idents(bednucleus)

# Reassign selected clusters
new.cluster.ids <- as.character(Idents(bednucleus))


FeaturePlot(bednucleus, features = c(
  "Snap25",       # neurons
  "Slc17a7",      # excitatory neurons
  "Gad1",         # inhibitory neurons
  "Slc1a2",       # astrocytes
  "Plp1",         # mature oligos
  "Pdgfra",       # OPCs
  "Pecam1",       # endothelial
  "Cx3cr1",       # microglia
  "Tmem119"
))
FeaturePlot(bednucleus, features = c(
  "Ecel1",
  "Cbln2",
  "Cdh7",
  "Prlr",
  "Col25a1"))

DotPlot(bednucleus, features = c("Gad2", "Slc1a2", "Mbp", "Pdgfra", "Slc17a6", 
                                 "Cd37", "Cdh5",  "Gata2"
                                ))
DotPlot(bednucleus, features = c( "Snap25",       # neurons
                                  "Camk2a",      # excitatory neurons
                                  "Gad1",         # inhibitory neurons
                                  "Slc1a2",       # astrocytes
                                  "Plp1",         # mature oligodendrocytes
                                  "Pdgfra",       # OPCs
                                  "Flt1",       # endothelial
                                  "Cx3cr1"
))


# Set new identities
bednucleus <- SetIdent(bednucleus, value = new.cluster.ids)

# Verify
DimPlot(bednucleus, reduction = "umap", label = TRUE)


new.cluster.ids <- c("Neurons_1", 
                     "Astrocytes_1", 
                     "Oligodendrocytes_1", #2
                     "Oligodendrocytes_2", 
                     "Neurons_2", 
                     "Inhibitory neurons_1", #5
                     "Excitatory neurons_1", 
                     "Oligodendrocyte Progenitor Cells", 
                     "Inhibitory neurons_2", #8
                     "Inhibitory neurons_3", 
                     "Excitatory neurons_2", 
                     "Inhibitory neurons_4", #11
                     "Microglia", #OK
                     "Inhibitory neurons_5", 
                     "Oligodendrocytes_3", #14
                     "Excitatory neurons_3", 
                     "Neurons_3", 
                     "Endothelial cells", #17
                     "Astrocytes_2", 
                     "Endothelial cells")  # OK
names(new.cluster.ids) <- levels(bednucleus)
bednucleus <- RenameIdents(bednucleus, new.cluster.ids)
DimPlot(bednucleus, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()




