#Import libraries
# Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. 
# tidyverse is an R programming package assists with data import, tidying, manipulation, and data visualization. 
#ggplot2 is an R package for producing statistical, or data, graphics
library(Seurat)
library(tidyverse)
library(ggplot2)

# Reading Data:
data <- read.csv("Liver-counts.csv", row.names = 1)
meta <- read.csv("metadata_FACS.csv")
annotations <- read.csv("annotations_FACS.csv")

# Filteration of the original annotations and meta datasets to include only the entries where the tissue is "Liver",
#The filtered datasets are saved into new variables: annotations_filtered and meta_filtered.
annotations_filtered <- filter(annotations, tissue == "Liver")
meta_filtered <- filter(meta, tissue == "Liver")



#Creating a Seurat object, which is a specialized data structure that holds and
#organizes the single-cell RNA-seq data and metadata.
seurat_obj <- CreateSeuratObject(data)
#calculating the percentage of mitochondrial gene expression for each cell and stores it as a new variable (percent_mito) in the Seurat object (seurate_obj).
seurat_obj[["percent_mito"]] = PercentageFeatureSet(seurat_obj, pattern = "ˆMT-")

seurat_obj@meta.data$IDs <- rownames(seurat_obj@meta.data)

#In a Seurat object, the row names of meta.data typically represent the cell identifiers (e.g., cell barcodes or cell names).

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, )
# Visualization of the distribution of a continuous variable,
#which allows to visually inspect the distribution of these quality control metrics across all the cells in the dataset


##creating a scatter plot that shows the relationship between two features in the dataset.
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

genes <- rownames(seurat_obj)
seurate_object <- seurat_obj %>%  #applying a series of functions sequentially to the Seurat object (seurat_obj), modifying it step by step.
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 3000) %>% #common QC step to remove low-quality cells
  NormalizeData() %>% #Ensuring that the data is comparable across cells.
  # Identifying the top 2000 most variable genes in the dataset
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = genes) %>% #making the data mean-centered and variance-stabilized
  RunPCA() #Reducing the dimensionality and extracting the most important features (principal components
rm(genes)

top10 <- head(VariableFeatures(seurate_object), 10)#This selects the top 10 most variable features.
plot1 <- VariableFeaturePlot(seurate_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2

seurate_object <- RunPCA(seurate_object, features = VariableFeatures(object = seurate_object))
# Examine and visualize PCA results a few different ways
print(seurate_object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurate_object, dims = 1:2, reduction = "pca")
DimPlot(seurate_object, reduction = "pca") + NoLegend()
#DimHeatmap isDimHeatmap is used to observe how gene expression varies across the different dimensions of the data.
DimHeatmap(seurate_object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurate_object, dims = 1:20, cells = 500, balanced = TRUE)
#Elbow plot visualizes the variance explained by each principal component (PC) in a scree plot format.
ElbowPlot(seurate_object, ndims = 50)

seurate_object <- seurate_object %>%
  RunUMAP(dims = 1:24) %>% #visualize high-dimensional data in lower dimensions
  FindNeighbors(dims = 1:24) %>% #calculates the pairwise distances between cells in the high-dimensional space
  FindClusters(resolution = 0.25)
# Look at cluster IDs of the first 5 cells
head(Idents(seurate_object), 5)
Idents(seurate_object)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(seurate_object, ident.1 = 2)
head(cluster2.markers, n = 5)
cluster2.markers<- FindMarkers(seurate_object, ident.1= 1:10)
head(cluster.markers, n = 5)
cluster.markers

#4 UMAP-visualization with separation by gender

seurate_object@meta.data <- seurate_object@meta.data  %>% mutate(gender = sub(pattern = "\\w\\d+.\\w+\\d+.\\d+_\\d+_", "",x = IDs)) %>% #It replaces part of a string in the IDs column with an empty string ("").
  mutate(gender = sub(pattern = ".1.1", "", x = gender))# processes the gender column by using sub
rownames(seurate_object@meta.data) <- seurate_object@meta.data$IDs # to ensure that each row of meta.data is indexed by the unique IDs of the cells
seurate_object@meta.data
meta_filtered %>% select(c(1))#select the first column of the meta_filtered data frame.
cols = c("#F73C7D", "#1597F7", "red", "green", "orange", "brown", "blue", "darkgray", "black", "pink","yellow")
D1 <- DimPlot(seurate_object, reduction = "umap") +
  scale_color_manual(values = cols)
D2 <- DimPlot(seurate_object, reduction = "umap", label = F, group.by = "gender")
D1+D2
# number of male and female
Malescount <- sum(seurate_object$gender == "M")
Femalescount <- sum(seurate_object$gender == "F")
Malescount
Femalescount
# number of all cells
Percentage_of_Males <- (Malescount / 639) * 100 #percentage of male 
Percentage_of_Males
Percentage_of_Females <- (Femalescount / 639) * 100 #percentage of female
Percentage_of_Females

#5
library(Seurat)
library(dplyr)

# Extracting Gender and Cluster data from metadata
clustergender_data <- seurate_object@meta.data %>%
  dplyr::select(seurat_clusters, gender)

# calculating the number of cells for each gender in each cluster
gender_counts <- clustergender_data %>%
  group_by(seurat_clusters, gender) %>%
  tally() %>%
  ungroup()

# calculating the total number of cells in each cluster
total_counts <- gender_counts %>%
  group_by(seurat_clusters) %>%
  summarize(total = sum(n))

# Merging Data and calculating the precentages
gender_percentage <- gender_counts %>%
  left_join(total_counts, by = "seurat_clusters") %>%
  mutate(percentage = (n / total) * 100)

head(gender_percentage)

#import ggplot2 library
library(ggplot2)

# stacked bar plot 

ggplot(gender_percentage, aes(x = factor(seurat_clusters), y = percentage, fill = gender)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Clusters", y = "Percentage of Cells (%)", fill = "Gender") +
  scale_fill_manual(values = c("blue", "pink")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#6

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
seurate_object.markers <- FindAllMarkers(seurate_object, only.pos = TRUE)
seurate_object.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
#print the first 20 rows of the data frame
head(gender_percentage, n = 20)

#Identifying markers with FindMarkers
cluster0.markers <- FindMarkers(seurate_object, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#Plotting gene expression with VlnPlot
VlnPlot(seurate_object, features = c("0610011F06Rik", "ERCC-00002"))


#Selecting top markers for heatmap visualization
seurate_object.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Showing with Heatmap that genes-markers separate 2 chosen sub tissues.

ungroup() -> top10
DoHeatmap(seurate_object, features = top10$gene) + NoLegend()
DoHeatmap(seurate_object, features = c("0610011F06Rik", "ERCC-00002")) + NoLegend()

