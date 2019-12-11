#install Seurat
install.packages("Seurat")
library(Seurat)

# Read in SPLiT-Seq datasets and metadata
data_split <- readRDS("~/data/BUKMAP_20190131_SPL-R_UMI_counts_dTn6_EmptyBC_Filter.rds")
meta_split <- read.table("~/data/BUKMAP_20190131_SPLiT-Seq_Meta_Table.txt",sep = "\t",header = T,row.names = 1)

#parse the SPLiT-Seq metadata
split_library <- sapply(Cells(data_split),ExtractField,1)
meta_split_parsed <- meta_split[split_library,]
rownames(meta_split_parsed) <- colnames(data_split)

# Read in 10X datasets. There is no metadata for the 10X datasets
data_10x_A <- readRDS("~/data/KPMP_20190529_K1900174_10X_UMI_counts_dTn6_EmptyBC_Filter.rds")
data_10x_B <- readRDS("~/data/KPMP_20190607L_K1900174_10X_UMI_counts_dTn6_EmptyBC_Filter.rds")

# create a list of Seurat objects, one for each datasets
objects <- list()
objects[[1]] <- CreateSeuratObject(counts = data_split, meta.data = meta_split_parsed)
objects[[2]] <- CreateSeuratObject(counts = data_10x_A)
objects[[3]] <- CreateSeuratObject(counts = data_10x_B,min.features = 600)

# For the first object, we will take only the cortex cells (the KPMP datasets are cortex)
objects[[1]] <- subset(objects[[1]], Region == "Cortex")

# Since the SPLiT-Seq dataset has more cells, we will treat it as a 'reference' for this example
# We process it through the Seurat SCTransform workflow (https://satijalab.org/seurat/v3.0/sctransform_vignette.html)
# We use 50PCs for all datasets, and learn technical noise models based on 3,000 cells. 
# These are standard parameters for large datasets, and have not been tuned here

objects[[1]] <- SCTransform(object = objects[[1]], ncells = 3000) %>% RunPCA() %>% RunUMAP(dims = 1:50) %>% FindNeighbors(dims = 1:50) %>% FindClusters(resolution = 1)

# We also normalize the other datasets as well and for visualization purposes, run dimensional reduction as well
objects[[2]] <- SCTransform(object = objects[[2]], ncells = 3000) %>% RunPCA() %>% RunUMAP(dims = 1:50)
objects[[3]] <- SCTransform(object = objects[[3]], ncells = 3000) %>% RunPCA() %>% RunUMAP(dims = 1:50)

# as a demonstration, save the objects to load back in later
#saveRDS(objects, file = "~/data/hubmap_kidney_objects.rds")

#readRDS("~/data/hubmap_kidney_objects.rds")
objects[[1]]@meta.data$tech <- 'SPLiT-Seq'
objects[[2]]@meta.data$tech <- '10X_Genomics'
objects[[3]]@meta.data$tech <- '10X_Genomics'

#Workflow  1 - integrate the datasets into a common space
# Now integrate the datasets using reciprocal PCA integration (
# Vignette at : https://satijalab.org/seurat/v3.1/integration.html

features <- SelectIntegrationFeatures(objects[c(1,3)])
objects <- PrepSCTIntegration(object.list = objects,anchor.features = features)
for(i in 1:3) {
  objects[[i]] <- RunPCA(objects[[i]], features = features)
}

int_anchors <- FindIntegrationAnchors(object.list = objects[c(1,3)],reduction = 'cca',dims = 1:30,anchor.features = features,normalization.method = 'SCT',nn.method = 'annoy')
i1 <- IntegrateData(anchorset = int_anchors,normalization.method = 'SCT')
i1 <- RunPCA(i1) %>% RunUMAP(dims=1:30)
DimPlot(i1,split.by = 'tech')

# the i1 dataset is now integrated

# Workflow 2 - transfer annotations 
trans_anchors <- FindTransferAnchors(reference = objects[[1]],query = objects[[3]],dims = 1:30,features = features,normalization.method = 'SCT',nn.method = 'annoy')
predicted_labels <- TransferData(anchorset = trans_anchors,refdata = Idents(objects[[1]]),weight.reduction = objects[[3]][["pca"]])
objects[[3]] <- AddMetaData(objects[[3]],predicted_labels)

trans_anchors_pca <- FindTransferAnchors(reference = objects[[3]],query = objects[[1]],dims = 1:30,features = features,normalization.method = 'SCT',nn.method = 'annoy')
predicted_labels <- TransferData(anchorset = trans_anchors_pca,refdata = Idents(objects[[3]]),weight.reduction = objects[[1]][["pca"]])
objects[[1]] <- AddMetaData(objects[[1]],predicted_labels)


trans_anchors <- FindTransferAnchors(reference = objects[[1]],query = objects[[3]],reduction = 'cca',dims = 1:30,features = features,normalization.method = 'SCT',nn.method = 'annoy')
predicted_labels <- TransferData(anchorset = trans_anchors,refdata = Idents(objects[[1]]),weight.reduction = "cca")
objects[[3]] <- AddMetaData(objects[[3]],predicted_labels)

#now, the second dataset has been annotated
p1 <- DimPlot(objects[[1]],label = T)+NoLegend()
p2 <- DimPlot(objects[[3]],label = F,group.by = 'predicted.id')+NoLegend()
plot_grid(p1,p2)

