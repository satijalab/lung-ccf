library(Seurat)


d175.data <- read.csv("/home/stuartt/github/lung-ccf/data/D175R_dist.csv", row.names = 1)
d175.lm <- read.csv("/home/stuartt/github/lung-ccf/data/D175R_lm.csv", stringsAsFactors = FALSE)
d175.data <- t(d175.data)
rownames(d175.data) <- d175.lm$label
d205.data <- read.csv("/home/stuartt/github/lung-ccf/data/D205R_dist.csv")
d205.lm <- read.csv("/home/stuartt/github/lung-ccf/data/D205R_lm.csv", stringsAsFactors = FALSE)
d205.lm$label <- gsub(pattern = "Cardiac_notch_1", replacement = "Lower_cardiac_notch_1", x = d205.lm$label)
rownames(d205.data) <- d205.lm$label

obj1 <- CreateSeuratObject(counts = d205.data, project = "lung1", assay = "CT")
obj2 <- CreateSeuratObject(counts = d175.data, project = "lung2", assay = "CT")

anchors <- FindIntegrationAnchors(
  object.list = list(obj1, obj2),
  anchor.features = rownames(obj1),
  reduction = 'cca',
  dims = 1:5
)

integrated <- IntegrateData(
  anchorset = anchors, dims = 1:5
)

integrated <- ScaleData(integrated, features = rownames(integrated))
integrated <- RunPCA(integrated, features = rownames(integrated), approx = FALSE)
integrated <- RunUMAP(integrated, dims = 1:5)

library(RANN)

# obj1 = 205 <- Neighbors
# obj2 = 175 <- Query
voxels1 <- colnames(obj1)
voxels2 <- colnames(obj2)

d1_embeddings <- Embeddings(object = integrated[['pca']])[voxels1, ]
d2_embeddings <- Embeddings(object = integrated[['pca']])[voxels2, ]

# collect points
k <- 20
n <- 100

neighbors <- nn2(data = d1_embeddings, query = d2_embeddings, k = k)
query_voxels <- sample(1:length(voxels2), n, replace = FALSE)
n1 <- matrix(voxels1[neighbors$nn.idx[query_voxels, ]], ncol = k)
# fix formatting
n1 <- gsub(pattern = "X.", replacement = "", x = n1)
n1 <- gsub(pattern = "..", replacement = ",", x = n1, fixed = TRUE)
n1 <- gsub(pattern = '^\\.|\\.$', replacement = '', x= n1)
qvoxel <- voxels2[query_voxels]
qvoxel <- gsub(pattern = "[()]| |", replacement = "", x = qvoxel)

rownames(n1) <- qvoxel

# save
write.table(
  x = n1,
  file = "/home/stuartt/github/lung-ccf/data/neighbor_points_D205R.tsv",
  quote = FALSE,
  sep = "\t"
)