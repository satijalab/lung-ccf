library(Seurat)
library(data.table)
mask <- as.matrix(read.table("~/Downloads/m1_mask.csv",sep = ','))
inds <- which((mask)>0)
inds_arr <- which((mask)>0,arr.ind = T)

data = list()
for(i in 1:14) {
  file=paste0("~/Downloads/m1_landmark",i,".csv")
  data[[i]] <- as.matrix(read.table(file,sep = ","))
  data[[i]] <- as.vector(data[[i]])[inds]
  data[[i]] <- exp(data[[i]]/255)
  print(i)
}
df <- data.frame(matrix(unlist(data), nrow=14, byrow=T),stringsAsFactors=FALSE)
rownames(df) <- paste("L",1:nrow(df),sep="")
colnames(df) <- paste0("cell1_",1:ncol(df))
srt <- CreateSeuratObject(counts = df)
drembed <- as.matrix(inds_arr)
rownames(drembed) <- Cells(srt)
colnames(drembed) <- c("img_1","img_2")
srt[["img"]] <- CreateDimReducObject(embeddings = drembed,key = "img")
srt <- ScaleData(srt)
srt <- RunPCA(srt,features = rownames(srt))
ElbowPlot(srt)

srt <- FindNeighbors(srt,dims = 1:7,k.param = 30)
srt <- FindClusters(srt,n.start = 1,resolution = 0.1)
#srt <- RunUMAP(srt,dims = 1:7)
FeaturePlot(srt,"PC_1",reduction = 'img')



####
mask <- as.matrix(read.table("~/Downloads/m2_mask.csv",sep = ','))
inds <- which((mask)>0)
inds_arr <- which((mask)>0,arr.ind = T)

data = list()
for(i in 1:14) {
  file=paste0("~/Downloads/m2_landmark",i,".csv")
  data[[i]] <- as.matrix(read.table(file,sep = ","))
  data[[i]] <- as.vector(data[[i]])[inds]
  data[[i]] <- exp(data[[i]]/255)
  print(i)
}
df <- data.frame(matrix(unlist(data), nrow=14, byrow=T),stringsAsFactors=FALSE)
rownames(df) <- paste("L",1:nrow(df),sep="")
colnames(df) <- paste0("cell2_",1:ncol(df))
srt2 <- CreateSeuratObject(counts = df)
drembed <- as.matrix(inds_arr)
rownames(drembed) <- Cells(srt2)
colnames(drembed) <- c("img_1","img_2")
srt2[["img"]] <- CreateDimReducObject(embeddings = drembed,key = "img")
srt2 <- ScaleData(srt2)
srt2 <- RunPCA(srt2,features = rownames(srt2))

srt2 <- FindNeighbors(srt2,dims = 1:7,k.param = 30)
srt2 <- FindClusters(srt2,n.start = 1,resolution = 0.1)
#srt2 <- RunUMAP(srt2,dims = 1:7)
FeaturePlot(srt2,"PC_1",reduction = 'img')

#optional subsampling
d1 <- subset(srt,cells = sample(Cells(srt),5000))
d2 <- subset(srt2,cells = sample(Cells(srt2),5000))

#optonal: transfer cluster  label from d1 onto d2
ta <- FindTransferAnchors(reference = d1,query = d2,features = rownames(d1),npcs = NULL,dims = 1:7)
tb <- TransferData(anchorset = ta,refdata = Idents(d1),dims = 1:7)
d2 <- AddMetaData(d2,metadata = tb)
DimPlot(d2,reduction = 'img',group.by = 'predicted.id')

d1$orig.ident <- 'cell1'
d2$orig.ident <- 'cell2'
ti <- FindIntegrationAnchors(object.list = list(d1,d2),reduction = 'rpca',dims = 1:7)

d3 <- IntegrateData(ti,dims = 1:7)
d3@meta.data[Cells(d1),"orig.ident"] <- 'cells2'
d3 <- ScaleData(d3,features = rownames(d3))
d3 <- RunPCA(d3,features = rownames(d3))

library(RANN)
library(cowplot)

ref_pca <- Embeddings(d3[["pca"]])[Cells(d1),1:7]
q_pca <- Embeddings(d3[["pca"]])[Cells(d2),1:7]
nn_refquery <- nn2(query = ref_pca,data = q_pca,k = 10)
cell <- 1
for(cell in seq(1,5000,50)) {
  query_nns <- nn_refquery$nn.idx[cell,]
  p1 <- DimPlot(d1,reduction = 'img',cells.highlight = Cells(d1)[cell],sizes.highlight = 2)+ggtitle("Slice 1")
  p2 <- DimPlot(d3,reduction = 'umap',cells.highlight = list(Cells(d1)[cell],Cells(d2)[query_nns]),cols.highlight = c("blue","red"),sizes.highlight = 2)+ggtitle("Shared Space")
  p3 <- DimPlot(d2,reduction = 'img',cells.highlight = Cells(d2)[query_nns],cols.highlight = c("blue"),sizes.highlight = 2)+ggtitle("Slice 3")
  p4=plot_grid(p1,p2,p3,ncol = 3)
  ggsave(p4,filename = paste0("~/Downloads/",cell,".png"),height = 5,width = 20)
}
