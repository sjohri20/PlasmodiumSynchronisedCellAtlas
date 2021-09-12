library(Seurat)

setwd("C:/Users/HP/Desktop/sem 8/project/IG/plasmodium/simulation2_re")
data.c = Read10X(data.dir = "C:/Users/HP/Desktop/sem 8/project/IG/plasmodium/Plasmo/10x/Single cell data/Control/filtered_feature_bc_matrix/")
data.t = Read10X(data.dir = "C:/Users/HP/Desktop/sem 8/project/IG/plasmodium/Plasmo/10x/Single cell data/Test/filtered_feature_bc_matrix/")

for (i in c(51:100))
{
    dat.c = data.c
    dat.t = data.t[,sample(ncol(data.t),ncol(data.c))]
    
    sobj.c <- CreateSeuratObject(counts = dat.c, project = "Plasmodium-control", min.cells = 3, min.features = 100)
    sobj.c$condition = "Control"
    
    sobj.t <- CreateSeuratObject(counts = dat.t, project = "Plasmodium-temp.stress", min.cells = 3, min.features = 100)
    sobj.t$condition = "temp.stress"
    
    sobj.c[["percent.mt"]] <- PercentageFeatureSet(object = sobj.c, pattern = "^mal")
    sobj.t[["percent.mt"]] <- PercentageFeatureSet(object = sobj.t, pattern = "^mal")
    
    sobj.c <- subset(x = sobj.c, subset = nFeature_RNA > quantile(sobj.c@meta.data$nFeature_RNA,probs=0.05) & nFeature_RNA < quantile(sobj.c@meta.data$nFeature_RNA,probs=0.95) & percent.mt < quantile(sobj.c@meta.data$percent.mt,probs=0.95))
    sobj.t <- subset(x = sobj.t, subset = nFeature_RNA > quantile(sobj.t@meta.data$nFeature_RNA,probs=0.05) & nFeature_RNA < quantile(sobj.t@meta.data$nFeature_RNA,probs=0.95)  & percent.mt < quantile(sobj.c@meta.data$percent.mt,probs=0.95) )
    
    sobj.c <- NormalizeData(object = sobj.c, normalization.method = "LogNormalize", scale.factor = 1000)
    sobj.t <- NormalizeData(object = sobj.t, normalization.method = "LogNormalize", scale.factor = 1000)
    
    sobj.c <- FindVariableFeatures(object = sobj.c, selection.method = "vst", nfeatures = 1000)
    sobj.t <- FindVariableFeatures(object = sobj.t, selection.method = "vst", nfeatures = 1000)
    
    all.genes <- rownames(x = sobj.c)
    sobj.c <- ScaleData(object = sobj.c, features = all.genes)
    all.genes <- rownames(x = sobj.t)
    sobj.t <- ScaleData(object = sobj.t, features = all.genes)
    
    sobj.c<- RunPCA(object = sobj.c, features = VariableFeatures(object = sobj.c),verbose = F,nfeatures.print = F,ndims.print = F )
    sobj.t <- RunPCA(object = sobj.t, features = VariableFeatures(object = sobj.t),verbose = F,nfeatures.print = F,ndims.print = F)
    
    sobj.c$condition = "control"
    sobj.t$condition = "treatment"
    obj_vector_1=c()
    anchors <- FindIntegrationAnchors(object.list = list(sobj.c, sobj.t), dims = 1:10)
    combined <- IntegrateData(anchorset = anchors, dims = 1:20)
    
    DefaultAssay(object = combined) <- "integrated"
    
    # Run the standard workflow for visualization and clustering
    combined <- ScaleData(object = combined, verbose = FALSE)
    combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
    # t-SNE and Clustering
    
    combined <- RunTSNE(object = combined, reduction = "pca", dims = 1:20)
    combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:20)
    combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:10)
    combined <- FindClusters(combined, resolution = 0.5)
    
    saveRDS(combined@meta.data,paste(paste("run_",i, sep = ""), ".rds", sep = ""))

}