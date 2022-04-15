plot_pca_tsne<-function(npc, resolution, object, method){
  # method should be a string
  object<-as(object,"dgCMatrix")
  X<-CreateSeuratObject(counts=object)
  #find features
  X <- FindVariableFeatures(X, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(X)
  X <- ScaleData(X, features = all.genes)
  X <- RunPCA(X, features = VariableFeatures(object = X_down),ndims.print = 1:5,npcs=npc)
  X <- FindNeighbors(X, dims = 1:npc)
  X <- FindClusters(X, resolution = resolution)
  tcell_tsne<-RunTSNE(X,reduction = "pca",cells = NULL,dims = 1:npc,features = NULL,seed.use = 1,tsne.method = "Rtsne",
                      dim.embed = 2,distance.matrix = NULL,reduction.name = "tsne",reduction.key = "tSNE_")
  
  DimPlot(tcell_tsne,reduction="tsne",combine=TRUE,label=TRUE)+labs(title=method, x="Dim1", y="Dim2")
  
}

plot_pca_umap<-function(npc, resolution, object, method){
  object<-as(object,"dgCMatrix")
  X<-CreateSeuratObject(counts=object)
  #find features
  X <- FindVariableFeatures(X, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(X)
  X <- ScaleData(X, features = all.genes)
  X <- RunPCA(X, features = VariableFeatures(object = X),ndims.print = 1:5,npcs=npc)
  X <- FindNeighbors(X, dims = 1:npc)
  X <- FindClusters(X, resolution = resolution)
  X <- RunUMAP(X, dims = 1:npc)
  DimPlot(X, reduction = "umap",combine=TRUE,label=TRUE)+labs(title=method, x="Dim1", y="Dim2")
}

