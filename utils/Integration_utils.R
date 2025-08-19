source("utils/utils.R")
#library(reticulate)
library(SeuratWrappers)
library(Seurat)
library(harmony)
#library(rliger)
library(Rfast)
library(dplyr)

Integration = function(seurat.obj, IntegrateMethod, n.pcs = 50, features=rownames(seurat.obj), is.sctransform=F){
  if(!(identical(IntegrateMethod, scVI) | identical(IntegrateMethod, LIGERv2) | identical(IntegrateMethod, FastMNN))){
    if(!is.sctransform){
      seurat.obj = FindVariableFeatures(seurat.obj, nfeatures = nrow(seurat.obj))
      seurat.obj <- ScaleData(seurat.obj) 
    }
    seurat.obj <- RunPCA(seurat.obj, npcs = n.pcs)
    seurat.obj = IntegrateMethod(seurat.obj, is.sctransform, n.pcs, features)
  }else{
    seurat.obj = IntegrateMethod(seurat.obj, is.sctransform, n.pcs, features)
  }
  return(seurat.obj)
}

Harmony = function(seurat.obj, is.sctransform, n.pcs, features){
  print("Running Harmony")
  if(is.sctransform==F){
    seurat.obj = IntegrateLayers(
      object = seurat.obj, method = HarmonyIntegration,
      orig.reduction = "pca", new.reduction = "integrated",
      verbose = T,
      npcs = n.pcs,
      features = features
    ) 
  }else{
    seurat.obj = IntegrateLayers(
      object = seurat.obj, method = HarmonyIntegration,
      normalization.method = "SCT",
      orig.reduction = "pca", new.reduction = "integrated",
      verbose = T,
      npcs = n.pcs,
      features = features
    ) 
  }
  return(seurat.obj)
}


`Seurat-RPCA` = function(seurat.obj, is.sctransform, n.pcs, features){
  print("Running Seurat RPCA")
  if(is.sctransform==F){
    seurat.obj = IntegrateLayers(
      object = seurat.obj, method = RPCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated",
      verbose = T,
      dims = 1:n.pcs,
      features = features
    ) 
  }else{
    seurat.obj = IntegrateLayers(
      object = seurat.obj, method = RPCAIntegration,
      normalization.method = "SCT",
      orig.reduction = "pca", new.reduction = "integrated",
      verbose = T,
      dims = 1:n.pcs,
      features = features
    )
  }
  return(seurat.obj)
}

FastMNN = function(seurat.obj, is.sctransform, n.pcs, features){
  print("Running Fast MNN")
  if(is.sctransform==F){
    seurat.obj = IntegrateLayers(object = seurat.obj, method = FastMNNIntegration,
                                 new.reduction = 'integrated', verbose = T, orig.reduction = NULL,
                                 features = features,
                                 assay.type = "logcounts") #batch = seurat.obj$batch)
  }else{
    seurat.obj = IntegrateLayers(
      object = seurat.obj, method = FastMNNIntegration,
      new.reduction = "integrated",orig.reduction = NULL,
      verbose = T,
      batch = seurat.obj$batch,
      features = features,
      assay = "SCT",
      assay.type = "scaledata"
    )
  }
  return(seurat.obj)
}

`Seurat-CCA` = function(seurat.obj, is.sctransform, n.pcs, features){
  print("Running Seruat CCA")
  if(is.sctransform==F){
    seurat.obj = IntegrateLayers(
      object = seurat.obj, method = CCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated",
      verbose = T,
      dims = 1:n.pcs,
      features = features
    ) 
  }else{
    seurat.obj = IntegrateLayers(
      object = seurat.obj, method = CCAIntegration,
      normalization.method = "SCT",
      orig.reduction = "pca", new.reduction = "integrated",
      verbose = T,
      dims = 1:n.pcs,
      features = features
    )
  }
}

# scVI = function(seurat.obj, is.sctransform, n.pcs, features){
#   print("Running scVI")
#   seurat.obj = IntegrateLayers(
#     object = seurat.obj,
#     method = scVIIntegration,
#     new.reduction = "integrated",
#     conda_env = "Benchmark",
#     verbose = T,
#     features = features,
#     ndims = n.pcs,
#     layers = "counts",
#     orig.reduction = NULL,
#     scale.layer = NULL,
#     assay = "RNA"
#   ) 
#   return(seurat.obj)
# }

# LIGERv2 = function(seurat.obj, is.sctransform, n.pcs, features){
#   seurat.obj <- seurat.obj %>%
#     normalize() %>%
#     selectGenes(nGenes = nrow(seurat.obj)) %>%
#     scaleNotCenter()
#   seurat.obj <- seurat.obj %>%
#     runINMF(k = n.pcs) %>%
#     quantileNorm()
#   seurat.obj[["integrated"]] = seurat.obj[["inmfNorm"]]
#   return(seurat.obj)
# }
