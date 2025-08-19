library(reticulate)
library(Seurat)
library(distances)
library(dplyr)
#use_condaenv("Benchmark")
# MR: comment out python for now; TODO: build/integrate reticulate env
# sc = import("scanpy", convert=F)
# source_python("utils/normalization.py")

Normalization = function(seurat.obj, NormMethod){
  seurat.obj = NormMethod(seurat.obj)
  return(seurat.obj)
}


log1pCP10k = function(seurat.obj){
  print("Running log1pCP10k")
  seurat.obj = NormalizeData(seurat.obj,scale.factor = 10^4)
  return(seurat.obj)
}

log1pCPM = function(seurat.obj){
  print("Running log1pCPM")
  seurat.obj = NormalizeData(seurat.obj,scale.factor= 10^6)
  return(seurat.obj)
}

`log1pCPMedian` = function(seurat.obj){
  print("Running Scanpy default normalization")
  counts.layer = Layers(object = seurat.obj, search = "counts")
  save <- make.unique(names = gsub(
    pattern = "counts",
    replacement = "data",
    x = counts.layer
  ))
  for (i in seq_along(along.with = counts.layer)) {
    l <- counts.layer[i]
    count = t(LayerData(object = seurat.obj, layer = l, fast = NA))
    adata = sc$AnnData(count)
    norm = log1p(py_to_r(sc$pp$normalize_total(adata,inplace=F)$X))
    LayerData(
      object = seurat.obj,
      layer = save[i],
      features = Features(x = seurat.obj, layer = l),
      cells = Cells(x = seurat.obj, layer = l)
    ) <- t(norm)
  }
  return(seurat.obj)
}

sctransform = function(seurat.obj){
  print("Running sctransform")
  seurat.obj = SCTransform(seurat.obj, vst.flavor = "v2", verbose = FALSE, 
                           variable.features.n = nrow(seurat.obj), method = "glmGamPoi",min_cells=0)
  return(seurat.obj)
}

log1pPF = function(seurat.obj){
  print("Running log1pPF")
  counts.layer = Layers(object = seurat.obj, search = "counts")
  save <- make.unique(names = gsub(
    pattern = "counts",
    replacement = "data",
    x = counts.layer
  ))
  for (i in seq_along(along.with = counts.layer)) {
    l <- counts.layer[i]
    LayerData(
      object = seurat.obj,
      layer = save[i],
      features = Features(x = seurat.obj, layer = l),
      cells = Cells(x = seurat.obj, layer = l)
    ) <- logPF(t(LayerData(object = seurat.obj, layer = l, fast = NA)))
  }
  return(seurat.obj)
}

PFlog1pPF = function(seurat.obj){
  print("Running PFlog1pPF")
  counts.layer = Layers(object = seurat.obj, search = "counts")
  save <- make.unique(names = gsub(
    pattern = "counts",
    replacement = "data",
    x = counts.layer
  ))
  for (i in seq_along(along.with = counts.layer)) {
    l <- counts.layer[i]
    LayerData(
      object = seurat.obj,
      layer = save[i],
      features = Features(x = seurat.obj, layer = l),
      cells = Cells(x = seurat.obj, layer = l)
    ) <- PFlogPF(t(LayerData(object = seurat.obj, layer = l, fast = NA)))
  }
  return(seurat.obj)
}













