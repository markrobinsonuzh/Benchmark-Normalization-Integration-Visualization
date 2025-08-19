library(reticulate)
library(SeuratWrappers)
library(Seurat)
library(densvis)
library(phateR)
#library(irlba)
library(distances)
library(dplyr)
#use_condaenv("Benchmark")
# MR: commenting out Python stuff for now
# sc = import("scanpy", convert=F)

# install FFTW for FIt-SNE
#system(c("cd /home/zhiqian/FIt-SNE", "g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member"))
#system("cd /home/zhiqian/FIt-SNE")
#source("/home/zhiqian/FIt-SNE/fast_tsne.R",chdir=T) # needs changing 
#FAST_TSNE_SCRIPT_DIR = "/home/zhiqian/FIt-SNE" # needs changing

Visualization = function(seurat.obj, VisualMethod, n.pcs=50, n.cores=10){
  reduction = VisualMethod(seurat.obj, n.pcs=n.pcs, n.cores = n.cores) 
  return(reduction)
}

`BH-tsne` = function(seurat.obj, n.pcs=50, n.cores = 10){
  print("Running BH-tSNE")
  seurat.obj = RunTSNE(seurat.obj, dims = 1:n.pcs, 
                       reduction="integrated", tsne.method = "Rtsne",check_duplicates = FALSE,
                       num_threads = n.cores)
  return(Embeddings(seurat.obj,"tsne"))
}

`Seurat-UMAP` = function(seurat.obj, n.pcs=50, n.cores = 10){
  print("Running UMAP")
  seurat.obj = RunUMAP(seurat.obj, dims = 1:n.pcs, 
                       reduction="integrated",
                       umap.method = "umap-learn",
                       n_threads = n.cores
  )
  return(Embeddings(seurat.obj,"umap"))
}

densMAP = function(seurat.obj, n.pcs=50, n.cores = 10){
  print("Running densMAP")
  seurat.obj = RunUMAP(seurat.obj, umap.method = "umap-learn", 
                       dims = 1:n.pcs, reduction="integrated",
                       densmap = T, 
                       n_jobs = n.cores)
  return(Embeddings(seurat.obj,"umap"))
}

# `FIT-SNE` =  function(seurat.obj, n.pcs=50, n.cores = 10){
#   print("Running FIT-SNE")
#   source("/home/zhiqian/FIt-SNE/fast_tsne.R",chdir=T) 
#   FAST_TSNE_SCRIPT_DIR = "/home/zhiqian/FIt-SNE" # needs to change for different desktop
#   latent = seurat.obj@reductions$integrated@cell.embeddings[,1:n.pcs]
#   pcaInit = prcomp(latent, rank=2)$x 
#   pcaInit = pcaInit / (sd(pcaInit[,1])* (nrow(pcaInit)-1) / nrow(pcaInit) ) * 0.0001
#   vis = fftRtsne(latent, nthreads = n.cores, 
#                  #perplexity = 30,
#                  perplexity_list = c(30, nrow(latent)/100),  # memory issue, set n/1000 instead of n/100
#                  initialization = pcaInit, 
#                  learning_rate = nrow(latent)/12)
#   rownames(vis) = colnames(seurat.obj)
#   return(vis)
# }

# graphFA = function(seurat.obj, n.pcs = 50, n.cores = 10){
#   print("Running graphFA")
#   temp_count = matrix(0, nrow = ncol(seurat.obj), ncol = nrow(seurat.obj))
#   adata = sc$AnnData(temp_count)
#   latent.method.key = "X_integrated"
#   adata$obsm[latent.method.key] = seurat.obj@reductions$integrated@cell.embeddings[,1:n.pcs]
#   sc$pp$neighbors(adata,use_rep=latent.method.key, n_pcs=as.integer(n.pcs))
#   sc$tl$draw_graph(adata,layout="fa")
#   vis = as.matrix(adata$obsm$get('X_draw_graph_fa'))
#   rownames(vis) = colnames(seurat.obj)
#   return(vis)
# }

# `scanpy-UMAP` = function(seurat.obj, n.pcs=50, n.cores = 10){
#   print("Running scanpy UMAP")
#   temp_count = matrix(0, nrow = ncol(seurat.obj), ncol = nrow(seurat.obj))
#   adata = sc$AnnData(temp_count)
#   latent.method.key = "X_integrated"
#   adata$obsm[latent.method.key] = seurat.obj@reductions$integrated@cell.embeddings[,1:n.pcs]
#   sc$pp$neighbors(adata, n_pcs=as.integer(n.pcs), use_rep=latent.method.key)
#   sc$tl$umap(adata)
#   vis = as.matrix(adata$obsm$get('X_umap'))
#   rownames(vis) = colnames(seurat.obj)
#   return(vis)
# }

denSNE = function(seurat.obj, n.pcs=50, n.cores=10){
  cat("Running denSNE")
  latent = seurat.obj@reductions$integrated@cell.embeddings[,1:n.pcs]
  vis = densne(latent, num_threads = n.cores)
  rownames(vis) = colnames(seurat.obj)
  return(vis)
}

PHATE = function(seurat.obj, n.pcs=50, n.cores=5){
  cat("Running PHATE")
  latent = seurat.obj@reductions$integrated@cell.embeddings[,1:n.pcs]
  vis = phate(latent, n.jobs = n.cores)$embedding
  rownames(vis) = colnames(seurat.obj)
  return(vis)
}

