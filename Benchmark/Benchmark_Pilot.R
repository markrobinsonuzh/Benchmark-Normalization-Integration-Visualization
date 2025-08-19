
# build to HTML with
# knitr::spin("Benchmark/Benchmark_Pilot.R") 

utils_dir <- here::here("utils")

suppressPackageStartupMessages({
  source(file.path(utils_dir,"utils.R"))
  source("utils/Normalization_utils.R")
  source("utils/Integration_utils.R")
  source("utils/Visualization_utils.R")
})

data_dir <- here::here("Benchmark")
    
data = readRDS(file.path(data_dir,"Data/Simulation/Pilot/MousePancreas_simulation.rds"))
dataset = "Pilot"

# NormalizationMethods = c("log1pCP10k","log1pPF","log1pCPMedian",
#                          "sctransform","log1pCPM","PFlog1pPF")
# IntegrationMethods = c("Seurat-CCA","Seurat-RPCA",
#                        "FastMNN","Harmony")
# VisualizationMethods = c("BH-tsne","FIt-SNE","Seurat-UMAP",
#                          "densMAP","graphFA","scanpy-UMAP")

# MR: reduce the space a bit
NormalizationMethods = c("log1pCP10k","sctransform")
IntegrationMethods = c("Seurat-CCA","Harmony")
VisualizationMethods = c("Seurat-UMAP","densMAP")


Time = expand.grid(Normalization = NormalizationMethods,
                   Integration = IntegrationMethods, 
                   Visualization = VisualizationMethods)
# Time = rbind(Time, expand.grid(Normalization = "None", Integration = c("scVI"), Visualization = VisualizationMethods))
# Time = rbind(Time, expand.grid(Normalization = "None", Integration = c("LIGER"), Visualization = VisualizationMethods))
Time$Norm.time = 0
Time$Integrate.time = 0
Time$Visual.time = 0
Time = Time[order(Time$Integration,Time$Normalization),]
rownames(Time) = 1:nrow(Time)
Time

counts = data$counts
meta = data$meta[,c("celltype","batch")]
rm(data)
rownames(meta) = colnames(counts)
seurat.obj = CreateSeuratObject(counts, meta.data=meta)
rm(counts, meta); gc()

seurat.obj[["RNA"]] <- split(seurat.obj[["RNA"]], f = seurat.obj$batch)
old_preprocessing = c("NA","NA")
for(i in 1:nrow(Time)){
  set.seed(42)
  print(paste0(i, "-th combination"))
  NormalizeMethod = as.character(Time$Normalization[i])
  IntegrateMethod = as.character(Time$Integration[i])
  VisualizeMethod =  as.character(Time$Visualization[i])
  new_preprocessing = c(NormalizeMethod, IntegrateMethod)
  if( !all(new_preprocessing == old_preprocessing) ){
    seurat = seurat.obj
    
    if( NormalizeMethod=="None" ){
      Time[i,]$Norm.time = 0
    }else{  
      start = Sys.time()
      seurat = Normalization(seurat, get(NormalizeMethod))
      end = Sys.time()
      print(end-start)
      Time[i,]$Norm.time = difftime(end,start,unit = "secs")
    }
    
    if( NormalizeMethod == "sctransform" ){
      start = Sys.time()
      seurat = Integration(seurat, get(IntegrateMethod), n.pcs = 50, 
                           features = row.names(seurat), is.sctransform=TRUE)
      end = Sys.time()
      print(end-start)
      Time[i,]$Integrate.time= difftime(end,start,unit = "secs")
    }else{
      start = Sys.time()
      seurat = Integration(seurat, get(IntegrateMethod), n.pcs = 50, 
                           features = row.names(seurat), is.sctransform=FALSE)
      end = Sys.time()
      print(end-start)
      Time[i,]$Integrate.time= difftime(end,start,unit = "secs")
    }
  }else{
    Time[i,]$Norm.time = Time[i-1,]$Norm.time
    Time[i,]$Integrate.time = Time[i-1,]$Integrate.time
  }
  gc()
  start = Sys.time()
  Visual = Visualization(seurat, get(VisualizeMethod))
  end = Sys.time()
  print(end-start)
  Time[i,]$Visual.time= difftime(end,start,unit = "secs")
  print(dim(Visual))
  saveRDS(Time, file.path(data_dir, "Results/Pilot/Runtime/Pilot_TimeComplexity.rds"))
  saveRDS(Visual, paste0(data_dir, "Results/Pilot/Embeddings/",dataset,"_",NormalizeMethod,"+",IntegrateMethod,"+",VisualizeMethod,".rds"))
  old_preprocessing = new_preprocessing
}
