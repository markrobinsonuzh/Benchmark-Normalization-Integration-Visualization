# Edit the FastMMNIntegration method to accomodate the sctransform results in Seurat
FastMNNIntegration <- function(
    object,
    assay = NULL,
    orig = NULL,
    groups = NULL,
    layers = NULL,
    scale.layer = NULL,
    features = 2000,
    new.reduction = "integrated.mnn",
    reduction.key = "mnn_",
    reconstructed.assay = "mnn.reconstructed",
    verbose = TRUE,
    ...
) {
  print("new fastmnn")
  object <- CreateSeuratObject(object)
  if (is.numeric(x = features)) {
    if (verbose) {
      message(paste("Computing", features, "integration features"))
    }
    features <- SelectIntegrationFeatures5(object = object, features = features)
  }
  if(assay == "SCT"){
    layers <- layers %||% Layers(object, search = "scale.data")  # change to scale.data for sctransform
  }else{
    layers <- layers %||% Layers(object, search = "data")
  }
  if (verbose) {
    message("Converting layers to SingleCellExperiment")
  }
  objects.sce <- lapply(
    X = layers,
    FUN = function(x, f) {
      return(as.SingleCellExperiment(
        x = subset(x = object,
                   features = f,
                   cells = colnames(LayerData(object, layer = x)))
      )
      )
    },
    f = features
  )
  if (verbose) {
    message("Running fastMNN")
  }
  out <- do.call(
    what = batchelor::fastMNN,
    args = c(
      objects.sce,
      list(...)
    )
  )
  colnames(x = SingleCellExperiment::reducedDim(x = out)) <- paste0(reduction.key, 1:ncol(x = SingleCellExperiment::reducedDim(x = out)))
  reduction <- CreateDimReducObject(
    embeddings = SingleCellExperiment::reducedDim(x = out),
    loadings = as.matrix(SingleCellExperiment::rowData(x = out)),
    assay = DefaultAssay(object = object),
    key = reduction.key
  )
  # Add reconstructed matrix (gene x cell)
  reconstructed_assay <- CreateAssayObject(
    data = as.matrix(SummarizedExperiment::assay(x = out)),
  )
  # Add variable features
  VariableFeatures(object = reconstructed_assay) <- features
  #Tool(object = object) <- S4Vectors::metadata(x = out)
  #object <- LogSeuratCommand(object = object)
  output.list <- list(reduction, reconstructed_assay)
  names(output.list) <- c(new.reduction, reconstructed.assay)
  return(output.list)
}

# Edit the labels_new in Seurat to speed up
labels_new <- function(object, values, select = c("first", "last", "common", "all"), simplify = TRUE, ...) {
  print("Using new labels")
  select <- select[1L]
  select <- match.arg(arg = select)
  values <- intersect(values, rownames(object))
  
  if (length(values) == 0) {
    return(character(0))
  }
  
  # Precompute logical matrix and column names
  cmap_data <- as.matrix(object)
  colnames_object <- colnames(object)
  rownames_object <- rownames(object)
  
  # Get row indices for values
  row_indices <- match(values, rownames_object)
  
  # Initialize the list to store results
  obs <- vector("list", length(values))
  names(obs) <- values
  
  # Direct indexing to replace sapply
  for (i in seq_along(row_indices)) {
    row_idx <- row_indices[i]
    vals <- colnames_object[cmap_data[row_idx, , drop = FALSE]]
    if (length(vals) > 0) {
      obs[[i]] <- vals
    }
  }
  
  obs <- Filter(length, obs)
  
  obs <- switch(select, 
                first = lapply(obs, `[[`, 1L), 
                last = lapply(obs, function(x) x[[length(x)]]), 
                common = {
                  counts <- table(unlist(obs))
                  tmp <- obs
                  obs <- vector("character", length(tmp))
                  names(obs) <- names(tmp)
                  for (i in seq_along(obs)) {
                    obs[i] <- names(which.max(counts[names(counts) %in% tmp[[i]]]))
                  }
                  obs
                }, 
                obs)
  
  if (isTRUE(simplify)) {
    tmp <- obs
    obs <- unlist(tmp)
    names(obs) <- make.unique(rep(names(tmp), times = lengths(tmp)))
  }
  
  return(obs)
}
godmode:::assignAnywhere("labels", labels_new)
godmode:::assignAnywhere("FastMNNIntegration", FastMNNIntegration)
