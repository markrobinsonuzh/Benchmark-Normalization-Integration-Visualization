options(repos="https://stat.ethz.ch/CRAN",
	BioC_mirror="https://ftp.gwdg.de/pub/misc/bioconductor")

install.packages("BiocManager")

pkgs <- c("scDesign3","SingleCellExperiment","readr","rmarkdown",
	  "Seurat","dplyr","GEOquery","R.utils", "here", "markdown",
          "remotes","distances","harmony","rliger")

BiocManager::install(pkgs, Ncpus=10)

BiocManager::install("miraisolutions/godmode", Ncpus=10)
BiocManager::install("satijalab/seurat-wrappers", Ncpus=10)


