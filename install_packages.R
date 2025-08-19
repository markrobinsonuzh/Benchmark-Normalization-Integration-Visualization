options(repos="https://stat.ethz.ch/CRAN",
	BioC_mirror="https://ftp.gwdg.de/pub/misc/bioconductor")

# notes: 
# - currently running it on Ubuntu 22.04
# - will need cmake >= 3.24
# - not installing Rliger yet

install.packages("BiocManager")

pkgs <- c("scDesign3","SingleCellExperiment","readr","rmarkdown",
	  "Seurat","dplyr","GEOquery","R.utils", "here", "markdown",
          "remotes","distances","harmony","Rfast")

BiocManager::install(pkgs, Ncpus=10)

BiocManager::install("miraisolutions/godmode", Ncpus=10)
BiocManager::install("satijalab/seurat-wrappers", Ncpus=10)

