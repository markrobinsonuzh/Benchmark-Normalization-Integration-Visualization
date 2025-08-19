options(repos="https://stat.ethz.ch/CRAN",
<<<<<<< HEAD
	BioC_mirror="https://ftp.gwdg.de/pub/misc/bioconductor")
=======
        BioC_mirror="https://ftp.gwdg.de/pub/misc/bioconductor")

# notes: 
# - currently running it on Ubuntu 22.04
# - needs cmake >= 3.24
>>>>>>> 515a61481b12a0382559e6527091259d7318410e

install.packages("BiocManager")

pkgs <- c("scDesign3","SingleCellExperiment","readr","rmarkdown",
<<<<<<< HEAD
	  "Seurat","dplyr","GEOquery","R.utils", "here", "markdown",
=======
          "Seurat","dplyr","GEOquery","R.utils", "here", "markdown",
>>>>>>> 515a61481b12a0382559e6527091259d7318410e
          "remotes","distances","harmony","rliger")

BiocManager::install(pkgs, Ncpus=10)

BiocManager::install("miraisolutions/godmode", Ncpus=10)
BiocManager::install("satijalab/seurat-wrappers", Ncpus=10)
<<<<<<< HEAD


=======
>>>>>>> 515a61481b12a0382559e6527091259d7318410e
