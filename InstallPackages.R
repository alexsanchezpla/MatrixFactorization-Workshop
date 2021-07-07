if(!require(BiocManager)){
  install.packages("BiocManager", dep=TRUE)
}

installifnot <- function (pckgName, BioC=TRUE){
  if(BioC){
    if(!require(pckgName, character.only=TRUE, quietly = TRUE)){
      BiocManager::install(pckgName)
    }
  }else{
    if(!require(pckgName, character.only=TRUE, quietly=TRUE)){
      install.packages(pckgName, dep=TRUE)
    }
  }
}
installifnot("Biobase")
installifnot("NMF", BioC=FALSE)
installifnot("devtools", BioC=FALSE)
installifnot("knitr", BioC=FALSE)
installifnot("kableExtra", BioC=FALSE)
installifnot("ggplot2", BioC=FALSE)
installifnot("dplyr", BioC=FALSE)
installifnot("fastICA", BioC=FALSE)
installifnot("Rtsne", BioC=FALSE)

installifnot("matrixStats", BioC=FALSE)
installifnot("FactoMineR", BioC=FALSE)

installifnot("leukemiasEset")
installifnot("SummarizedExperiment")
installifnot("pheatmap")
installifnot("iClusterPlus")

installifnot("enrichR")

devtools::install_github("compgenomr/compGenomRData")
# If this does not work try to clone/download the repo into you computer
# Once you have it in your disk install it using:
# devtools::install_local("PATH_TO_compGenomRData")
