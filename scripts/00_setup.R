packages <- c(
  "BiocManager",
  "Biostrings",
  "DECIPHER",
  "ape",
  "seqinr",
  "tidyverse",
  "ggplot2",
  "readr",
  "patchwork",
  "knitr",
  "rmarkdown"
)

install_if_missing <- function(pkg){

  if(!requireNamespace(pkg, quietly = TRUE)){

    if(pkg %in% c("Biostrings","DECIPHER")){

      BiocManager::install(pkg)

    } else {

      install.packages(pkg)

    }

  }

}

lapply(packages, install_if_missing)