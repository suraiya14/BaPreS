pkgTest_decipher <- function(x)
{if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
  if (!require(x,character.only = TRUE))
  {
    BiocManager::install(x)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}