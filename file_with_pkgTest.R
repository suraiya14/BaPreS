pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE,repos="http://cran.us.r-project.org")
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}