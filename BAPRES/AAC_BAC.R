#' Amino Acid Composition Descriptor
#'
#' This function calculates the Amino Acid Composition descriptor (dim: 20).
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @return A length 20 named vector
#'


extractAAC_BAC <- function(x) {
  
  AADict <- c(
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V" 
  )
  
  if (all(strsplit(x, split = "")[[1]] %in% AADict) == FALSE) {
    stop("x has unrecognized amino acid type")
  }

  # 20 Amino Acid Abbrevation Dictionary from
  # https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties

  

  AAC <- summary(
    factor(strsplit(x, split = "")[[1]], levels = AADict),
    maxsum = 21
  ) / nchar(x)

  return (AAC)
}



