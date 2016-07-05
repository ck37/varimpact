#' Get missingness for each column
#'
#' Function for getting total number missing values for vector
#'
#' @param x Vector, matrix, or dataframe
sum_na = function(x) {
  sum(is.na(x))
}
