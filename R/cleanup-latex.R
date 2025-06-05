#' Clean up LaTeX files created by exportLatex
#'
#' This function removes LaTeX files that are typically created by the exportLatex() function.
#' It's designed to be used after exportLatex() calls to clean up temporary files.
#'
#' @param dir Directory where LaTeX files are located (default: current directory)
#' @param outname Prefix for the LaTeX files (default: empty string)
#' @param verbose If TRUE, print messages about which files were removed
#'
#' @return Invisibly returns a logical vector indicating which files were successfully removed
#'
#' @examples
#' \dontrun{
#' # After calling exportLatex()
#' exportLatex(vim)
#' cleanup_latex_files()
#' 
#' # With custom directory and prefix
#' exportLatex(vim, outname = "myresults_", dir = "output/")
#' cleanup_latex_files(dir = "output/", outname = "myresults_")
#' }
#'
#' @export
cleanup_latex_files <- function(dir = ".", outname = "", verbose = FALSE) {
  
  # Define the standard LaTeX file names that exportLatex() creates
  latex_files <- c(
    paste0(dir, "/", outname, "varimpByFold.tex"),
    paste0(dir, "/", outname, "varimpAll.tex"),
    paste0(dir, "/", outname, "varimpConsistent.tex")
  )
  
  # Check which files exist
  existing_files <- latex_files[file.exists(latex_files)]
  
  if (length(existing_files) == 0) {
    if (verbose) {
      cat("No LaTeX files found to clean up.\n")
    }
    return(invisible(logical(0)))
  }
  
  if (verbose) {
    cat("Cleaning up LaTeX files:\n")
    cat(paste("  -", basename(existing_files), collapse = "\n"), "\n")
  }
  
  # Remove the files
  removal_success <- suppressWarnings({
    file.remove(existing_files)
  })
  
  if (verbose) {
    successful_removals <- sum(removal_success)
    cat("Successfully removed", successful_removals, "of", length(existing_files), "files.\n")
  }
  
  return(invisible(removal_success))
}