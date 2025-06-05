# Global variable declarations to avoid R CMD check NOTEs
# These variables are used in dplyr operations and ggplot2

# Variables used in dplyr operations
utils::globalVariables(c(
  "name", "level", "level_label", "test_msg", "train_msg", 
  "cv_fold", "train_cell_size", "test_cell_size",
  "rawp", "BH", "AvePsi", "Consistent",
  "test_theta_tmle", "color"
))

# Function used in dplyr operations
utils::globalVariables("desc")