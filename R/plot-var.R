#' Plot the varimpact results for a given variable
#'
#' Displays the adjusted treatment-specific means and the impact estimate.
#' @param var_name String name of the variable
#' @param vim Varimpact result object that contains the variable
#' @param digits Number of digits for rounding purposes.
#' @param verbose If true, display extra output.
# TODO: @example
#' @import ggplot2
#' @export
plot_var =
  function(var_name, vim, digits = 2L, verbose = FALSE) {

  # Confirm that we can plot this variable.

  # Can only plot numerics currently, need to expand to factors.
  numeric_vars = unique(vim$numeric_vims$results_by_level$name)
  factor_vars = unique(vim$factor_vims$results_by_level$name)
  if (var_name %in% numeric_vars) {
    results = vim$numeric_vims$results_by_level
  } else if (var_name %in% factor_vars) {
    results = vim$factor_vims$results_by_level
  } else {
    stop("There is no variable called", var_name, ".\n")
  }


  # Create plot dataframe.
  plot_data = results[results$name == var_name,
                      c("level", "level_label", "test_theta_tmle",
                        "test_var_tmle")]

  # Create color column based on min, max, and other.
  plot_data$color = "other"
  plot_data$color[which.min(plot_data$test_theta_tmle)] = "min"
  plot_data$color[which.max(plot_data$test_theta_tmle)] = "max"

  # Add "Impact" row to dataframe.
  result_row = rownames(vim$results_all) == var_name
  plot_data = rbind(plot_data,
                    list(level = NA,
                         level_label = "Impact",
                         test_theta_tmle = vim$results_all$Estimate[result_row],
                         # Don't have this yet.
                         test_var_tmle = NA,
                         color = "impact"
                         ))

  # red, blue, orange, gray
  # impact, max, min, other
  plot_palette = c("#d09996", "#95b4df", "#f2c197", "#d9d9d9")

  # Plot TSMs from $numeric_vims$results_by_level and the varimpact estimate.
  p = ggplot(data = plot_data,
             aes(x = factor(level_label),
                 y = test_theta_tmle,
                 #label = round(test_theta_tmle, 1L),
                 fill = factor(color))) +
    geom_col(width = 0.5) +
    geom_label(aes(label = round(test_theta_tmle, digits)),
               size = 6,
               fill = "white",
               color = "gray10",
               # label.color = "gray50",
               hjust = -0.2) +
    # Make the bar char horizontal rather than vertical
    coord_flip() +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 18),
          plot.title = element_text(size = 20),
          panel.background = element_rect(fill = "white", color = "gray50"),
          #plot.background = element_rect(fill = "gray95"),
          plot.background = element_rect(fill = "#f2f2f2"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    geom_hline(yintercept = 0, color = "gray90") +
    scale_x_discrete(limits = c("Impact", rev(setdiff(unique(plot_data$level_label), "Impact")))) +
    scale_y_continuous(expand = c(0.05, 0.08)) +
    scale_fill_manual(values = plot_palette) +
    labs(title = paste("Impact of", var_name))

  p
}
