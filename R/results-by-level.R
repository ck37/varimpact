#' Aggregate the results_by_fold_and_level df into a results_by_level df.
#'
#' @param results_by_fold_and_level Dataframe containing the VIM results for
#' all levels of each variable across all CV folds.
#' @param verbose If true, display extra output.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize_all select funs mutate first
#' @importFrom modeest mlv
results_by_level =
  function(results_by_fold_and_level,
           verbose = FALSE) {

  results = results_by_fold_and_level %>%
    # NOTE: because we group by level_label, we are assuming that any histogram
    # penalization happened outside of the CV to ensure that the levels are
    # the same across training folds.
    group_by(name, level) %>%
    # TEMP: restrict to the most common label.
    #mutate(level_label = as.character(mlv(as.factor(level_label), method = "mfv")$M)) %>%
    # TEMP: restrict to the first label
    mutate(level_label = first(level_label)) %>%
    # Now we can also group by level_label because they will be the same for a given level.
    group_by(name, level, level_label) %>%
    # Remove test_msg for now.
    # TODO: take mode of test_msg or first value, rather than mean.
    select(-c(test_msg, train_msg)) %>%
    # this generates a warning in mean() because test_msg is a character not a numeric.
    summarize_all(funs(mean)) %>%
    select(-c(cv_fold, train_cell_size, test_cell_size))

  # Don't keep this as a tibble.
  as.data.frame(results)
}
