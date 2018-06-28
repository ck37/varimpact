#' Aggregate the results_by_fold_and_level df into a results_by_level df.
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
    summarize_all(funs(mean)) %>%
    select(-c(cv_fold, train_msg, test_msg, train_cell_size, test_cell_size))

  # Don't keep this as a tibble.
  as.data.frame(results)
}
