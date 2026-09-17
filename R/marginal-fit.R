#' A fit that predicts a constant probability
#'
#' Represents a binary regression that was not estimated from covariates,
#' because the outcome had too few observations in its rarer class to support
#' one. Standing in for the model rather than returning \code{NULL} lets the
#' same value be reused when the fit is applied to a validation fold, so the
#' training and validation folds agree.
#'
#' @param p Probability to predict, typically the marginal proportion of the
#'   outcome in the training data.
#'
#' @return An object of class \code{varimpact_marginal}.
#'
#' @seealso [predict.varimpact_marginal()]
#'
#' @keywords internal
marginal_fit = function(p) {
  stopifnot(length(p) == 1L, is.finite(p), p >= 0, p <= 1)
  structure(list(p = p), class = "varimpact_marginal")
}

#' Predict from a constant-probability fit
#'
#' @param object A \code{varimpact_marginal} object, from [marginal_fit()].
#' @param newdata Data to predict on. Only its number of rows is used.
#' @param ... Ignored. Present so that the method can be called with the
#'   arguments a \code{glm} fit would accept, such as \code{type = "response"}.
#'
#' @return A numeric vector of \code{object$p}, one element per row of
#'   \code{newdata}.
#'
#' @export
predict.varimpact_marginal = function(object, newdata, ...) {
  rep(object$p, NROW(newdata))
}
