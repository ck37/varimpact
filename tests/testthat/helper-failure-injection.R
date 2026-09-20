# Helpers for reaching varimpact's recovery paths on purpose.
#
# The package wraps its per-bin and per-fold estimation steps in try() and
# carries on when one fails. Ordinary data does not trip those paths reliably
# - estimate_tmle2() and apply_tmle_to_validation() are hardened enough that
# degenerate inputs are handled before they can fail - so these helpers replace
# one internal function for the duration of an expression.
#
# testthat's local_mocked_bindings() does the namespace surgery (unlocking the
# binding, restoring it afterwards) that an earlier version of this helper did
# by hand with unlockBinding().

# Run expr with `fun_name` replaced by `make_replacement(original)`, where
# `original` is the real function, so the replacement can wrap it.
with_replaced = function(fun_name, make_replacement, expr) {
  original = get(fun_name, envir = asNamespace("varimpact"))
  args = list(make_replacement(original), .package = "varimpact", .env = environment())
  names(args)[1] = fun_name
  do.call(testthat::local_mocked_bindings, args)
  force(expr)
}

# Run expr with `fun_name` failing on the calls for which `should_fail(i)` is
# TRUE, counting calls from 1, and behaving normally otherwise.
with_failing = function(fun_name, should_fail, expr) {
  with_replaced(fun_name, function(original) {
    call_i = 0L
    function(...) {
      call_i <<- call_i + 1L
      if (should_fail(call_i)) {
        stop("synthetic failure in ", fun_name, "(), call ", call_i)
      }
      original(...)
    }
  }, expr)
}
