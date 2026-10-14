# Development tools for varimpact. renv records them in renv.lock through
# these library() calls (the package itself does not depend on them, and dev/
# is in .Rbuildignore). roxygen2 is pinned to the version DESCRIPTION records
# (RoxygenNote), and devtools to 2.4.5, the last release that accepts it. See
# "Development" in the README.
library(roxygen2)
library(devtools)
library(testthat)
