# Small helpers.

# Base R has this from 4.4; defined here for older R.
`%||%` <- function(a, b) if (is.null(a)) b else a
