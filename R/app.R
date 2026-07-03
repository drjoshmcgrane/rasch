# rmt :: Shiny launcher

#' Launch the rmt graphical interface
#'
#' Opens the Shiny application: data upload with ID, person-factor, and item
#' column nomination; the full analysis; interactive tables and
#' plots; and one-click export of every table and plot.
#'
#' @param ... Passed to \code{shiny::runApp}.
#' @return Called for its side effect of launching the app.
#' @examples
#' if (interactive()) run_app()
#' @export
run_app <- function(...) {
  for (pkg in c("shiny", "bslib", "DT")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("the rmt app needs the '", pkg, "' package: install.packages(\"", pkg, "\")")
  }
  dir <- system.file("shiny", package = "rmt")
  if (dir == "") stop("app not found: reinstall rmt")
  shiny::runApp(dir, ...)
}
