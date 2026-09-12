# rasch :: Shiny launcher

#' Launch the rasch point-and-click graphical interface
#'
#' Opens the Shiny application for fitting models and examining their tables,
#' plots and diagnostics. The R code for each result is available in the app.
#' Analyses can be saved and reopened, or exported as HTML, Word or PDF
#' reports.
#'
#' The app's interface packages ('shiny', 'bslib', 'DT', 'bsicons', and
#' 'callr' for cancellable EFRM and fit-bootstrap estimation) are suggested rather than
#' required by the package. If any are missing,
#' \code{run_app} lists them all and, in an interactive session, offers to
#' install them before launching.
#'
#' Older saved analyses are checked against the current person-scoring
#' algorithm. If their scores differ, refit the analysis before reopening it;
#' the original file is left unchanged. Its source data can be recovered with
#' \code{readRDS(file)$data}.
#' Saved EFRM and CJ frame fits without a current likelihood-check record
#' also require refitting, as does an extended frame fit whose stored score
#' curves predate the shared design enumeration. Their settings remain in
#' \code{readRDS(file)$settings}.
#' Superseded DIF results, or CJ DIF without verified judge-role alignment,
#' are omitted with a warning. Rerun those analyses before reporting them.
#'
#' @param ... Passed to \code{shiny::runApp}.
#' @return Called for its side effect of launching the app.
#' @examples
#' if (interactive()) run_app()
#' @export
run_app <- function(...) {
  .app_require(c("shiny", "bslib", "DT", "bsicons", "callr"))
  dir <- system.file("shiny", package = "rasch")
  if (dir == "") stop("app not found: reinstall rasch")
  old_limit <- getOption("shiny.maxRequestSize")
  if (is.null(old_limit) || old_limit < 100 * 1024^2) {
    options(shiny.maxRequestSize = 100 * 1024^2)
    on.exit(options(shiny.maxRequestSize = old_limit), add = TRUE)
  }
  shiny::runApp(dir, ...)
}

# Reconstruct one of the datasets bundled with the graphical interface. This
# is internal: it exists so the analysis code shown for an example run is
# executable without embedding hundreds of lines of generated data.
.app_example_data <- function(name) {
  name <- match.arg(name, c("pcm", "dich", "rsm", "mfrm", "efrm", "btl"))
  path <- system.file("shiny", "examples.R", package = "rasch")
  if (!nzchar(path)) stop("bundled app examples not found: reinstall rasch")
  env <- new.env(parent = asNamespace("rasch"))
  sys.source(path, envir = env)
  switch(name,
    pcm = env$.demo_data(), dich = env$.demo_dich(),
    rsm = env$.demo_rsm(), mfrm = env$.demo_mfrm(),
    efrm = env$.demo_efrm(), btl = env$.demo_btl())
}

# Check the app's display packages in one pass. Reports everything missing
# at once, and in an interactive session offers to install the whole set
# (with explicit consent) rather than sending the user through one
# stop-install-retry loop per package.
.app_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1),
                          quietly = TRUE)]
  if (!length(missing)) return(invisible(TRUE))
  cmd <- sprintf("install.packages(c(%s))",
                 paste(sprintf('\"%s\"', missing), collapse = ", "))
  if (interactive()) {
    ans <- utils::askYesNo(sprintf(
      "The rasch app needs %d more package(s): %s. Install now?",
      length(missing), paste(missing, collapse = ", ")), default = TRUE)
    if (isTRUE(ans)) {
      utils::install.packages(missing)
      still <- missing[!vapply(missing, requireNamespace, logical(1),
                               quietly = TRUE)]
      if (!length(still)) return(invisible(TRUE))
      stop("installation did not complete for: ",
           paste(still, collapse = ", "), "; try ", cmd, call. = FALSE)
    }
  }
  stop("the rasch app needs ", paste(missing, collapse = ", "),
       "; run ", cmd, call. = FALSE)
}
