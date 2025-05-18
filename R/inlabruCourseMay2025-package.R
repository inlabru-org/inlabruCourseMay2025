#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
#' @import methods
#' @import stats
NULL

get_exports_pkg <- function(pkg) {
  nms <- getNamespaceExports(pkg)
  paste0(
    paste0("#' @import ", pkg, "\n"),
    paste0("#' @export\n", pkg, "::`", nms, "`", collapse = "\n"),
    "\n"
  )
}

get_exports <- function() {
  pkg <- c("inlabru", "fmesher")
  names(pkg) <- pkg
  lapply(pkg, get_exports_pkg)
}

update_exports <- function() {
  a <- get_exports()
  for (pkg in names(a)) {
    cat(a[[pkg]], file = paste0("R/export-", pkg, ".R"), append = FALSE)
  }
}
