#' @rdname terraTCGAdata-defunct
#' @aliases findTCGAworkspaces
#'
#' @title findTCGAworkspaces is defunct in terraTCGAdata
#'
#' @description The function has been replaced by `selectTCGAworkspace` which
#'   provides an interactive workspace selection method.
#'
#' @return Defunct functions return an error
#'
#' @export
findTCGAworkspaces <- function() {
    .Defunct(
        "selectTCGAworkspace", "terraTCGAdata",
        c(
            "'findTCGAworkspaces' is defunct.\n",
            "Use 'selectTCGAworkspace' instead.\n",
            "See help(\"terraTCGAdata-defunct\")"
        )
    )
}
