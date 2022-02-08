#' Obtain or set the Terra Workspace Project Dataset
#'
#' Terra allows access to about 71 open access TCGA datasets. A dataset
#' workspace can be set using the `terraWorkspace` function with a `projectCode`
#' input. Use the `findTCGAworkspaces` function to list all of the available
#' open access TCGA data workspaces.
#'
#' @aliases findTCGAworkspaces
#'
#' @param projectCode character(1) A project code usually in the form of
#' `TCGA_CODE_OpenAccess_V1-0_DATA`. See `findTCGAworkspaces` for a list of
#' project codes.
#'
#' @return A Terra TCGA Workspace name
#'
#' @md
#'
#' @examples
#'
#' findTCGAworkspaces()
#'
#' @export
terraWorkspace <- function(projectCode = NULL) {
    getOption(
        "terraTCGAdata.workspace",
        setTerraWorkspace(projectCode = projectCode)
    )
}

#' @describeIn terraWorkspace Function to enumerate the available TCGA data
#'     workspaces in Terra
#'
#' @export
findTCGAworkspaces <- function(project = "TCGA", cancerCode = ".*") {
    avs <- avworkspaces()
    project_code <- paste(project, cancerCode, sep = "_")
    grep(project_code, avs[["name"]], value = TRUE)
}

setTerraWorkspace <- function(projectCode) {
    ws <- getOption("terraTCGAdata.workspace")
    if (!nzchar(ws)) {
        tcga_choices <- findTCGAworkspaces()
        if (!is.null(projectCode)) {
            validPC <- projectCode %in% tcga_choices
            if (!validPC)
                stop("'projectCode' not in the 'findTCGAworkspaces()' list ")
        } else {
            ws <- menu(
                tcga_choices,
                title = "Select a TCGA terra Workspace: "
            )
            ws <- tcga_choices[ws]
        }
        options("terraTCGAdata.workspace" = ws)
    }
    ws
}
