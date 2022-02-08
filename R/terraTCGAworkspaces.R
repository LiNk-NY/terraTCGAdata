#' Obtain or set the Terra Workspace Project Dataset
#'
#' Terra allows access to about 71 open access TCGA datasets. A dataset
#' workspace can be set using the `terraWorkspace` function with a `projectName`
#' input. Use the `findTCGAworkspaces` function to list all of the available
#' open access TCGA data workspaces.
#'
#' @aliases findTCGAworkspaces
#'
#' @param projectName character(1) A project code usually in the form of
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
terraWorkspace <- function(projectName = NULL) {
    getOption(
        "terraTCGAdata.workspace",
        setTerraWorkspace(projectName = projectName)
    )
}

#' @describeIn terraWorkspace Function to enumerate the available TCGA data
#'     workspaces in Terra
#'
#' @export
findTCGAworkspaces <- function(project = "^TCGA", cancerCode = ".*") {
    avs <- avworkspaces()
    project_code <- paste(project, cancerCode, sep = "_")
    ind <- grep(project_code, avs[["name"]])
    avs[ind, ]
}

setTerraWorkspace <-
    function(projectName, namespace = "broad-firecloud-tcga")
{
    ws <- getOption("terraTCGAdata.workspace")
    if (!nzchar(ws) || is.null(ws)) {
        tcga_choices <- findTCGAworkspaces()[["name"]]
        if (!is.null(projectName)) {
            validPC <- projectName %in% tcga_choices
            if (!validPC)
                stop("'projectName' not in the 'findTCGAworkspaces()' list ")
        } else {
            wsi <- menu(
                tcga_choices,
                title = "Select a TCGA terra Workspace: "
            )
            ws <- tcga_choices[wsi]
        }
        options("terraTCGAdata.workspace" = ws)
    }
    c(namespace = namespace, name = ws)
}
