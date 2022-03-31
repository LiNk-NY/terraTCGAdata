.DEFAULT_NAMESPACE <- "broad-firecloud-tcga"

#' Obtain or set the Terra Workspace Project Dataset
#'
#' Terra allows access to about 71 open access TCGA datasets. A dataset
#' workspace can be set using the `terraTCGAworkspace` function with a `projectName`
#' input. Use the `findTCGAworkspaces` function to list all of the available
#' open access TCGA data workspaces.
#'
#' @details 
#'     Note that GDC workspaces are not supported and are excluded
#'     from the search results. GDC workspaces use a Terra workflow to download
#'     TCGA data rather than providing Google Bucket storage locations for easy
#'     data retrieval.
#'
#' @aliases findTCGAworkspaces
#'
#' @param projectName character(1) A project code usually in the form of
#' `TCGA_CODE_OpenAccess_V1-0_DATA`. See `findTCGAworkspaces` for a list of
#' project codes.
#'
#' @param project character(1) A prefix for the regex search across all public
#' projects on the terra platform (default: `"^TCGA"`). Usually, this does not
#' change.
#'
#' @param cancerCode character(1) Corresponds to the TCGA cancer code (e.g,
#'   "ACC" for AdrenoCortical Carcinoma) of interest. The default value of
#'   (`.*`) provides all available cancer datasets.
#'
#' @return A Terra TCGA Workspace name
#'
#' @md
#'
#' @examples
#' if (AnVIL::gcloud_exists())
#'   findTCGAworkspaces()
#'
#' @export
terraTCGAworkspace <- function(projectName = NULL) {
    getOption(
        "terraTCGAdata.workspace",
        setTerraWorkspace(projectName = projectName)
    )
}

#' @describeIn terraTCGAworkspace Function to enumerate the available TCGA data
#'     workspaces in Terra
#'
#' @export
findTCGAworkspaces <- function(project = "^TCGA", cancerCode = ".*") {
    avs <- avworkspaces()
    gdcind <- -grep("GDCDR", avs[["name"]])
    avs <- avs[gdcind, ]
    project_code <- paste(project, cancerCode, sep = "_")
    ind <- grep(project_code, avs[["name"]])
    avs[ind, ]
}

setTerraWorkspace <-
    function(projectName, namespace = .DEFAULT_NAMESPACE)
{
    ws <- getOption("terraTCGAdata.workspace")
    if ((!nzchar(ws) || is.null(ws)) && missing(projectName)) {
        tcga_choices <- findTCGAworkspaces()[["name"]]
        wsi <- utils::menu(
            tcga_choices,
            title = "Select a TCGA terra Workspace: "
        )
        ws <- tcga_choices[wsi]
        options("terraTCGAdata.workspace" = ws)
    } else if (!missing(projectName)) {
        .isSingleChar(projectName)
        tcga_choices <- findTCGAworkspaces()[["name"]]
        validPC <- projectName %in% tcga_choices
        if (!validPC)
            stop("'projectName' not in the 'findTCGAworkspaces()' list ")
        if (!identical(ws, projectName))
            warning("Replacing 'terraTCGAData.workspace' with ", projectName)
        ws <- projectName
        options("terraTCGAdata.workspace" = ws)
    }
    ws
}

.isSingleChar <- function(ch) {
    stopifnot(is.character(ch), !is.na(ch), identical(length(ch), 1L))
}