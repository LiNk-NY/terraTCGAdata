findTCGAworkspaces <- function(project = "TCGA", cancerCode = ".*") {
    avs <- avworkspaces()
    project_code <- paste(project, cancerCode, sep = "_")
    grep(project_code, avs[["name"]], value = TRUE)
}

setTCGAworkspace <- function(projectCode = NULL) {
    ## TODO: set package-wide option and/or ENV for namespace / name
    getOption("terraTCGAdata.workspace", terraSetWorkspace())
}

terraSetWorkspace <- function(projectCode) {
    ws <- getOption("terraTCGAdata.workspace")
    if (!nzchar(ws)) {
        tcga_choices <- findTCGAworkspaces()
        if (!missing(projectCode)) {
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
