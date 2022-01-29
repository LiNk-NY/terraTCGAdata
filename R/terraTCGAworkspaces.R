findTCGAworkspaces <- function(project = "TCGA", cancerCode = ".*") {
    avs <- avworkspaces()
    project_code <- paste(project, cancerCode, sep = "_")
    grep(project_code, avs[["name"]], value = TRUE)
}

setTCGAworkspace <- function(projectCode = NULL) {
    ## TODO: set package-wide option and/or ENV for namespace / name
    getOption("terraTCGAdata.workspace", projectCode)
}