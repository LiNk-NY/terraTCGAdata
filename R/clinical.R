#' @import MultiAssayExperiment AnVIL
NULL

.PARTICIPANT_METADATA_COLS <- c("sample_id", "sample_type", "participant",
    "tcga_sample_id", "submitter_id", "participant_id", "project_id")

.isGDC <- function(workspacename) {
    if (grepl("GDC", workspacename, fixed = TRUE))
        stop("GDC Workspaces are not supported.")
}

#' Obtain the reference table for clinical data
#'
#' The column names in the output table can be used in the `getClinical`
#' function.
#'
#' @inheritParams getClinical
#'
#' @md
#'
#' @return A tibble of Google Storage resource locations e.g.,
#'     `gs://firecloud...`
#'
#' @export
getClinicalTable <-
    function(
        tablename = .DEFAULT_TABLENAME, metacols = .PARTICIPANT_METADATA_COLS,
        workspace = terraTCGAworkspace(), namespace = .DEFAULT_NAMESPACE,
        verbose = TRUE
    )
{
    if (verbose)
        message(
            "Using namespace/workspace: ", paste0(namespace, "/", workspace)
        )
    .isGDC(workspace)
    
    avtab <- avtable(
        tablename, namespace = namespace, name = workspace
    )
    
    metadata <- avtab[, names(avtab) %in% metacols]
    avtab <- avtab[, !names(avtab) %in% metacols]
        
    avtab[, grep("clin", names(avtab))]
}

#' Obtain clinical data
#'
#' The participant table may contain curated demographic information e.g.,
#' sex, age, etc.
#'
#' @param columnName The name of the column to extract files, see
#'     `getClinicalTable` table. If not provided, the first column in the table
#'     will be used to obtain the clinical information.
#'
#' @param participants logical(1) Whether to merge the participant table
#'     from `avtable("participant")` to the clinical data
#'
#' @param tablename The Terra data model table from which to extract the
#'     clinical data (default: "sample")
#'
#' @param workspace character(1) The Terra Data Resources workspace from which
#'     to pull TCGA data (default: see `terraTCGAworkspace()`). This is set to a
#'     package-wide option.
#'
#' @param namespace character(1) The Terra Workspace Namespace that
#'     defaults to "broad-firecloud-tcga" and rarely needs to be changed.
#'
#' @param verbose logical(1) Whether to output additional information regarding
#'     the workspace and namespace (default: `TRUE`).
#'
#' @param metacols The set of columns that comprise of the metadata columns.
#'     See the `.PARTICIPANT_METADATA_COLS` global variable
#'
#' @return A `DataFrame` with clinical information from TCGA. The metadata i.e.,
#'     `metadata(object)` includes the `columnName` used to obtain the data.
#'
#' @export
#'
#' @md
#'
#' @examples
#'
#' getClinical(workspace = "TCGA_ACC_OpenAccess_V1-0_DATA")
#'
getClinical <-
    function(columnName, participants = TRUE, tablename = .DEFAULT_TABLENAME,
        workspace = terraTCGAworkspace(), namespace = .DEFAULT_NAMESPACE,
        verbose = TRUE, metacols = .PARTICIPANT_METADATA_COLS
    )
{
    allclin <- getClinicalTable(
        tablename = tablename, workspace = workspace, namespace = namespace,
        metacols = metacols
    )
    if (missing(columnName))
        columnName <- names(allclin)[1L]
    clinfiles <- unlist(unique(allclin[, columnName]))
    bfc <- BiocFileCache::BiocFileCache()
    rpath <- BiocFileCache::bfcquery(bfc, columnName, exact = TRUE)[["rpath"]]
    if (!length(rpath)) {
        rpath <- BiocFileCache::bfcnew(bfc, rname = columnName, rtype = "local")
        dir.create(rpath)
        gsutil_cp(clinfiles, rpath)
    }

    clinrows <- lapply(
        list.files(rpath, full.names = TRUE),
        readr::read_tsv,
        show_col_types = FALSE
    )
    allclins <- lapply(clinrows,
        tidyr::pivot_wider, names_from = "node_name", values_from = "node_value"
    )
    clinical <- dplyr::bind_rows(allclins)
    readr::type_convert(clinical)
}

#' Get an overview of the samples available in the workspace
#'
#' The function provides an overview of samples from the `avtables("sample")`
#' table for the current workspace. Along with the sample codes and frequencies,
#' the output provides a description for each code and the short letter codes.
#'
#' @inheritParams getClinical
#'
#' @return A `tibble` of sample codes and frequency along with their
#'     definition and short letter code
#'
#' @md
#'
#' @examples
#'
#' sampleTypesTable(workspace = "TCGA_COAD_OpenAccess_V1-0_DATA")
#'
#' @export
sampleTypesTable <-
    function(
        workspace = terraTCGAworkspace(),
        namespace = .DEFAULT_NAMESPACE,
        tablename = .DEFAULT_TABLENAME,
        verbose = TRUE
    )
{
    if (verbose)
        message(
            "Using namespace/workspace: ", paste0(namespace, "/", workspace)
        )
    avtab <- avtable(tablename, namespace = namespace, name = workspace)
    stable <- table(avtab[["sample_type"]])
    dataenv <- new.env(parent = emptyenv())
    utils::data("sampleTypes", envir = dataenv, package = "TCGAutils")
    stypes <- dataenv[["sampleTypes"]]
    sidx <- match(names(stable), stypes[["Short.Letter.Code"]])
    rest <- cbind(stypes[sidx, ], Frequency = as.numeric(stable))
    rownames(rest) <- NULL
    dplyr::as_tibble(rest)
}
