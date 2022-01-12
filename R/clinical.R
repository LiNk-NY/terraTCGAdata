#' @import MultiAssayExperiment AnVIL
NULL

.PARTICIPANT_METADATA_COLS <-
    c("sample_id", "sample_type", "participant", "tcga_sample_id")

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
    function(tablename = "sample", metacols = .PARTICIPANT_METADATA_COLS)
{
    samples <- avtable(tablename)
    metadata <- samples[, metacols]
    samples <- samples[, !names(samples) %in% metacols]
    samples[, grep("clin", names(samples))]
}

#' Obtain clinical data
#'
#' The participant table may contain curated demographic information e.g.,
#' sex, age, etc.
#'
#' @param tablename The terra data model table from which to extract the
#'     clinical data (default: "sample")
#'
#' @param participants logical(1) Whether to merge the participant table
#'     from `avtable("participant")` to the clinical data
#'
#' @param columnName The name of the column to extract files, see
#'     `getClinicalTable` table. If not provided, the first column in the table
#'     will be used to obtain the clinical information.
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
#' getClinical()
#'
getClinical <-
    function(columnName, participants = TRUE, tablename = "sample",
        metacols = .PARTICIPANT_METADATA_COLS)
{
    allclin <- getClinicalTable(tablename = tablename, metacols = metacols)
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
#' @return A `tibble` of sample codes and frequency along with their
#'     definition and short letter code
#'
#' @md
#'
#' @examples
#'
#' sampleTypesTable()
#'
#' @export
sampleTypesTable <- function() {
    stable <- table(avtable("sample")[["sample_type"]])
    dataenv <- new.env(parent = emptyenv())
    utils::data("sampleTypes", envir = dataenv, package = "TCGAutils")
    stypes <- dataenv[["sampleTypes"]]
    sidx <- match(names(stable), stypes[["Short.Letter.Code"]])
    rest <- cbind(stypes[sidx, ], Frequency = as.numeric(stable))
    rownames(rest) <- NULL
    dplyr::as_tibble(rest)
}
