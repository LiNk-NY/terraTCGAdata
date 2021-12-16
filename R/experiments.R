#' Obtain a reference table for assay data in the Terra data model
#'
#' The column names in the output can be used in the `getAssayData` function.
#'
#' @inheritParams getAssayData
#'
#' @return A tibble of pointers to resources within the Terra data model
#'
#' @export
getAssayTable <-
    function(tablename = "sample", metacols = .PARTICIPANT_METADATA_COLS)
{
    samples <- avtable(tablename)
    metadata <- samples[, metacols]
    samples <- samples[, !names(samples) %in% participant_meta]
    samples[, grep("clin", names(samples), invert = TRUE)]
}

# assayName <- "protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data"

#' Obtain assay datasets from Terra
#'
#' @inheritParams getClinical
#'
#' @param assayName character() The name of the assay dataset column from
#'     `getAssayTable` to import into the current workspace.
#'
#' @param sampleCode character(1) The sample code used to filtering samples
#'     e.g., "01" for Primary Solid Tumors, see
#'     `data("sampleTypes", package = "TCGAutils")` for reference
#'
#' @return An `ExperimentList` of assays selected
#'
#' @export
getAssayData <-
    function(assayName, sampleCode = "01", tablename = "sample",
        metacols = .PARTICIPANT_METADATA_COLS)
{
    if (missing(assayName))
        stop("Select an assay name from 'getAssayTable'")

    assayTable <- getAssayTable(tablename = tablename, metacols = metacols)
    assayfiles <- unlist(unique(na.omit(assayTable[, assayName])))

    bfc <- BiocFileCache()
    rpath <- bfcquery(bfc, assayName)[["rpath"]]
    if (!length(rpath)) {
        rpath <- bfcnew(bfc, rname = assayName, rtype = "local")
        dir.create(rpath)
        gsutil_cp(assayfiles, rpath)
    }

    assaycopied <- list.files(rpath, pattern = "\\.txt", full.names = TRUE)
    tcgaids <- gsub(".data.txt", "", basename(assaycopied), fixed = TRUE)

    datarows <- lapply(
        list.files(rpath, pattern = "TCGA.*data.txt", full.names = TRUE),
        readr::read_tsv
    )
    cnames <- unlist(datarows[[1]][1, ], use.names = FALSE)
    datarows <- lapply(datarows, function(x) x[-1, ])
    rowlist <- lapply(datarows, `[[`, 1L)
}
