.PARTICIPANT_METADATA_COLS <-
    c("sample_id", "sample_type", "participant", "tcga_sample_id")

#' Obtain the reference table for clinical data
#'
#' The column names in the output table can be used in the `getClinical`
#' function.
#'
#' @inheritParams getClinical
#'
#' @return A tibble with reference links to data resources on Terra
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
#' @param tablename The terra data model table from which to extract the
#'     clinical data (default: "sample")
#'
#' @param participants logical(1) Whether to merge the participant table
#'     from `avtable("participant")` to the clinical data
#'
#' @param columnName The name of the column to extract files, see
#'     `getClinicalTable` table.
#'
#' @param metacols The set of columns that comprise of the metadata columns.
#'     See the `.PARTICIPANT_METADATA_COLS` global variable
#'
#' @return A tibble of Google Storage resource locations e.g.,
#'     `gs://firecloud...`
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
    bfc <- BiocFileCache()
    rpath <- bfcquery(bfc, columnName, exact = TRUE)[["rpath"]]
    if (!length(rpath)) {
        rpath <- bfcnew(bfc, rname = columnName, rtype = "local")
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
    clinical <- readr::type_convert(clinical)
    coldata <- as(clinical, "DataFrame")
    if (!is.null(coldata[["patient.bcr_patient_barcode"]])) {
        rownames(coldata) <- coldata[["participant_id"]] <-
            toupper(coldata[["patient.bcr_patient_barcode"]])
        if (participants) {
            parts <- avtable("participant")
            parts <- as(parts, "DataFrame")
            rownames(parts) <- parts[["participant_id"]]
            ## replace munged COAD with TCGA
            parts[["participant_id"]] <-
                gsub("^[A-Z]{4}", "TCGA", parts[["participant_id"]])
            coldata <- merge(
                parts, coldata, by = "participant_id",
            )
            rownames(coldata) <- coldata[["participant_id"]]
        }
    }
    coldata
}

tablecolumn <- function(tablename = "sample", column = "sample_type") {
    table(avtable(tablename)[[column]])
}

#' Get an overview of the samples available in the workspace
#'
#' The function provides an overview of samples from the `avtables("sample")`
#' table for the current workspace. Along with the sample codes and frequencies,
#' the output provides a description for each code and the short letter codes.
#'
#' @return A `data.frame` of sample codes and frequency along with their
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
    data("sampleTypes", envir = dataenv, package = "TCGAutils")
    stypes <- dataenv[["sampleTypes"]]
    sidx <- match(names(stable), stypes[["Short.Letter.Code"]])
    rest <- cbind(stypes[sidx, ], Frequency = as.numeric(stable))
    rownames(rest) <- NULL
    rest
}
