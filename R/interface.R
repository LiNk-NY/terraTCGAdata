.PARTICIPANT_METADATA_COLS <-
    c("sample_id", "sample_type", "participant", "tcga_sample_id")

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
#' @examples
#'
#' getClinical()
#'
getClinical <-
    function(tablename = "sample", participants = TRUE, columnName = NULL,
        metacols = .PARTICIPANT_METADATA_COLS)
{
    allclin <- getClinicalTable(tablename = tablename, metacols = metacols)
    if (is.null(columnName))
        columnName <- names(allclin)[1L]
    clinfiles <- unlist(unique(allclin[, columnName]))
    bfc <- BiocFileCache()
    rpath <- bfcquery(bfc, columnName)[["rpath"]]
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
        rownames(coldata) <- toupper(coldata[["patient.bcr_patient_barcode"]])
        if (participants) {
            parts <- avtable("participant")
            coldata <- merge(
                parts, coldata, by.x = "participant_id",
                by.y = "patient.bcr_patient_barcode"
            )
        }
    }
    coldata
}

#' Obtain the reference table for clinical data
#'
#' The column names in the output table can be used in the `getClinical`
#' function.
#'
#' @inheritParams getClinical
#'
#' @return A table with reference links to data resources on Terra
#'
#' @export
getClinicalTable <-
    function(tablename = "sample", metacols = .PARTICIPANT_METADATA_COLS)
{
    samples <- avtable(tablename)
    metadata <- samples[, metacols]
    samples <- samples[, !names(samples) %in% participant_meta]
    samples[, grep("clin", names(samples))]
}
