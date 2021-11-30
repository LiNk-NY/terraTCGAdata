.PARTICIPANT_METADATA_COLS <-
    c("sample_id", "sample_type", "participant", "tcga_sample_id")


#' Get the list of clinical datasets
#'
#' @param tablename The terra data model table from which to extract the
#'     clinical data (default: "sample")
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
getClinical <- function(tablename = "sample", metacols = .PARTICIPANT_METADATA_COLS) {
    samples <- avtable(tablename)
    metadata <- samples[, metacols]
    samples <- samples[, !names(samples) %in% participant_meta]
    allclin <- samples[, grep("clin", names(samples))]
}

getParticipants <- function() {
    avtable("participant")
}

getSamples <- function() {
    avtable("sample")
}
