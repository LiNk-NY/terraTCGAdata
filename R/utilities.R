.is_character <- function(x, na.ok = FALSE, zchar = FALSE)
{
    is.character(x) &&
        (na.ok || all(!is.na(x))) &&
        (zchar || all(nzchar(x)))
}

.is_scalar_character <- function(x, na.ok = FALSE, zchar = FALSE)
    length(x) == 1L && .is_character(x, na.ok, zchar)


.is_scalar_logical <- function(x, na.ok = FALSE) {
    is.logical(x) && length(x) == 1L && (na.ok || !is.na(x))
}

.mergeRangedSamples <- function(sample_list) {
    sample_names <- vapply(
        sample_list, function(x) { unique(x[["Sample"]]) }, character(1L)
    )
    samples <- lapply(sample_list, function(x) x[, names(x) != "Sample"])
    names(samples) <- sample_names
    listgranges <- lapply(
        samples, GenomicRanges::makeGRangesFromDataFrame,
        keep.extra.columns = TRUE
    )
    grls <- as(listgranges, "GRangesList")
    RaggedExperiment::RaggedExperiment(grls)
}
