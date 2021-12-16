getAssayTable <-
    function(tablename = "sample", metacols = .PARTICIPANT_METADATA_COLS)
{
    samples <- avtable(tablename)
    metadata <- samples[, metacols]
    samples <- samples[, !names(samples) %in% participant_meta]
    samples[, grep("clin", names(samples), invert = TRUE)]
}

# assayName <- "protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data"

getAssayData <-
    function(tablename = "sample", metacols = .PARTICIPANT_METADATA_COLS,
        assayName, sampleCode = "01")
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
