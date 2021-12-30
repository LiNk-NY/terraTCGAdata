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
    samples <- samples[, !names(samples) %in% metacols]
    samples[, grep("clin", names(samples), invert = TRUE)]
}

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
#' @return Either a matrix or RaggedExperiment depending on the assay selected
#'
#' @examples
#'
#' getAssayData(
#'     assayName = "protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
#'     sampleCode = c("01", "10"),
#'     tablename = "sample",
#'     metacols = .PARTICIPANT_METADATA_COLS
#' )
#'
#' getAssayData(
#'     assayName = "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg",
#'     sampleCode = c("01", "10"),
#' )
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
    rpath <- bfcquery(bfc, assayName, exact = TRUE)[["rpath"]]
    if (!length(rpath)) {
        rpath <- bfcnew(bfc, rname = assayName, rtype = "local")
        dir.create(rpath)
        if (length(assayfiles) > 800)
            lapply(
                split(assayfiles,
                      cut(seq_along(assayfiles), 3, labels = letters[1:3])),
                gsutil_cp,
                destination = rpath
            )
        else
            gsutil_cp(assayfiles, rpath)
    }

    assaycopied <- list.files(rpath, pattern = "\\.txt", full.names = TRUE)
    suffix <-
        if (grepl("^cna|snp", assayName)) { ".seg.txt" } else { ".data.txt" }
    tcgaids <- gsub(suffix, "", basename(assaycopied), fixed = TRUE)
    if (!missing(sampleCode) && .is_character(sampleCode))
        tcgaids <- tcgaids[
            TCGAutils::TCGAsampleSelect(tcgaids, sampleCodes = sampleCode)
        ]
    assayselect <- file.path(rpath, paste0(tcgaids, suffix))
    datarows <- lapply(
        assayselect,
        readr::read_tsv,
        show_col_types = FALSE
    )
    if (grepl("seg", suffix)) {
        .mergeRangedSamples(datarows)
    } else {
        cnames <- unlist(datarows[[1]][1, ], use.names = FALSE)
        datarows <- lapply(datarows, function(x) x[-1, ])
        rowlist <- lapply(datarows, `[[`, 1L)
        dataonly <- lapply(datarows, function(x) {
            suppressMessages(readr::type_convert(x[2L]))
        })
        dups <- rowlist[!duplicated(rowlist)]
        logilist <- vector("list", length(dups))
        bindlist <- vector("list", length(dups))
        for (i in seq_along(dups)) {
            logilist[[i]] <-
                vapply(
                    rowlist,
                    function(rl) identical(rl, dups[[i]]),
                    logical(1L)
                )
            bindlist[[i]] <-
                cbind(rownames = dups[[i]],
                    dplyr::bind_cols(dataonly[logilist[[i]]]))
        }
        df <- Reduce(function(...) {
            merge(..., all = TRUE, by = "rownames")
        }, bindlist)
        rownames(df) <- df[["rownames"]]
        data.matrix(df[, -1])
    }
}

#' Import Terra TCGA data as an ExperimentList
#'
#' @inheritParams getAssayData
#'
#' @param assayNames character() A vector of assays selected from the colnames
#'     of `getAssayTable`.
#'
#' @return An `ExperimentList` of assays
#'
#' @examples
#'
#' ExperimentListData(
#'     assayNames = c("protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
#'     "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg"),
#'     sampleCode = c("01", "10")
#' )
#'
#' @export
ExperimentListData <- function(assayNames, sampleCode) {
    exps <- lapply(
        stats::setNames(nm = assayNames), getAssayData, sampleCode = sampleCode
    )
    ExperimentList(exps)
}

#' @examples
#'
#' terraTCGAdata(
#'     assays = c("protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
#'     "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data"),
#'     sampleCode = c("01", "10")
#' )
#'
#' @export
terraTCGAdata <- function(assays, sampleCode = NULL, split = TRUE) {
    el <- ExperimentListData(assayNames = assays, sampleCode = sampleCode)
    coldata <- getClinical()
    samap <- TCGAutils::generateMap(el, coldata, TCGAutils::TCGAbarcode)
    mae <- MultiAssayExperiment(
        experiments = el, colData = coldata, sampleMap = samap
    )
    if (split && (length(sampleCode) || !is.null(sampleCode))) {
        TCGAsplitAssays(mae, sampleCodes = sampleCode, exclusive = TRUE)
    } else {
        mae
    }
}
