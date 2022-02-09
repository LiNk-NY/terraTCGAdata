.DEFAULT_TABLENAME <- "sample"

#' Obtain a reference table for assay data in the Terra data model
#'
#' The column names in the output can be used in the `getAssayData` function.
#'
#' @inheritParams getAssayData
#'
#' @md
#'
#' @return A tibble of pointers to resources within the Terra data model
#'
#' @export
getAssayTable <-
    function(
        tablename = .DEFAULT_TABLENAME, metacols = .PARTICIPANT_METADATA_COLS,
        workspace = terraWorkspace(), namespace = .DEFAULT_NAMESPACE
    )
{
    samples <- avtable(
        table = tablename, namespace = namespace, name = workspace
    )
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
#' @seealso [getAssayTable()]
#'
#' @md
#'
#' @examples
#'
#' getAssayData(
#'     assayName = "protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
#'     sampleCode = c("01", "10")
#' )
#'
#' getAssayData(
#'     assayName = "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg",
#'     sampleCode = c("01", "10")
#' )
#'
#' getAssayData(
#'     assayName = "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg"
#' )
#'
#' @export
getAssayData <-
    function(assayName, sampleCode = "01", tablename = .DEFAULT_TABLENAME,
        workspace = terraWorkspace(), namespace = .DEFAULT_NAMESPACE,
        metacols = .PARTICIPANT_METADATA_COLS
    )
{
    if (missing(assayName))
        stop("Select an assay name from 'getAssayTable'")

    assayTable <- getAssayTable(
        tablename = tablename, metacols = metacols,
        workspace = workspace, namespace = namespace
    )
    assayfiles <- unlist(unique(stats::na.omit(assayTable[, assayName])))

    bfc <- BiocFileCache::BiocFileCache()
    rpath <- BiocFileCache::bfcquery(bfc, assayName, exact = TRUE)[["rpath"]]
    if (!length(rpath)) {
        rpath <- BiocFileCache::bfcnew(bfc, rname = assayName, rtype = "local")
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
    if (!missing(sampleCode) && .is_character(sampleCode) && !is.null(sampleCode))
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
#' @md
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
ExperimentListData <-
    function(
        assayNames, sampleCode, workspace = terraWorkspace(),
        namespace = .DEFAULT_NAMESPACE, tablename = .DEFAULT_TABLENAME,
        verbose = TRUE
    )
{
    if (verbose)
        message(
            "Using namespaece/workspace: ", paste0(namespace, "/", workspace)
        )
    elist <- structure(vector("list", length(assayNames)), .Names = assayNames)
    for (assay in assayNames) {
        elist[[assay]] <-
            getAssayData(
                assay, sampleCode = sampleCode, tablename = tablename,
                workspace = workspace, namespace = namespace
            )
    }
    ExperimentList(elist)
}

#' Obtain a MultiAssayExperiment from the Terra workspace
#'
#' Workspaces on Terra come pre-loaded with TCGA Data. The examples in the
#' documentation correspond to the TCGA_COAD_OpenAccess_V1 workspace that
#' can be found on \url{app.terra.bio}.
#'
#' @inheritParams getClinical
#'
#' @param clinicalName character(1) The column name taken from
#'     `getClinicalTable()` and downloaded to be included as the `colData`.
#'
#' @param assays character() A character vector of assay names taken from
#'     `getAssayTable()`
#'
#' @param sampleCode character() A character vector of sample codes from
#'     `sampleTypesTable()`. By default, (NULL) all samples are downloaded and
#'     kept in the data.
#'
#' @param split logical(1L) Whether or not to split the `MultiAssayExperiment`
#'     by sample types using `splitAssays` helper function (default FALSE).
#'
#' @return A `MultiAssayExperiment` object with n number of assays corresponding
#'     to the `assays` argument.
#'
#' @md
#'
#' @examples
#'
#' terraTCGAdata(
#'     clinicalName = "clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin",
#'     assays = c("protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
#'     "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data"),
#'     sampleCode = NULL,
#'     split = FALSE
#' )
#'
#' @export
terraTCGAdata <-
    function(
        clinicalName, assays, participants = TRUE,
        sampleCode = NULL, split = FALSE,
        workspace = terraWorkspace(), namespace = .DEFAULT_NAMESPACE,
        tablename = .DEFAULT_TABLENAME, verbose = TRUE
    )
{
    if (verbose)
        message(
            "Using namespaece/workspace: ", paste0(namespace, "/", workspace)
        )
    el <- ExperimentListData(
        assayNames = assays, sampleCode = sampleCode, workspace = workspace,
        namespace = namespace, tablename = tablename, verbose = FALSE
    )
    coldata <- getClinical(
        columnName = clinicalName, workspace = workspace,
        namespace = namespace, verbose = FALSE
    )
    coldata <- .transform_clinical_to_coldata(
        clinical_data = coldata,
        columnName = clinicalName,
        participants = participants
    )
    samap <- TCGAutils::generateMap(el, coldata, TCGAutils::TCGAbarcode)
    mae <- MultiAssayExperiment(
        experiments = el, colData = coldata, sampleMap = samap
    )
    if (split && (length(sampleCode) || !is.null(sampleCode))) {
        TCGAutils::TCGAsplitAssays(
            mae, sampleCodes = sampleCode, exclusive = TRUE
        )
    } else {
        mae
    }
}

.transform_clinical_to_coldata <-
    function(clinical_data, columnName, participants)
{
    coldata <- as.data.frame(clinical_data)
    if (!is.null(coldata[["patient.bcr_patient_barcode"]])) {
        rownames(coldata) <- coldata[["participant_id"]] <-
            toupper(coldata[["patient.bcr_patient_barcode"]])
        if (participants) {
            parts <- avtable("participant")
            parts <- methods::as(parts, "DataFrame")
            rownames(parts) <- parts[["participant_id"]]
            ## replace munged COAD with TCGA
            parts[["participant_id"]] <-
                gsub("^[A-Z]{4}", "TCGA", parts[["participant_id"]])
            coldata <- merge(
                parts, coldata, by = "participant_id",
            )
            rownames(coldata) <- coldata[["participant_id"]]
        }
        coldata <- methods::as(coldata, "DataFrame")
    }
    S4Vectors::metadata(coldata)[["columnName"]] <- columnName
    coldata
}
