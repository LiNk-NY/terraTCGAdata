
# terraTCGAData

## Installation

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("terraTCGAdata")
```

# Overview

The `terraTCGAdata` R package aims to import TCGA datasets, as
[MultiAssayExperiment](http://bioconductor.org/packages/MultiAssayExperiment/),
available on the Terra platform. The package provides a set of functions
that allow the discovery of relevant datasets. It provides one main
function and two helper functions:

1.  `terraTCGAdata` allows the creation of the `MultiAssayExperiment`
    object from the different indicated resources.

2.  The `getClinicalTable` and `getAssayTable` functions allow for the
    discovery of datasets within the Terra data model. The column names
    from these tables can be provided as inputs to the `terraTCGAdata`
    function.

## Data

Some public Terra workspaces come pre-packaged with TCGA data (i.e.,
cloud data resources are linked within the data model). Particularly the
workspaces that are labelled `OpenAccess_V1-0`. Datasets harmonized to
the hg38 genome, such as those from the Genomic Data Commons data
repository, use a different data model / workflow and are not compatible
with the functions in this package. For those that are, we make use of
the Terra data model and represent the data as `MultiAssayExperiment`.

For more information on `MultiAssayExperiment`, please see the vignette
in that package.

# Requirements

## Loading packages

``` r
library(AnVIL)
library(terraTCGAdata)
```

## gcloud sdk installation

A valid GCloud SDK installation is required to use the package. To get
set up, see the Bioconductor tutorials for running RStudio on Terra. Use
the `gcloud_exists()` function from the
*[AnVIL](https://bioconductor.org/packages/3.16/AnVIL)* package to
identify whether it is installed in your system.

``` r
gcloud_exists()
#> [1] TRUE
```

You can also use the `gcloud_project` to set a project name by
specifying the project argument:

``` r
gcloud_project()
#> [1] "bioconductor-rpci-anvil"
```

# Default Data Workspace

To get a table of available TCGA workspaces, use the
`selectTCGAworkspace()` function:

``` r
selectTCGAworkspace()
#> [1] "TCGA_COAD_OpenAccess_V1-0_DATA"
```

You can also set the package-wide option with the `terraTCGAworkspace`
function and check the setting with
`getOption('terraTCGAdata.workspace')` or by running
`terraTCGAworkspace` function.

``` r
terraTCGAworkspace("TCGA_COAD_OpenAccess_V1-0_DATA")
#> [1] "TCGA_COAD_OpenAccess_V1-0_DATA"
getOption("terraTCGAdata.workspace")
#> [1] "TCGA_COAD_OpenAccess_V1-0_DATA"
```

# Clinical data resources

In order to determine what datasets to download, use the
`getClinicalTable` function to list all of the columns that correspond
to clinical data from the different collection centers.

``` r
ct <- getClinicalTable(workspace = "TCGA_COAD_OpenAccess_V1-0_DATA")
#> Using namespace/workspace: broad-firecloud-tcga/TCGA_COAD_OpenAccess_V1-0_DATA
ct
#> # A tibble: 960 × 6
#>    clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin                                                                                            clin_…¹ clin_…² clin_…³ clin_…⁴ clin_…⁵
#>    <chr>                                                                                                                                                     <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 gs://firecloud-tcga-open-access/tcga/dcc/coad/clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin/nationwidechildrens.org_COAD.bio.Level_1.38… gs://f… gs://f… <NA>    <NA>    <NA>   
#>  2 gs://firecloud-tcga-open-access/tcga/dcc/coad/clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin/nationwidechildrens.org_COAD.bio.Level_1.38… gs://f… gs://f… <NA>    <NA>    <NA>   
#>  3 gs://firecloud-tcga-open-access/tcga/dcc/coad/clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin/nationwidechildrens.org_COAD.bio.Level_1.38… gs://f… gs://f… <NA>    <NA>    <NA>   
#>  4 gs://firecloud-tcga-open-access/tcga/dcc/coad/clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin/nationwidechildrens.org_COAD.bio.Level_1.38… gs://f… gs://f… <NA>    <NA>    <NA>   
#>  5 gs://firecloud-tcga-open-access/tcga/dcc/coad/clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin/nationwidechildrens.org_COAD.bio.Level_1.42… gs://f… gs://f… <NA>    <NA>    <NA>   
#>  6 gs://firecloud-tcga-open-access/tcga/dcc/coad/clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin/nationwidechildrens.org_COAD.bio.Level_1.42… gs://f… gs://f… <NA>    <NA>    <NA>   
#>  7 gs://firecloud-tcga-open-access/tcga/dcc/coad/clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin/nationwidechildrens.org_COAD.bio.Level_1.42… gs://f… gs://f… <NA>    <NA>    <NA>   
#>  8 gs://firecloud-tcga-open-access/tcga/dcc/coad/clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin/nationwidechildrens.org_COAD.bio.Level_1.42… gs://f… gs://f… <NA>    <NA>    <NA>   
#>  9 gs://firecloud-tcga-open-access/tcga/dcc/coad/clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin/nationwidechildrens.org_COAD.bio.Level_1.42… gs://f… <NA>    <NA>    <NA>    <NA>   
#> 10 gs://firecloud-tcga-open-access/tcga/dcc/coad/clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin/nationwidechildrens.org_COAD.bio.Level_1.42… gs://f… <NA>    <NA>    <NA>    <NA>   
#> # … with 950 more rows, and abbreviated variable names ¹​clin__bio__nationwidechildrens_org__Level_1__auxiliary__clin, ²​clin__bio__nationwidechildrens_org__Level_1__clinical__clin,
#> #   ³​clin__bio__intgen_org__Level_1__auxiliary__clin, ⁴​clin__bio__intgen_org__Level_1__clinical__clin, ⁵​clin__bio__intgen_org__Level_1__biospecimen__clin
names(ct)
#> [1] "clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin" "clin__bio__nationwidechildrens_org__Level_1__auxiliary__clin"  
#> [3] "clin__bio__nationwidechildrens_org__Level_1__clinical__clin"    "clin__bio__intgen_org__Level_1__auxiliary__clin"               
#> [5] "clin__bio__intgen_org__Level_1__clinical__clin"                 "clin__bio__intgen_org__Level_1__biospecimen__clin"
```

# Clinical data download

After picking the column in the `getClinicalTable` output, use the
column name as input to the `getClinical` function to obtain the data:

``` r
column_name <- "clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin"
clin <- getClinical(
    columnName = column_name,
    participants = TRUE,
    workspace = "TCGA_COAD_OpenAccess_V1-0_DATA"
)
#> Using namespace/workspace: broad-firecloud-tcga/TCGA_COAD_OpenAccess_V1-0_DATA
#> 
#> ── Column specification ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> cols(
#>   .default = col_character(),
#>   admin.day_of_dcc_upload = col_double(),
#>   admin.month_of_dcc_upload = col_double(),
#>   admin.year_of_dcc_upload = col_double(),
#>   patient.additional_studies = col_logical(),
#>   patient.days_to_index = col_double(),
#>   patient.samples.sample.additional_studies = col_logical(),
#>   patient.samples.sample.biospecimen_sequence = col_logical(),
#>   patient.samples.sample.longest_dimension = col_double(),
#>   patient.samples.sample.intermediate_dimension = col_double(),
#>   patient.samples.sample.shortest_dimension = col_double(),
#>   patient.samples.sample.initial_weight = col_double(),
#>   patient.samples.sample.current_weight = col_logical(),
#>   patient.samples.sample.freezing_method = col_logical(),
#>   patient.samples.sample.oct_embedded = col_logical(),
#>   patient.samples.sample.preservation_method = col_logical(),
#>   patient.samples.sample.tissue_type = col_logical(),
#>   patient.samples.sample.composition = col_logical(),
#>   patient.samples.sample.tumor_descriptor = col_logical(),
#>   patient.samples.sample.days_to_collection = col_double(),
#>   patient.samples.sample.time_between_clamping_and_freezing = col_logical()
#>   # ... with 1225 more columns
#> )
#> ℹ Use `spec()` for the full column specifications.
clin[, 1:6]
#> # A tibble: 460 × 6
#>    admin.bcr                      admin.file_uuid                      admin.batch_number admin.project_code admin.disease_code admin.day_of_dcc_upload
#>    <chr>                          <chr>                                <chr>              <chr>              <chr>                                <dbl>
#>  1 nationwide children's hospital a93e6bbe-80de-41a1-9cc6-41fd0f56a4e9 385.38.0           tcga               coad                                     1
#>  2 nationwide children's hospital 8b055cbc-b2ff-4c62-a07c-ccfa44964937 385.38.0           tcga               coad                                     1
#>  3 nationwide children's hospital 61f5baab-8b35-45f4-a188-7d4f3d1a2a8b 422.33.0           tcga               coad                                     1
#>  4 nationwide children's hospital fbad35cb-8be3-4b36-a05d-e93aee1c3975 422.33.0           tcga               coad                                     1
#>  5 nationwide children's hospital 5620a991-2a62-446a-a26e-41ad5c1a92c7 422.33.0           tcga               coad                                     1
#>  6 nationwide children's hospital d7563bda-caea-473f-82fd-905c2bee66ea 422.33.0           tcga               coad                                     1
#>  7 nationwide children's hospital ef41a4ba-feb2-47c2-9292-a0a0680cf9f6 422.33.0           tcga               coad                                     1
#>  8 nationwide children's hospital 96b2bc07-30bf-4e67-b776-371a791249c0 422.33.0           tcga               coad                                     1
#>  9 nationwide children's hospital d1fedef8-53a4-42ff-9cf7-194fd92c004b 76.73.0            tcga               coad                                     1
#> 10 nationwide children's hospital 51c274ce-f952-45da-a0b3-285559d5c361 29.77.0            tcga               coad                                     1
#> # … with 450 more rows
dim(clin)
#> [1]  460 2376
```

# Assay data resources

We use the same approach for assay data. We first produce a list of
assays from the `getAssayTable` and then we select one along with any
sample codes of interest.

``` r
at <- getAssayTable(workspace = "TCGA_COAD_OpenAccess_V1-0_DATA")
at
#> # A tibble: 960 × 29
#>    snp__ge…¹ snp__…² snp__…³ snp__…⁴ rnase…⁵ rnase…⁶ rnase…⁷ prote…⁸ rnase…⁹ rnase…˟ methy…˟ rnase…˟ cna__…˟ trans…˟ rnase…˟ mirna…˟ rnase…˟ mirna…˟ rnase…˟ rnase…˟ rnase…˟ rnase…˟ methy…˟ rnase…˟
#>    <chr>     <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 gs://fir… gs://f… gs://f… gs://f… <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>   
#>  2 gs://fir… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>   
#>  3 gs://fir… gs://f… gs://f… gs://f… <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>   
#>  4 gs://fir… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>   
#>  5 gs://fir… gs://f… gs://f… gs://f… <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>   
#>  6 gs://fir… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>   
#>  7 gs://fir… gs://f… gs://f… gs://f… <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>   
#>  8 gs://fir… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… <NA>    gs://f… gs://f… gs://f… gs://f… <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>   
#>  9 gs://fir… gs://f… gs://f… gs://f… <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>   
#> 10 gs://fir… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… gs://f… <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>    <NA>   
#> # … with 950 more rows, 5 more variables: mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_isoform_expression__data <chr>,
#> #   mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data <chr>, rnaseq__illuminaga_rnaseq__unc_edu__Level_3__gene_expression__data <chr>,
#> #   rnaseq__illuminaga_rnaseq__unc_edu__Level_3__exon_expression__data <chr>, rnaseq__illuminaga_rnaseq__unc_edu__Level_3__splice_junction_expression__data <chr>, and abbreviated variable names
#> #   ¹​snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg, ²​snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg18__seg,
#> #   ³​snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg, ⁴​snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg,
#> #   ⁵​rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data, ⁶​rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data,
#> #   ⁷​rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data, ⁸​protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data, …
names(at)
#>  [1] "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg"    
#>  [2] "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg18__seg"                       
#>  [3] "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg"                       
#>  [4] "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg"    
#>  [5] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data"                           
#>  [6] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data"                
#>  [7] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data"             
#>  [8] "protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data"               
#>  [9] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__exon_quantification__data"                  
#> [10] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__junction_quantification__data"              
#> [11] "methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data"
#> [12] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms__data"                        
#> [13] "cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg"                        
#> [14] "transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data" 
#> [15] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data"                
#> [16] "mirnaseq__illuminaga_mirnaseq__bcgsc_ca__Level_3__miR_isoform_expression__data"                 
#> [17] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data"                   
#> [18] "mirnaseq__illuminaga_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data"                    
#> [19] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__junction_quantification__data"                 
#> [20] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes__data"                              
#> [21] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_isoforms__data"                           
#> [22] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__exon_quantification__data"                     
#> [23] "methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data" 
#> [24] "rnaseq__illuminaga_rnaseq__unc_edu__Level_3__coverage__data"                                    
#> [25] "mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_isoform_expression__data"              
#> [26] "mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data"                 
#> [27] "rnaseq__illuminaga_rnaseq__unc_edu__Level_3__gene_expression__data"                             
#> [28] "rnaseq__illuminaga_rnaseq__unc_edu__Level_3__exon_expression__data"                             
#> [29] "rnaseq__illuminaga_rnaseq__unc_edu__Level_3__splice_junction_expression__data"
```

# Summary of sample types in the data

You can get a summary table of all the samples in the adata by using the
`sampleTypesTable`:

``` r
sampleTypesTable(workspace = "TCGA_COAD_OpenAccess_V1-0_DATA")
#> Using namespace/workspace: broad-firecloud-tcga/TCGA_COAD_OpenAccess_V1-0_DATA
#> # A tibble: 5 × 4
#>   Code  Definition            Short.Letter.Code Frequency
#>   <chr> <chr>                 <chr>                 <dbl>
#> 1 10    Blood Derived Normal  NB                      406
#> 2 11    Solid Tissue Normal   NT                       92
#> 3 06    Metastatic            TM                        1
#> 4 01    Primary Solid Tumor   TP                      460
#> 5 02    Recurrent Solid Tumor TR                        1
```

# Intermediate function for obtaining only the data

Note that if you have the package-wide option set, the workspace
argument is not needed in the function call.

``` r
prot <- getAssayData(
    assayName = "protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
    sampleCode = c("01", "10"),
    workspace = "TCGA_COAD_OpenAccess_V1-0_DATA",
    sampleIdx = 1:4
)
head(prot)
#>                     TCGA-3L-AA1B-01A-21-A45F-20 TCGA-4N-A93T-01A-21-A45F-20 TCGA-4T-AA8H-01A-21-A45F-20 TCGA-5M-AAT5-01A-11-A45F-20
#> 14-3-3_beta-R-V                    -0.080527936                 -0.15754027                  -0.3840605                 -0.08742583
#> 14-3-3_epsilon-M-C                  0.055408025                  0.05978939                   0.1628557                 -0.15276783
#> 14-3-3_zeta-R-V                    -0.002073837                 -0.13374613                   0.2685011                 -0.09958612
#> 4E-BP1-R-V                         -0.026154748                 -0.35821838                   0.3263404                 -0.15502503
#> 4E-BP1_pS65-R-V                    -0.110226155                 -0.15277484                  -0.1381699                 -0.09373361
#> 4E-BP1_pT37_T46-R-V                -0.202870876                 -0.17585007                  -0.1931612                  0.34677646
```

# MultiAssayExperiment

Finally, once you have collected all the relevant column names, these
can be inputs to the main `terraTCGAdata` function:

``` r
mae <- terraTCGAdata(
    clinicalName = "clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin",
    assays =
        c("protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data",
        "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data"),
    sampleCode = NULL,
    split = FALSE,
    sampleIdx = 1:4,
    workspace = "TCGA_COAD_OpenAccess_V1-0_DATA"
)
#> Using namespace/workspace: broad-firecloud-tcga/TCGA_COAD_OpenAccess_V1-0_DATA
#> Using namespace/workspace: broad-firecloud-tcga/TCGA_COAD_OpenAccess_V1-0_DATA
#> Warning in .checkBarcodes(barcodes): Inconsistent barcode lengths: 27, 28
#> Using namespace/workspace: broad-firecloud-tcga/TCGA_COAD_OpenAccess_V1-0_DATA
#> 
#> ── Column specification ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> cols(
#>   .default = col_character(),
#>   admin.day_of_dcc_upload = col_double(),
#>   admin.month_of_dcc_upload = col_double(),
#>   admin.year_of_dcc_upload = col_double(),
#>   patient.additional_studies = col_logical(),
#>   patient.days_to_index = col_double(),
#>   patient.samples.sample.additional_studies = col_logical(),
#>   patient.samples.sample.biospecimen_sequence = col_logical(),
#>   patient.samples.sample.longest_dimension = col_double(),
#>   patient.samples.sample.intermediate_dimension = col_double(),
#>   patient.samples.sample.shortest_dimension = col_double(),
#>   patient.samples.sample.initial_weight = col_double(),
#>   patient.samples.sample.current_weight = col_logical(),
#>   patient.samples.sample.freezing_method = col_logical(),
#>   patient.samples.sample.oct_embedded = col_logical(),
#>   patient.samples.sample.preservation_method = col_logical(),
#>   patient.samples.sample.tissue_type = col_logical(),
#>   patient.samples.sample.composition = col_logical(),
#>   patient.samples.sample.tumor_descriptor = col_logical(),
#>   patient.samples.sample.days_to_collection = col_double(),
#>   patient.samples.sample.time_between_clamping_and_freezing = col_logical()
#>   # ... with 1225 more columns
#> )
#> ℹ Use `spec()` for the full column specifications.
#> Warning in .checkBarcodes(barcodes): Inconsistent barcode lengths: 27, 28
#> harmonizing input:
#>   removing 455 colData rownames not in sampleMap 'primary'
mae
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data: matrix with 200 rows and 4 columns
#>  [2] rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data: matrix with 20531 rows and 4 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```

We expect that most `OpenAccess_V1-0` cancer datasets follow this data
model. If you encounter any errors, please provide a minimally
reproducible example at <https://github.com/waldronlab/terraTCGAdata>.

# Session Info

``` r
sessionInfo()
#> R version 4.2.1 Patched (2022-09-26 r82921)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 22.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
#> LAPACK: /home/mr148/src/svn/r-4-2/R/lib/R/lib/libRlapack.so
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
#>  [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] terraTCGAdata_1.1.1         shiny_1.7.2                 testthat_3.1.4              MultiAssayExperiment_1.23.9 SummarizedExperiment_1.27.3 Biobase_2.57.1             
#>  [7] GenomicRanges_1.49.1        GenomeInfoDb_1.33.7         IRanges_2.31.2              S4Vectors_0.35.4            BiocGenerics_0.43.4         MatrixGenerics_1.9.1       
#> [13] matrixStats_0.62.0          AnVIL_1.9.9                 dplyr_1.0.10               
#> 
#> loaded via a namespace (and not attached):
#>  [1] rjson_0.2.21              ellipsis_0.3.2            futile.logger_1.4.3       XVector_0.37.1            rstudioapi_0.14           DT_0.25                   bit64_4.0.5              
#>  [8] AnnotationDbi_1.59.1      fansi_1.0.3               xml2_1.3.3                codetools_0.2-18          cachem_1.0.6              knitr_1.40                pkgload_1.3.0            
#> [15] jsonlite_1.8.0            Rsamtools_2.13.4          dbplyr_2.2.1              png_0.1-7                 BiocManager_1.30.18       readr_2.1.2               compiler_4.2.1           
#> [22] httr_1.4.4                assertthat_0.2.1          Matrix_1.5-1              fastmap_1.1.0             cli_3.4.1                 later_1.3.0               formatR_1.12             
#> [29] htmltools_0.5.3           prettyunits_1.1.1         tools_4.2.1               glue_1.6.2                GenomeInfoDbData_1.2.9    rappdirs_0.3.3            Rcpp_1.0.9               
#> [36] rapiclient_0.1.3          vctrs_0.4.2               Biostrings_2.65.6         rtracklayer_1.57.0        xfun_0.33                 stringr_1.4.1             brio_1.1.3               
#> [43] rvest_1.0.3               mime_0.12                 miniUI_0.1.1.1            lifecycle_1.0.2           restfulr_0.0.15           XML_3.99-0.10             zlibbioc_1.43.0          
#> [50] BiocStyle_2.25.0          vroom_1.5.7               hms_1.1.2                 promises_1.2.0.1          parallel_4.2.1            lambda.r_1.2.4            yaml_2.3.5               
#> [57] curl_4.3.2                memoise_2.0.1             biomaRt_2.53.2            stringi_1.7.8             RSQLite_2.2.17            BiocIO_1.7.1              GenomicDataCommons_1.21.4
#> [64] GenomicFeatures_1.49.7    filelock_1.0.2            BiocParallel_1.31.12      rlang_1.0.6               pkgconfig_2.0.3           bitops_1.0-7              evaluate_0.16            
#> [71] lattice_0.20-45           purrr_0.3.4               GenomicAlignments_1.33.1  htmlwidgets_1.5.4         bit_4.0.4                 tidyselect_1.1.2          magrittr_2.0.3           
#> [78] R6_2.5.1                  generics_0.1.3            DelayedArray_0.23.2       DBI_1.1.3                 pillar_1.8.1              KEGGREST_1.37.3           RCurl_1.98-1.8           
#> [85] tibble_3.1.8              crayon_1.5.1              futile.options_1.0.1      utf8_1.2.2                BiocFileCache_2.5.0       tzdb_0.3.0                rmarkdown_2.16           
#> [92] progress_1.2.2            grid_4.2.1                blob_1.2.3                digest_0.6.29             xtable_1.8-4              tidyr_1.2.1               httpuv_1.6.6             
#> [99] TCGAutils_1.17.3
```
