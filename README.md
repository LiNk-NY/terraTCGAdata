
# terraTCGAdata

## cBioPortal data and MultiAssayExperiment

### Overview

The `terraTCGAdata` R package aims to import TCGA datasets, as
[MultiAssayExperiment](http://bioconductor.org/packages/MultiAssayExperiment/),
that are available on the terra platform. The package provides a set of
functions that allow the discovery of relevant datasets. There is one
central function and two helper functions.

1.  The main function `terraTCGAdata` allows the creation of the
    `MultiAssayExperiment` object from the different indicated
    resources.
2.  The `getClinicalTable` and `getAssayTable` functions allow for the
    discovery of datasets within the terra data model. The column names
    from these tables can be provided as inputs to the `terraTCGAdata`
    function.

## Quick Start

### Installation

To install from Bioconductor (recommended for most users, this will
install the release or development version corresponding to your version
of Bioconductor):

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

To install from GitHub (for bleeding-edge, generally not necessary
because changes here are also pushed to
[bioc-devel](https://contributions.bioconductor.org/use-devel.html?q=devel#mid-october)):

``` r
if (!require("cBioPortalData", quietly = TRUE))
    BiocManager::install("waldronlab/terraTCGAdata")
```

Note. The above generally requires the installation of [Bioconductor
devel](https://contributions.bioconductor.org/use-devel.html?q=devel#mid-october)
and the most recent R version.

In order to make full use of this package, a local installation of
GCloud SDK is required. The common use case for this package lies within
the Terra environment. Please see the Bioconductor tutorials for running
RStudio on Terra.

To load the package:

``` r
library(AnVIL)
library(terraTCGAdata)
```

## gcloud sdk installation

A valid GCloud SDK installation is required to use the package. Use the
`gcloud_exists()` function from the `AnVIL` package to identify whether
it is installed in your system.

``` r
gcloud_exists()
#> [1] TRUE
```

You can also use the `gcloud_project` to set a project name by
specifying the project argument:

``` r
gcloud_project()
#> [1] "landmarkanvil2"
```

# Default Data Workspace

To get a list of available TCGA workspaces, use the
`findTCGAworkspaces()` function:

``` r
findTCGAworkspaces()
#> # A tibble: 36 × 5
#>   name                        lastModified createdBy namespace accessLevel
#>   <chr>                       <date>       <chr>     <chr>     <chr>      
#> 1 TCGA_ACC_OpenAccess_V1-0_D… 2017-01-31   birger@b… broad-fi… READER     
#> 2 TCGA_BLCA_OpenAccess_V1-0_… 2017-01-31   birger@b… broad-fi… READER     
#> 3 TCGA_BRCA_OpenAccess_V1-0_… 2017-01-31   birger@b… broad-fi… READER     
#> 4 TCGA_CESC_OpenAccess_V1-0_… 2017-01-31   birger@b… broad-fi… READER     
#> 5 TCGA_CHOL_OpenAccess_V1-0_… 2017-01-31   birger@b… broad-fi… READER     
#> # … with 31 more rows
```

You can then set a package-wide option with the `terraTCGAworkspace`
function and check the setting with the
`getOption('terraTCGAdata.workspace')` option.

``` r
terraTCGAworkspace()
#> [1] "TCGA_ACC_OpenAccess_V1-0_DATA"
getOption("terraTCGAdata.workspace")
#> [1] "TCGA_ACC_OpenAccess_V1-0_DATA"
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
#>   clin__bio__nationwid… clin__bio__nati… clin__bio__nati… clin__bio__intg…
#>   <chr>                 <chr>            <chr>            <chr>           
#> 1 gs://firecloud-tcga-… gs://firecloud-… gs://firecloud-… <NA>            
#> 2 gs://firecloud-tcga-… gs://firecloud-… gs://firecloud-… <NA>            
#> 3 gs://firecloud-tcga-… gs://firecloud-… gs://firecloud-… <NA>            
#> 4 gs://firecloud-tcga-… gs://firecloud-… gs://firecloud-… <NA>            
#> 5 gs://firecloud-tcga-… gs://firecloud-… gs://firecloud-… <NA>            
#> # … with 955 more rows, and 2 more variables:
#> #   clin__bio__intgen_org__Level_1__clinical__clin <chr>,
#> #   clin__bio__intgen_org__Level_1__biospecimen__clin <chr>
names(ct)
#> [1] "clin__bio__nationwidechildrens_org__Level_1__biospecimen__clin"
#> [2] "clin__bio__nationwidechildrens_org__Level_1__auxiliary__clin"  
#> [3] "clin__bio__nationwidechildrens_org__Level_1__clinical__clin"   
#> [4] "clin__bio__intgen_org__Level_1__auxiliary__clin"               
#> [5] "clin__bio__intgen_org__Level_1__clinical__clin"                
#> [6] "clin__bio__intgen_org__Level_1__biospecimen__clin"
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
#> ── Column specification ────────────────────────────────────────────────────────
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
#>   admin.bcr              admin.file_uuid admin.batch_num… admin.project_c…
#>   <chr>                  <chr>           <chr>            <chr>           
#> 1 nationwide children's… a93e6bbe-80de-… 385.38.0         tcga            
#> 2 nationwide children's… 8b055cbc-b2ff-… 385.38.0         tcga            
#> 3 nationwide children's… 61f5baab-8b35-… 422.33.0         tcga            
#> 4 nationwide children's… fbad35cb-8be3-… 422.33.0         tcga            
#> 5 nationwide children's… 5620a991-2a62-… 422.33.0         tcga            
#> # … with 455 more rows, and 2 more variables: admin.disease_code <chr>,
#> #   admin.day_of_dcc_upload <dbl>
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
#>   snp__genome_wide_snp… snp__genome_wid… snp__genome_wid… snp__genome_wid…
#>   <chr>                 <chr>            <chr>            <chr>           
#> 1 gs://firecloud-tcga-… gs://firecloud-… gs://firecloud-… gs://firecloud-…
#> 2 gs://firecloud-tcga-… gs://firecloud-… gs://firecloud-… gs://firecloud-…
#> 3 gs://firecloud-tcga-… gs://firecloud-… gs://firecloud-… gs://firecloud-…
#> 4 gs://firecloud-tcga-… gs://firecloud-… gs://firecloud-… gs://firecloud-…
#> 5 gs://firecloud-tcga-… gs://firecloud-… gs://firecloud-… gs://firecloud-…
#> # … with 955 more rows, and 25 more variables:
#> #   protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data <chr>,
#> #   rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data <chr>,
#> #   rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data <chr>,
#> #   rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__exon_quantification__data <chr>,
#> #   rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__junction_quantification__data <chr>,
#> #   rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data <chr>, …
names(at)
#>  [1] "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg"    
#>  [2] "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg18__seg"                       
#>  [3] "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg"                       
#>  [4] "snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg"    
#>  [5] "protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data"               
#>  [6] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data"                           
#>  [7] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data"                
#>  [8] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__exon_quantification__data"                  
#>  [9] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__junction_quantification__data"              
#> [10] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data"             
#> [11] "methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data"
#> [12] "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms__data"                        
#> [13] "cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg"                        
#> [14] "transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data" 
#> [15] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data"                
#> [16] "mirnaseq__illuminaga_mirnaseq__bcgsc_ca__Level_3__miR_isoform_expression__data"                 
#> [17] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data"                   
#> [18] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__junction_quantification__data"                 
#> [19] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes__data"                              
#> [20] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_isoforms__data"                           
#> [21] "rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__exon_quantification__data"                     
#> [22] "mirnaseq__illuminaga_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data"                    
#> [23] "methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data" 
#> [24] "rnaseq__illuminaga_rnaseq__unc_edu__Level_3__coverage__data"                                    
#> [25] "mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_isoform_expression__data"              
#> [26] "mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data"                 
#> [27] "rnaseq__illuminaga_rnaseq__unc_edu__Level_3__splice_junction_expression__data"                  
#> [28] "rnaseq__illuminaga_rnaseq__unc_edu__Level_3__gene_expression__data"                             
#> [29] "rnaseq__illuminaga_rnaseq__unc_edu__Level_3__exon_expression__data"
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
#>                     TCGA-3L-AA1B-01A-21-A45F-20
#> 14-3-3_beta-R-V                    -0.080527936
#> 14-3-3_epsilon-M-C                  0.055408025
#> 14-3-3_zeta-R-V                    -0.002073837
#> 4E-BP1-R-V                         -0.026154748
#> 4E-BP1_pS65-R-V                    -0.110226155
#> 4E-BP1_pT37_T46-R-V                -0.202870876
#>                     TCGA-4N-A93T-01A-21-A45F-20
#> 14-3-3_beta-R-V                     -0.15754027
#> 14-3-3_epsilon-M-C                   0.05978939
#> 14-3-3_zeta-R-V                     -0.13374613
#> 4E-BP1-R-V                          -0.35821838
#> 4E-BP1_pS65-R-V                     -0.15277484
#> 4E-BP1_pT37_T46-R-V                 -0.17585007
#>                     TCGA-4T-AA8H-01A-21-A45F-20
#> 14-3-3_beta-R-V                      -0.3840605
#> 14-3-3_epsilon-M-C                    0.1628557
#> 14-3-3_zeta-R-V                       0.2685011
#> 4E-BP1-R-V                            0.3263404
#> 4E-BP1_pS65-R-V                      -0.1381699
#> 4E-BP1_pT37_T46-R-V                  -0.1931612
#>                     TCGA-5M-AAT5-01A-11-A45F-20
#> 14-3-3_beta-R-V                     -0.08742583
#> 14-3-3_epsilon-M-C                  -0.15276783
#> 14-3-3_zeta-R-V                     -0.09958612
#> 4E-BP1-R-V                          -0.15502503
#> 4E-BP1_pS65-R-V                     -0.09373361
#> 4E-BP1_pT37_T46-R-V                  0.34677646
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
#> ── Column specification ────────────────────────────────────────────────────────
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
