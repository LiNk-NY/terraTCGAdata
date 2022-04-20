test_that("terraTCGAworkspace works", {
    withr::with_options(
        list("terraTCGAdata.workspace" = ""),
        expect_warning(
            terraTCGAworkspace()
        )
    )
    
    withr::with_options(
        list("terraTCGAdata.workspace" = NULL), 
        expect_warning(
            terraTCGAworkspace()
        )
    )
    
    withr::with_options(
        list("terraTCGAdata.workspace" = "TCGA_ACC_OpenAccess_V1-0_DATA"), 
        expect_identical(
            terraTCGAworkspace(),
            "TCGA_ACC_OpenAccess_V1-0_DATA"
        )
    )
    
    ## No OP when option is already set
    withr::with_options(
        list("terraTCGAdata.workspace" = "TCGA_ACC_OpenAccess_V1-0_DATA"), 
        expect_identical(
            terraTCGAworkspace("TCGA_TGCT_OpenAccess_V1-0_DATA"),
            "TCGA_ACC_OpenAccess_V1-0_DATA"
        )
    )
    
    skip_if_not(AnVIL::gcloud_exists()) 
    expect_error(
        terraTCGAworkspace("ABDC")
    )
})
