test_that("the bundled ENA inspection makes a complete offline dry-run plan", {
    local_mocked_bindings(
        .ena_perform_request = function(...) stop("Unexpected network access"),
        .package = "selexprepR"
    )
    path <- system.file("extdata", "ena-inspection.rds", package = "selexprepR")
    inspection <- readRDS(path)
    destination <- tempfile("selexprep-offline-plan-")
    result <- fetch_selexprep_reads(inspection, destination)
    expect_s3_class(inspection, "selexprep_inspection")
    expect_identical(inspection$accession, "PRJDB19098")
    expect_true(result$dry_run)
    expect_identical(sort(result$plan$round_number), 1:6)
    expect_true(all(result$plan$round_unambiguous))
    expect_true(all(nzchar(result$plan$expected_md5)))
    expect_length(result$downloaded_files, 0L)
    expect_false(dir.exists(destination))
})
