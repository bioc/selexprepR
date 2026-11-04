.real_example_paths <- function() {
    paths <- system.file("extdata", sprintf("cs_round_%02d.fastq.gz",
        c(3L, 5L, 7L)), package = "selexprepR", mustWork = TRUE)
    stats::setNames(paths, c("round_03", "round_05", "round_07"))
}

.real_example_arguments <- function() {
    list(primer_5p = "CCAGAATGGAGAACGCAG", primer_3p = "GTCGTCCCGAGAGCC",
        extraction_mode = "BOTH_PRIMERS_SINGLE_READ")
}

.real_capture_warnings <- function(code) {
    observed <- character()
    result <- withCallingHandlers(code, warning = function(condition) {
        observed <<- c(observed, conditionMessage(condition))
        invokeRestart("muffleWarning")
    })
    list(result = result, warnings = observed)
}

test_that("bundled real FASTQ records match their source-content manifest", {
    paths <- .real_example_paths()
    metadata <- jsonlite::fromJSON(system.file("extdata",
        "cs-source-metadata.json", package = "selexprepR", mustWork = TRUE),
        simplifyVector = FALSE)
    expect_identical(metadata$license, "CC-BY-4.0")
    pools <- lapply(paths, read_selexprep_fastq)
    expect_identical(vapply(pools, length, integer(1)),
        c(round_03 = 600L, round_05 = 600L, round_07 = 600L))
    for (index in seq_along(paths)) {
        expect_identical(digest::digest(file = paths[[index]], algo = "sha256"),
            metadata$subsets[[index]]$sha256)
        connection <- gzfile(paths[[index]], "rb")
        content <- readBin(connection, "raw", n = 1000000L)
        close(connection)
        expect_identical(digest::digest(content, algo = "sha256",
            serialize = FALSE), metadata$subsets[[index]]$uncompressed_sha256)
    }
    expect_false(identical(as.character(pools[[1L]]), as.character(pools[[2L]])))
    expect_false(identical(as.character(pools[[2L]]), as.character(pools[[3L]])))
})

test_that("real FASTQ workflow records safe inference and honest retention", {
    paths <- .real_example_paths()
    inputs <- read_selexprep_inputs(paths)
    report <- selexprep_detect(inputs$sequences_by_round)
    expect_identical(report$status, "UNABLE_TO_INFER")
    expect_identical(report$failure_reason,
        "No candidate 5' or 3' flank found in the earliest round")
    expect_error(run_selexprep_files(paths), "does not permit extraction")
    captured <- .real_capture_warnings(do.call(run_selexprep_files,
        c(list(x = paths), .real_example_arguments())))
    experiment <- captured$result
    expect_length(captured$warnings, 2L)
    expect_true(any(grepl("3 round(s)", captured$warnings, fixed = TRUE)))
    expect_true(any(grepl("3 input FASTQ(s)", captured$warnings, fixed = TRUE)))
    expect_identical(as.numeric(Matrix::colSums(
        SummarizedExperiment::assay(experiment, "counts"))), c(228, 278, 297))
    expect_identical(dim(experiment), c(689L, 3L))
    expect_identical(S4Vectors::metadata(experiment)$extraction$input_reads,
        c(round_03 = 600, round_05 = 600, round_07 = 600))
    expect_true("low_total_reads" %in%
        S4Vectors::metadata(experiment)$qc$flags$name)
    report <- S4Vectors::metadata(experiment)$library_report
    expect_equal(report$match_rate_5p, 236 / 600)
    expect_equal(report$match_rate_3p, 228 / 600)
    expect_identical(report$n_length_mode, 39)
    expect_equal(report$n_length_confidence, 219 / 228)
    expect_identical(S4Vectors::metadata(experiment)$manifest$parameters[[
        "manual_primer_override"]], "true")
})

test_that("teaching barcodes recover real pools and their exact counts", {
    inputs <- read_selexprep_inputs(.real_example_paths())
    barcodes <- c(AAAAA = 3L, CCCCC = 5L, GGGGG = 7L)
    tagged <- Map(function(pool, barcode) paste0(barcode, as.character(pool)),
        inputs$sequences_by_round, names(barcodes))
    demultiplexed <- selexprep_demultiplex(
        c(unlist(tagged, use.names = FALSE), "TTTTTACGT"), barcodes)
    expected <- lapply(inputs$sequences_by_round,
        function(pool) unname(as.character(pool)))
    observed <- lapply(demultiplexed$sequences_by_round, as.character)
    expect_identical(observed, expected)
    expect_identical(length(demultiplexed$unassigned$r1), 1L)
    file_result <- .real_capture_warnings(do.call(run_selexprep_files,
        c(list(x = .real_example_paths()), .real_example_arguments())))$result
    barcode_result <- .real_capture_warnings(do.call(run_selexprep,
        c(list(sequences_by_round = demultiplexed), .real_example_arguments())))$result
    expect_identical(SummarizedExperiment::assay(barcode_result, "counts"),
        SummarizedExperiment::assay(file_result, "counts"))
})
