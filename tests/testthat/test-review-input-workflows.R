.review_download_inspection <- function(md5) {
    structure(list(
        accession = "PRJNA000001", library_strategy = "SELEX",
        runs = S4Vectors::DataFrame(
            run_accession = "SRR000001", sample_accession = "SRS000001",
            sample_title = "Round 1",
            fastq_urls = I(list(c("https://example.org/SRR000001_1.fastq",
                "https://example.org/SRR000001_2.fastq"))),
            fastq_md5 = I(list(rep(md5, 2L))),
            fastq_bytes = I(list(c(18, 18)))
        )
    ), class = c("selexprep_inspection", "list"))
}

test_that("ENA downloads verify content before exposing completed FASTQ paths", {
    payload <- charToRaw("@read\nACGT\n+\nIIII\n")
    checksum <- digest::digest(payload, algo = "md5", serialize = FALSE)
    inspection <- .review_download_inspection(toupper(checksum))
    requests <- character()
    timeouts <- numeric()
    testthat::local_mocked_bindings(req_perform = function(request, path) {
        requests <<- c(requests, request$url)
        timeouts <<- c(timeouts, request$options$timeout_ms)
        writeBin(payload, path)
        httr2::response(status_code = 200L)
    }, .package = "httr2")
    destination <- withr::local_tempdir()
    result <- fetch_selexprep_reads(inspection, destination, dry_run = FALSE,
        timeout_seconds = 2)
    expect_false(result$dry_run)
    expect_identical(result$downloaded_files, normalizePath(file.path(destination,
        c("SRR000001_1.fastq", "SRR000001_2.fastq")), winslash = "/"))
    expect_identical(requests, inspection$runs$fastq_urls[[1L]])
    expect_identical(timeouts, c(2000, 2000))
    expect_identical(unname(vapply(result$downloaded_files, digest::digest,
        character(1), algo = "md5", file = TRUE)), rep(checksum, 2L))
    expect_error(fetch_selexprep_reads(inspection, destination, dry_run = FALSE),
        "Refusing to overwrite")
    expect_length(requests, 2L)
    expect_identical(length(list.files(destination, all.files = TRUE,
        pattern = "^\\.selexprep-download-")), 0L)
})

test_that("failed ENA transfers retain diagnostics without publishing a FASTQ", {
    payload <- charToRaw("@read\nACGT\n+\nIIII\n")
    inspection <- .review_download_inspection(strrep("0", 32L))
    status <- 200L
    testthat::local_mocked_bindings(req_perform = function(request, path) {
        writeBin(payload, path)
        httr2::response(status_code = status)
    }, .package = "httr2")
    for (status in c(200L, 500L)) {
        destination <- withr::local_tempdir()
        failure <- tryCatch(fetch_selexprep_reads(inspection, destination,
            dry_run = FALSE), error = identity)
        expect_s3_class(failure, "error")
        expect_match(conditionMessage(failure), if (status == 200L) {
            "MD5 verification failed.*temporary file retained"
        } else {
            "Download failed.*HTTP status 500"
        })
        expect_false(file.exists(file.path(destination, "SRR000001_1.fastq")))
        retained <- list.files(destination, pattern = "^\\.selexprep-download-",
            all.files = TRUE, full.names = TRUE)
        expect_length(retained, 1L)
        expect_identical(readBin(retained[[1L]], "raw", n = length(payload)),
            payload)
    }
})

test_that("FASTQ pooling preserves lane order while sorting file provenance", {
    destination <- withr::local_tempdir()
    sequences <- list("CCCC", c("TTTT", "GGGG"), c("AAAA", "TTTT"),
        "GGGG", c("ACGT", "TGCA"))
    paths <- file.path(destination, paste0("file_", seq_along(sequences), ".fastq"))
    for (index in seq_along(paths)) {
        reads <- sequences[[index]]
        writeLines(as.vector(rbind(paste0("@read_", seq_along(reads)), reads,
            "+", strrep("I", nchar(reads)))), paths[[index]])
    }
    table <- data.frame(path = paths, round_number = c(0, 2, 0, 0, 2),
        mate = c("R1", "R2", "R2", "R1", "R1"))
    result <- read_selexprep_inputs(table)
    expect_identical(result$read_source, "R1_AND_R2")
    expect_identical(names(result$sequences_by_round), c("round_00", "round_02"))
    expect_identical(unname(as.character(result$sequences_by_round$round_00)),
        c("CCCC", "GGGG"))
    expect_identical(unname(as.character(result$paired_mate_streams$round_00)),
        c("AAAA", "TTTT"))
    expect_identical(unname(as.character(result$sequences_by_round$round_02)),
        c("ACGT", "TGCA"))
    expect_identical(as.character(result$files$path),
        normalizePath(paths[c(1, 4, 3, 5, 2)], winslash = "/"))
    expect_identical(result$files$reads_loaded, c(1L, 1L, 2L, 2L, 2L))
})
