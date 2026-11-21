# Guards for a user running selexprepR on their own reads: each test is a
# place where a result could be wrong, or a call could stop, without a message
# that says what to do.

.own_p5 <- "ATCCAGAGTGACGCAGCATTGC"
.own_p3 <- "TGGACACGGTGGCTTAGTCAA"
.own_truseq <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

.own_inserts <- function(n, length = 30L, seed = 1L) {
    withr::with_seed(seed, vapply(seq_len(n), function(i) {
        paste(sample(c("A", "C", "G", "T"), length, replace = TRUE),
            collapse = "")
    }, character(1)))
}

.own_report <- function(mode = "BOTH_PRIMERS_SINGLE_READ", primer_5p = .own_p5,
    primer_3p = .own_p3) {
    reads <- paste0(.own_p5, .own_inserts(600L), .own_p3)
    fields <- unclass(suppressMessages(
        selexprep_detect(list(round_00 = reads))
    ))
    fields["primer_5p"] <- list(primer_5p)
    fields["primer_3p"] <- list(primer_3p)
    fields$extraction_mode <- mode
    fields$full_insert_recovered <- identical(mode, "BOTH_PRIMERS_SINGLE_READ")
    .new_library_report(fields)
}

test_that("a known adapter is found anywhere in a called primer", {
    truseq <- "AGATCGGAAGAGC"
    within <- .known_adapter_within
    expect_identical(within(paste0(truseq, "ACAC")), "TRUSEQ_R1")
    expect_identical(within(paste0("AA", truseq, "ACAC")), "TRUSEQ_R1")
    expect_identical(within(paste0(.own_p3, truseq)), "TRUSEQ_R1")
    expect_identical(within(paste0(reverse_complement(truseq), "GG")),
        "TRUSEQ_R1")
    expect_null(within(.own_p3))
    expect_null(within(NULL))
})

test_that("detection warns when a called primer contains an adapter", {
    report <- .own_report("THREE_PRIME_ONLY", primer_5p = NULL,
        primer_3p = "AAAGATCGGAAGAGCACAC")
    expect_warning(
        expect_message(
            .warn_about_inferred_library(report),
            "Only one constant was found.*THREE_PRIME_ONLY"
        ),
        "3' primer \\(AAAGATCGGAAGAGCACAC\\) contains the TRUSEQ_R1"
    )
    expect_silent(.warn_about_inferred_library(.own_report()))
})

test_that("detection warns on reads that run into the adapter", {
    pools <- lapply(c(round_00 = 1L, round_03 = 2L), function(seed) {
        paste0(.own_p5, .own_inserts(1500L, seed = seed), .own_p3, .own_truseq)
    })
    expect_warning(
        report <- suppressMessages(selexprep_detect(pools)),
        "sequencing adapter"
    )
    expect_false(identical(report$status, "UNABLE_TO_INFER"))
})

test_that("both primers with a one-sided mode warn that one is not trimmed", {
    reads <- list(round_00 = paste0(.own_p5, .own_inserts(60L), .own_p3,
        .own_truseq))
    report <- .own_report("FIVE_PRIME_ONLY", primer_3p = NULL)
    expect_warning(
        one_sided <- selexprep_extract(reads, report, primer_5p = .own_p5,
            primer_3p = .own_p3),
        "Both primers are set.*FIVE_PRIME_ONLY.*BOTH_PRIMERS_SINGLE_READ"
    )
    expect_true(all(Biostrings::width(one_sided$sequences_by_round$round_00) >
        30L))
    expect_no_warning(
        both <- selexprep_extract(reads, report, primer_5p = .own_p5,
            primer_3p = .own_p3, extraction_mode = "BOTH_PRIMERS_SINGLE_READ")
    )
    expect_identical(unique(Biostrings::width(both$sequences_by_round$round_00)),
        30L)
    expect_true(both$full_insert_recovered)
})

test_that("a reverse primer typed as ordered is named in the warning", {
    reads <- list(round_00 = paste0(.own_p5, .own_inserts(60L), .own_p3))
    as_ordered <- reverse_complement(.own_p3)
    warnings <- character()
    withCallingHandlers(
        selexprep_extract(reads, .own_report(), primer_3p = as_ordered),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    hit <- grep("supplied 3' primer", warnings, value = TRUE)
    expect_length(hit, 1L)
    expect_match(hit, sprintf("reverse complement, %s, is in 100.0%%", .own_p3),
        fixed = TRUE)
})

test_that("a supplied primer may carry IUPAC codes", {
    reads <- list(round_00 = paste0(.own_p5, .own_inserts(60L), .own_p3))
    degenerate <- "ATCCAGRGTGACGCAGCATYGC"
    expect_no_warning(
        extraction <- selexprep_extract(reads, .own_report(),
            primer_5p = degenerate, primer_3p = .own_p3)
    )
    expect_identical(unname(extraction$output_reads), 60)
    expect_error(
        selexprep_extract(reads, .own_report(), primer_5p = "T7-fwd"),
        "IUPAC ambiguity codes are allowed"
    )
})

test_that("an invalid extraction mode lists the valid ones", {
    reads <- list(round_00 = paste0(.own_p5, .own_inserts(60L), .own_p3))
    expect_error(
        selexprep_extract(reads, .own_report(), extraction_mode = "both"),
        paste("BOTH_PRIMERS_SINGLE_READ, FIVE_PRIME_ONLY, THREE_PRIME_ONLY,",
            "PAIRED_END_SPLIT_PRIMERS"),
        fixed = TRUE
    )
})

test_that("file inputs name an empty or a missing FASTQ", {
    full <- tempfile(fileext = ".fastq")
    empty <- tempfile(fileext = ".fastq")
    on.exit(unlink(c(full, empty)), add = TRUE)
    writeLines(c("@read_1", "ACGT", "+", "IIII"), full)
    file.create(empty)
    expect_warning(
        inputs <- read_selexprep_inputs(c(full, empty), round = c(0L, 1L)),
        "with no reads"
    )
    expect_identical(inputs$files$reads_loaded, c(1L, 0L))
    expect_error(
        read_selexprep_inputs(c(full, "no_such_file.fastq"), round = c(0L, 1L)),
        "not found: no_such_file.fastq",
        fixed = TRUE
    )
})
