.cutadapt_reference <- function(name) {
    utils::read.delim(gzfile(test_path("fixtures", name)),
        colClasses = "character", na.strings = ".", quote = "")
}

.trimmed_or_na <- function(reads, mode, primer_5p, primer_3p) {
    kept <- selexprepR:::.trim_round(reads, mode, primer_5p, primer_3p,
        keep_mask = TRUE)
    trimmed <- rep(NA_character_, length(reads))
    trimmed[attr(kept, "kept")] <- unname(kept)
    trimmed
}

test_that("primer trimming matches cutadapt 5.2 read by read", {
    reference <- .cutadapt_reference("cutadapt-5.2-synthetic-reference.tsv.gz")
    reads <- reference$read
    primer_5p <- "CCAGAATGGAGAACGCAG"
    primer_3p <- "GTCGTCCCGAGAGCC"
    expect_identical(.trimmed_or_na(reads, "BOTH_PRIMERS_SINGLE_READ",
        primer_5p, primer_3p), reference$both)
    expect_identical(.trimmed_or_na(reads, "FIVE_PRIME_ONLY", primer_5p, NULL),
        reference$five_prime_only)
    expect_identical(.trimmed_or_na(reads, "THREE_PRIME_ONLY", NULL, primer_3p),
        reference$three_prime_only)
    # Degenerate primers: IUPAC codes in the primer are wildcards.
    expect_identical(.trimmed_or_na(reads, "BOTH_PRIMERS_SINGLE_READ",
        "CCAGAATGGNGAACRCAG", "GTCGTCCCGAGAGCY"), reference$both_iupac)
})

test_that("the aligner alone reproduces cutadapt without the shortcut", {
    reference <- .cutadapt_reference("cutadapt-5.2-synthetic-reference.tsv.gz")
    reads <- reference$read
    located_5p <- selexprepR:::.cutadapt_locate(reads, "CCAGAATGGAGAACGCAG",
        front = TRUE)
    trimmed_5p <- ifelse(is.na(located_5p[, "stop"]), NA_character_,
        substring(reads, located_5p[, "stop"] + 1L))
    expect_identical(trimmed_5p, reference$five_prime_only)
    located_3p <- selexprepR:::.cutadapt_locate(reads, "GTCGTCCCGAGAGCC",
        front = FALSE)
    trimmed_3p <- ifelse(is.na(located_3p[, "start"]), NA_character_,
        substr(reads, 1L, located_3p[, "start"]))
    expect_identical(trimmed_3p, reference$three_prime_only)
})

test_that("bundled real reads are trimmed exactly as cutadapt 5.2 trims them", {
    # Insert coordinates cutadapt reported for every read, so the fixture
    # holds no sequence from the CC BY 4.0 subsets.
    reference <- .cutadapt_reference("cutadapt-5.2-real-coordinates.tsv.gz")
    for (round in unique(reference$round)) {
        path <- system.file("extdata", sprintf("cs_%s.fastq.gz", round),
            package = "selexprepR", mustWork = TRUE)
        reads <- unname(as.character(read_selexprep_fastq(path)))
        rows <- reference[reference$round == round, ]
        expect_identical(as.integer(rows$index), seq_along(reads))
        expected <- ifelse(is.na(rows$start), NA_character_,
            substr(reads, as.integer(rows$start), as.integer(rows$end)))
        expect_identical(.trimmed_or_na(reads, "BOTH_PRIMERS_SINGLE_READ",
            "CCAGAATGGAGAACGCAG", "GTCGTCCCGAGAGCC"), expected, info = round)
    }
})

test_that("cutadapt choices: partial primers, first occurrence, substitution", {
    primer <- "CCAGAATGGAGAACGCAG"
    trim_5p <- function(reads) {
        .trimmed_or_na(reads, "FIVE_PRIME_ONLY", primer, NULL)
    }
    insert <- "ACGTTGCAACGTTGCAACGT"
    # A 5' primer cut short at the read start, down to 3 bases.
    expect_identical(trim_5p(paste0(substring(primer, 9L), insert)), insert)
    expect_identical(trim_5p(paste0("CAG", insert)), insert)
    # A full later copy replaces a partial one; a one-error copy at the read
    # start is kept over a later exact copy that does not overlap it.
    expect_identical(trim_5p(paste0("GAACGCAGTTTT", primer, insert)), insert)
    one_error <- paste0(substr(primer, 1L, 9L), "T", substring(primer, 11L))
    expect_identical(trim_5p(paste0(one_error, "TT", primer, insert)),
        paste0("TT", primer, insert))
    # An error on the last primer base is a substitution, not a deletion.
    expect_identical(trim_5p(paste0(substr(primer, 1L, 17L), "T", insert)),
        insert)
    # A 3' primer cut short at the read end, with at least 3 bases.
    trim_3p <- function(reads) {
        .trimmed_or_na(reads, "THREE_PRIME_ONLY", NULL, "GTCGTCCCGAGAGCC")
    }
    expect_identical(trim_3p(paste0(insert, "GTC")), insert)
    expect_identical(trim_3p(paste0(insert, "GT")), NA_character_)
})
