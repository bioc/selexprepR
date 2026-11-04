.review_extraction_report <- function(mode = 'BOTH_PRIMERS_SINGLE_READ',
    orientation = 'FORWARD') {
    fields <- unclass(read_library_report(test_path('fixtures',
        'library-report-both-primers.json')))
    fields$primer_5p <- 'ACGTCAGTACGATCGTACGA'
    fields$primer_3p <- 'TGCATCGACTGACGTAGCTA'
    fields$extraction_mode <- mode
    fields$orientation <- orientation
    fields$full_insert_recovered <- identical(mode, 'BOTH_PRIMERS_SINGLE_READ')
    fields$required_action <- if (identical(mode, 'PAIRED_END_SPLIT_PRIMERS')) {
        'READ_MERGING_RECOMMENDED'
    } else {
        'NONE'
    }
    selexprepR:::.new_library_report(fields)
}

test_that('reverse complement keeps names and case, complements IUPAC codes', {
    expect_identical(reverse_complement(c(rna = 'aCGUn', ambiguity = 'Ryn?-')),
        c(rna = 'nACGt', ambiguity = '-?nrY'))
    expect_identical(reverse_complement('RYKMSWBVDHN'), 'NDHBVWSKMRY')
    codes <- 'ACGTRYKMSWBVDHNacgtrykmswbvdhn'
    expect_identical(reverse_complement(reverse_complement(codes)), codes)
    expect_identical(reverse_complement(character()), character())
    expect_identical(reverse_complement(c('', 'ACGT')), c('', 'ACGT'))
    expect_identical(reverse_complement('\u00e9Ac'), 'gT\u00e9')
})

test_that('adapter ranking chooses exact matches and the earliest exact tie', {
    report <- .review_extraction_report()
    p5 <- report$primer_5p
    p3 <- report$primer_3p
    mutated <- paste0('T', substring(p5, 2L))
    reads <- c(
        paste0(mutated, 'AAAA', p5, 'cGNu', p3),
        paste0(p5, 'ACGT', p5, 'TGCA', p3),
        paste0(p5, 'GATT', p3, 'CCCC', p3)
    )
    extracted <- selexprep_extract(list(round_00 = reads), report)
    expect_identical(as.character(extracted$sequences_by_round$round_00),
        c('cGNu', paste0('ACGT', p5, 'TGCA'), 'GATT'))
})

test_that('batch extraction keeps input order and exact paired retention masks', {
    report <- .review_extraction_report('PAIRED_END_SPLIT_PRIMERS')
    p5 <- report$primer_5p
    r2_primer <- reverse_complement(report$primer_3p)
    r1 <- c(paste0(p5, 'AAAA'), paste0(p5, 'CCCC'),
        'NNNN', paste0(p5, 'TTTT'))
    r2 <- c(paste0(r2_primer, 'GGGG'), 'NNNN',
        paste0(r2_primer, 'CCCC'), paste0(r2_primer, 'AAAA'))
    extracted <- selexprepR:::.extract_selexprep(list(round_00 = r1), report,
        list(round_00 = r2), keep_masks = TRUE)
    expect_identical(as.character(extracted$sequences_by_round$round_00$r1),
        c('AAAA', 'TTTT'))
    expect_identical(as.character(extracted$sequences_by_round$round_00$r2),
        c('GGGG', 'AAAA'))
    expect_identical(attr(extracted, 'kept_by_round')$round_00,
        c(TRUE, FALSE, FALSE, TRUE))
    expect_identical(extracted$output_reads, c(round_00 = 2))
})

test_that('reverse-oriented reads use the same ambiguity policy as counting', {
    report <- .review_extraction_report(orientation = 'REVERSE')
    forward <- paste0(report$primer_5p, 'rY-N?', report$primer_3p)
    extracted <- selexprep_extract(
        list(round_00 = reverse_complement(forward)), report)
    expect_identical(as.character(extracted$sequences_by_round$round_00), 'rY-N?')
})

test_that('strand assignment requires an anchored prefix within one mismatch', {
    observed <- selexprepR:::.strand_distribution(
        c('CGTC', 'ACGT', 'ACGG', 'AC', ''), 'ACGT', NULL)
    expect_identical(observed, c(forward = 2L, reverse = 0L, ambiguous = 3L))
})

test_that('paired extraction accepts an empty round without inventing a mate', {
    report <- .review_extraction_report('PAIRED_END_SPLIT_PRIMERS')
    extraction <- selexprep_extract(list(round_00 = character()), report,
        paired_mate_streams = list(round_00 = character()))
    expect_length(extraction$sequences_by_round$round_00$r1, 0L)
    expect_length(extraction$sequences_by_round$round_00$r2, 0L)
    expect_identical(extraction$output_reads, c(round_00 = 0))
})
