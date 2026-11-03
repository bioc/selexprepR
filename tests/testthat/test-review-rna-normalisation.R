test_that('lowercase RNA barcodes route paired reads without changing inserts', {
    r1 <- c('uuuuuacGuN', 'AaAaAgCuA', 'cccccTGca')
    r2 <- c('UgCa', 'aCgU', 'nNuU')
    result <- selexprep_demultiplex(r1, c(uuuuu = 3L, aaaaa = 7L),
        paired_mates = r2, max_mismatches = 0L)
    expect_identical(result$barcodes, c(TTTTT = 3L, AAAAA = 7L))
    expect_identical(as.character(result$sequences_by_round$round_03), 'acGuN')
    expect_identical(as.character(result$sequences_by_round$round_07), 'gCuA')
    expect_identical(as.character(result$paired_mate_streams$round_03), 'UgCa')
    expect_identical(as.character(result$paired_mate_streams$round_07), 'aCgU')
    expect_identical(as.character(result$unassigned$r1), 'cccccTGca')
    expect_identical(as.character(result$unassigned$r2), 'nNuU')
    expect_identical(result$summary$assigned_reads, c(1L, 1L))
})

test_that('DNA and RNA barcode case produce equivalent assignments', {
    r2 <- c('UgCa', 'aCgU', 'nNuU')
    expected <- selexprep_demultiplex(
        c('TTTTTacGuN', 'AAAAAgCuA', 'cccccTGca'),
        c(TTTTT = 3L, AAAAA = 7L), paired_mates = r2, max_mismatches = 0L)
    for (barcode in c('TTTTT', 'ttttt', 'UUUUU', 'uuuuu', 'uUtTu')) {
        observed <- selexprep_demultiplex(
            c(paste0(barcode, 'acGuN'), 'AaAaAgCuA', 'cccccTGca'),
            stats::setNames(c(3L, 7L), c(barcode, 'aaaaa')),
            paired_mates = r2, max_mismatches = 0L)
        expect_identical(observed, expected, info = barcode)
    }
})

test_that('matching normalisation does not alter an untrimmed RNA read', {
    result <- selexprep_demultiplex('uuuuuacGu', c(TTTTT = 0L),
        max_mismatches = 0L, trim_barcode = FALSE)
    expect_identical(as.character(result$sequences_by_round$round_00),
        'uuuuuacGu')
})

test_that('RNA and DNA aliases cannot define conflicting barcode assignments', {
    expect_error(selexprep_demultiplex('uuuuuACGT', c(uuuuu = 0L, TTTTT = 1L),
        max_mismatches = 0L), 'at least 1')
})
