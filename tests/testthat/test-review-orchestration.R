test_that("detection subsampling restores the caller's random state", {
    reads <- withr::with_seed(91L, {
        inserts <- replicate(600L, paste0(
            sample(c("A", "C", "G", "T"), 24L, replace = TRUE), collapse = ""
        ))
        paste0("GGTAATACGACTCACTATAGGG", inserts, "CCATGCATGCATGCATGCAT")
    })
    pools <- list(round_01 = reads, round_00 = reads)
    withr::with_seed(903L, {
        before <- .Random.seed
        first <- selexprep_detect(pools, sampling_seed = 55L,
            max_reads_per_round = 500L)
        expect_identical(.Random.seed, before)
        expect_identical(first, selexprep_detect(pools, sampling_seed = 55L,
            max_reads_per_round = 500L))
        expect_identical(.Random.seed, before)
        expect_identical(first$sampling_seed, 55)
        expect_equal(first$read_fraction_used_for_inference, 5 / 6)
        failed <- selexprep_detect(pools, max_reads_per_round = 100L)
        expect_identical(failed$status, "UNABLE_TO_INFER")
        expect_identical(.Random.seed, before)
    })
})

test_that("QC preserves random state with custom rarefaction settings", {
    experiment <- selexprepR:::.as_selexprep_experiment(list(
        round_01 = selexprep_count(rep(c("AAAAAA", "CCCCCC", "GGGGGG"),
            c(12L, 8L, 4L))),
        round_00 = selexprep_count(rep(c("AAAAAA", "ACGTAC"), c(20L, 4L)))
    ))
    withr::with_seed(738L, {
        before <- .Random.seed
        first <- selexprep_qc(experiment, low_total_reads = 0L,
            rarefaction_depth = 8L, rarefaction_seed = 9L,
            kmer_size = 3L, kmer_top_sequences = 2L)
        expect_identical(.Random.seed, before)
        expect_identical(first, selexprep_qc(experiment, low_total_reads = 0L,
            rarefaction_depth = 8L, rarefaction_seed = 9L,
            kmer_size = 3L, kmer_top_sequences = 2L))
        expect_identical(.Random.seed, before)
        expect_identical(as.character(first$per_round$round),
            c("round_00", "round_01"))
        expect_identical(first$settings, list(
            rarefaction_depth = 8L, effective_rarefaction_depth = 8,
            rarefaction_seed = 9L, kmer_size = 3L, kmer_top_sequences = 2L
        ))
    })
})

test_that("QC reports empty rounds without inventing rarefaction depth", {
    experiment <- selexprepR:::.as_selexprep_experiment(list(
        round_10 = selexprep_count(character()),
        round_02 = selexprep_count(character()),
        round_00 = selexprep_count(character())
    ))
    qc <- selexprep_qc(experiment, low_total_reads = 0)
    expect_identical(as.character(qc$per_round$round),
        c("round_00", "round_02", "round_10"))
    expect_equal(qc$per_round$n_reads, c(0, 0, 0))
    expect_true(all(is.na(qc$per_round$rarefied_unique)))
    expect_true(is.na(qc$settings$effective_rarefaction_depth))
    expect_equal(qc$round_similarity$jaccard_distance, c(0, 0, 0))
    expect_identical(nrow(qc$flags), 0L)
})
