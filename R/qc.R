.shannon_entropy <- function(reads) {
    total <- sum(reads)
    if (!total) {
        return(0)
    }
    probabilities <- reads[reads > 0]/total
    -sum(probabilities * log2(probabilities))
}
.top_n_coverage <- function(reads, n = 100L) {
    total <- sum(reads)
    if (n <= 0L || !total) {
        return(0)
    }
    sum(utils::head(sort(reads, decreasing = TRUE), n))/total
}
.singleton_fraction <- function(reads) {
    if (!length(reads)) {
        return(0)
    }
    sum(reads == 1)/length(reads)
}
.weighted_gc_content <- function(sequences, reads) {
    if (!length(sequences) || !sum(reads)) {
        return(c(mean = 0, sd = 0))
    }
    widths <- nchar(sequences)
    usable <- widths > 0L
    if (!any(usable)) {
        return(c(mean = 0, sd = 0))
    }
    gc <- (nchar(gsub("[^GCgc]", "", sequences[usable]))/widths[usable])
    weights <- reads[usable]
    mean_gc <- sum(gc * weights)/sum(weights)
    variance_gc <- sum(weights * (gc - mean_gc)^2)/sum(weights)
    c(mean = mean_gc, sd = sqrt(variance_gc))
}
.weighted_bad_alphabet_fraction <- function(sequences, reads) {
    total <- sum(reads)
    if (!length(sequences) || !total) {
        return(0)
    }
    bad <- grepl("[^ACGTUN]", toupper(sequences))
    sum(reads[bad])/total
}
.sequence_kmers <- function(sequence, k) {
    n_windows <- nchar(sequence) - k + 1L
    if (n_windows < 1L) {
        return(character())
    }
    starts <- seq_len(n_windows)
    substring(sequence, starts, starts + k - 1L)
}
.top_kmer_fraction <- function(sequences, reads, k = 6L, top_sequences = 10000L) {
    if (!length(sequences) || k < 1L || top_sequences < 1L) {
        return(0)
    }
    selected <- utils::head(order(reads, decreasing = TRUE, method = "radix"), top_sequences)
    kmer_counts <- numeric()
    for (index in selected) {
        kmers <- .sequence_kmers(toupper(sequences[[index]]), k)
        if (!length(kmers)) {
            next
        }
        unique_kmers <- unique(kmers)
        occurrences <- tabulate(match(kmers, unique_kmers), nbins = length(unique_kmers))
        current <- kmer_counts[unique_kmers]
        current[is.na(current)] <- 0
        kmer_counts[unique_kmers] <- current + occurrences * reads[[index]]
    }
    total <- sum(kmer_counts)
    if (!total)
        0 else max(kmer_counts)/total
}
.with_qc_seed <- function(seed, code) {
    withr::with_seed(as.integer(seed), code)
}
# Sample a count table without replacement via sequential hypergeometric
# draws.
.rarefy_counts <- function(counts, depth, seed = 42L) {
    if (!is.numeric(depth) || length(depth) != 1L || is.na(depth) || depth < 1L) {
        stop("`depth` must be one positive integer.", call. = FALSE)
    }
    counts <- counts[counts > 0]
    total <- sum(counts)
    if (!length(counts) || depth >= total) {
        return(counts)
    }
    depth <- as.integer(depth)
    .with_qc_seed(seed, {
        sampled <- numeric(length(counts))
        remaining_population <- total
        remaining_draws <- depth
        if (length(counts) > 1L) {
            for (index in seq_len(length(counts) - 1L)) {
                successes <- counts[[index]]
                sampled[[index]] <- stats::rhyper(1L, m = successes, n = remaining_population -
                  successes, k = remaining_draws)
                remaining_population <- remaining_population - successes
                remaining_draws <- remaining_draws - sampled[[index]]
            }
        }
        sampled[[length(counts)]] <- remaining_draws
        keep <- sampled > 0
        stats::setNames(sampled[keep], names(counts)[keep])
    })
}
.round_kmer_set <- function(sequences, reads, k = 6L, top_sequences = 10000L) {
    if (!length(sequences) || k < 1L) {
        return(character())
    }
    selected <- utils::head(order(reads, decreasing = TRUE, method = "radix"),
        top_sequences)
    sequences <- .normalise_adapter_alphabet(sequences[selected])
    kmers <- unique(unlist(lapply(sequences, .sequence_kmers, k = k),
        use.names = FALSE))
    if (!length(kmers)) {
        return(kmers)
    }
    reverse <- reverse_complement(kmers)
    unique(ifelse(kmers < reverse, kmers, reverse))
}
.jaccard_distance <- function(left, right) {
    if (!length(left) && !length(right)) {
        return(0)
    }
    1 - length(intersect(left, right))/length(union(left, right))
}
.qc_round_order <- function(round_names) {
    suffix <- suppressWarnings(as.numeric(sub(".*?([0-9]+)$", "\\1", round_names)))
    order(is.na(suffix), suffix, round_names, method = "radix")
}
.counts_by_round <- function(experiment) {
    counts <- SummarizedExperiment::assay(experiment, "counts")
    sequences <- as.character(SummarizedExperiment::rowData(experiment)$sequence)
    result <- stats::setNames(lapply(seq_len(ncol(counts)), function(index) {
        if (inherits(counts, "Matrix")) {
            entries <- Matrix::summary(counts[, index, drop = FALSE])
            return(stats::setNames(as.numeric(entries$x), sequences[entries$i]))
        }
        values <- as.numeric(counts[, index])
        keep <- values > 0
        stats::setNames(values[keep], sequences[keep])
    }), colnames(counts))
    result[.qc_round_order(names(result))]
}
.modal_sequence_length <- function(sequences, reads) {
    if (!length(sequences)) {
        return(NA_real_)
    }
    lengths <- nchar(sequences)
    unique_lengths <- unique(lengths)
    weighted <- vapply(unique_lengths, function(length) sum(reads[lengths == length]),
        numeric(1))
    unique_lengths[[which.max(weighted)]]
}
.round_kmer_similarity <- function(counts_by_round, k, top_sequences) {
    names_by_round <- names(counts_by_round)
    if (length(names_by_round) < 2L) {
        return(S4Vectors::DataFrame(from_round = character(), to_round = character(),
            jaccard_distance = numeric()))
    }
    kmer_sets <- lapply(counts_by_round, function(counts) {
        .round_kmer_set(names(counts), counts, k = k, top_sequences = top_sequences)
    })
    pairs <- utils::combn(seq_along(kmer_sets), 2L)
    S4Vectors::DataFrame(from_round = names_by_round[pairs[1L, ]], to_round = names_by_round[pairs[2L,
        ]], jaccard_distance = vapply(seq_len(ncol(pairs)), function(index) {
        .jaccard_distance(kmer_sets[[pairs[1L, index]]], kmer_sets[[pairs[2L, index]]])
    }, numeric(1)))
}
.round_kmer_monotonicity_flag <- function(round_similarity) {
    round_names <- unique(c(as.character(round_similarity$from_round), as.character(round_similarity$to_round)))
    if (length(round_names) < 3L) {
        return(NULL)
    }
    lookup <- stats::setNames(round_similarity$jaccard_distance, paste(round_similarity$from_round,
        round_similarity$to_round, sep = "\r"))
    distance <- function(left, right) lookup[[paste(left, right, sep = "\r")]]
    violations <- list()
    for (start in seq_len(length(round_names) - 2L)) {
        neighbouring <- distance(round_names[[start]], round_names[[start + 1L]])
        for (end in seq.int(start + 2L, length(round_names))) {
            distant <- distance(round_names[[start]], round_names[[end]])
            if (!is.null(neighbouring) && !is.null(distant) && neighbouring > distant *
                1.1) {
                violations <- append(violations, list(list(from_round = round_names[[start]],
                  neighbour_round = round_names[[start + 1L]], distant_round = round_names[[end]],
                  consecutive_distance = neighbouring, distant_distance = distant)))
            }
        }
    }
    if (!length(violations)) {
        return(NULL)
    }
    .qc_flag("round_kmer_nonmonotonic", "warn", list(tolerance = 1.1, violations = violations))
}
.adapter_contamination_flag <- function(library_report, experiment) {
    hits <- library_report$known_adapter_hits
    if (!length(hits) || !any(hits > 0)) {
        return(NULL)
    }
    extraction <- S4Vectors::metadata(experiment)$extraction
    input_reads <- if (inherits(extraction, "selexprep_extraction"))
        extraction$input_reads else NULL
    denominator <- if (is.null(input_reads) || !length(input_reads)) {
        0
    } else {
        as.numeric(input_reads[[1L]]) * library_report$read_fraction_used_for_inference
    }
    if (!is.finite(denominator) || denominator <= 0) {
        return(.qc_flag("adapter_contamination_high", "info", list(reason = "denominator_unavailable",
            known_adapter_hits = as.list(hits))))
    }
    fractions <- hits/denominator
    above <- which(fractions > 0.05)
    if (!length(above)) {
        return(NULL)
    }
    .qc_flag("adapter_contamination_high", "warn", list(max_fraction = 0.05, denominator = denominator,
        adapters_above_threshold = lapply(above, function(index) list(adapter = names(hits)[[index]],
            hits = as.numeric(hits[[index]]), fraction_of_reads = as.numeric(fractions[[index]])))))
}
.strand_mix_flag <- function(experiment) {
    extraction <- S4Vectors::metadata(experiment)$extraction
    distributions <- if (inherits(extraction, "selexprep_extraction")) {
        extraction$strand_distribution
    } else {
        NULL
    }
    if (is.null(distributions) || !length(distributions)) {
        return(NULL)
    }
    flagged <- lapply(names(distributions), function(round) {
        distribution <- distributions[[round]]
        total <- sum(distribution)
        if (!total || distribution[["reverse"]]/total <= 0.2) {
            return(NULL)
        }
        list(round = round, reverse_fraction = as.numeric(distribution[["reverse"]]/total),
            forward = as.integer(distribution[["forward"]]), reverse = as.integer(distribution[["reverse"]]),
            ambiguous = as.integer(distribution[["ambiguous"]]))
    })
    flagged <- Filter(Negate(is.null), flagged)
    if (!length(flagged)) {
        return(NULL)
    }
    .qc_flag("strand_mix", "warn", list(max_reverse_fraction = 0.2, rounds_above_threshold = flagged))
}
.qc_flag <- function(name, severity, evidence) {
    list(name = name, severity = severity, evidence = evidence)
}
.flags_dataframe <- function(flags) {
    if (!length(flags)) {
        return(S4Vectors::DataFrame(name = character(), severity = character(),
            evidence = I(list())))
    }
    S4Vectors::DataFrame(name = vapply(flags, `[[`, character(1), "name"), severity = vapply(flags,
        `[[`, character(1), "severity"), evidence = I(lapply(flags, `[[`, "evidence")))
}
#' Quality-control summary for a SELEX count experiment
#'
#' Calculates depth, diversity, enrichment and sequence-length summaries for
#' each round of a `SummarizedExperiment` returned by `run_selexprep()`. It
#' also returns actionable flags for low depth, unstable random-region length,
#' weak primer evidence and paired-end libraries that need read merging.
#'
#' @param experiment A `SummarizedExperiment` with an assay named `'counts'`
#'   and a `sequence` column in `rowData`.
#' @param library_report Optional `selexprep_library_report`. When omitted, a
#'   report stored in the experiment metadata is used when available.
#' @param low_total_reads Minimum number of reads expected in every round.
#' @param rarefaction_depth Maximum common depth used when comparing observed
#'   diversity across rounds.
#' @param rarefaction_seed Integer seed for deterministic rarefaction.
#' @param kmer_size Length of canonical k-mers used for between-round
#'   consistency diagnostics.
#' @param kmer_top_sequences Maximum number of abundant sequences considered
#'   from each round for k-mer diagnostics.
#'
#' @return A `selexprep_qc` list with `per_round`, `round_similarity`, and
#'   `flags` DataFrames plus the applied `settings`.
#' @export
#' @examples
#' p5 <- 'GGTAATACGACTCACTATAGGG'
#' p3 <- 'CCATGCATGCATGCATGCAT'
#' inserts <- c('ACGTACGTACGTACGT', 'TGCATGCATGCATGCA')
#' pools <- list(round_00 = rep(paste0(p5, inserts, p3), each = 250))
#' result <- run_selexprep(pools, low_total_reads = 0)
#' selexprep_qc(result)
selexprep_qc <- function(experiment, library_report = NULL,
    low_total_reads = 10000L, rarefaction_depth = 10000L,
    rarefaction_seed = 42L, kmer_size = 6L, kmer_top_sequences = 10000L) {
    settings <- list(
        rarefaction_depth = rarefaction_depth,
        rarefaction_seed = rarefaction_seed,
        kmer_size = kmer_size, kmer_top_sequences = kmer_top_sequences
    )
    library_report <- .validate_qc_inputs(
        experiment, library_report, low_total_reads, settings
    )
    counts <- .counts_by_round(experiment)
    totals <- vapply(counts, sum, numeric(1))
    positive_totals <- totals[totals > 0]
    depth <- if (length(positive_totals)) {
        min(as.numeric(rarefaction_depth), min(positive_totals))
    } else {
        NA_real_
    }
    per_round <- .qc_per_round(counts, totals, depth, settings)
    similarity <- .round_kmer_similarity(
        counts, k = as.integer(kmer_size),
        top_sequences = as.integer(kmer_top_sequences)
    )
    flags <- .qc_flags(
        per_round, similarity, experiment, library_report,
        low_total_reads, rarefaction_depth, depth
    )
    structure(list(
        per_round = per_round, round_similarity = similarity,
        flags = .flags_dataframe(flags), settings = list(
            rarefaction_depth = as.integer(rarefaction_depth),
            effective_rarefaction_depth = depth,
            rarefaction_seed = as.integer(rarefaction_seed),
            kmer_size = as.integer(kmer_size),
            kmer_top_sequences = as.integer(kmer_top_sequences)
        )
    ), class = c("selexprep_qc", "list"))
}

.validate_qc_inputs <- function(experiment, report, low_total_reads, settings) {
    if (!inherits(experiment, "SummarizedExperiment")) {
        stop("`experiment` must be a SummarizedExperiment.", call. = FALSE)
    }
    if (!"counts" %in% SummarizedExperiment::assayNames(experiment) ||
        !"sequence" %in% colnames(SummarizedExperiment::rowData(experiment))) {
        stop("`experiment` must contain a counts assay and rowData$sequence.",
            call. = FALSE)
    }
    if (is.null(report)) {
        report <- S4Vectors::metadata(experiment)$library_report
    }
    if (!is.null(report) && !inherits(report, "selexprep_library_report")) {
        stop("`library_report` must be NULL or a selexprep_library_report.",
            call. = FALSE)
    }
    if (!is.numeric(low_total_reads) || length(low_total_reads) != 1L ||
        is.na(low_total_reads) || low_total_reads < 0L) {
        stop("`low_total_reads` must be one non-negative number.",
            call. = FALSE)
    }
    for (argument in names(settings)) {
        value <- settings[[argument]]
        if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            value < 1L || value != floor(value)) {
            stop(sprintf("`%s` must be one positive integer.", argument),
                call. = FALSE)
        }
    }
    report
}

.qc_rarefied_unique <- function(counts, depth, seed) {
    if (is.na(depth) || !sum(counts)) {
        return(NA_integer_)
    }
    length(.rarefy_counts(counts, depth, seed))
}

.qc_truseq_fraction <- function(counts) {
    if (!sum(counts)) {
        return(0)
    }
    hits <- grepl("AGATCGGAAGAGC", toupper(names(counts)), fixed = TRUE)
    sum(counts[hits]) / sum(counts)
}

.qc_per_round <- function(counts, totals, depth, settings) {
    gc <- lapply(counts, function(values) {
        .weighted_gc_content(names(values), values)
    })
    S4Vectors::DataFrame(
        round = names(counts), n_reads = unname(totals),
        n_unique = unname(vapply(counts, length, integer(1))),
        shannon_entropy_bits = unname(vapply(
            counts, .shannon_entropy, numeric(1)
        )),
        rarefied_unique = unname(vapply(
            counts, .qc_rarefied_unique, integer(1), depth = depth,
            seed = settings$rarefaction_seed
        )),
        singleton_fraction = unname(vapply(
            counts, .singleton_fraction, numeric(1)
        )),
        top_1_coverage = unname(vapply(
            counts, .top_n_coverage, numeric(1), n = 1L
        )),
        top_100_coverage = unname(vapply(counts, .top_n_coverage, numeric(1))),
        modal_sequence_length = unname(vapply(counts, function(values) {
            .modal_sequence_length(names(values), values)
        }, numeric(1))),
        mean_gc = unname(vapply(gc, `[[`, numeric(1), "mean")),
        gc_sd = unname(vapply(gc, `[[`, numeric(1), "sd")),
        nonstandard_alphabet_fraction = unname(vapply(counts, function(values) {
            .weighted_bad_alphabet_fraction(names(values), values)
        }, numeric(1))),
        truseq_reads_fraction = unname(vapply(
            counts, .qc_truseq_fraction, numeric(1)
        )),
        top_kmer_fraction = unname(vapply(counts, function(values) {
            .top_kmer_fraction(
                names(values), values, k = settings$kmer_size,
                top_sequences = settings$kmer_top_sequences
            )
        }, numeric(1)))
    )
}

.qc_threshold_flag <- function(per_round, column, threshold, name,
    evidence_name = "fractions", above = FALSE) {
    values <- per_round[[column]]
    indices <- which(if (above) values > threshold else values < threshold)
    if (!length(indices)) {
        return(NULL)
    }
    evidence <- list(
        threshold = threshold, rounds = as.list(per_round$round[indices])
    )
    evidence[[evidence_name]] <- as.list(values[indices])
    .qc_flag(name, "warn", evidence)
}

.qc_length_flag <- function(per_round) {
    lengths <- unique(per_round$modal_sequence_length[
        !is.na(per_round$modal_sequence_length)
    ])
    if (length(lengths) <= 2L) {
        return(NULL)
    }
    .qc_flag("n_length_variation_across_rounds", "warn", list(
        modal_lengths = as.list(per_round$modal_sequence_length),
        distinct_modes = as.list(sort(lengths))
    ))
}

.qc_diversity_flag <- function(per_round, configured_depth, depth) {
    if (nrow(per_round) < 2L || is.na(depth)) {
        return(NULL)
    }
    increases <- which(diff(per_round$rarefied_unique) > 0)
    if (!length(increases)) {
        return(NULL)
    }
    .qc_flag("unexpected_rarefied_diversity_increase", "warn", list(
        configured_depth = as.integer(configured_depth),
        effective_depth = depth,
        rarefied_unique_per_round = as.list(stats::setNames(
            per_round$rarefied_unique, per_round$round
        )),
        increases = lapply(increases, function(index) list(
            from_round = per_round$round[[index]],
            to_round = per_round$round[[index + 1L]],
            rarefied_unique_prev = per_round$rarefied_unique[[index]],
            rarefied_unique_curr = per_round$rarefied_unique[[index + 1L]]
        ))
    ))
}

.qc_library_flags <- function(report, experiment) {
    if (is.null(report)) {
        return(list())
    }
    sides <- character()
    if (report$match_rate_5p < 0.4) sides <- c(sides, "5p")
    if (report$match_rate_3p < 0.4) sides <- c(sides, "3p")
    list(
        if (length(sides)) .qc_flag("low_primer_match", "warn", list(
            threshold = 0.4, sides = as.list(sides)
        )),
        if (identical(report$required_action, "READ_MERGING_RECOMMENDED")) {
            .qc_flag("requires_read_merging_for_full_insert", "info", list(
                extraction_mode = report$extraction_mode
            ))
        },
        .adapter_contamination_flag(report, experiment)
    )
}

.qc_flags <- function(per_round, similarity, experiment, report,
    low_total_reads, configured_depth, depth) {
    flags <- list(
        .qc_threshold_flag(
            per_round, "n_reads", low_total_reads, "low_total_reads", "reads"
        ),
        .qc_length_flag(per_round),
        .qc_diversity_flag(per_round, configured_depth, depth),
        .qc_threshold_flag(
            per_round, "nonstandard_alphabet_fraction", 0.001,
            "nonstandard_alphabet", above = TRUE
        ),
        .qc_threshold_flag(
            per_round, "truseq_reads_fraction", 0.01, "truseq_residual",
            above = TRUE
        ),
        .round_kmer_monotonicity_flag(similarity),
        .strand_mix_flag(experiment)
    )
    Filter(Negate(is.null), c(flags, .qc_library_flags(report, experiment)))
}
