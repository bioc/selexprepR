.adapter_match <- function(sequence, primer) {
    sequence <- .normalise_adapter_alphabet(sequence)
    primer <- .normalise_adapter_alphabet(primer)
    max_errors <- floor(nchar(primer) * 0.1)
    subject <- Biostrings::BString(sequence)
    pattern <- Biostrings::BString(primer)
    hits <- Biostrings::matchPattern(pattern, subject, max.mismatch = max_errors,
        with.indels = TRUE, fixed = TRUE)
    if (!length(hits)) {
        return(NULL)
    }
    starts <- Biostrings::start(hits)
    ends <- Biostrings::end(hits)
    errors <- Biostrings::neditAt(pattern, subject, at = starts, with.indels = TRUE,
        fixed = TRUE)
    widths <- ends - starts + 1L
    scores <- pmin(nchar(primer), widths) - 2L * errors
    best <- order(-scores, errors, starts, -ends, method = "radix")[[1L]]
    c(start = starts[[best]], end = ends[[best]])
}
# Exact matches have the maximum possible adapter score. Find those in one
# native Biostrings call; only distinct residual reads need indel alignment.
.adapter_matches <- function(sequences, primer) {
    coordinates <- matrix(NA_integer_, nrow = length(sequences), ncol = 2L,
        dimnames = list(NULL, c("start", "end")))
    if (!length(sequences)) {
        return(coordinates)
    }
    normalised <- .normalise_adapter_alphabet(sequences)
    primer <- .normalise_adapter_alphabet(primer)
    hits <- Biostrings::vmatchPattern(Biostrings::BString(primer),
        Biostrings::BStringSet(normalised), fixed = TRUE)
    counts <- S4Vectors::elementNROWS(hits)
    matched <- counts > 0L
    first <- cumsum(counts) - counts + 1L
    starts <- unlist(Biostrings::start(hits), use.names = FALSE)
    coordinates[matched, "start"] <- starts[first[matched]]
    coordinates[matched, "end"] <- coordinates[matched, "start"] +
        nchar(primer) - 1L
    residual <- which(!matched)
    if (length(residual)) {
        coordinates[residual, ] <- .approximate_adapter_matches(
            normalised[residual], primer)
    }
    coordinates
}

.approximate_adapter_matches <- function(sequences, primer) {
    distinct <- unique(sequences)
    coordinates <- matrix(NA_integer_, nrow = length(distinct), ncol = 2L,
        dimnames = list(NULL, c("start", "end")))
    candidates <- Biostrings::vcountPattern(Biostrings::BString(primer),
        Biostrings::BStringSet(distinct),
        max.mismatch = floor(nchar(primer) * 0.1),
        with.indels = TRUE, fixed = TRUE) > 0L
    # Biostrings vectorizes indel-aware counts but not their coordinates.
    # Resolve only matching, distinct residual reads with the original ranking.
    if (any(candidates)) {
        coordinates[candidates, ] <- t(vapply(distinct[candidates],
            function(sequence) .adapter_match(sequence, primer), integer(2)))
    }
    coordinates[match(sequences, distinct), , drop = FALSE]
}

.trim_adapter <- function(sequences, primer, side) {
    coordinates <- .adapter_matches(sequences, primer)
    kept <- rep(NA_character_, length(sequences))
    matched <- !is.na(coordinates[, "start"])
    kept[matched] <- if (identical(side, "after")) {
        substring(sequences[matched], coordinates[matched, "end"] + 1L)
    } else {
        substr(sequences[matched], 1L, coordinates[matched, "start"] - 1L)
    }
    kept
}

.trim_round <- function(sequences, mode, primer_5p, primer_3p,
    keep_mask = FALSE) {
    kept <- switch(mode,
        BOTH_PRIMERS_SINGLE_READ = {
            first <- .trim_adapter(sequences, primer_5p, "after")
            matched <- !is.na(first)
            first[matched] <- .trim_adapter(first[matched], primer_3p, "before")
            first
        },
        FIVE_PRIME_ONLY = .trim_adapter(sequences, primer_5p, "after"),
        THREE_PRIME_ONLY = .trim_adapter(sequences, primer_3p, "before"),
        stop("Unsupported single-read extraction mode.", call. = FALSE)
    )
    names(kept) <- if (is.null(names(sequences))) sequences else names(sequences)
    keep <- !is.na(kept)
    result <- kept[keep]
    if (keep_mask) {
        attr(result, "kept") <- keep
    }
    result
}

.prefix_matches <- function(sequences, primer) {
    matched <- rep(FALSE, length(sequences))
    if (is.null(primer)) {
        return(matched)
    }
    sequences <- .normalise_adapter_alphabet(sequences)
    primer <- .normalise_adapter_alphabet(primer)
    long_enough <- nchar(sequences) >= nchar(primer)
    matched[long_enough] <- as.vector(Biostrings::neditStartingAt(
        Biostrings::BString(primer),
        Biostrings::BStringSet(sequences[long_enough]),
        starting.at = 1L, with.indels = FALSE, fixed = TRUE
    )) <= 1L
    matched
}

.strand_distribution <- function(sequences, primer_5p, primer_3p) {
    reverse_3p <- if (is.null(primer_3p)) NULL else reverse_complement(primer_3p)
    forward <- .prefix_matches(sequences, primer_5p)
    reverse <- !forward & .prefix_matches(sequences, reverse_3p)
    c(forward = as.integer(sum(forward)), reverse = as.integer(sum(reverse)),
        ambiguous = as.integer(sum(!forward & !reverse)))
}
.low_yield_rounds <- function(input_reads, output_reads, label = "round") {
    nonempty <- input_reads > 0
    fractions <- output_reads[nonempty] / input_reads[nonempty]
    if (!length(fractions)) {
        return(character())
    }
    best <- max(fractions)
    low <- fractions < 0.5 | fractions < 0.5 * best
    rounds <- names(fractions)[low]
    descriptions <- sprintf(
        "%s: kept %s of %s reads (%.1f%%; best %s %.1f%%)",
        rounds, output_reads[rounds], input_reads[rounds],
        100 * fractions[low], label, 100 * best
    )
    stats::setNames(descriptions, rounds)
}

.new_extraction <- function(mode, full_insert_recovered, sequences_by_round, input_reads,
    output_reads, strand_distribution = list(), kept_by_round = NULL) {
    low_yield <- .low_yield_rounds(input_reads, output_reads)
    if (length(low_yield)) {
        warning(paste0(
            "Extraction: ", length(low_yield),
            " round(s) kept under half of their reads or under half the best ",
            "round's yield; the primers may not fit their reads ",
            "(run-specific outer sequence or a different construct): ",
            paste(low_yield, collapse = "; ")
        ), call. = FALSE)
    }
    result <- structure(list(extraction_mode = mode, full_insert_recovered = full_insert_recovered,
        sequences_by_round = sequences_by_round, input_reads = input_reads, output_reads = output_reads,
        strand_distribution = strand_distribution, low_yield_rounds = low_yield), class = c("selexprep_extraction",
        "list"))
    if (!is.null(kept_by_round)) {
        attr(result, "kept_by_round") <- kept_by_round
    }
    result
}

.validate_primer_override <- function(primer, argument) {
    if (is.null(primer)) {
        return(NULL)
    }
    valid <- is.character(primer) && length(primer) == 1L &&
        !is.na(primer) && nzchar(primer)
    if (!valid || grepl("[^ACGTU]", toupper(primer))) {
        stop(
            sprintf(
                "`%s` must be one non-empty DNA or RNA sequence.",
                argument
            ),
            call. = FALSE
        )
    }
    toupper(primer)
}

.override_extraction_mode <- function(fields, extraction_mode) {
    if (is.null(extraction_mode)) {
        extraction_mode <- fields$extraction_mode
        if (identical(extraction_mode, "UNABLE_TO_EXTRACT")) {
            extraction_mode <- if (!is.null(fields$primer_5p) &&
                !is.null(fields$primer_3p)) {
                "BOTH_PRIMERS_SINGLE_READ"
            } else if (!is.null(fields$primer_5p)) {
                "FIVE_PRIME_ONLY"
            } else if (!is.null(fields$primer_3p)) {
                "THREE_PRIME_ONLY"
            } else {
                "UNABLE_TO_EXTRACT"
            }
        }
    }
    valid <- is.character(extraction_mode) && length(extraction_mode) == 1L &&
        !is.na(extraction_mode) &&
        extraction_mode %in% setdiff(.extraction_modes, "UNABLE_TO_EXTRACT")
    if (!valid) {
        stop("`extraction_mode` must name an extractable mode.", call. = FALSE)
    }
    extraction_mode
}

.validate_override_primers <- function(fields, extraction_mode) {
    required <- switch(extraction_mode,
        BOTH_PRIMERS_SINGLE_READ = c("primer_5p", "primer_3p"),
        FIVE_PRIME_ONLY = "primer_5p",
        THREE_PRIME_ONLY = "primer_3p",
        PAIRED_END_SPLIT_PRIMERS = c("primer_5p", "primer_3p")
    )
    if (any(vapply(fields[required], is.null, logical(1)))) {
        stop("The selected extraction mode requires additional primers.",
            call. = FALSE)
    }
}

.apply_primer_overrides <- function(report, primer_5p = NULL,
    primer_3p = NULL, extraction_mode = NULL) {
    requested <- !is.null(primer_5p) || !is.null(primer_3p) ||
        !is.null(extraction_mode)
    if (!requested) {
        return(report)
    }
    primer_5p <- .validate_primer_override(primer_5p, "primer_5p")
    primer_3p <- .validate_primer_override(primer_3p, "primer_3p")
    fields <- unclass(report)
    if (!is.null(primer_5p)) {
        fields$primer_5p <- primer_5p
        fields$variants_5p <- list()
    }
    if (!is.null(primer_3p)) {
        fields$primer_3p <- primer_3p
        fields$variants_3p <- list()
    }
    extraction_mode <- .override_extraction_mode(fields, extraction_mode)
    .validate_override_primers(fields, extraction_mode)
    fields$extraction_mode <- extraction_mode
    fields$full_insert_recovered <- identical(extraction_mode,
        "BOTH_PRIMERS_SINGLE_READ")
    fields$required_action <- if (identical(extraction_mode,
        "PAIRED_END_SPLIT_PRIMERS")) "READ_MERGING_RECOMMENDED" else "NONE"
    fields$status <- "LOW"
    fields$confidence <- 0
    fields["failure_reason"] <- list(NULL)
    .new_library_report(fields)
}

#' Extract random regions from a SELEX library
#'
#' Uses a `selexprep_library_report` to remove the inferred constant regions.
#' The implementation operates directly on R sequence objects and preserves the
#' explicit safe-failure behavior of the original pipeline: reports marked as
#' unable to extract stop with an error unless the caller supplies a corrected
#' report. In split-primer paired-end mode, the two trimmed sides remain
#' separate; they are never silently merged.
#'
#' Full-primer matching follows the cutadapt defaults used by the Python
#' pipeline: adapters may occur away from the read boundary, substitutions
#' and indels are allowed up to a 10 percent error rate, and RNA `U` is
#' matched as DNA `T`. The returned insert retains the alphabet and case of
#' the input read.
#'
#' @param sequences_by_round A named list of character vectors or
#'   `Biostrings::XStringSet` objects containing R1 sequences.
#' @param library_report A `selexprep_library_report` from `selexprep_detect()`
#'   or `read_library_report()`.
#' @param paired_mate_streams Optional named R2 sequence pools, required for
#'   `PAIRED_END_SPLIT_PRIMERS` reports.
#' @param primer_5p Optional manually reviewed 5-prime primer override.
#' @param primer_3p Optional manually reviewed 3-prime primer override.
#' @param extraction_mode Optional manually reviewed extraction-mode override.
#'
#' @return A `selexprep_extraction`. Full and single-primer modes store one
#'   `Biostrings::BStringSet` per round; paired split mode stores an `r1` and
#'   `r2` `BStringSet` for every round. `low_yield_rounds` is a named character
#'   vector describing nonempty rounds retaining less than half their reads or
#'   less than half the best round's yield; these rounds also emit a warning.
#' @export
#'
#' @examples
#' primer_5p <- 'GGTAATACGACTCACTATAGGG'
#' primer_3p <- 'CCATGCATGCATGCATGCAT'
#' bases <- c('A', 'C', 'G', 'T')
#' random_regions <- vapply(0:499, function(i) {
#'     paste0(bases[(i * 7 + 0:15 * 13) %% 4 + 1], collapse = '')
#' }, character(1))
#' reads <- paste0(primer_5p, random_regions, primer_3p)
#' report <- selexprep_detect(list(round_00 = reads))
#' selexprep_extract(list(round_00 = reads[seq_len(3L)]), report)
selexprep_extract <- function(sequences_by_round, library_report,
    paired_mate_streams = NULL, primer_5p = NULL, primer_3p = NULL,
    extraction_mode = NULL) {
    .extract_selexprep(
        sequences_by_round, library_report,
        paired_mate_streams = paired_mate_streams,
        primer_5p = primer_5p, primer_3p = primer_3p,
        extraction_mode = extraction_mode
    )
}

.validated_extraction_report <- function(report, primer_5p, primer_3p,
    extraction_mode) {
    if (!inherits(report, "selexprep_library_report")) {
        stop("`library_report` must be a selexprep_library_report.", call. = FALSE)
    }
    report <- .apply_primer_overrides(report, primer_5p, primer_3p,
        extraction_mode)
    if (identical(report$status, "UNABLE_TO_INFER") ||
        identical(report$extraction_mode, "UNABLE_TO_EXTRACT")) {
        stop(
            "LibraryReport does not permit extraction; provide explicit corrected primers.",
            call. = FALSE
        )
    }
    report
}

.validate_single_primers <- function(report) {
    mode <- report$extraction_mode
    if (identical(mode, "BOTH_PRIMERS_SINGLE_READ") &&
        (is.null(report$primer_5p) || is.null(report$primer_3p))) {
        stop("Both-primer extraction requires both primer sequences.",
            call. = FALSE)
    }
    if (identical(mode, "FIVE_PRIME_ONLY") && is.null(report$primer_5p)) {
        stop("Five-prime extraction requires a 5' primer.", call. = FALSE)
    }
    if (identical(mode, "THREE_PRIME_ONLY") && is.null(report$primer_3p)) {
        stop("Three-prime extraction requires a 3' primer.", call. = FALSE)
    }
}

.single_extraction <- function(pools, report, strand_distribution, keep_masks) {
    .validate_single_primers(report)
    extracted <- lapply(pools, .trim_round, mode = report$extraction_mode,
        primer_5p = report$primer_5p, primer_3p = report$primer_3p,
        keep_mask = keep_masks)
    kept_by_round <- if (keep_masks) {
        lapply(extracted, attr, which = "kept")
    } else {
        NULL
    }
    extracted <- lapply(extracted, function(sequences) {
        Biostrings::BStringSet(unname(sequences))
    })
    .new_extraction(report$extraction_mode, report$full_insert_recovered,
        extracted, vapply(pools, length, numeric(1)),
        vapply(extracted, length, numeric(1)), strand_distribution,
        kept_by_round = kept_by_round)
}

.validated_mate_pools <- function(pools, paired_mate_streams, report) {
    if (is.null(paired_mate_streams)) {
        stop("PAIRED_END_SPLIT_PRIMERS extraction requires `paired_mate_streams`.",
            call. = FALSE)
    }
    mates <- .as_round_pools(paired_mate_streams, "paired_mate_streams")
    if (!identical(names(pools), names(mates))) {
        stop(
            "`paired_mate_streams` must have the same round names as `sequences_by_round`.",
            call. = FALSE
        )
    }
    if (is.null(report$primer_5p) || is.null(report$primer_3p)) {
        stop("Split-primer extraction requires both primer sequences.",
            call. = FALSE)
    }
    sizes <- vapply(pools, length, integer(1))
    mate_sizes <- vapply(mates, length, integer(1))
    if (any(sizes != mate_sizes)) {
        first <- names(pools)[which(sizes != mate_sizes)[[1L]]]
        stop(sprintf("Paired streams have different read counts in %s.", first),
            call. = FALSE)
    }
    mates
}

.trim_paired_round <- function(r1, r2, primer_5p, primer_r2, keep_mask) {
    trimmed_r1 <- .trim_adapter(r1, primer_5p, "after")
    trimmed_r2 <- .trim_adapter(r2, primer_r2, "after")
    keep <- !is.na(trimmed_r1) & !is.na(trimmed_r2)
    result <- list(r1 = Biostrings::BStringSet(unname(trimmed_r1[keep])),
        r2 = Biostrings::BStringSet(unname(trimmed_r2[keep])))
    if (keep_mask) {
        attr(result, "kept") <- keep
    }
    result
}

.paired_extraction <- function(pools, mates, report, strand_distribution,
    keep_masks) {
    primer_r2 <- reverse_complement(report$primer_3p)
    extracted <- Map(function(r1, r2) {
        .trim_paired_round(r1, r2, report$primer_5p, primer_r2, keep_masks)
    }, pools, mates)
    kept_by_round <- NULL
    if (keep_masks) {
        kept_by_round <- lapply(extracted, attr, which = "kept")
        extracted <- lapply(extracted, function(round) {
            attr(round, "kept") <- NULL
            round
        })
    }
    .new_extraction(report$extraction_mode, FALSE, extracted,
        vapply(pools, length, numeric(1)),
        vapply(extracted, function(round) length(round$r1), numeric(1)),
        strand_distribution, kept_by_round = kept_by_round)
}

.extract_selexprep <- function(sequences_by_round, library_report,
    paired_mate_streams = NULL, primer_5p = NULL, primer_3p = NULL,
    extraction_mode = NULL, keep_masks = FALSE) {
    report <- .validated_extraction_report(library_report, primer_5p, primer_3p,
        extraction_mode)
    pools <- .as_round_pools(sequences_by_round, "sequences_by_round")
    strand_distribution <- lapply(pools, .strand_distribution,
        primer_5p = report$primer_5p, primer_3p = report$primer_3p)
    if (identical(report$orientation, "REVERSE")) {
        pools <- lapply(pools, function(pool) {
            stats::setNames(reverse_complement(pool), pool)
        })
    }
    if (identical(report$extraction_mode, "PAIRED_END_SPLIT_PRIMERS")) {
        mates <- .validated_mate_pools(pools, paired_mate_streams, report)
        return(.paired_extraction(pools, mates, report, strand_distribution,
            keep_masks))
    }
    .single_extraction(pools, report, strand_distribution, keep_masks)
}
