.validate_demultiplex_barcodes <- function(barcodes, max_mismatches) {
    valid_rounds <- is.numeric(barcodes) && length(barcodes) > 0L &&
        !is.null(names(barcodes)) && all(nzchar(names(barcodes))) &&
        !anyNA(barcodes) && all(is.finite(barcodes)) &&
        all(barcodes >= 0) && all(barcodes == floor(barcodes))
    if (!valid_rounds || anyDuplicated(names(barcodes))) {
        stop(
            "`barcodes` must be a named non-negative integer vector.",
            call. = FALSE
        )
    }
    valid_mismatches <- is.numeric(max_mismatches) &&
        length(max_mismatches) == 1L && !is.na(max_mismatches) &&
        max_mismatches >= 0 && max_mismatches == floor(max_mismatches)
    if (!valid_mismatches) {
        stop(
            "`max_mismatches` must be one non-negative integer.",
            call. = FALSE
        )
    }
    barcode_sequences <- .normalise_adapter_alphabet(names(barcodes))
    widths <- nchar(barcode_sequences)
    if (length(unique(widths)) != 1L || any(widths == 0L) ||
        any(grepl("[^ACGT]", barcode_sequences))) {
        stop(
            "Barcode names must be equal-length DNA or RNA sequences.",
            call. = FALSE
        )
    }
    .validate_barcode_distance(barcode_sequences, names(barcodes), max_mismatches)
    stats::setNames(as.integer(barcodes), barcode_sequences)
}

.validate_barcode_distance <- function(barcode_sequences, display_names,
    max_mismatches) {
    minimum_distance <- 2L * as.integer(max_mismatches) + 1L
    barcode_bases <- strsplit(barcode_sequences, "", fixed = TRUE)
    if (length(barcode_sequences) > 1L) {
        for (left in seq_len(length(barcode_sequences) - 1L)) {
            for (right in seq.int(left + 1L, length(barcode_sequences))) {
                distance <- sum(
                    barcode_bases[[left]] != barcode_bases[[right]]
                )
                if (distance < minimum_distance) {
                    stop(
                        sprintf(
                            "Barcodes %s and %s have distance %d; at least %d is required.",
                            display_names[[left]],
                            display_names[[right]],
                            distance,
                            minimum_distance
                        ),
                        call. = FALSE
                    )
                }
            }
        }
    }
}

.demultiplex_assignments <- function(sequences, barcodes, max_mismatches) {
    barcode_width <- nchar(names(barcodes)[[1L]])
    barcode_bases <- stats::setNames(
        strsplit(names(barcodes), "", fixed = TRUE), names(barcodes)
    )
    vapply(sequences, function(sequence) {
        sequence <- .normalise_adapter_alphabet(sequence)
        if (nchar(sequence) < barcode_width) {
            return(NA_integer_)
        }
        prefix <- strsplit(substr(sequence, 1L, barcode_width), "")[[1L]]
        distances <- vapply(barcode_bases, function(barcode) {
            sum(prefix != barcode)
        }, integer(1))
        matched <- which(distances <= max_mismatches)
        if (length(matched) == 1L) barcodes[[matched]] else NA_integer_
    }, integer(1))
}

#' Demultiplex SELEX reads with user-supplied barcodes
#'
#' Assigns reads to selection rounds from an equal-length barcode at the R1
#' 5-prime end. Barcode sequences are supplied explicitly; they are never
#' inferred. R2 mates follow their R1 assignment and remain untrimmed.
#'
#' Barcodes must have pairwise Hamming distance of at least
#' `2 * max_mismatches + 1`, preventing ambiguous assignments. RNA `U` and DNA
#' `T` are treated as equivalent during matching.
#'
#' @param sequences R1 sequences as a character vector or
#'   `Biostrings::XStringSet`.
#' @param barcodes Named integer vector mapping barcode sequences to
#'   non-negative round numbers, for example `c(AAAAA = 0L, TTTTT = 1L)`.
#' @param paired_mates Optional R2 sequences in the same order as `sequences`.
#' @param max_mismatches Maximum Hamming distance allowed from a barcode.
#' @param trim_barcode Whether to remove the matched barcode from assigned R1
#'   sequences.
#'
#' @return A `selexprep_demultiplex` object. Its first four fields follow the
#'   `selexprep_read_inputs` contract and can be supplied directly to
#'   [run_selexprep()]. Additional fields contain unassigned reads, a summary,
#'   and the validated barcode map.
#' @export
#'
#' @examples
#' reads <- c("AAAAACCCCC", "AAAATGGGGG", "TTTTTACGTA", "GGGGGAAAAA")
#' demultiplexed <- selexprep_demultiplex(
#'     reads,
#'     c(AAAAA = 0L, TTTTT = 1L)
#' )
#' demultiplexed$summary
selexprep_demultiplex <- function(sequences, barcodes, paired_mates = NULL,
    max_mismatches = 1L, trim_barcode = TRUE) {
    sequences <- .as_sequence_character(sequences)
    paired_mates <- .prepare_demultiplex_mates(
        sequences, paired_mates, trim_barcode
    )
    barcodes <- .validate_demultiplex_barcodes(barcodes, max_mismatches)
    assignments <- .demultiplex_assignments(
        sequences, barcodes, as.integer(max_mismatches)
    )
    rounds <- sort(unique(unname(barcodes)), method = "radix")
    r1 <- .demultiplex_round_pools(
        sequences, assignments, rounds,
        trim_width = if (trim_barcode) nchar(names(barcodes)[[1L]]) else 0L
    )
    r2 <- if (is.null(paired_mates)) {
        NULL
    } else {
        .demultiplex_round_pools(paired_mates, assignments, rounds)
    }
    summary <- .demultiplex_summary(rounds, assignments, length(sequences))
    unassigned <- is.na(assignments)
    unassigned_reads <- list(
        r1 = Biostrings::BStringSet(unname(sequences[unassigned])),
        r2 = if (is.null(paired_mates)) {
            NULL
        } else {
            Biostrings::BStringSet(unname(paired_mates[unassigned]))
        }
    )
    structure(list(
        sequences_by_round = r1, paired_mate_streams = r2,
        read_source = if (is.null(r2)) "R1" else "R1_AND_R2",
        files = S4Vectors::DataFrame(), unassigned = unassigned_reads,
        summary = summary, barcodes = barcodes,
        max_mismatches = as.integer(max_mismatches), trim_barcode = trim_barcode
    ), class = c("selexprep_demultiplex", "selexprep_read_inputs", "list"))
}

.prepare_demultiplex_mates <- function(sequences, paired_mates, trim_barcode) {
    if (!is.null(paired_mates)) {
        paired_mates <- .as_sequence_character(paired_mates)
        if (length(sequences) != length(paired_mates)) {
            stop("`paired_mates` must contain one mate per R1 sequence.",
                call. = FALSE)
        }
    }
    if (!is.logical(trim_barcode) || length(trim_barcode) != 1L ||
        is.na(trim_barcode)) {
        stop("`trim_barcode` must be one logical value.", call. = FALSE)
    }
    paired_mates
}

.demultiplex_round_pools <- function(sequences, assignments, rounds,
    trim_width = 0L) {
    pools <- lapply(rounds, function(round) {
        selected <- sequences[!is.na(assignments) & assignments == round]
        if (trim_width) selected <- substring(selected, trim_width + 1L)
        Biostrings::BStringSet(unname(selected))
    })
    stats::setNames(pools, sprintf("round_%02d", rounds))
}

.demultiplex_summary <- function(rounds, assignments, total) {
    assigned <- vapply(rounds, function(round) {
        sum(assignments == round, na.rm = TRUE)
    }, integer(1))
    S4Vectors::DataFrame(
        round = sprintf("round_%02d", rounds), round_number = rounds,
        assigned_reads = assigned,
        assigned_fraction = if (total) assigned / total else numeric(length(rounds))
    )
}
