# This file ports to R the aligner of cutadapt 5.2
# (https://github.com/marcelm/cutadapt, src/cutadapt/_align.pyx,
# Aligner.locate).
# cutadapt: Copyright (c) 2010 Marcel Martin <marcel.martin@scilifelab.se>
# and contributors. MIT license; the full notice is in
# inst/licenses/cutadapt-MIT.txt (see also inst/COPYRIGHTS).
#
# Python selexprep trims primers with cutadapt. Here a 5' primer is a cutadapt
# regular 5' adapter (-g) and a 3' primer a regular 3' adapter (-a): up to 10%
# errors (substitutions or indels) over the aligned part of the primer, at
# least 3 aligned primer bases, so a primer may be cut short at the read start
# (5') or end (3'). Among several occurrences cutadapt keeps the first
# acceptable one and replaces it only by one that overlaps it or covers more
# of the primer, and has a higher score (+1 match, -1 mismatch, -2 indel).

.PRIMER_ERROR_RATE <- 0.1
.PRIMER_MIN_OVERLAP <- 3L
.PRIMER_DP_CHUNK <- 5000L

.read_code_table <- local({
    table <- rep(128L, 256L)
    bases <- c(A = 1L, C = 2L, G = 4L, T = 8L, U = 8L)
    for (base in names(bases)) {
        table[utf8ToInt(base) + 1L] <- bases[[base]]
        table[utf8ToInt(tolower(base)) + 1L] <- bases[[base]]
    }
    table
})

.primer_code_table <- local({
    table <- integer(256L)
    codes <- c(A = 1L, C = 2L, G = 4L, T = 8L, U = 8L, R = 5L, Y = 10L,
        S = 6L, W = 9L, K = 12L, M = 3L, B = 14L, D = 13L, H = 11L, V = 7L,
        N = 15L + 128L)
    for (base in names(codes)) {
        table[utf8ToInt(base) + 1L] <- codes[[base]]
        table[utf8ToInt(tolower(base)) + 1L] <- codes[[base]]
    }
    table
})

.encode_reads <- function(sequences) {
    widths <- nchar(sequences, type = "bytes")
    codes <- matrix(0L, nrow = length(sequences), ncol = max(c(widths, 0L)))
    if (sum(widths) > 0L) {
        bytes <- as.integer(charToRaw(paste(sequences, collapse = "")))
        codes[cbind(rep.int(seq_along(sequences), widths), sequence(widths))] <-
            .read_code_table[bytes + 1L]
    }
    list(codes = codes, widths = widths)
}

# DP state of Aligner.locate for many reads, kept in an environment so the
# steps below update it in place. Rows of the matrices are reads; column i + 1
# holds DP row i (primer prefix of length i).
.aligner_state <- function(sequences, primer, front) {
    state <- new.env(parent = emptyenv())
    encoded <- .encode_reads(sequences)
    state$codes <- encoded$codes
    state$widths <- encoded$widths
    state$front <- front
    state$reference <- .primer_code_table[utf8ToInt(primer) + 1L]
    m <- state$m <- length(state$reference)
    state$k <- as.integer(.PRIMER_ERROR_RATE * m)
    state$wildcard <- grepl("[^ACGT]", primer)
    is_n <- strsplit(primer, "", fixed = TRUE)[[1L]] == "N"
    state$n_counts <- c(0L, cumsum(is_n))
    state$effective_length <- m - if (state$wildcard) sum(is_n) else 0L
    n <- length(sequences)
    rows <- seq_len(m + 1L) - 1L
    state$cost <- matrix(if (front) 0L else rows, n, m + 1L, byrow = TRUE)
    state$score <- matrix(if (front) 0L else -2L * rows, n, m + 1L,
        byrow = TRUE)
    state$origin <- matrix(if (front) -rows else 0L, n, m + 1L, byrow = TRUE)
    # Ukkonen's trick: rows beyond `last` cannot be within k errors.
    state$last <- rep(if (front) m else min(m, state$k + 1L), n)
    state$unset <- state$widths + m + 1L
    state$best_cost <- state$unset
    state$best_score <- integer(n)
    state$best_origin <- integer(n)
    state$best_ref_stop <- rep(m, n)
    state$best_query_stop <- state$widths
    # cutadapt's search of the last column compares the `origin` last assigned
    # in the column loop, not that of the cell being tested; kept as is.
    state$origin_variable <- integer(n)
    state$last_filled <- integer(n)
    state$searching <- state$widths > 0L
    state
}

# Offer cells as the best match; returns the reads whose best match changed.
.aligner_offer <- function(state, index, origin, cost, score, ref_stop,
    query_stop, compared_origin = origin) {
    length <- ref_stop + pmin(origin, 0L)
    effective <- length
    if (state$wildcard) {
        ref_start <- -pmin(origin, 0L)
        effective <- ifelse(length < state$m, length -
            (state$n_counts[ref_stop + 1L] - state$n_counts[ref_start + 1L]),
            state$effective_length)
    }
    best_origin <- state$best_origin[index]
    best_length <- state$best_ref_stop[index] + pmin(best_origin, 0L)
    higher <- score > state$best_score[index]
    update <- length >= min(.PRIMER_MIN_OVERLAP, state$m) &
        cost <= effective * .PRIMER_ERROR_RATE &
        (state$best_cost[index] == state$unset[index] |
            (compared_origin <= best_origin + state$m %/% 2L & higher) |
            (length > best_length & higher))
    chosen <- index[update]
    state$best_score[chosen] <- score[update]
    state$best_cost[chosen] <- cost[update]
    state$best_origin[chosen] <- origin[update]
    state$best_ref_stop[chosen] <- ref_stop
    state$best_query_stop[chosen] <- query_stop[update]
    update
}

.aligner_fill_column <- function(state, j, active) {
    query <- state$codes[, j]
    diag_cost <- state$cost[, 1L]
    diag_score <- state$score[, 1L]
    diag_origin <- state$origin[, 1L]
    state$origin[active, 1L] <- state$origin[active, 1L] + 1L
    for (i in seq_len(state$m)) {
        index <- which(active & state$last >= i)
        if (!length(index)) {
            break
        }
        equal <- bitwAnd(state$reference[[i]], query[index]) != 0L
        current_cost <- state$cost[index, i + 1L]
        current_score <- state$score[index, i + 1L]
        current_origin <- state$origin[index, i + 1L]
        cost_diag <- diag_cost[index] + 1L
        cost_insertion <- current_cost + 1L
        cost_deletion <- state$cost[index, i] + 1L
        mismatch <- cost_diag <= cost_deletion & cost_diag <= cost_insertion
        deletion <- !mismatch & cost_deletion <= cost_insertion
        state$cost[index, i + 1L] <- ifelse(equal, diag_cost[index],
            ifelse(mismatch, cost_diag,
                ifelse(deletion, cost_deletion, cost_insertion)))
        state$score[index, i + 1L] <- ifelse(equal, diag_score[index] + 1L,
            ifelse(mismatch, diag_score[index] - 1L, ifelse(deletion,
                state$score[index, i] - 2L, current_score - 2L)))
        state$origin[index, i + 1L] <- ifelse(equal | mismatch,
            diag_origin[index],
            ifelse(deletion, state$origin[index, i], current_origin))
        diag_cost[index] <- current_cost
        diag_score[index] <- current_score
        diag_origin[index] <- current_origin
    }
    computed <- which(active & state$last >= 1L)
    state$origin_variable[computed] <- state$origin[cbind(computed,
        state$last[computed] + 1L)]
    state$last_filled[active] <- state$last[active]
}

# Move `last` to the last row within k errors; where that is the full primer,
# offer the match ending in column j. An exact match ends the search.
.aligner_close_column <- function(state, j, active) {
    m <- state$m
    within <- rep(-1L, length(active))
    for (i in seq.int(0L, m)) {
        within[state$last >= i & state$cost[, i + 1L] <= state$k] <- i
    }
    state$last[active] <- within[active]
    grow <- active & state$last < m
    state$last[grow] <- state$last[grow] + 1L
    at_end <- which(active & !grow)
    if (length(at_end)) {
        cost <- state$cost[at_end, m + 1L]
        origin <- state$origin[at_end, m + 1L]
        state$origin_variable[at_end] <- origin
        update <- .aligner_offer(state, at_end, origin, cost,
            state$score[at_end, m + 1L], m, rep(j, length(at_end)))
        state$searching[at_end[update & cost == 0L & origin >= 0L]] <- FALSE
    }
}

# After the last column of a read: partial 3' primers at the read end.
.aligner_search_last_column <- function(state, j) {
    ending <- which(state$searching & state$widths == j)
    for (i in rev(seq.int(if (state$front) state$m else 0L, state$m))) {
        index <- ending[state$last_filled[ending] >= i]
        if (length(index)) {
            .aligner_offer(state, index, state$origin[index, i + 1L],
                state$cost[index, i + 1L], state$score[index, i + 1L], i,
                state$widths[index],
                compared_origin = state$origin_variable[index])
        }
    }
    state$searching[ending] <- FALSE
}

# Port of Aligner.locate for many reads at once. Returns the 0-based start
# and the end (exclusive) of each match in the read, NA without a match.
.cutadapt_locate <- function(sequences, primer, front) {
    state <- .aligner_state(sequences, primer, front)
    for (j in seq_len(ncol(state$codes))) {
        active <- state$searching & state$widths >= j
        if (!any(active)) {
            break
        }
        .aligner_fill_column(state, j, active)
        .aligner_close_column(state, j, active)
        .aligner_search_last_column(state, j)
    }
    found <- state$best_cost != state$unset
    cbind(start = ifelse(found, pmax(state$best_origin, 0L), NA_integer_),
        stop = ifelse(found, state$best_query_stop, NA_integer_))
}

# Reads that may hold an acceptable alignment other than their leftmost exact
# occurrence. The tests are supersets of cutadapt's conditions, so a read
# left out has no other alignment the aligner could choose.
.primer_needs_alignment <- function(sequences, primer, front, exact_start) {
    m <- nchar(primer)
    k <- as.integer(.PRIMER_ERROR_RATE * m)
    pattern <- Biostrings::BString(primer)
    approximate <- function(subjects, pattern, errors) {
        Biostrings::vcountPattern(pattern, Biostrings::BStringSet(subjects),
            max.mismatch = errors, with.indels = TRUE, fixed = TRUE) > 0L
    }
    needs <- rep(FALSE, length(sequences))
    exact <- !is.na(exact_start)
    # A full-length alignment starting more than m %/% 2 bases before the exact
    # one is recorded first and is not replaced by it.
    prefix_length <- exact_start - 2L - m %/% 2L + m + k
    check <- which(exact & exact_start - 2L - m %/% 2L >= 0L)
    if (length(check)) {
        needs[check] <- approximate(substr(sequences[check], 1L,
            prefix_length[check]), pattern, k)
    }
    residual <- which(!exact)
    if (!length(residual)) {
        return(needs)
    }
    reads <- sequences[residual]
    candidate <- approximate(reads, pattern, k)
    # Partial primers: the end of a 5' primer at the read start, the start of
    # a 3' primer at the read end.
    overlaps <- seq_len(m - 1L)
    for (overlap in overlaps[overlaps >= min(.PRIMER_MIN_OVERLAP, m)]) {
        errors <- floor(overlap * .PRIMER_ERROR_RATE)
        part <- if (front) {
            substring(primer, m - overlap + 1L)
        } else {
            substr(primer, 1L, overlap)
        }
        candidate <- candidate | if (errors == 0) {
            if (front) startsWith(reads, part) else endsWith(reads, part)
        } else {
            window <- if (front) {
                substr(reads, 1L, overlap + errors)
            } else {
                substring(reads, pmax(nchar(reads) - overlap - errors + 1L, 1L))
            }
            approximate(window, Biostrings::BString(part), errors)
        }
    }
    needs[residual] <- candidate
    needs
}

# 1-based start and end of the primer in each read, NA without a match.
.primer_matches <- function(sequences, primer, front) {
    coordinates <- matrix(NA_integer_, nrow = length(sequences), ncol = 2L,
        dimnames = list(NULL, c("start", "end")))
    if (!length(sequences)) {
        return(coordinates)
    }
    normalised <- .normalise_adapter_alphabet(sequences)
    primer <- chartr("I", "N", .normalise_adapter_alphabet(primer))
    m <- nchar(primer)
    align <- rep(TRUE, length(normalised))
    if (!grepl("[^ACGT]", primer)) {
        hits <- Biostrings::vmatchPattern(Biostrings::BString(primer),
            Biostrings::BStringSet(normalised), fixed = TRUE)
        counts <- S4Vectors::elementNROWS(hits)
        exact_start <- rep(NA_integer_, length(normalised))
        first <- cumsum(counts) - counts + 1L
        starts <- unlist(Biostrings::start(hits), use.names = FALSE)
        exact_start[counts > 0L] <- starts[first[counts > 0L]]
        align <- .primer_needs_alignment(normalised, primer, front, exact_start)
        resolved <- !is.na(exact_start) & !align
        coordinates[resolved, "start"] <- exact_start[resolved]
        coordinates[resolved, "end"] <- exact_start[resolved] + m - 1L
    }
    to_align <- which(align)
    if (length(to_align)) {
        distinct <- unique(normalised[to_align])
        located <- do.call(rbind, lapply(
            split(distinct, ceiling(seq_along(distinct) / .PRIMER_DP_CHUNK)),
            .cutadapt_locate, primer = primer, front = front))
        position <- match(normalised[to_align], distinct)
        coordinates[to_align, "start"] <- located[position, "start"] + 1L
        coordinates[to_align, "end"] <- located[position, "stop"]
    }
    coordinates
}
