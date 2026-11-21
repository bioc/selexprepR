#' Run the in-memory SELEX preprocessing workflow
#'
#' Composes primer detection, extraction, sequence counting and quality control
#' into a single Bioconductor result. The primary output is a
#' `SummarizedExperiment` with one column per round and a sparse `counts`
#' assay. The LibraryReport, extraction result and QC report are retained in
#' `metadata()` for reproducibility and inspection.
#'
#' Split-primer paired-end inputs return a valid zero-row experiment, with a
#' warning, because full inserts cannot be counted without read merging. Their
#' separate trimmed sides remain available in `metadata(result)$extraction`,
#' and `metadata(result)$qc` is `NULL`.
#'
#' Library status rates read support, not a verified constant/random boundary.
#' Check the published library design when diverse early rounds are absent.
#' When reads are longer than the construct they run into the sequencing
#' adapter: detection can then call the adapter as a primer, or find one
#' constant only. Both cases are reported, the first with a warning; supply
#' `primer_5p`, `primer_3p` and `extraction_mode` to extract the random region.
#' With supplied primers the report has `status = "LOW"` and `confidence = 0`:
#' nothing was inferred, so there is no inference to rate. The evidence for the
#' supplied primers is in `match_rate_5p` and `match_rate_3p`.
#' Extraction warns about low read retention; diagnostics remain in
#' `metadata(result)$extraction$low_yield_rounds`. With FASTQ input bundles,
#' per-file counts and warnings also remain in `input_yields` and
#' `low_yield_inputs`, so pooling files cannot hide a poorly matching input.
#'
#' @param sequences_by_round A named list of R1 pools, as accepted by
#'   `selexprep_detect()`, or a `selexprep_read_inputs` object.
#' @param read_source Input read layout.
#' @param paired_mate_streams Optional named R2 pools.
#' @param accession Optional public-study accession stored in metadata.
#' @param sampling_seed Seed forwarded to `selexprep_detect()`.
#' @param max_reads_per_round Optional detection subsampling cap.
#' @param low_total_reads QC threshold for per-round read depth.
#' @param primer_5p Optional manually reviewed 5-prime primer override,
#'   written as it reads at the start of the reads.
#' @param primer_3p Optional manually reviewed 3-prime primer override,
#'   written as it reads at the end of the reads: a reverse primer as ordered
#'   is its reverse complement.
#' @param extraction_mode Optional manually reviewed extraction-mode override:
#'   `"BOTH_PRIMERS_SINGLE_READ"`, `"FIVE_PRIME_ONLY"`, `"THREE_PRIME_ONLY"` or
#'   `"PAIRED_END_SPLIT_PRIMERS"`. A primer override keeps the inferred mode,
#'   so give `"BOTH_PRIMERS_SINGLE_READ"` with both primers when detection
#'   reported a one-sided mode.
#'
#' @return A `SummarizedExperiment` with the LibraryReport, extraction and QC
#'   result in its metadata.
#' @export
#' @examples
#' p5 <- 'GGTAATACGACTCACTATAGGG'
#' p3 <- 'CCATGCATGCATGCATGCAT'
#' inserts <- c('ACGTACGTACGTACGT', 'TGCATGCATGCATGCA')
#' pools <- list(round_00 = rep(paste0(p5, inserts, p3), each = 250))
#' run_selexprep(pools, low_total_reads = 0)
run_selexprep <- function(sequences_by_round,
    read_source = c("R1", "R2", "R1_AND_R2", "INTERLEAVED", "UNKNOWN"),
    paired_mate_streams = NULL, accession = NULL, sampling_seed = 42L,
    max_reads_per_round = NULL, low_total_reads = 10000L,
    primer_5p = NULL, primer_3p = NULL, extraction_mode = NULL) {
    inputs <- .resolve_run_inputs(
        sequences_by_round, read_source, paired_mate_streams,
        read_source_supplied = !missing(read_source)
    )
    detected <- system.time(report <- .detect_library(
        inputs$sequences_by_round, inputs$read_source,
        inputs$paired_mate_streams, sampling_seed, max_reads_per_round
    ))
    manual_override <- !is.null(primer_5p) || !is.null(primer_3p) ||
        !is.null(extraction_mode)
    report <- .apply_primer_overrides(
        report, primer_5p = primer_5p, primer_3p = primer_3p,
        extraction_mode = extraction_mode,
        sequences_by_round = inputs$sequences_by_round,
        paired_mate_streams = inputs$paired_mate_streams
    )
    .warn_about_inferred_library(
        report,
        supplied = c("5'"[!is.null(primer_5p)], "3'"[!is.null(primer_3p)]),
        mode_supplied = !is.null(extraction_mode)
    )
    extracted <- system.time(extraction <- .extract_run_inputs(inputs, report))
    counted <- system.time(counts <- .run_count_tables(extraction, report))
    experiment <- .as_selexprep_experiment(counts, accession = accession)
    metadata <- S4Vectors::metadata(experiment)
    metadata$library_report <- report
    metadata$extraction <- extraction
    S4Vectors::metadata(experiment) <- metadata
    quality_controlled <- system.time(qc <- if (.run_needs_merging(report)) {
        NULL
    } else {
        selexprep_qc(experiment, report, low_total_reads = low_total_reads)
    })
    manifest <- .build_run_manifest(
        report, accession, inputs$read_source, sampling_seed,
        max_reads_per_round, low_total_reads, manual_override, qc,
        c(detect = unname(detected[["elapsed"]]),
            extract = unname(extracted[["elapsed"]]),
            count = unname(counted[["elapsed"]]),
            qc = unname(quality_controlled[["elapsed"]]))
    )
    .finalise_run_metadata(experiment, qc, manifest, inputs$bundle)
}

.resolve_run_inputs <- function(sequences, read_source, mates,
    read_source_supplied) {
    bundle <- NULL
    if (inherits(sequences, "selexprep_read_inputs")) {
        if (read_source_supplied || !is.null(mates)) {
            stop("Sequence layout is already defined by the input bundle.",
                call. = FALSE)
        }
        bundle <- sequences
        mates <- bundle$paired_mate_streams
        read_source <- bundle$read_source
        sequences <- bundle$sequences_by_round
    } else {
        read_source <- match.arg(read_source, .read_sources)
    }
    list(sequences_by_round = sequences, paired_mate_streams = mates,
        read_source = read_source, bundle = bundle)
}

.extract_run_inputs <- function(inputs, report) {
    files <- inputs$bundle$files
    file_inputs <- !is.null(files) && nrow(files) > 0L
    extraction <- .extract_selexprep(
        inputs$sequences_by_round, report,
        paired_mate_streams = inputs$paired_mate_streams,
        keep_masks = file_inputs
    )
    if (file_inputs) {
        extraction <- .add_input_yields(extraction, files)
    }
    extraction
}

.run_needs_merging <- function(report) {
    identical(report$extraction_mode, "PAIRED_END_SPLIT_PRIMERS")
}

.run_count_tables <- function(extraction, report) {
    if (.run_needs_merging(report)) {
        warning(paste(
            "The primers are split between R1 and R2",
            "(PAIRED_END_SPLIT_PRIMERS): no sequence was counted and the",
            "experiment has zero rows. The trimmed sides are in",
            "metadata(result)$extraction. To count full inserts, merge the",
            "read pairs with a read merger and run the merged reads as",
            "single-end input."
        ), call. = FALSE)
        round_names <- names(extraction$sequences_by_round)
        empty <- lapply(round_names, function(round) {
            selexprep_count(character())
        })
        return(stats::setNames(empty, round_names))
    }
    lapply(extraction$sequences_by_round, selexprep_count)
}

.run_qc_flags <- function(qc, report) {
    if (!is.null(qc)) {
        return(as.character(qc$flags$name))
    }
    if (identical(report$required_action, "READ_MERGING_RECOMMENDED")) {
        "requires_read_merging_for_full_insert"
    } else {
        character()
    }
}

.build_run_manifest <- function(report, accession, read_source, seed,
    max_reads, low_total_reads, manual_override, qc, timings) {
    build_selexprep_manifest(
        report, accession = accession,
        parameters = c(
            read_source = read_source, sampling_seed = as.character(seed),
            max_reads_per_round = if (is.null(max_reads)) {
                "all"
            } else {
                as.character(max_reads)
            },
            low_total_reads = as.character(low_total_reads),
            manual_primer_override = if (manual_override) "true" else "false"
        ),
        runtime_seconds_per_stage = timings,
        flags = .run_qc_flags(qc, report)
    )
}

.finalise_run_metadata <- function(experiment, qc, manifest, bundle) {
    metadata <- S4Vectors::metadata(experiment)
    metadata$qc <- qc
    metadata$manifest <- manifest
    metadata$counting_skipped <- .run_needs_merging(metadata$library_report)
    if (!is.null(bundle) && inherits(bundle, "selexprep_demultiplex")) {
        metadata$demultiplex <- list(
            summary = bundle$summary, barcodes = bundle$barcodes,
            max_mismatches = bundle$max_mismatches,
            trim_barcode = bundle$trim_barcode,
            unassigned_reads = length(bundle$unassigned$r1)
        )
    }
    S4Vectors::metadata(experiment) <- metadata
    experiment
}
