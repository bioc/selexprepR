#' Reproducible preprocessing of public HT-SELEX data
#'
#' selexprepR connects raw sequencing reads to a sparse
#' `SummarizedExperiment`, with explicit primer inference, extraction,
#' quality-control and provenance records.
#'
#' @section Workflows:
#' Use [run_selexprep_files()] for local FASTQ files, or [run_selexprep()]
#' for named lists of sequence pools. [selexprep_demultiplex()] assigns
#' inline-barcoded reads when a reviewed barcode map is supplied.
#' Inspect `S4Vectors::metadata(result)` for the library report, extraction
#' yields, QC and manifest. Plot the QC report with [plot.selexprep_qc()].
#'
#' @section Individual stages:
#' [selexprep_detect()] infers constant regions. [selexprep_extract()]
#' removes them, and [selexprep_count()] counts recovered sequences.
#' [selexprep_qc()] summarizes depth, diversity, composition and flags.
#' Reviewed primer overrides are available when inference is unsafe.
#'
#' @section Public data:
#' [selexprep_catalog()] searches a frozen discovery catalog.
#' [inspect_selexprep_accession()] queries current ENA metadata, and
#' [fetch_selexprep_reads()] creates a download plan before transferring
#' reads. Network availability and a public accession are required for
#' live ENA inspection; bundled examples also work offline.
#'
#' @section Further information:
#' The vignette, `vignette("selexprep", package = "selexprepR")`, walks
#' through a complete FASTQ workflow with small real data and a QC plot.
#' Data sources and redistribution terms are recorded in the installed
#' `DATA_LICENSES.md` and the corresponding provenance files.
#'
#' @name selexprepR-package
#' @aliases selexprepR
#' @docType package
#' @keywords package
#' @examples
#' selexprep_catalog("FGF-9")
#' selexprep_count(c("ACGT", "ACGT", "GGGG"))
NULL
