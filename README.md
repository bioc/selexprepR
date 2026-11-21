# selexprepR

`selexprepR` is an R/Bioconductor package for reproducible preprocessing of
high-throughput SELEX reads. It infers constant regions, extracts variable
regions, creates sparse multi-round count experiments, records provenance, and
reports quality-control signals.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("selexprepR")
```

The main workflow is deliberately small:

```r
library(selexprepR)

experiment <- run_selexprep_files(c(
    round_00 = "round_00.fastq.gz",
    round_01 = "round_01.fastq.gz"
))
```

Inline-barcoded pools can be split with an explicit barcode-to-round map:

```r
inputs <- selexprep_demultiplex(
    read_selexprep_fastq("multiplexed.fastq.gz"),
    c(AAAAA = 0L, TTTTT = 1L)
)
experiment <- run_selexprep(inputs)
```

`experiment` is a `SummarizedExperiment` with a sparse `counts` assay. Its
metadata retains the validated library report, extraction summary, QC report,
and a versioned reproducibility manifest. See the package vignette for the full
workflow.

The bundled study catalog is both a package dataset and a queryable table:

```r
data("selexprep_public_catalog", package = "selexprepR")
selexprep_catalog("FGF-9")
```

It works offline and carries snapshot provenance in
`S4Vectors::metadata(selexprep_public_catalog)`. ENA requests are explicit: use
`inspect_selexprep_accession()` before calling `fetch_selexprep_reads()`; the
latter defaults to a checksum-aware dry-run plan.

The R implementation is aligned to Python `selexprep v0.4.8`
([doi:10.5281/zenodo.24427025](https://doi.org/10.5281/zenodo.24427025)), with
independent R versioning. Shared constants are checked across rounds;
extraction warns when fewer than half the reads survive, with separate FASTQ
diagnostics for file workflows. `HIGH` and `MEDIUM` describe read support
rather than a verified constant/random boundary. Consult the study's library
design when a diverse early round is missing.

## With your own reads

Look at what detection called before counting anything:

```r
report <- S4Vectors::metadata(experiment)$library_report
report[c("status", "extraction_mode", "orientation", "primer_5p", "primer_3p",
    "n_length_mode")]
```

`primer_5p`, `primer_3p` and `n_length_mode` (the most common random-region
length) should match your library design. If your reads are longer than the
construct they run into the sequencing adapter: detection can then call the
adapter as a primer, with a warning, or find one constant only. Give your own
primers, written as they read in the FASTQ from left to right, and the mode:

```r
experiment <- run_selexprep_files(
    files,
    primer_5p = "ATCCAGAGTGACGCAGCATTGC",
    primer_3p = "TGGACACGGTGGCTTAGTCAA",
    extraction_mode = "BOTH_PRIMERS_SINGLE_READ"
)
```

The vignette section "Start with your own reads" covers paired-end files,
reads in both orientations, replicates and export.

## Citation

`citation("selexprepR")` gives the references for the package and for the
selexprep method. Primer location is a port of cutadapt, which should be cited
with it.

## Development status

selexprepR is in the Bioconductor development branch. Development of the R
package is hosted at <https://github.com/marcorotanegroni/selexprepR>; the
original Python implementation remains available at
<https://github.com/marcorotanegroni/selexprep>.

## Example data and redistribution

The vignette uses small real sequencing subsets and includes executable
FASTQ, demultiplexing and QC examples. File-level sources, licenses and
subsetting details are distributed with the data. See
[the data-license inventory](inst/DATA_LICENSES.md) and
[catalog provenance](inst/CATALOG_PROVENANCE.md).

The catalog retains its upstream curation statuses. Human spot checks are
recorded separately; agreement between the two extraction models is not
presented as human validation.
