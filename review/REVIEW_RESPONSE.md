# Draft response to Bioconductor review

For the existing issue:
https://github.com/Bioconductor/BiocContributions/issues/90#issuecomment-6596200402

Maintainer publication checklist: publish the local 0.99.7 candidate only after
reviewing the changes; add supported-environment CI/build links below; update
the ten-record human audit if completed. This file has not been posted.

---

Thank you for the detailed review. The 0.99.7 candidate addresses the
requested documentation, implementation and data-provenance changes.

## Vignette and examples

The vignette now uses 600 complete R1 reads from each of three real
hydrocortisone DL-SELEX rounds (3, 5 and 7), from Zenodo 14272647. The source
reads, identifiers and quality strings are preserved. The source license,
attribution, subsetting/recompression method and checksums are documented in
`inst/REAL_DATA_PROVENANCE.md`.

It shows FASTQ import, round-aware inputs, the library report, explicit primer
overrides based on the published Supplementary Table S8, sparse counts,
file-level provenance, QC tables and a QC plot. The sample has mixed
orientations and adapter read-through: its safe inference failure, retained
forward-read subset and QC warnings are shown and explained. The small
convenience subset is not presented as evidence of biological enrichment.

Demultiplexing runs offline on these real sequences with explicitly artificial
teaching barcodes; removing the barcodes recovers the original round pools
and count matrix. No vignette chunks are disabled. Installation and live ENA
operations are described in prose with their external requirements.

Package-level help has been added. All 25 help pages contain executable
offline examples, including `run_selexprep_files()`, which creates its input
FASTQ. A bundled ENA metadata snapshot makes download-planning examples
runnable offline. A fresh `inspect_selexprep_accession()` call is the sole
optional `donttest` block, explicitly documented as requiring live ENA access.

## R implementation

Counting and extraction now share a vectorized Biostrings-based reverse
complement. The public Python-compatible case, U/T and ambiguity policy is
preserved; `reverseComplement(DNAStringSet)` would change case and complement
IUPAC codes under that contract. The former internal duplicate has been
removed, including its failure on other nucleotide symbols.

Adapter matching now uses native vectorized exact matching and indel-aware
candidate filtering. Only matching distinct residual reads need the original
coordinate/ranking calculation, because Biostrings does not vectorize
indel-aware match coordinates. Strand evidence is also vectorized.

The long workflows were simplified by sharing repeated flank-scoring and QC
threshold logic, mate pooling, barcode preprocessing and fetch-result
construction. Validation, extraction and assembly helpers have explicit
responsibilities. Error messages, RNG behavior, ordering, schemas and paired
retention were checked against 0.99.6. A separate regression fixes the
pre-existing lowercase-RNA barcode normalization inconsistency.

A local illustrative benchmark used 3,000 unique reads per scenario and
checked complete extraction-object equality. Exact/mixed/no-match scenarios
improved by 7.63x, 4.33x and 3.89x respectively. These are single-run local
measurements, not performance guarantees; the script and environment are
recorded under `review/` and in `inst/REVIEW_VALIDATION.md`.

## Bundled-data redistribution

`inst/DATA_LICENSES.md` inventories the data and their terms. The upstream
MIT grant for the catalog compilation and annotations is retained verbatim;
the CSV remains pinned to Python v0.4.6. The real FASTQ subsets retain the
source CC BY 4.0 license and attribution. ENA example metadata and synthetic
regression fixtures are documented separately. No source publication,
figure, supplement or complete upstream FASTQ file is bundled.

## Human catalog checks and AI assistance

The maintainer manually spot-checked five randomly selected catalog records
against their original sources and found no discrepancies in the fields
checked. Those earlier identifiers were not recorded. This was a limited
spot check. The original field-level curation status and source provenance
remain available; model agreement is not described as human verification.

A further reproducibly drawn ten-record sample and its source locations are
prepared in `review/catalog-spot-check-2026-11-03.md`. Its human results are
currently pending and must not be reported as completed. Codex assistance,
including parallel implementation and audit agents, is disclosed in
`inst/AI_USAGE.md`.

## Validation and build reports

Local source-mode verification passed 115 test blocks and 771 assertions,
all 25 offline help examples and all vignette chunks. The package retains
R >= 4.6.0; supported-environment limitations are recorded in
`inst/REVIEW_VALIDATION.md`.

Before posting, add links to the 0.99.7 personal CI and BiocStaging build
reports after they have run; do not describe the previous 0.99.6 report as a
check of this candidate.
