# selexprepR 0.99.10

* Gabriele Sales is listed as a contributor. The affiliation,
  Department of Biology, University of Padova, is recorded in `DESCRIPTION`.
  The code does not change.

# selexprepR 0.99.9

* Aligned to Python selexprep v0.4.8
  ([doi:10.5281/zenodo.24427025](https://doi.org/10.5281/zenodo.24427025)).
  No call, status or extracted read changes. This version is about running
  the package on your own reads.
* `selexprep_detect()` warns when a called primer contains a sequencing
  adapter, which happens when the reads run past the construct, and prints a
  message when it finds one constant only, saying what the sequences will
  keep.
* Extraction warns when both primers are set but the mode is one-sided, so
  that one primer is left untrimmed, and names
  `extraction_mode = "BOTH_PRIMERS_SINGLE_READ"`. It warns when a supplied
  primer is in under 5% of the earliest round's reads, and says so when its
  reverse complement is there instead (a reverse primer typed as ordered).
* `run_selexprep()` warns when the primers are split between R1 and R2 and no
  sequence is counted. `read_selexprep_inputs()` warns about FASTQ files with
  no reads and names a path that does not exist.
* Supplied primers may carry IUPAC ambiguity codes, as the help of
  `selexprep_extract()` already said. An invalid `extraction_mode` lists the
  valid values.
* The help pages give the values of `extraction_mode`, `status` and
  `required_action`, say how to write the primers, and explain why a report
  with supplied primers has `status = "LOW"`. The vignette and the README
  have a guide for your own reads: reads longer than the construct, supplied
  primers, replicates, run time, and counting both read orientations.
* Catalog provenance links point to the archived release
  ([doi:10.5281/zenodo.23104868](https://doi.org/10.5281/zenodo.23104868))
  and to Software Heritage permalinks of the source commit. The dataset and
  its metadata are unchanged.
* Add `inst/CITATION`: the package, the selexprep method and cutadapt.

# selexprepR 0.99.8

* Credit cutadapt for the aligner ported in `R/primer-match.R`: its MIT
  license and copyright notice are in `inst/licenses/cutadapt-MIT.txt`,
  summarized in `inst/COPYRIGHTS`, and Marcel Martin is listed as copyright
  holder in `DESCRIPTION`.

# selexprepR 0.99.7

* Replace the identical synthetic-round vignette with licensed, attributed
  600-read R1 subsets of three real hydrocortisone DL-SELEX rounds. Show
  safe inference failure, published Table S8 flanks, explicit primer
  overrides, extraction yields, sparse counts, QC plots and demultiplexing.
* Add a vignette guide for starting from your own reads. Explain the chosen
  read orientation, paired-end inputs and each QC flag in the example, and
  export one count table per round with `write_selexprep_counts()`.
* Show the QC retention, length and diversity plots as separate figures.
  Length fractions omit lengths not observed in a round, and the retention
  plot leaves room for its legend.
* Add package-level help and runnable offline examples to all help pages.
  FASTQ examples create their inputs; a cached ENA inspection supports
  offline download planning. Only fresh ENA inspection requires a network.
* Locate primers as cutadapt 5.2 does in the Python pipeline. A primer cut
  short at the read start (5') or end (3') is accepted when at least 3 bases
  align; among several occurrences the first acceptable one is kept unless a
  later one overlaps it or covers more of the primer with a higher score; an
  error at a primer edge is a substitution rather than a deletion. Synthetic
  edge cases and the bundled real reads are now trimmed exactly as cutadapt
  trims them; the example retains 228, 278 and 297 reads instead of 225, 277
  and 294.
* Vectorize exact primer matching and strand evidence with Biostrings; only
  reads with a possible alternative match are aligned. Batch the
  reverse-complement canonicalization of distinct QC k-mers.
* Use one reverse-complement implementation for counting and extraction,
  fixing the internal failure on other nucleotide symbols.
  `reverse_complement()` now complements IUPAC ambiguity codes (R/Y, K/M,
  B/V, D/H) while keeping case, RNA input and other symbols.
* Measure match rates and insert lengths for manually supplied primers on
  the earliest round, so the library report and the `low_primer_match` QC
  flag describe the primers actually used for extraction.
* Report sides without a candidate flank as such in detection failures,
  instead of a 0.00 match rate.
* Name reads in the opposite orientation as a possible cause in low-yield
  extraction warnings.
* Describe the QC summaries and flags in `selexprep_qc()` as computed.
* Simplify shared flank evidence, QC thresholds, mate pooling, barcode
  comparisons and verified downloads; preserve public result schemas, RNG
  behavior and failure diagnostics while shortening workflow functions.
* Normalize lowercase RNA barcodes and reads consistently with uppercase
  RNA and DNA, reusing the common alphabet-normalization helper.
* Document redistribution terms separately for the MIT catalog compilation,
  CC BY 4.0 FASTQ subsets, public metadata and synthetic regression fixtures.
  Preserve the Python v0.4.6 catalog values and original curation statuses;
  record human spot checks separately from model agreement.

# selexprepR 0.99.6

* Simplify ENA condition messages to address the two new BiocCheck coding
  practice notes, preserving four-attempt retry and sample-metadata fallback.
* Keep diagnostic accessions, underlying sample failures and retry advice
  explicit; remove redundant failure wording from the filereport message.

# selexprepR 0.99.5

* Align library inference with Python selexprep v0.4.6: retain shared library
  constants across rounds using per-round offsets, preserve the insert boundary,
  and avoid ambiguous repeated-core anchors.
* Warn about extraction pools retaining fewer than half their reads or half
  the best pool's yield, retaining diagnostics in `low_yield_rounds`. File
  workflows also retain `input_yields` and `low_yield_inputs` per FASTQ, so
  aggregation cannot mask a poorly matching input.
* Normalize lowercase RNA consistently during detection and primer matching.
* Read ENA BioSample attributes for round inference; handle conflicting and
  weak PCR-cycle evidence, and stop interpreting ArrayExpress RV sample labels
  as selection rounds.
* Retry ENA filereport service-error bodies, retain run-level library strategy
  and sample evidence, and warn when a round pools distinct BioSamples.
* Pin the 240-record catalog and provenance to Python v0.4.6, including the
  corrected PFOA target class. Add offline Python/R regression references.

# selexprepR 0.99.4

* Refuse inferred libraries with a zero-length random region instead of
  reporting a valid full-insert extraction.
* Rescue well-supported primer cores when noisy outer flanks lower the
  full-length match rate, without changing the extraction boundary.
* Update the bundled public catalog to the adjudicated Python v0.4.1
  snapshot, retaining immutable source and evidence provenance.

# selexprepR 0.99.3

* Write LibraryReport and manifest JSON with UTF-8 LF line endings on every
  platform, preserving deterministic hashes and Python golden-file parity on
  Windows.

# selexprepR 0.99.1

* Initial Bioconductor submission candidate.
* Provides primer inference, extraction, sparse multi-round counting, quality
  control, reproducibility manifests, FASTQ input, and ENA data discovery.
* Adds `selexprep_public_catalog` as a documented, reproducibly generated
  package dataset with immutable snapshot provenance.
* Adds R-native barcode demultiplexing and direct FASTQ/fetch-result input
  bridges with file hashes retained in the run manifest.
