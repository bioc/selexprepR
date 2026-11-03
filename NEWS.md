# selexprepR 0.99.7

* Replace the identical synthetic-round vignette with licensed, attributed
  600-read R1 subsets of three real hydrocortisone DL-SELEX rounds. Show
  safe inference failure, published Table S8 flanks, explicit primer
  overrides, extraction yields, sparse counts, QC plots and demultiplexing.
* Add package-level help and runnable offline examples to all help pages.
  FASTQ examples create their inputs; a cached ENA inspection supports
  offline download planning. Only fresh ENA inspection requires a network.
* Vectorize exact adapter matching, indel-aware candidate filtering and
  strand evidence with Biostrings. Retain original approximate-match ranking
  for matching residual reads, including paired retention masks. Batch the
  reverse-complement canonicalization of distinct QC k-mers.
* Use one reverse-complement implementation for counting and extraction,
  preserving the public DNA/RNA, case and ambiguity policy and fixing the
  internal failure on other nucleotide symbols.
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
