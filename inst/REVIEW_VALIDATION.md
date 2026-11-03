# November 2026 review candidate: selexprepR 0.99.7

Baseline: R 0.99.6, commit c2f7abdc07a12278bb2c459f064f119f0cb7f7e1.
The behavioral reference remains Python selexprep v0.4.6, commit
63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5; later Python development is excluded.
The Python working tree was not modified. The author name and R >= 4.6.0
requirement are preserved.

## Verification scope

Validation was performed in source mode with local R 4.3.1 and Biostrings
2.68.0, below the package's supported R requirement. Source-mode tests and
knitting do not establish a supported-environment package installation or
HTML vignette build. The installed mixed Bioconductor library emits the
pre-existing S4Arrays/DelayedArray import-replacement warning at namespace
load, outside test assertions.

Component checks preceding integration:

- Core detection/QC/run refactor: 37 direct object/warning/RNG comparisons
  against 0.99.6; only measured manifest timings were normalized.
- Shared-constant helper: 350 additional direct equivalence cases.
- Input/demultiplex/download refactor: 207 direct comparisons, including
  error/warning behavior, file content, lane ordering and mocked transfer
  failures. Lowercase-RNA normalization is a separately documented bug fix.
- Extraction: 7,200 adapter-coordinate comparisons; 7,200 trimming/mask
  comparisons; 36 strand-distribution comparisons; 12 complete objects
  across four extraction modes, three orientations and an empty round.
- QC k-mers: 120 direct set/fraction comparisons against 0.99.6, including
  mixed DNA/RNA case, ambiguity and empty inputs. On 1,000 randomly generated
  39-nt sequences, canonical sets were identical (0.464 s before, 0.013 s
  after in a single local measurement). See `review/validate-kmer-qc.R`.
- Existing exact-tag Python report, extraction and round-inference fixtures
  remain part of the full R suite. Ambiguous-base extraction now follows
  the public reverse-complement policy instead of failing in its duplicate.
- All 25 help pages parse and have successfully executed offline examples;
  the sole optional fresh ENA call is excluded because it requires a live
  service. A cached metadata snapshot exercises offline planning.
- All vignette chunks knit offline with errors stopping execution. The QC
  figure was visually inspected. Full BiocStyle HTML rendering still requires
  the destination build environment.
- Three real FASTQ subsets: 600 unchanged R1 records each, source-prefix and
  subset hashes verified. Published flanks verified against Supplementary
  Table S8; the source document's MD5 matches the current PMC Cloud record.
  Full upstream FASTQ MD5 values are explicitly labelled unverified.

## Local extraction timing

Single elapsed-time measurements on 3,000 distinct reads per scenario:

| 5-prime matches: exact / approximate / none | 0.99.6 | Revised | Ratio |
| --- | ---: | ---: | ---: |
| 100% / 0% / 0% | 6.276 s | 0.823 s | 7.63x |
| 80% / 20% / 0% | 6.088 s | 1.405 s | 4.33x |
| 60% / 20% / 20% | 5.287 s | 1.358 s | 3.89x |

Complete extraction objects were identical. These illustrative measurements
are not throughput guarantees. `review/benchmark-extraction.R` reproduces
them from the immutable Git baseline; `review/extraction-differential.R`
reproduces the differential cases. Developer artifacts under `review/` are
available in the Git repository and excluded from the source tarball.

## Final integrated checks

The final source-mode suite was run on the actual repository after integrating
all code changes: **115 test blocks, 771 assertions passed**, with zero
failures, skips, errors or test warnings. The 25 offline help examples and
all vignette chunks (including the QC plot) also passed on that checkout.
The static code-usage scan found no diagnostics apart from the pre-existing
namespace-load warning described above.

| Package check | Local result |
| --- | --- |
| `R CMD build --no-build-vignettes` | Source tarball built; three FASTQ subsets, license inventory and package help included; developer review artifacts excluded. |
| Full `R CMD build` | Stops at vignette-build installation: R 4.3.1 does not meet R >= 4.6.0. |
| `R CMD check --no-manual --no-build-vignettes` | Stops at dependency check: `arrow` and `BiocStyle` unavailable locally. No complete check result. |
| BiocCheck / BiocCheckGitClone | Not installed in this local environment. Run through the configured CI. |

No package requirement was lowered. The 0.99.7 GitHub Actions and BiocStaging
reports are pending publication by the maintainer; the older 0.99.6 report is
not a validation of this candidate. The current destination setup uses
R 4.6.1/Bioconductor devel. The existing submission and staging branch are
retained.

## Human catalog checks

The five manual source checks reported by the maintainer are distinct from
all automated checks above. Ten further randomly drawn records have source
locations and pending outcome fields in `review/`. No new human result is
asserted until the maintainer records it.
