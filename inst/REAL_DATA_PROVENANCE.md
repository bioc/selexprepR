# Real FASTQ example: hydrocortisone DL-SELEX

## Source and attribution

Zhao, Zibin (2025). *Structure-enhanced Deep Learning Accelerates Aptamer
Selection for Small Molecule Families like Steroids*. Zenodo,
published 2025-10-18. DOI: <https://doi.org/10.5281/zenodo.14272647>.
Creator ORCID: <https://orcid.org/0000-0002-3121-9131>.
Data collectors credited by the deposit: Zibin Zhao and Haosi Lin.

Record: <https://zenodo.org/records/14272647>.
Metadata API: <https://zenodo.org/api/records/14272647>.
Associated publication: <https://doi.org/10.1093/bib/bbaf680>.
Retrieved 2026-11-03. The record explicitly describes paired FASTQ data from
hydrocortisone DL-SELEX selection rounds 3, 5 and 7.

## License and modifications

The source record metadata declares `license.id = "cc-by-4.0"`.
The three FASTQ subsets remain under Creative Commons Attribution 4.0
International: <https://creativecommons.org/licenses/by/4.0/>.
Legal code: <https://creativecommons.org/licenses/by/4.0/legalcode>.
This license is distinct from the package software license.

Modification by Marco Rota Negroni for this package: retain only the first
600 complete R1 records from each source file and recompress those records.
No read identifier, sequence, separator or quality line was edited. R2 is
not bundled, and the subsets must not be represented as complete paired data.
The illustrative inline barcodes in the vignette are added in memory for
demonstration; they are not part of the deposited experiment or these files.

## Acquisition and integrity

Only the first 1,048,576 compressed bytes of each source file were downloaded,
using HTTP Range bytes 0-1048575. Zenodo returned HTTP 206 with the expected
Content-Range and original full-file size. These prefixes contain more than
600 complete records. The first 600 records were read before the end of the
prefix; incomplete trailing compressed data were never included in a subset.

The full upstream MD5 values below are supplied by Zenodo and were **not
locally verified**, since the full files were not downloaded. SHA-256 values
for the retrieved prefixes and both compressed and uncompressed subsets were
computed locally. `extdata/cs-source-metadata.json` records the exact values,
file URLs, sizes, license and creator. The checks below establish that the
bundled data match the acquired source prefixes; they do not validate the
original experiment.

| Round | Original R1 file | Original bytes | Zenodo full-file MD5 (not locally verified) |
| --- | --- | ---: | --- |
| 3 | `3th-cs-LFL5718_L3_1.fq.gz` | 909798034 | `ce663275300a74b03b6d18d3aa90c487` |
| 5 | `5th-cs-LFL5722_L3_1.fq.gz` | 872260986 | `e31888a485180c23e79cf99db8538fb6` |
| 7 | `7th-cs-LFL5726_L3_1.fq.gz` | 768565752 | `fe5dc9eb7158e56fe0b90c18b8f4efd4` |

| Bundled subset | Records | Compressed SHA-256 |
| --- | ---: | --- |
| `cs_round_03.fastq.gz` | 600 | `1a8725235b55866e40c178dfce7de82f521c76706a97916fc8cf244d871bb4f3` |
| `cs_round_05.fastq.gz` | 600 | `b93a365fd5cd4e52e3311a7e66d246050a34cb6e5d47c9cd7c29dd025447cc29` |
| `cs_round_07.fastq.gz` | 600 | `daa8f6e8f3d76392605b2c90b88768295a0a4c09545a17b2ef69434300a83099` |

## Rebuilding

Run the following from a package source checkout after acquiring the original
files or the documented compressed prefixes into a local directory:

```sh
Rscript inst/scripts/make-real-selex-example.R /path/to/raw /path/to/selexprepR
```

The builder makes no network requests. It verifies either the full source
MD5 or the prefix SHA-256, validates complete FASTQ records, and checks the
uncompressed subset SHA-256 before writing gzip files. Gzip timestamps are
zero. Compression bytes can depend on the R/zlib version; the uncompressed
bytes and their hashes are the reproducible source-content reference.

## Interpretation

The first records are a deterministic convenience subset, not a random or
representative sample of a selection population. No claim about enrichment
or affinity follows from this example. The R1 reads include mixed strand
orientations and sequencing-adapter read-through. The automated detector
returns a safe failure on the bundled raw sample. The vignette shows this
failure, then supplies source-checked primer overrides and displays the
retained fractions and QC flags. It does not silently lower detection
thresholds or filter the source files to make the example appear successful.

## Primer source and chosen orientation

The published Supplementary Table S8 contains aptamer sequences with the
constant flanks marked in lowercase. For example, AI-CSR3CC and AI-CSR5CC
have 5-prime flank `GGCTCTCGGGACGAC` and 3-prime flank
`CTGCGTTCTCCATTCTGG`. These strings were checked directly in the supplementary
DOCX against the original source (document MD5
`80fd8959e0d6d85f854e106a2929b698`).

Source: <https://doi.org/10.1093/bib/bbaf680>, Supplementary Table S8.
PMC supplementary object:
<https://pmc-oa-opendata.s3.amazonaws.com/PMC12713628.1/dl-selex-si-final_bbaf680.docx?md5=80fd8959e0d6d85f854e106a2929b698>.
The complete supplementary document is not distributed with this package.

The example targets the reverse-complement orientation of those printed
aptamers: `primer_5p = reverse_complement(published_3p)` gives
`CCAGAATGGAGAACGCAG`, and
`primer_3p = reverse_complement(published_5p)` gives `GTCGTCCCGAGAGCC`.
This is an explicit choice to demonstrate preprocessing one orientation of
the raw R1 pool. It does not recover or reorient every read in the deposit.
The extraction counts, retention warnings and strand-mix flag retain that
limitation in the result.
