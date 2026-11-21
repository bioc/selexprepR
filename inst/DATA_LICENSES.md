# Data licenses and attributions

This inventory distinguishes the license of package-authored data and
annotations from the terms of materials retrieved from other repositories.
Source policies were consulted on 2026-11-03.

## Frozen public catalog

Files: `extdata/selexprep_catalog.csv` and
`data/selexprep_public_catalog.rda` in the source package.

The catalog compilation and curation annotations are distributed under the
upstream MIT grant, copyright 2026 Marco Rota Negroni. The R dataset is a
format conversion of the CSV with package-specific provenance metadata. The
CSV remains byte-identical to Python selexprep v0.4.6.

The upstream grant is retained verbatim in `licenses/selexprep-MIT.txt`:
<https://archive.softwareheritage.org/browse/revision/63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5/?path=LICENSE>.

The distributed flat table contains identifiers, study titles and short
factual summaries of experimental parameters, together with curation status.
It does not include the underlying papers, figures, supporting-information
files, evidence quotations or sequencing reads. Its record-level evidence
sources and curation process are described in `CATALOG_PROVENANCE.md`.

Policies for the source repository metadata:

- Figshare.com metadata: CC0 1.0, as specified in its
  [copyright and license policy](https://help.figshare.com/article/copyright-and-license-policy).
- Zenodo metadata: CC0, except email addresses, under its
  [general policies](https://about.zenodo.org/policies/).
  No email addresses are included in this catalog.
- ENA/INSDC records: public open-access records under the
  [INSDC principles](https://www.insdc.org/about-insdc/) and
  [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/).
  EMBL-EBI adds no restrictions beyond those of the contributors; this is
  not represented as a blanket CC0 license for all contributed materials.

The metadata policies do not replace the licenses of files downloadable from
those services. This package does not redistribute such files through the
catalog. Consult the originating deposit's file license before obtaining or
redistributing its source material.

## Offline ENA metadata example

File: `extdata/ena-inspection.rds`.

This contains the public run and sample metadata of PRJDB19098 retrieved
from ENA on 2026-11-03. No reads are included. Its use and redistribution
follow the ENA/INSDC policies linked above. The study attribution, endpoints,
retrieval date, checksum and rebuild script are in
`ENA_EXAMPLE_PROVENANCE.md`.

## Real FASTQ examples

Files: `extdata/cs_round_03.fastq.gz`, `cs_round_05.fastq.gz`, and
`cs_round_07.fastq.gz`; source details in `extdata/cs-source-metadata.json`.

These are 600-record R1 subsets of Zibin Zhao’s hydrocortisone DL-SELEX data,
Zenodo 14272647, under CC BY 4.0. They retain that license, separately from
the MIT software/catalog license. The attribution, source URLs, original
file sizes, subsetting/recompression changes and checksums are documented in
`REAL_DATA_PROVENANCE.md`, together with the Creative Commons license and
legal-code links. Only file prefixes were acquired: full upstream MD5 values
are labelled unverified, while prefix and subset hashes were checked.
Repository fields in the bundled metadata summary derive from Zenodo
metadata under CC0; primer-source annotations also cite Table S8 of the
associated publication, as detailed in REAL_DATA_PROVENANCE.md.

## Regression fixtures and illustrative data

Files under `tests/testthat/fixtures/` contain synthetic regression inputs
and reference outputs authored for selexprep/selexprepR. These are covered by
the MIT grant, copyright 2026 Marco Rota Negroni. They do not contain copied
publication text or experimental sequencing reads. Commit-pinned provenance
is documented in their README files.

Any barcodes added to real sequences solely for the demultiplexing example
are illustrative package-authored data; the original sequencing subset
continues to carry its source license. The vignette labels this transformation.
