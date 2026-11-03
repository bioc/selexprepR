# Catalog provenance

`selexprep_public_catalog` is a frozen flat view of 240 public HT-SELEX
deposits. It was migrated without modification from
`src/selexprep/catalog/data/curated_metadata.csv` in the original Python
repository at commit `63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5`.

- Snapshot version: `selexprep-v0.4.6`
- Snapshot date: 2026-10-02
- CSV SHA-256:
  `121b599786a718d877bd950660a3a4e0ec178db2f1b54636815d575185c3e0fe`
- Immutable CSV: <https://github.com/marcorotanegroni/selexprep/blob/63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5/src/selexprep/catalog/data/curated_metadata.csv>
- Canonical provenance JSON: <https://github.com/marcorotanegroni/selexprep/blob/63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5/src/selexprep/catalog/data/curated_metadata.json>
- Curation materials: <https://github.com/marcorotanegroni/selexprep/tree/63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5/benchmarks/dual_extraction>

Relative to the v0.4.1 flat catalog, v0.4.6 changes only PRJEB114397
(PFOA)'s target class to `small molecule`; the 240 identifiers and curation
statuses remain unchanged. Its canonical JSON also updates evidence and
concordance notes for PRJEB114397 and PRJNA1481083.

The eight experimental fields were extracted independently by Claude and
Codex/GPT, then reconciled. The paired `_curation` columns preserve whether a
cell was concordant, adjudicated after disagreement, absent, independently
verified, or supplied by one extraction arm. The 47 adjudicated cells have a
single resolved flat value; the canonical JSON preserves both extraction arms,
the resolution rationale, and record-level evidence sources and locations.

## Redistribution terms

The catalog compilation and curation annotations are distributed under the
MIT license of the upstream selexprep repository, copyright 2026 Marco Rota
Negroni. The full upstream notice is retained in `licenses/selexprep-MIT.txt`;
the immutable source license is
<https://github.com/marcorotanegroni/selexprep/blob/63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5/LICENSE>.
This applies to both the CSV and the R DataFrame derived from it.

The flat catalog contains identifiers, study titles and short factual
summaries of experimental parameters. It does not bundle source articles,
figures, full text, supporting-information files, evidence quotations, or
sequencing reads. Source metadata policies and the separate terms applying
to downloaded source materials are listed in `DATA_LICENSES.md`.

## Human spot checks

On 2026-11-03 the maintainer reported having manually checked five randomly
selected records against the original sources, finding no discrepancies in
the fields checked. The record identifiers and checked fields were not
recorded in this package. This is a limited spot check, not validation of all
240 records. The original curation statuses are preserved.

A further reproducible sample of ten records was drawn on 2026-11-03 for
manual checking. The sample, source locations and pending result fields are
in the source repository's `review/` directory. Its human results are pending;
AI cross-checks and software regression tests are not counted as human review.

Inclusion in the frozen catalog does not assert current read availability.
