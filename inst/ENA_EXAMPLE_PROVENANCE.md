# Offline ENA inspection example

`extdata/ena-inspection.rds` is the result of
`inspect_selexprep_accession("PRJDB19098")`, retrieved on 2026-11-03.
It contains the six public FGF-9 HT-SELEX run records DRR615208–DRR615213,
BioSample SAMD00829157, the custom sample attributes returned by ENA, and
FASTQ download locations/checksums. It contains no sequencing reads.

- Study: https://www.ebi.ac.uk/ena/browser/view/PRJDB19098
- Run metadata endpoint: https://www.ebi.ac.uk/ena/portal/api/filereport
- Sample metadata: https://www.ebi.ac.uk/ena/browser/api/xml/SAMD00829157
- SHA-256: `2aef011c06cea7f179a75165c034bc95739b3f63d9decbd8f224f9e09372405a`
- Rebuild script: `scripts/make-ena-example.R` (requires network access).

The snapshot supports runnable download-planning examples without contacting
ENA or transferring FASTQ files. It is a dated example, not an assertion of
current read availability. Query ENA again before downloading a study.

The submitted factual run metadata is made available through ENA/INSDC under
its open-access policy and EMBL-EBI terms. It is not labelled CC0 by this
package. See DATA_LICENSES.md for the policy links and the distinction between
metadata and third-party sequencing data.
