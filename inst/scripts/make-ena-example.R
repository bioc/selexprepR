# Run with an installed selexprepR and network access from the source root.
# This script refreshes metadata only; it does not download sequencing reads.
inspection <- selexprepR::inspect_selexprep_accession("PRJDB19098")
saveRDS(
    inspection, file = "inst/extdata/ena-inspection.rds",
    compress = "xz", version = 3L
)
# Record the retrieval date and new checksum in ENA_EXAMPLE_PROVENANCE.md.
digest::digest("inst/extdata/ena-inspection.rds", algo = "sha256", file = TRUE)
