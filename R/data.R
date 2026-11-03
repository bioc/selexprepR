#' Curated public HT-SELEX catalog
#'
#' A frozen, offline snapshot of 240 public HT-SELEX deposits. Eight
#' experimental fields were extracted independently by two language models,
#' reconciled, and retained with a per-field curation status. The flat dataset
#' intentionally omits evidence quotations; immutable links to the canonical
#' provenance-rich JSON and its generation materials are stored in
#' `S4Vectors::metadata(selexprep_public_catalog)`.
#'
#' The catalog compilation and curation annotations are MIT-licensed under
#' the pinned upstream grant. Source papers and sequencing reads are not
#' distributed as part of this dataset; their terms are recorded separately
#' in the installed `DATA_LICENSES.md` file.
#'
#' Missing values mean that a field was not stated in the curated sources. A
#' `adjudicated` means that the two extraction arms disagreed and the value
#' was resolved using the source evidence. The canonical provenance JSON
#' retains both original values and the adjudication rationale.
#'
#' @format An `S4Vectors::DataFrame` with 240 rows and 19 columns:
#' \describe{
#'   \item{`bioproject_id`}{ENA BioProject accession or a platform-prefixed
#'     Figshare or Zenodo deposit identifier.}
#'   \item{`source`}{`ena`, or the platform-prefixed Figshare or Zenodo
#'     deposit identifier used as the source key.}
#'   \item{`study_title`}{Title recorded for the public study or deposit.}
#'   \item{`study_type`, `study_type_curation`}{Curated experiment type and
#'     its reconciliation status.}
#'   \item{`target`, `target_curation`}{Selection target and its
#'     reconciliation status.}
#'   \item{`target_class`, `target_class_curation`}{Target class and its
#'     reconciliation status.}
#'   \item{`chemistry`, `chemistry_curation`}{DNA or RNA library chemistry
#'     and its reconciliation status.}
#'   \item{`n_random`, `n_random_curation`}{Reported random-region length and
#'     its reconciliation status.}
#'   \item{`n_rounds`, `n_rounds_curation`}{Reported selection-round count and
#'     its reconciliation status.}
#'   \item{`selection_format`, `selection_format_curation`}{Reported selection
#'     format and its reconciliation status.}
#'   \item{`counter_selection`, `counter_selection_curation`}{Reported
#'     counter-selection conditions and their reconciliation status.}
#' }
#'
#' Curation statuses are `adjudicated`, `concordant`, `not_stated`, `verified`,
#' `single_source:claude`, or `single_source:codex`. These describe upstream
#' curation; a status alone does not establish human validation. Maintainer
#' spot checks are documented separately in `CATALOG_PROVENANCE.md`.
#'
#' @source The snapshot was migrated without modification from the original
#'   selexprep v0.4.6 catalog at commit
#'   `63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5`. See the object metadata and
#'   `inst/CATALOG_PROVENANCE.md` for immutable source links and terms.
#' @examples
#' data("selexprep_public_catalog", package = "selexprepR")
#' head(selexprep_public_catalog[, c("bioproject_id", "study_type", "target")])
#' S4Vectors::metadata(selexprep_public_catalog)$snapshot_version
#' @usage data(selexprep_public_catalog)
#' @keywords datasets
"selexprep_public_catalog"
