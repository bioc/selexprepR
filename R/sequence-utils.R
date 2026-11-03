.normalise_adapter_alphabet <- function(sequence) {
    chartr("U", "T", toupper(sequence))
}

#' Reverse-complement nucleotide sequences
#'
#' This function reproduces the sequence-counting behavior of the Python
#' implementation. It accepts DNA or RNA bases, keeps input case, maps U/u to
#' A/a, and preserves all other symbols, including IUPAC ambiguity codes.
#' It is vectorized over `sequence` using Biostrings. A `BStringSet` reversal
#' follows the explicit complement map because `reverseComplement()` on a
#' `DNAStringSet` would complement IUPAC codes and discard input case, changing
#' the Python-compatible counting contract.
#'
#' @param sequence A character vector of nucleotide sequences without missing
#'   values.
#'
#' @return A character vector of reverse-complemented sequences.
#' @export
#'
#' @examples
#' reverse_complement(c('ACGT', 'ACGU', 'ANNT'))
reverse_complement <- function(sequence) {
    if (!is.character(sequence)) {
        stop("`sequence` must be a character vector.", call. = FALSE)
    }
    if (anyNA(sequence)) {
        stop("`sequence` must not contain missing values.", call. = FALSE)
    }
    complements <- chartr("ACGTUNacgtun", "TGCAANtgcaan", sequence)
    bytewise <- nchar(complements, type = "bytes") == nchar(complements)
    result <- complements
    result[bytewise] <- as.character(Biostrings::reverse(
        Biostrings::BStringSet(complements[bytewise])))
    if (any(!bytewise)) {
        result[!bytewise] <- vapply(strsplit(complements[!bytewise], "",
            fixed = TRUE), function(bases) paste0(rev(bases), collapse = ""),
            character(1))
    }
    result
}
