.normalise_adapter_alphabet <- function(sequence) {
    chartr("U", "T", toupper(sequence))
}

#' Reverse-complement nucleotide sequences
#'
#' Accepts DNA or RNA bases and keeps input case. A/T, C/G and the IUPAC
#' ambiguity codes are complemented (R/Y, K/M, B/V, D/H; S, W and N are their
#' own complements), and U/u maps to A/a. Other symbols are kept unchanged.
#' It is vectorized over `sequence` using Biostrings. Unlike
#' `Biostrings::reverseComplement()` on a `DNAStringSet`, it keeps lowercase
#' letters and accepts RNA and non-IUPAC symbols, so read sequences are
#' returned in their original alphabet.
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
    complements <- chartr("ACGTUNRYKMSWBVDHacgtunrykmswbvdh",
        "TGCAANYRMKSWVBHDtgcaanyrmkswvbhd", sequence)
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
