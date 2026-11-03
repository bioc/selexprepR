# Run from the source repository (or pass the candidate and baseline paths).
Sys.unsetenv("LC_ALL")
arguments <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(arguments)) arguments[[1L]] else ".")
baseline_repo <- if (length(arguments) > 1L) arguments[[2L]] else repo
pkgload::load_all(repo, quiet = TRUE)
namespace <- asNamespace("selexprepR")
old <- new.env(parent = namespace)
reference <- "c2f7abdc07a12278bb2c459f064f119f0cb7f7e1"
for (path in c("R/sequence-utils.R", "R/qc.R")) {
    code <- system2("git", c("-C", shQuote(baseline_repo), "show",
        shQuote(paste0(reference, ":", path))), stdout = TRUE)
    stopifnot(is.null(attr(code, "status")))
    eval(parse(text = code), old)
}
new_set <- get(".round_kmer_set", namespace)
new_fraction <- get(".top_kmer_fraction", namespace)
set.seed(738L)
alphabet <- strsplit("ACGTUNacgtunRYSW", "", fixed = TRUE)[[1L]]
comparisons <- 0L
for (n in c(0L, 1L, 5L, 40L)) {
    for (k in c(1L, 3L, 6L, 20L, 100L)) {
        for (top in c(1L, 10L, 10000L)) {
            sequences <- vapply(seq_len(n), function(i) {
                paste0(sample(alphabet, sample(0:50, 1L), TRUE), collapse = "")
            }, character(1))
            reads <- sample(1:10, n, TRUE)
            stopifnot(identical(new_set(sequences, reads, k, top),
                old$.round_kmer_set(sequences, reads, k, top)))
            stopifnot(identical(new_fraction(sequences, reads, k, top),
                old$.top_kmer_fraction(sequences, reads, k, top)))
            comparisons <- comparisons + 2L
        }
    }
}
sequences <- replicate(1000L, paste0(sample(c("A", "C", "G", "T"),
    39L, TRUE), collapse = ""))
reads <- rep(1L, length(sequences))
before <- system.time(baseline <- old$.round_kmer_set(sequences, reads))
after <- system.time(candidate <- new_set(sequences, reads))
stopifnot(identical(baseline, candidate))
cat(sprintf("PASS: %d identical k-mer set/fraction comparisons.\n", comparisons))
cat(sprintf("1,000 39-nt sequences: 0.99.6 %.3f s, vectorized %.3f s; identical sets.\n",
    before[["elapsed"]], after[["elapsed"]]))
