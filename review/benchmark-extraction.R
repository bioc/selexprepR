# Run from a selexprepR checkout, or pass its path as the first argument.
# A second argument can supply the baseline Git checkout for a source-only copy.
# Source-mode developer verification; no package dependency is relaxed.
Sys.unsetenv('LC_ALL')
arguments <- commandArgs(trailingOnly = TRUE)
repository <- normalizePath(if (length(arguments)) arguments[[1L]] else '.')
pkgload::load_all(repository, quiet = TRUE)
baseline_repository <- normalizePath(if (length(arguments) >= 2L) {
    arguments[[2L]]
} else {
    repository
})
namespace <- asNamespace('selexprepR')
baseline <- new.env(parent = namespace)
reference <- 'c2f7abdc07a12278bb2c459f064f119f0cb7f7e1'
for (path in c('R/sequence-utils.R', 'R/extract.R')) {
    original <- system2('git', c('-C', shQuote(baseline_repository), 'show',
        shQuote(paste0(reference, ':', path))), stdout = TRUE)
    stopifnot(is.null(attr(original, 'status')))
    eval(parse(text = original), envir = baseline)
}
current <- get('.extract_selexprep', namespace)
report <- read_library_report(file.path(repository, 'tests', 'testthat',
    'fixtures', 'library-report-both-primers.json'))
fields <- unclass(report)
fields$primer_5p <- 'ACGTCAGTACGATCGTACGA'
fields$primer_3p <- 'TGCATCGACTGACGTAGCTA'
report <- get('.new_library_report', namespace)(fields)
set.seed(2506L)
number <- 3000L
inserts <- replicate(number, paste(sample(c('A', 'C', 'G', 'T'), 32L, TRUE),
    collapse = ''))
exact <- paste0(report$primer_5p, inserts, report$primer_3p)
mutated <- exact
substitutions <- seq.int(1L, number, by = 5L)
substr(mutated[substitutions], 5L, 5L) <- 'G'
missing <- mutated
unmatched <- seq.int(2L, number, by = 5L)
missing[unmatched] <- paste0(inserts[unmatched], report$primer_3p)
scenarios <- list(exact = exact, mixed_20pct_substitution = mutated,
    mixed_20pct_substitution_20pct_unmatched = missing)
cat('Reference:', reference, '\n')
cat('R:', R.version.string, '\n')
cat('Biostrings:', as.character(utils::packageVersion('Biostrings')), '\n')
cat('Seed: 2506; single elapsed-time measurement per scenario/implementation.\n')
cat('Counts, inserts and all extraction metadata must match exactly.\n')
cat('Rates classify the 5-prime adapter; the 3-prime adapter is unchanged.\n')
for (scenario in names(scenarios)) {
    reads <- scenarios[[scenario]]
    stopifnot(length(unique(reads)) == number)
    subjects <- Biostrings::BStringSet(reads)
    pattern <- Biostrings::BString(report$primer_5p)
    exact_hit <- Biostrings::vcountPattern(pattern, subjects, fixed = TRUE) > 0L
    any_hit <- Biostrings::vcountPattern(pattern, subjects, max.mismatch = 2L,
        with.indels = TRUE, fixed = TRUE) > 0L
    before <- system.time(old <- baseline$.extract_selexprep(
        list(round_00 = reads), report))[['elapsed']]
    after <- system.time(new <- current(list(round_00 = reads), report))[['elapsed']]
    stopifnot(identical(old, new))
    cat(sprintf(paste0(
        'scenario=%s unique_reads=%d exact_pct=%.1f approximate_pct=%.1f ',
        'unmatched_pct=%.1f baseline_seconds=%.3f revised_seconds=%.3f ',
        'baseline_reads_per_second=%.1f revised_reads_per_second=%.1f ',
        'speedup=%.2f identical=TRUE\n'),
        scenario, number, 100 * mean(exact_hit),
        100 * mean(any_hit & !exact_hit), 100 * mean(!any_hit),
        before, after, number / before, number / after, before / after))
}
cat('PASS: all three complete extraction objects agree with 0.99.6.\n')
