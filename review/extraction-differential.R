Sys.unsetenv('LC_ALL')
arguments <- commandArgs(trailingOnly = TRUE)
repository <- normalizePath(if (length(arguments)) arguments[[1L]] else '.')
baseline_repository <- normalizePath(if (length(arguments) >= 2L) {
    arguments[[2L]]
} else {
    repository
})
pkgload::load_all(repository, quiet = TRUE)
ns <- asNamespace('selexprepR')
baseline <- new.env(parent = ns)
reference <- 'c2f7abdc07a12278bb2c459f064f119f0cb7f7e1'
for (path in c('R/sequence-utils.R', 'R/extract.R')) {
    original <- system2('git', c('-C', shQuote(baseline_repository), 'show',
        shQuote(paste0(reference, ':', path))), stdout = TRUE)
    stopifnot(is.null(attr(original, 'status')))
    eval(parse(text = original), envir = baseline)
}
cat('Reference:', reference, '\n')
current <- function(name) get(name, ns)
p5 <- 'ACGTCAGTACGATCGTACGA'
p3 <- 'TGCATCGACTGACGTAGCTA'
make_reads <- function(n) {
    inserts <- replicate(n, paste(sample(c('A','C','G','T'), 32L, TRUE), collapse=''))
    reads <- paste0(p5, inserts, p3)
    for (i in seq_along(reads)) {
        kind <- (i - 1L) %% 16L
        reads[[i]] <- switch(as.character(kind),
            '0' = reads[[i]],
            '1' = paste0('T', substring(p5, 2L), inserts[[i]], p3),
            '2' = paste0(substr(p5,1L,6L), substring(p5,8L), inserts[[i]],p3),
            '3' = paste0(substr(p5,1L,7L),'T',substring(p5,8L),inserts[[i]],p3),
            '4' = paste0(p5,inserts[[i]],'A',substring(p3,2L)),
            '5' = paste0(p5,inserts[[i]],substr(p3,1L,8L),substring(p3,10L)),
            '6' = paste0(p5,inserts[[i]],substr(p3,1L,9L),'A',substring(p3,10L)),
            '7' = paste0('TG',substring(p5,3L),inserts[[i]],p3),
            '8' = paste0('TGT',substring(p5,4L),inserts[[i]],p3),
            '9' = paste0('NN',p5,inserts[[i]],p3,'NN'),
            '10' = paste0(p5,'ACGT',p5,inserts[[i]],p3,'ACGT',p3),
            '11' = paste0(inserts[[i]],p3),
            '12' = paste0(p5,inserts[[i]]),
            '13' = tolower(chartr('T','U',reads[[i]])),
            '14' = '',
            '15' = substr(p5,1L,18L))
    }
    reads
}
set.seed(934157)
reads <- make_reads(2400L)
for (primer in c(p5,p3,'ACGTACGTACGTACGTACGT')) {
    expected <- t(vapply(reads, function(read) {
        match <- baseline$.adapter_match(read, primer)
        if (is.null(match)) c(NA_integer_, NA_integer_) else match
    }, integer(2L)))
    colnames(expected) <- c('start','end')
    rownames(expected) <- NULL
    actual_coords <- current('.adapter_matches')(reads,primer)
    if (!identical(actual_coords,expected)) {
        print(all.equal(actual_coords,expected)); print(head(which(rowSums(is.na(actual_coords) != is.na(expected) | (!is.na(actual_coords) & !is.na(expected) & actual_coords != expected)) > 0L),20L)); stop('coords mismatch')
    }
}
cat('PASS: 7,200 exact/approximate adapter coordinate comparisons.\n')
for (mode in c('BOTH_PRIMERS_SINGLE_READ','FIVE_PRIME_ONLY','THREE_PRIME_ONLY')) {
    expected <- baseline$.trim_round(reads,mode,p5,p3,TRUE)
    actual <- current('.trim_round')(reads,mode,p5,p3,TRUE)
    stopifnot(identical(actual,expected))
}
cat('PASS: 7,200 per-read extraction results/masks across all single-read modes.\n')
for (inputs in list(reads,character(),c('ACGT','A','',tolower(p5)),
    c('CGTC','ACGT','ACGG'))) {
    for (left in list(p5,NULL,'ACGT')) {
        for (right in list(p3,NULL,'ACGT')) {
            stopifnot(identical(current('.strand_distribution')(inputs,left,right),
                baseline$.strand_distribution(inputs,left,right)))
        }
    }
}
cat('PASS: 36 strand-distribution comparisons, including empty/short/shifted reads.\n')
report <- read_library_report(file.path(repository, 'tests', 'testthat',
    'fixtures', 'library-report-both-primers.json'))
fields <- unclass(report)
fields$primer_5p <- p5
fields$primer_3p <- p3
pools <- list(round_02=reads[seq_len(300L)], round_00=reads[seq_len(400L)], round_01=character())
for (mode in c('BOTH_PRIMERS_SINGLE_READ','FIVE_PRIME_ONLY','THREE_PRIME_ONLY',
    'PAIRED_END_SPLIT_PRIMERS')) {
    for (orientation in c('FORWARD','REVERSE','MIXED')) {
        fields$extraction_mode <- mode
        fields$orientation <- orientation
        fields$full_insert_recovered <- identical(mode,'BOTH_PRIMERS_SINGLE_READ')
        fields$required_action <- if (identical(mode,'PAIRED_END_SPLIT_PRIMERS'))
            'READ_MERGING_RECOMMENDED' else 'NONE'
        r <- current('.new_library_report')(fields)
        r2 <- lapply(pools, function(pool) {
            if (!length(pool)) return(character())
            paste0(reverse_complement(p3),pool)
        })
        expected <- suppressWarnings(baseline$.extract_selexprep(pools,r,r2,keep_masks=TRUE))
        actual <- suppressWarnings(current('.extract_selexprep')(pools,r,r2,keep_masks=TRUE))
        if (!identical(actual,expected)) {
            print(c(mode=mode,orientation=orientation)); print(all.equal(actual,expected)); stop('whole extraction mismatch')
        }
    }
}
cat('PASS: 12 complete extraction-object comparisons, four modes/three orientations/three rounds, including an empty round.\n')
set.seed(2506)
benchmark_inserts <- replicate(3000L,paste(sample(c('A','C','G','T'),32L,TRUE),collapse=''))
exact_reads <- paste0(p5,benchmark_inserts,p3)
stopifnot(length(unique(exact_reads))==length(exact_reads))
mixed_reads <- exact_reads
for (i in seq(1L,length(mixed_reads),by=5L)) substr(mixed_reads[[i]],5L,5L) <- 'G'
for (scenario in c('exact','20pct_substitution')) {
    input <- if (scenario=='exact') exact_reads else mixed_reads
    before <- system.time(old <- baseline$.trim_round(input,'BOTH_PRIMERS_SINGLE_READ',p5,p3))[['elapsed']]
    after <- system.time(new <- current('.trim_round')(input,'BOTH_PRIMERS_SINGLE_READ',p5,p3))[['elapsed']]
    stopifnot(identical(old,new))
    cat(sprintf('BENCHMARK scenario=%s unique_reads=%d baseline_seconds=%.3f vector_seconds=%.3f speedup=%.2f baseline_reads_per_second=%.1f vector_reads_per_second=%.1f\n',
        scenario,length(unique(input)),before,after,before/after,length(input)/before,length(input)/after))
}
