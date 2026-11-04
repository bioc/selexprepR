# cutadapt 5.2 primer-trimming references

Python selexprep trims primers with cutadapt; selexprepR ports its aligner
(`R/primer-match.R`). These references record what cutadapt 5.2 (Python 3.12)
returned, with the options Python selexprep passes. In both files `.` marks a
read that cutadapt discarded.

## `cutadapt-5.2-synthetic-reference.tsv.gz`

2,000 synthetic reads, drawn with Python's `random.Random(7)`, and the output
of four runs over them:

```text
cutadapt --discard-untrimmed -g CCAGAATGGAGAACGCAG...GTCGTCCCGAGAGCC   both
cutadapt --discard-untrimmed -g CCAGAATGGAGAACGCAG                     five_prime_only
cutadapt --discard-untrimmed -a GTCGTCCCGAGAGCC                        three_prime_only
cutadapt --discard-untrimmed -g CCAGAATGGNGAACRCAG...GTCGTCCCGAGAGCY   both_iupac
```

Each read joins random bases with copies of the two primers. A primer copy is
intact, carries one to three substitutions, insertions or deletions (some
substitutions are `N`), or is absent. Some 5' primers are cut short at the
read start, some reads are cut short at the end, some carry a second pair of
primer copies, and some are random reads. Insert lengths include 0, 39 and
40 nucleotides. The reads are stored in the file, so the generator is not
needed to run the tests.

## `cutadapt-5.2-real-coordinates.tsv.gz`

For every read of the three bundled R1 subsets (`inst/extdata`), the start and
end (1-based, inclusive) of the insert that
`cutadapt --discard-untrimmed -g CCAGAATGGAGAACGCAG...GTCGTCCCGAGAGCC` kept,
taken from its `--info-file`. Only coordinates are stored, not sequences.
