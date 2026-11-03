# Intentional RNA normalisation correction

The shared alphabet normaliser now also serves demultiplexing. It uppercases
first and then converts U to T for matching, while preserving original read
text in trimmed inserts, untrimmed reads and R2 mates. The helper was moved
from R/extract.R to R/sequence-utils.R without changing its implementation.

This intentionally differs from 0.99.6 only for lowercase RNA barcode matching:
lowercase RNA barcode names, previously rejected, are accepted, and lowercase
RNA prefixes can now be assigned to their equivalent DNA/RNA codes. This
corrects the already documented U/T equivalence. RNA/DNA aliases still fail
barcode distance validation if they define conflicting assignments.

The remaining baseline comparison covers 72 complete demultiplex objects:
DNA upper/lowercase and uppercase RNA, barcode/read encodings combined,
paired/unpaired, trimming enabled/disabled and mismatch limits 0/1. All agree
exactly with immutable baseline commit c2f7abdc07a12278bb2c459f064f119f0cb7f7e1.
The new four regression blocks add 15 assertions. Focused tests: 30 blocks,
94 assertions passed; no failed/skipped/error/test warning. Full suite belongs
to the root integration.
