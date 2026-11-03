# Rebuild the real-data subsets from local original files or their HTTP prefixes.
# Usage: Rscript inst/scripts/make-real-selex-example.R RAW_DIRECTORY PACKAGE_ROOT
# No network request is made by this script.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
    stop("Supply RAW_DIRECTORY and PACKAGE_ROOT.", call. = FALSE)
}
raw_directory <- normalizePath(arguments[[1L]], mustWork = TRUE)
package_root <- normalizePath(arguments[[2L]], mustWork = TRUE)
metadata_path <- file.path(package_root, "inst", "extdata",
    "cs-source-metadata.json")
metadata <- jsonlite::fromJSON(metadata_path, simplifyVector = FALSE)

.read_example_records <- function(path, record_count) {
    connection <- gzfile(path, "rb")
    on.exit(close(connection))
    lines <- readLines(connection, n = 4L * record_count, warn = FALSE)
    stopifnot(length(lines) == 4L * record_count)
    offsets <- seq.int(1L, length(lines), 4L)
    stopifnot(all(startsWith(lines[offsets], "@")))
    stopifnot(all(startsWith(lines[offsets + 2L], "+")))
    stopifnot(all(nchar(lines[offsets + 1L], type = "bytes") ==
        nchar(lines[offsets + 3L], type = "bytes")))
    charToRaw(paste0(paste(lines, collapse = "\n"), "\n"))
}

for (subset in metadata$subsets) {
    source_path <- file.path(raw_directory, subset$source_file)
    if (file.info(source_path)$size == subset$source_bytes) {
        stopifnot(identical(digest::digest(file = source_path, algo = "md5"),
            subset$source_md5))
    } else {
        stopifnot(file.info(source_path)$size == subset$prefix_bytes)
        stopifnot(identical(digest::digest(file = source_path, algo = "sha256"),
            subset$prefix_sha256))
    }
    records <- .read_example_records(source_path, subset$n_records)
    stopifnot(identical(digest::digest(records, algo = "sha256", serialize = FALSE),
        subset$uncompressed_sha256))
    target_path <- file.path(package_root, "inst", "extdata", subset$subset_file)
    connection <- gzfile(target_path, "wb", compression = 9L)
    writeBin(records, connection)
    close(connection)
    message(subset$subset_file, ": verified ", subset$n_records,
        " original records; compressed SHA-256 ",
        digest::digest(file = target_path, algo = "sha256"))
}
