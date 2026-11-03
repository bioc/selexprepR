# Developer check: run documented offline examples from a source checkout.
# Live ENA calls inside donttest are deliberately excluded, as in R CMD check.
Sys.unsetenv("LC_ALL")
arguments <- commandArgs(trailingOnly = TRUE)
repository <- normalizePath(if (length(arguments)) arguments[[1L]] else ".")
pkgload::load_all(repository, quiet = TRUE)
run_examples <- function() {
    output <- tempfile("selexprep-examples-")
    dir.create(output)
    grDevices::pdf(file.path(output, "example-plots.pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    files <- list.files(file.path(repository, "man"), "[.]Rd$", full.names = TRUE)
    optional_network <- 0L
    for (path in files) {
        rd <- tools::parse_Rd(path)
        stopifnot(is.null(tools::checkRd(rd)) ||
            length(tools::checkRd(rd)) == 0L)
        text <- readLines(path, warn = FALSE)
        stopifnot(any(grepl("\\examples{", text, fixed = TRUE)))
        optional_network <- optional_network + sum(grepl(
            "\\donttest{", text, fixed = TRUE
        ))
        script <- file.path(output, paste0(basename(path), ".R"))
        tools::Rd2ex(rd, out = script, commentDontrun = TRUE,
            commentDonttest = TRUE)
        cat("EXAMPLES", basename(path), "\n")
        sys.source(script, envir = new.env(parent = globalenv()))
    }
    cat("PASS:", length(files), "help pages with runnable offline examples;",
        optional_network, "optional live-network block excluded.\n")
}
run_examples()
