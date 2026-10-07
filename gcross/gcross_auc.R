#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
#
# UM-MAP, G-cross step: G-cross area-under-the-curve (Gfx) summaries for every slide in a
# folder of per-cell tables.
#
# Dependencies (tested versions): R 4.4.0, spatstat 3.0-8 (spatstat.geom 3.2-9,
# spatstat.explore 3.2-7, spatstat.random 3.2-3), pracma 2.4.6, and openxlsx 4.2.9 (needed
# only for --xlsx). Install with: install.packages(c("spatstat", "pracma", "openxlsx"))
#
# Usage:
#   Rscript gcross/gcross_auc.R --input <folder> --output <results.csv> [options]
#
#   --input DIR       folder with one per-cell table per slide, named *.tsv or *.tsv.gz
#                     (sub-folders are not searched); required
#   --output FILE     results table (CSV), one row per slide; required
#   --radii LIST      comma-separated whole-number radii from 1 to 300; default 10,20,40
#   --xlsx FILE       also write the results table as an Excel workbook
#   --classes LIST    cell classes that form the pairs, in order;
#                     default Fibroblast,Tumor,Lymphoid
#   --class-col NAME  name of the cell-class column; default Class
#   --overwrite       replace existing output files (without it, the script stops if they exist)
#   --help            print the usage
#
# Input: tab-separated text with a header line and one row per detected cell, as written by
# qupath/classify_export.groovy, with columns Class, Centroid X µm and Centroid Y µm. Other
# columns are ignored. The centroid columns are found by name ("Centroid X ...",
# "Centroid Y ..."), so the result does not depend on the R locale. The per-cell tables of
# the UMICH cohorts are not distributed; demo/example_input holds one public TCGA table.
#
# Method, for each slide:
#   1. Spaces are removed from class names. Cells without a class are counted as
#      "Unclassified".
#   2. Window: the convex hull of the centroids of all cells on the slide, whatever their
#      class. Classes that are in no pair (for example Myeloid or Garbage) still shape the
#      window. The window is not a tissue outline.
#   3. G-cross for each ordered pair (x = reference class, y = target class):
#        G_xy(r) = probability that the nearest type-y cell of a type-x cell is within
#                  distance r,
#      with Kaplan-Meier edge correction, on the grid r = 0, 1, ..., 300:
#      spatstat Gcross(pp, x, y, r = 0:300)$km. For x = y, spatstat uses the nearest other
#      cell of the same class (the nearest-neighbour function G of that class).
#   4. AUC up to each radius R: trapezoidal rule (pracma::trapz) over r = 0, 1, ..., R.
#      The AUC up to R is at most R.
#   Pairs for the default classes: Fibroblast->Tumor, Fibroblast->Lymphoid, Tumor->Lymphoid,
#   Fibroblast->Fibroblast, Tumor->Tumor, Lymphoid->Lymphoid (every pair of different
#   classes in --classes order, then each class with itself).
#
# Missing-class rule: if either class of a pair has no cells on a slide, its curve is set to 0
# at every r, so its AUC is 0, not NA. A class with a single cell has no other cell of its
# class, so its self-pair AUC is 0 as well. A value of 0 can therefore mean that a class is
# absent or nearly absent; check the n_<class> columns. A slide whose cells do not span an
# area (fewer than 3 cells, or all on one line) has no window: its AUCs are NA and a warning
# is printed.
#
# Units: coordinates are expected in micrometres. Radii, the r grid and the AUC are in the
# units of the coordinates, so a table with pixel coordinates gives radii in pixels and AUCs
# in pixel units, which cannot be compared with micrometre values (at 0.25 um per pixel,
# 10 um is 40 px and the AUC comes out about 4 times larger). The script warns when the
# centroid column names give "px" or no unit.
#
# Output, one row per slide (the optional workbook holds the same table):
#   slide            file name without .tsv or .tsv.gz
#   file             file name
#   n_cells          number of cells (rows)
#   n_<class>        cells per class: the --classes first, then every other class found in
#                    any input file (0 where absent), e.g. n_Myeloid, n_Garbage
#   hull_area        area of the convex-hull window, in squared coordinate units
#   Gfx<R>_<x>_<y>   AUC of the x->y G-cross up to radius R, e.g. Gfx10_Tumor_Lymphoid
#   Numbers in the CSV have 15 significant digits.
#
# References:
#   G-cross AUC method: Barua S, Fang P, Sharma A, Fujimoto J, Wistuba I, Rao AUK, Lin SH.
#   Spatial interaction of tumor cells and regulatory T cells correlates with survival in
#   non-small cell lung cancer. Lung Cancer 2018;117:73-79.
#   https://doi.org/10.1016/j.lungcan.2018.01.022
#   spatstat: Baddeley A, Rubak E, Turner R. Spatial Point Patterns: Methodology and
#   Applications with R. Chapman and Hall/CRC Press, London, 2015.

USAGE <- paste(
  "Usage: Rscript gcross/gcross_auc.R --input <folder> --output <results.csv> [options]",
  "  --input DIR       folder of per-cell tables, one per slide (*.tsv or *.tsv.gz)",
  "  --output FILE     results table (CSV), one row per slide",
  "  --radii LIST      comma-separated whole-number radii, 1 to 300 (default 10,20,40)",
  "  --xlsx FILE       also write the results as an Excel workbook (needs openxlsx)",
  "  --classes LIST    classes that form the pairs, in order (default Fibroblast,Tumor,Lymphoid)",
  "  --class-col NAME  cell-class column (default Class)",
  "  --overwrite       replace existing output files",
  "  --help            print this usage",
  sep = "\n")

R_GRID <- 0:300  # evaluation grid of the G-cross curve, in coordinate units (um)

parse_args <- function(argv) {
  opts <- list(input = NULL, output = NULL, radii = "10,20,40", xlsx = NULL,
               classes = "Fibroblast,Tumor,Lymphoid", class_col = "Class",
               overwrite = FALSE, help = FALSE)
  keys <- c("--input" = "input", "--output" = "output", "--radii" = "radii",
            "--xlsx" = "xlsx", "--classes" = "classes", "--class-col" = "class_col")
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (a %in% c("-h", "--help")) {
      opts$help <- TRUE
      i <- i + 1L
      next
    }
    if (a == "--overwrite") {
      opts$overwrite <- TRUE
      i <- i + 1L
      next
    }
    key <- sub("=.*$", "", a)
    if (!key %in% names(keys)) stop("unknown argument '", a, "'\n", USAGE, call. = FALSE)
    if (grepl("=", a, fixed = TRUE)) {
      val <- sub("^[^=]*=", "", a)
      i <- i + 1L
    } else {
      if (i == length(argv)) stop("missing value after ", key, "\n", USAGE, call. = FALSE)
      val <- argv[[i + 1L]]
      i <- i + 2L
    }
    opts[[keys[[key]]]] <- val
  }
  opts
}

split_list <- function(s) {
  v <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
  v[nzchar(v)]
}

parse_radii <- function(s) {
  v <- split_list(s)
  r <- suppressWarnings(as.numeric(v))
  if (!length(r) || anyNA(r) || any(r != round(r)) || any(r < 1) || any(r > max(R_GRID)))
    stop("--radii must be comma-separated whole numbers from 1 to ", max(R_GRID),
         ", got '", s, "'", call. = FALSE)
  unique(as.integer(r))
}

parse_classes <- function(s) {
  # Spaces are removed from class names in the data, so they are removed here as well.
  v <- gsub(" ", "", split_list(s), fixed = TRUE)
  v <- v[nzchar(v)]
  if (!length(v)) stop("--classes is empty", call. = FALSE)
  if (anyDuplicated(v)) stop("--classes has a repeated class: '", s, "'", call. = FALSE)
  v
}

# Ordered pairs (reference, target): each pair of different classes in the given order,
# then each class with itself. For Fibroblast,Tumor,Lymphoid this gives F-T, F-L, T-L,
# F-F, T-T, L-L.
make_pairs <- function(classes) {
  diff_pairs <- if (length(classes) > 1L) utils::combn(classes, 2L, simplify = FALSE) else list()
  c(diff_pairs, lapply(classes, function(k) c(k, k)))
}

auc_names <- function(pairs, radii) {
  unlist(lapply(radii, function(rad)
    vapply(pairs, function(p) sprintf("Gfx%d_%s_%s", rad, p[[1]], p[[2]]), character(1))))
}

# Unit named in a centroid column header: "um", "px" or "unknown". Matched on bytes so that
# the micro sign (UTF-8 C2 B5) and the Greek mu (CE BC) work in every locale.
coord_unit <- function(nm) {
  if (grepl("\xc2\xb5m|\xce\xbcm|[Mm]icron|(^|[^A-Za-z])um([^A-Za-z]|$)", nm, useBytes = TRUE))
    return("um")
  if (grepl("(^|[^A-Za-z])px([^A-Za-z]|$)|[Pp]ixel", nm, useBytes = TRUE)) return("px")
  "unknown"
}

pick_centroid <- function(nms, axis, file) {
  pat <- sprintf("^[[:space:]]*[Cc]entroid[ ._]?[%s%s]([^A-Za-z0-9]|$)", axis, tolower(axis))
  hit <- nms[grepl(pat, nms, useBytes = TRUE)]
  if (length(hit) > 1L) {
    in_um <- hit[vapply(hit, coord_unit, character(1)) == "um"]
    if (length(in_um) == 1L) hit <- in_um
  }
  if (length(hit) != 1L)
    stop(sprintf("%s: expected one 'Centroid %s' column, found %d (columns: %s)",
                 file, axis, length(hit), paste(nms, collapse = ", ")), call. = FALSE)
  hit
}

read_cells <- function(path, class_col) {
  file <- basename(path)
  # read.table, tab-separated, with a header line; column names are kept as written
  # (check.names = FALSE) instead of being converted by make.names.
  d <- tryCatch(
    utils::read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                      stringsAsFactors = FALSE, encoding = "UTF-8"),
    error = function(e) stop(file, ": cannot read the table: ", conditionMessage(e), call. = FALSE))
  nms <- names(d)
  if (!class_col %in% nms)
    stop(sprintf("%s: no '%s' column (columns: %s); set --class-col", file, class_col,
                 paste(nms, collapse = ", ")), call. = FALSE)
  xcol <- pick_centroid(nms, "X", file)
  ycol <- pick_centroid(nms, "Y", file)
  units <- unique(c(coord_unit(xcol), coord_unit(ycol)))
  if (nrow(d) == 0L)
    return(list(x = numeric(0), y = numeric(0), cls = character(0), units = units))
  x <- d[[xcol]]
  y <- d[[ycol]]
  if (!is.numeric(x) || !is.numeric(y) || any(!is.finite(x)) || any(!is.finite(y)))
    stop(file, ": centroid columns must hold finite numbers", call. = FALSE)
  cls <- gsub(" ", "", as.character(d[[class_col]]), fixed = TRUE)
  cls[is.na(cls) | !nzchar(cls)] <- "Unclassified"
  list(x = x, y = y, cls = cls, units = units)
}

# G-cross AUCs of one slide. The bounding rectangle below (0 to max + 100 for non-negative
# coordinates) is only used to build the pattern before its window is replaced by the
# convex hull.
gcross_slide <- function(x, y, cls, pairs, radii) {
  auc <- stats::setNames(rep(NA_real_, length(pairs) * length(radii)), auc_names(pairs, radii))
  win <- NULL
  if (length(x) >= 3L) {
    # convexhull.xy() returns NULL when the points do not span an area (all on one line);
    # owin(poly = NULL) would then silently give the unit square, so the slide gets no
    # window and NA values instead.
    cw <- tryCatch(convexhull.xy(x, y), error = function(e) NULL)
    if (!is.null(cw)) win <- owin(poly = cw$bdry[[1]])
  }
  if (is.null(win)) return(list(hull_area = NA_real_, auc = auc))
  xmin <- min(0, min(x))
  ymin <- min(0, min(y))
  xmax <- max(x) + 100
  ymax <- max(y) + 100
  pp <- as.ppp(cbind(x, y), W = owin(c(xmin, xmax), c(ymin, ymax)))
  pp <- pp %mark% factor(cls)
  pp$window <- win
  present <- levels(pp$marks)
  for (p in pairs) {
    if (all(c(p[[1]], p[[2]]) %in% present)) {
      g <- Gcross(pp, p[[1]], p[[2]], r = R_GRID)$km
    } else {
      g <- rep(0, length(R_GRID))  # missing-class rule
    }
    for (rad in radii) {
      sel <- R_GRID <= rad
      auc[[sprintf("Gfx%d_%s_%s", rad, p[[1]], p[[2]])]] <- pracma::trapz(R_GRID[sel], g[sel])
    }
  }
  list(hull_area = area(win), auc = auc)
}

pkg_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "not installed")
}

main <- function(argv) {
  options(warn = 1)
  opts <- parse_args(argv)
  if (isTRUE(opts$help)) {
    cat(USAGE, "\n", sep = "")
    return(invisible(NULL))
  }
  if (is.null(opts$input) || is.null(opts$output))
    stop("--input and --output are required\n", USAGE, call. = FALSE)
  if (!dir.exists(opts$input)) stop("input folder not found: ", opts$input, call. = FALSE)
  if (grepl("\\.xlsx?$", opts$output, ignore.case = TRUE))
    stop("--output is written as CSV; use --xlsx for a workbook", call. = FALSE)
  radii <- parse_radii(opts$radii)
  classes <- parse_classes(opts$classes)
  pairs <- make_pairs(classes)

  outputs <- c(opts$output, opts$xlsx)
  existing <- outputs[file.exists(outputs)]
  if (length(existing) && !opts$overwrite)
    stop("output file exists: ", paste(existing, collapse = ", "),
         " (use --overwrite to replace it)", call. = FALSE)

  files <- list.files(opts$input, pattern = "\\.tsv(\\.gz)?$", ignore.case = TRUE,
                      full.names = TRUE)
  files <- files[!dir.exists(files)]
  if (!length(files)) stop("no .tsv or .tsv.gz files in ", opts$input, call. = FALSE)
  slides <- sub("\\.tsv(\\.gz)?$", "", basename(files), ignore.case = TRUE)
  if (anyDuplicated(slides))
    stop("more than one file for slide(s): ",
         paste(unique(slides[duplicated(slides)]), collapse = ", "), call. = FALSE)

  for (pkg in c("spatstat", "pracma", if (!is.null(opts$xlsx)) "openxlsx"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("R package '", pkg, "' is not installed: install.packages(\"", pkg, "\")",
           call. = FALSE)
  suppressPackageStartupMessages(library(spatstat))

  cat(sprintf("G-cross AUC: %d file(s) in %s; radii %s; pairs %s\n", length(files), opts$input,
              paste(radii, collapse = ","),
              paste(vapply(pairs, paste, character(1), collapse = "->"), collapse = ", ")))

  rows <- vector("list", length(files))
  for (k in seq_along(files)) {
    t0 <- proc.time()[["elapsed"]]
    cells <- read_cells(files[[k]], opts$class_col)
    if (!identical(cells$units, "um"))
      warning(sprintf(paste("%s: centroid columns are in '%s' units, not um; radii and AUCs",
                            "are in the units of the coordinates"),
                      basename(files[[k]]), paste(cells$units, collapse = "/")), call. = FALSE)
    res <- tryCatch(gcross_slide(cells$x, cells$y, cells$cls, pairs, radii),
                    error = function(e) stop(basename(files[[k]]), ": ", conditionMessage(e),
                                             call. = FALSE))
    if (is.na(res$hull_area))
      warning(sprintf(paste("%s: %d cell(s) do not span an area (fewer than 3 cells, or all on",
                            "one line); AUCs set to NA"), basename(files[[k]]),
                      length(cells$x)), call. = FALSE)
    counts <- table(cells$cls)
    count_names <- as.character(names(counts))
    rows[[k]] <- list(counts = stats::setNames(as.integer(counts), count_names),
                      n_cells = length(cells$x), hull_area = res$hull_area, auc = res$auc)
    shown <- c(classes[classes %in% count_names],
               sort(setdiff(count_names, classes), method = "radix"))
    cat(sprintf("[%d/%d] %s: %d cells (%s), %.1f s\n", k, length(files), basename(files[[k]]),
                length(cells$x), paste(shown, as.integer(counts[shown]), collapse = ", "),
                proc.time()[["elapsed"]] - t0))
  }

  all_classes <- as.character(unique(unlist(lapply(rows, function(r) names(r$counts)))))
  count_classes <- c(classes, sort(setdiff(all_classes, classes), method = "radix"))
  out <- data.frame(slide = slides, file = basename(files),
                    n_cells = vapply(rows, function(r) r$n_cells, integer(1)),
                    stringsAsFactors = FALSE, check.names = FALSE)
  for (cl in count_classes)
    out[[paste0("n_", cl)]] <- vapply(rows, function(r)
      if (cl %in% names(r$counts)) r$counts[[cl]] else 0L, integer(1))
  out$hull_area <- vapply(rows, function(r) r$hull_area, numeric(1))
  auc <- do.call(rbind, lapply(rows, function(r) r$auc))
  for (nm in colnames(auc)) out[[nm]] <- unname(auc[, nm])

  out_dir <- dirname(opts$output)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  utils::write.csv(out, opts$output, row.names = FALSE)
  cat(sprintf("Wrote %d row(s) to %s\n", nrow(out), opts$output))
  if (!is.null(opts$xlsx)) {
    xlsx_dir <- dirname(opts$xlsx)
    if (!dir.exists(xlsx_dir)) dir.create(xlsx_dir, recursive = TRUE)
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "gcross_auc")
    openxlsx::writeData(wb, "gcross_auc", out)
    openxlsx::saveWorkbook(wb, opts$xlsx, overwrite = opts$overwrite)
    cat(sprintf("Wrote %s\n", opts$xlsx))
  }
  cat(sprintf("%s; spatstat %s (spatstat.geom %s, spatstat.explore %s); pracma %s\n",
              R.version.string, pkg_version("spatstat"), pkg_version("spatstat.geom"),
              pkg_version("spatstat.explore"), pkg_version("pracma")))
  invisible(out)
}

# Run main() only when the file is executed with Rscript, not when it is source()d (the
# functions above can then be reused, e.g. for testing).
if (sys.nframe() == 0L) {
  status <- tryCatch({
    main(commandArgs(trailingOnly = TRUE))
    0L
  }, error = function(e) {
    message("Error: ", conditionMessage(e))
    1L
  })
  quit(save = "no", status = status)
}
