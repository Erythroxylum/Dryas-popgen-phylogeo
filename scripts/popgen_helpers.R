required_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Install the conda environment from environment.yml."
    )
  }
}

read_vcf_ids <- function(path) {
  bcftools <- Sys.which("bcftools")
  if (nzchar(bcftools)) {
    ids <- system2(bcftools, c("query", "-l", path), stdout = TRUE, stderr = TRUE)
    status <- attr(ids, "status")
    if (is.null(status) || status == 0) return(ids)
  }
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) stop("VCF has no #CHROM header: ", path)
    if (startsWith(line, "#CHROM")) return(strsplit(line, "\t", fixed = TRUE)[[1]][-(1:9)])
  }
}

read_selected_map <- function(path, vcf_ids) {
  map <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c(
    "raw_id", "sampleID", "phyloID", "pop_gen", "pop_gen_population",
    "spID", "geoID", "Latitude", "Longitude", "Plate", "selected"
  )
  absent <- setdiff(required, names(map))
  if (length(absent)) stop("Crosswalk columns not found: ", paste(absent, collapse = ", "))
  is_selected <- map$selected %in% c(TRUE, "TRUE", "True", "true", 1, "1")
  map <- map[is_selected, , drop = FALSE]
  if (anyDuplicated(map$raw_id)) stop("Duplicated raw_id in selected crosswalk")
  if (!setequal(vcf_ids, map$raw_id)) {
    stop(
      "VCF and selected crosswalk contain different samples. Only in VCF: ",
      paste(setdiff(vcf_ids, map$raw_id), collapse = ", "),
      "; only in crosswalk: ", paste(setdiff(map$raw_id, vcf_ids), collapse = ", ")
    )
  }
  map[match(vcf_ids, map$raw_id), , drop = FALSE]
}

load_dosage <- function(vcf_file, map_file) {
  required_packages(c("dartR", "adegenet", "vcfR"))
  suppressPackageStartupMessages(library(adegenet))
  vcf_ids <- read_vcf_ids(vcf_file)
  map <- read_selected_map(map_file, vcf_ids)
  dartr_metadata <- data.frame(
    id = map$raw_id,
    pop = map$pop_gen_population,
    lat = suppressWarnings(as.numeric(map$Latitude)),
    lon = suppressWarnings(as.numeric(map$Longitude)),
    sampleID = map$sampleID,
    phyloID = map$phyloID,
    pop_gen = map$pop_gen,
    stringsAsFactors = FALSE
  )
  metadata_temp <- tempfile(fileext = ".csv")
  on.exit(unlink(metadata_temp), add = TRUE)
  write.csv(dartr_metadata, metadata_temp, row.names = FALSE, na = "")
  gl <- dartR::gl.read.vcf(vcf_file, ind.metafile = metadata_temp, verbose = 2)
  observed_ids <- adegenet::indNames(gl)
  if (!identical(observed_ids, vcf_ids)) {
    stop("dartR changed VCF sample order or identifiers")
  }
  gl@pop <- factor(map$pop_gen_population)
  dos <- as.matrix(gl)
  storage.mode(dos) <- "integer"
  rownames(dos) <- vcf_ids
  invalid <- !is.na(dos) & !(dos %in% 0:2)
  if (any(invalid)) stop("Genotype matrix contains values other than 0, 1, 2, or NA")
  list(gl = gl, dos = dos, map = map)
}

make_hierfstat <- function(dos, population) {
  pop_factor <- factor(population)
  geno <- matrix(NA_integer_, nrow = nrow(dos), ncol = ncol(dos),
                 dimnames = dimnames(dos))
  geno[dos == 0] <- 11L
  geno[dos == 1] <- 12L
  geno[dos == 2] <- 22L
  hf <- data.frame(pop = as.integer(pop_factor), geno, check.names = FALSE)
  list(data = hf, levels = levels(pop_factor))
}

rename_population_matrix <- function(x, population_levels) {
  if (is.null(dim(x))) return(x)
  old <- colnames(x)
  if (!is.null(old) && all(old %in% as.character(seq_along(population_levels)))) {
    colnames(x) <- population_levels[as.integer(old)]
  } else if (ncol(x) == length(population_levels)) {
    colnames(x) <- population_levels
  }
  old_rows <- rownames(x)
  if (!is.null(old_rows) && all(old_rows %in% as.character(seq_along(population_levels)))) {
    rownames(x) <- population_levels[as.integer(old_rows)]
  } else if (nrow(x) == length(population_levels) && ncol(x) == length(population_levels)) {
    rownames(x) <- population_levels
  }
  x
}

matrix_to_long <- function(x, value_name) {
  out <- as.data.frame(as.table(x), stringsAsFactors = FALSE)
  names(out) <- c("population_1", "population_2", value_name)
  out
}

mean_na <- function(x) {
  result <- colMeans(x, na.rm = TRUE)
  result[!is.finite(result)] <- NA_real_
  result
}

