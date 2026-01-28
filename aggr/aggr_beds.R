#!/usr/bin/env Rscript

library(data.table)
library(yaml)

cfg <- yaml::read_yaml("input.yml")

Sys.setenv(
  GTF_PATH = cfg$reference$gtf,
  COMMON_PY = file.path(cfg$reference$pipeline_dir, "common/common.py")
)

args = commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0)
if (length(args) == 1) {
  InDir = OutDir = args[1]
} else {
  InDir = args[1]
  OutDir = args[2]
}

# Use environment variables or default paths
GTF = Sys.getenv("GTF_PATH")
commonPy = Sys.getenv("COMMON_PY")

stopifnot(file.exists(commonPy), file.exists(GTF), dir.exists(InDir))
if (!dir.exists(OutDir)) {
  dir.create(OutDir, recursive = TRUE)
}

message("Reading GTF annotation data...")
gtf_data = fread(GTF, header = FALSE, sep = "\t", na.strings = c(".", "NA"))
gtf_data = gtf_data[V3 == "gene"]

message("Extracting gene and region information...")
gtf_data[, c("gene_id", "gene") := {
  attrs <- tstrsplit(V9, "; ")
  gene_id <- gsub('"', '', attrs[grep("^gene_id", attrs)])
  gene <- gsub('"', '', attrs[grep("^gene_name", attrs)])
  list(gene_id = gene_id, gene = gene)
}, by = 1:nrow(gtf_data)]

gtf_data = unique(gtf_data[, list(chr = V1, start = V4, end = V5, gene_id = substr(sub("gene_id ", "", gene_id), 1, 15), gene = sub("gene_name ", "", gene))])

message("Checking gtf data...")
print(head(gtf_data))

message("Scanning result bed files and extracting sample IDs...")
methods = c("circexplorer2", "circRNA_finder", "CIRI", "find_circ")
fileList = list.files(InDir, pattern = paste0(paste0("(", paste(methods, collapse = ")|("), ")"), ".bed"))
fileListAll = list.files(InDir, pattern = ".bed", full.names = TRUE)
fileList = fileList[file.info(fileListAll)$size > 0]

# Remove .bed suffix
sample_ids <- gsub("\\.bed$", "", fileList)

# Remove method suffixes
for (suffix in methods) {
  sample_ids <- gsub(paste0("\\.", suffix, "$"), "", sample_ids)
}

sample_ids = unique(sample_ids)
stopifnot(length(sample_ids) > 0)
message(length(sample_ids), " sample(s) detected")

# Functions ------------------------------------------
overlaps <- function(x, y) {
  if (!is.data.frame(x)) {
    stop("x must be a data.frame")
  }
  if (!is.data.frame(y)) {
    stop("y must be a data.frame")
  }

  x <- data.table::as.data.table(x)
  y <- data.table::as.data.table(y)

  colnames(x)[1:3] <- colnames(y)[1:3] <- c("chr", "start", "end")

  data.table::setkey(y, chr, start, end)
  out <- data.table::foverlaps(x, y)[!is.na(start)]
  out
}

aggr_circRNA_beds = function(sample, methods) {
  bed_files = file.path(InDir, paste0(sample, ".", methods, ".bed"))
  nonexists = !file.exists(bed_files)
  len_nonexists = sum(nonexists)
  if (len_nonexists > 0) {
    warning(len_nonexists, " result bed file(s) doesn't exist for sample ", sample, immediate. = TRUE)
    cat(sample, "\t", len_nonexists, "\n", file = file.path(OutDir, "non_exist_results.report.txt"), append = TRUE)
  }
  if (sum(!nonexists) > 1) {
    # Get common regions
    cmd = paste("cat", paste(bed_files[!nonexists], collapse = " "), "|", commonPy, ">", file.path(OutDir, paste0(sample, ".common.txt")))
    system(cmd)
    bed_common = fread(file.path(OutDir, paste0(sample, ".common.txt")), header = FALSE, sep = "\t")
    file.remove(file.path(OutDir, paste0(sample, ".common.txt")))

    # Read all files and filter them
    bed_list = parallel::mclapply(methods[!nonexists], function(x) {
      if (x == "CIRI") {
        d = fread(file.path(InDir, paste0(sample, ".", x, ".bed")), select = c(1:4, 7), header = FALSE, sep = "\t")
      } else {
        d = fread(file.path(InDir, paste0(sample, ".", x, ".bed")), select = 1:5, header = FALSE, sep = "\t")
      }

      if (nrow(d) == 0) {
        warning(sprintf("void data detected for %s with method %s", sample, x), immediate. = TRUE)
        return(d)
      }

      if (is.character(d[[5]])) {
        warning(sprintf("invalid (character) count detected for %s with method %s", sample, x), immediate. = TRUE)
        d[[5]] = as.integer(d[[5]])
      }

      d$tool = x
      d
    }, mc.cores = getOption("mc.cores", 4L))
    bed_dt = rbindlist(bed_list, use.names = FALSE)
    colnames(bed_dt) = c("chr", "start", "end", "strand", "count", "tool")
    bed_dt = bed_dt[!is.na(bed_dt$count)]

    # Add sample-level circRNA ID output
    bed_dt_o = copy(bed_dt)
    bed_dt_o$sample = sample
    fwrite(x = bed_dt_o, file = file_to_all, sep = "\t", append = TRUE)

    # Get and analyze common circRNA
    if (nrow(bed_common) == 0) {
      warning("no common circRNA detected for sample ", sample, immediate. = TRUE)
      return(invisible(NULL))
    }
    colnames(bed_common) = c("chr", "start", "end")
    bed_common = bed_common[chr %in% paste0("chr", c(1:22, "X", "Y"))]
    if (nrow(bed_common) == 0) {
      warning("no common circRNA detected for sample in chr 1-22,X,Y", sample, immediate. = TRUE)
      return(invisible(NULL))
    }

    bed_dt[, id := paste(chr, start, end, sep = "-")]
    bed_common[, id := paste(chr, start, end, sep = "-")]

    annot = unique(overlaps(bed_common, gtf_data)[
      , .(id, gene_id, gene, ovp_len = fcase(
        i.start <= start, i.end - start + 1,
        i.end >= end, end - i.start + 1,
        i.start > start, i.end - i.start + 1
      ))])

    annot = annot[ , list(gene = gene[which.max(ovp_len)]), by = .(id)]

    # Final output for one sample
    bed_dt2 = merge(bed_dt, annot, by = "id", all.x = FALSE, all.y = TRUE)
    solid_ids = bed_dt2[, .(N = sum(count >= 2)), by = .(id)][N >= 1]$id

    if (length(solid_ids) > 0) {
      message("\t=>", length(solid_ids), " solid circRNAs detected")
      bed_dt2 = bed_dt2[id %in% solid_ids]
      bed_dt2$id = NULL
      rv = bed_dt2[, .(tool = paste(tool, collapse = ","),
                       count = mean(count, na.rm=TRUE)),
                   by = .(chr, start, end, strand, gene)]
      rv$sample = sample
      fwrite(x = rv, file = file.path(OutDir, paste0(sample, ".aggr.txt")), sep = "\t")
    } else {
      warning("no solid circRNA detected for sample ", sample, immediate. = TRUE)
    }
  } else {
    warning("no solid circRNA detected for sample ", sample, immediate. = TRUE)
  }
}

# Processing in batch
message("Reading, joining & annotating circRNAs...")

file_to_all = file.path(OutDir, paste0(basename(InDir), ".circRNA_all.txt"))
if (file.exists(file_to_all)) invisible(file.remove(file_to_all))

invisible(
  parallel::mclapply(sample_ids, function(sample) {
    message("\thandling ", sample)
    aggr_circRNA_beds(sample, methods)
  }, mc.cores = min(parallel::detectCores(), 20L, length(sample_ids)))
)

message("Done. Please check result *.aggr.txt in ", OutDir)
