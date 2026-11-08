library(data.table)
library(plinkbinr)

# Configuration
WORK_DIR <- "~/lj74/CD4T"
EXPOSURE_FILE <- "CD4all1.8.txt"
EUR_BFILE <- "/data/h004/lj74/1000G.EUR/1000G_EUR"
CLUMP_WORK_DIR <- "/data/h004/lj74/Rbao/plink_clump_sensitive3_simplified"

P_THRESHOLD <- 0.05
F_THRESHOLD <- 5
CLUMP_KB <- 10000
CLUMP_R2 <- 0.1
CLUMP_P <- 1
N_WORKERS <- 4L
PLINK_THREADS_PER_JOB <- 1L
PLINK_MEMORY_MB <- 4096L
STOP_ON_FIRST_FAILED_BATCH <- TRUE
RESET_CLUMP_WORK <- FALSE  # TRUE removes this script's checkpoint directory.

EXPOSURE_FILTERED_FILE <- "CD4T_exposure_after_P_F_filter.csv"
EXPOSURE_CLUMP_INDEX_FILE <- "CD4T_exposure_LD_clump_index.csv"
EXPOSURE_CLUMPED_FILE <- "CD4T_exposure_after_LD_clumping.csv"

# Validate inputs before changing checkpoints.
setwd(path.expand(WORK_DIR))
plink_path <- plinkbinr::get_plink_exe()
panel_files <- paste0(EUR_BFILE, c(".bed", ".bim", ".fam"))
required_files <- c(EXPOSURE_FILE, plink_path, panel_files)
if (any(!file.exists(required_files))) {
  stop("Missing files: ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}
stopifnot(P_THRESHOLD > 0, P_THRESHOLD <= 1, F_THRESHOLD >= 0,
          CLUMP_R2 >= 0, CLUMP_R2 <= 1, CLUMP_KB > 0,
          CLUMP_P > 0, CLUMP_P <= 1, N_WORKERS >= 1,
          PLINK_THREADS_PER_JOB >= 1, PLINK_MEMORY_MB >= 1)
cores <- parallel::detectCores()
if (is.na(cores)) cores <- 1L
N_WORKERS <- as.integer(min(N_WORKERS, cores))
if (.Platform$OS.type == "windows") N_WORKERS <- 1L

# Read, clean and filter exposure data.
b <- fread(EXPOSURE_FILE)
required_cols <- c("id", "SNP", "chr", "pos", "other allele", "effect allele",
                   "maf", "beta", "se", "pval", "samplesize", "source")
missing_cols <- setdiff(required_cols, names(b))
if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))
setnames(b, c("id", "other allele", "effect allele"), c("gene", "OA_exp", "EA_exp"))
b[, `:=`(
  gene = trimws(as.character(gene)), SNP = trimws(as.character(SNP)),
  source = trimws(as.character(source)), chr = as.character(chr),
  pos = as.integer(pos), OA_exp = toupper(trimws(as.character(OA_exp))),
  EA_exp = toupper(trimws(as.character(EA_exp))), maf = as.numeric(maf),
  beta = as.numeric(beta), se = as.numeric(se), pval = as.numeric(pval),
  samplesize = as.numeric(samplesize)
)]
b <- b[
  !is.na(gene) & gene != "" & !is.na(source) & source != "" &
  !is.na(SNP) & SNP != "" & SNP != "." &
  is.finite(beta) & is.finite(se) & se > 0 &
  is.finite(pval) & pval >= 0 & pval <= 1 & pval < P_THRESHOLD &
  !is.na(chr) & !is.na(pos) &
  !is.na(EA_exp) & EA_exp != "" & !is.na(OA_exp) & OA_exp != ""
]
b[, F_statistic := (beta / se)^2]
b <- b[is.finite(F_statistic) & F_statistic >= F_THRESHOLD]
b[, phenotype_id := paste(gene, source, sep = "__")]
b[, exposure_id := phenotype_id]
setorder(b, phenotype_id, SNP, pval)
b <- unique(b, by = c("phenotype_id", "SNP"))
if (!nrow(b)) stop("No SNPs passed the P/F filters.")

# Match SNPs to the reference panel.
panel <- fread(paste0(EUR_BFILE, ".bim"), header = FALSE, select = 2,
               colClasses = "character", col.names = "rsid")
panel <- unique(panel[!is.na(rsid) & rsid != "" & rsid != "."])
clump_input <- b[, .(rsid = SNP, pval, id = phenotype_id)]
clump_input <- clump_input[panel, on = "rsid", nomatch = 0L]
setorder(clump_input, id, pval, rsid)
if (!nrow(clump_input)) stop("No SNPs matched the reference panel.")
setindex(clump_input, id)
all_ids <- sort(unique(clump_input$id))
rm(panel)
gc(verbose = FALSE)

# One atomic RDS result per exposure replaces separate result/done files.
# Resume only when the input, parameters, panel metadata and ID mapping match.
manifest <- list(
  version = 1L,
  exposure_md5 = unname(tools::md5sum(EXPOSURE_FILE)),
  plink_md5 = unname(tools::md5sum(plink_path)),
  panel_paths = normalizePath(panel_files, winslash = "/"),
  panel_metadata = file.info(panel_files)[, c("size", "mtime")],
  parameters = c(P_THRESHOLD, F_THRESHOLD, CLUMP_KB, CLUMP_R2, CLUMP_P),
  ids = all_ids
)
if (RESET_CLUMP_WORK && dir.exists(CLUMP_WORK_DIR)) {
  # Refuse broad or unexpected deletion targets.
  target <- normalizePath(CLUMP_WORK_DIR, winslash = "/", mustWork = TRUE)
  if (basename(target) != "plink_clump_sensitive3_simplified") {
    stop("Reset refused: use a new checkpoint directory instead.")
  }
  unlink(target, recursive = TRUE)
  if (dir.exists(target)) stop("Could not reset checkpoint directory.")
}
dir.create(CLUMP_WORK_DIR, recursive = TRUE, showWarnings = FALSE)
manifest_file <- file.path(CLUMP_WORK_DIR, "manifest.rds")
if (file.exists(manifest_file)) {
  if (!identical(readRDS(manifest_file), manifest)) {
    stop("Checkpoint inputs or settings changed. Use a new directory or reset checkpoints.")
  }
} else {
  if (length(list.files(CLUMP_WORK_DIR, all.files = TRUE, no.. = TRUE))) {
    stop("Checkpoint directory is not empty and has no manifest. Use a new directory.")
  }
  saveRDS(manifest, manifest_file)
}
job_map <- data.table(id = all_ids, tag = sprintf("job_%06d", seq_along(all_ids)))
job_map[, result_file := file.path(CLUMP_WORK_DIR, paste0(tag, ".rds"))]
fwrite(job_map, file.path(CLUMP_WORK_DIR, "job_map.csv"))
fwrite(b, EXPOSURE_FILTERED_FILE)

read_result <- function(i) {
  path <- job_map$result_file[i]
  if (!file.exists(path)) return(NULL)
  x <- readRDS(path)
  if (!is.data.table(x) || !identical(names(x), c("rsid", "pval", "id")) ||
      !nrow(x) || anyNA(x) || !all(x$id == job_map$id[i]) ||
      anyDuplicated(x$rsid)) {
    stop("Invalid checkpoint: ", path)
  }
  x
}
completed <- vapply(seq_len(nrow(job_map)), function(i) !is.null(read_result(i)), logical(1))
todo <- which(!completed)
message(sprintf("Exposures: %d; completed: %d; pending: %d; workers: %d",
                nrow(job_map), sum(completed), length(todo), N_WORKERS))

# Run PLINK. Successful jobs keep only their RDS result; failed jobs keep logs.
run_one_clump <- function(i) {
  tryCatch({
    current_id <- job_map$id[i]
    dat <- clump_input[id == current_id, .(rsid, pval)]
    prefix <- file.path(CLUMP_WORK_DIR, job_map$tag[i])
    assoc <- paste0(prefix, ".assoc")
    extract <- paste0(prefix, ".extract")
    clumped <- paste0(prefix, ".clumped")
    console <- paste0(prefix, ".console.log")
    # Prevent a failed retry from reading an old PLINK output.
    unlink(clumped)
    if (file.exists(clumped)) stop("Could not remove stale PLINK output.")
    fwrite(dat[, .(SNP = rsid, P = pval)], assoc, sep = "\t", quote = FALSE)
    fwrite(dat[, .(rsid)], extract, col.names = FALSE, quote = FALSE)
    args <- c(
      "--bfile", shQuote(EUR_BFILE), "--extract", shQuote(extract),
      "--clump", shQuote(assoc), "--clump-snp-field", "SNP", "--clump-field", "P",
      "--clump-p1", CLUMP_P, "--clump-p2", CLUMP_P,
      "--clump-r2", CLUMP_R2, "--clump-kb", CLUMP_KB,
      "--threads", PLINK_THREADS_PER_JOB, "--memory", PLINK_MEMORY_MB,
      "--out", shQuote(prefix)
    )
    status <- system2(plink_path, args, stdout = console, stderr = console,
                      env = paste0(c("OMP_NUM_THREADS=", "OPENBLAS_NUM_THREADS=",
                                     "MKL_NUM_THREADS="), PLINK_THREADS_PER_JOB))
    if (status != 0L || !file.exists(clumped)) {
      stop("PLINK failed (status ", status, "). See ", console)
    }
    selected <- fread(clumped, fill = TRUE)
    if (!"SNP" %in% names(selected)) stop("Missing SNP column in PLINK output.")
    kept <- unique(as.character(selected$SNP))
    kept <- kept[!is.na(kept) & kept != ""]
    if (!length(kept) || !all(kept %chin% dat$rsid)) {
      stop("Empty or unexpected index SNPs in PLINK output.")
    }
    result <- dat[rsid %chin% kept, .(rsid, pval, id = current_id)]
    path <- job_map$result_file[i]
    saveRDS(result, paste0(path, ".tmp"))
    if (!file.rename(paste0(path, ".tmp"), path)) stop("Could not commit checkpoint.")
    unlink(c(assoc, extract, clumped, console, paste0(prefix, c(".log", ".nosex"))))
    list(ok = TRUE, error = NA_character_)
  }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
}

# Small batches provide progress updates and stop promptly on failures.
failed <- data.table(id = character(), error = character())
failure_file <- file.path(CLUMP_WORK_DIR, "failed_ids.csv")
fwrite(failed, failure_file)
batch_size <- N_WORKERS * 2L
batches <- split(todo, ceiling(seq_along(todo) / batch_size))
for (batch in batches) {
  results <- parallel::mclapply(batch, run_one_clump,
                               mc.cores = min(N_WORKERS, length(batch)),
                               mc.preschedule = FALSE, mc.set.seed = FALSE)
  for (j in seq_along(batch)) {
    x <- results[[j]]
    if (is.null(x) || inherits(x, "try-error")) {
      x <- list(ok = FALSE, error = if (is.null(x)) "Worker exited without a result." else as.character(x))
    }
    i <- batch[j]
    completed[i] <- isTRUE(x$ok)
    if (!completed[i]) {
      failed <- rbind(failed, data.table(id = job_map$id[i], error = x$error))
    }
  }
  fwrite(failed, failure_file)
  message(sprintf("Completed: %d/%d; failed this run: %d",
                  sum(completed), nrow(job_map), nrow(failed)))
  if (nrow(failed) && STOP_ON_FIRST_FAILED_BATCH) {
    stop("Clumping failed. See ", failure_file, ". Rerun with RESET_CLUMP_WORK = FALSE.")
  }
}

# Publish final outputs only when every current exposure has a valid result.
if (!all(completed)) stop("Some jobs failed. Final outputs were not updated. See ", failure_file)
clump_result <- rbindlist(lapply(seq_len(nrow(job_map)), read_result))
if (!setequal(clump_result$id, all_ids)) stop("Incomplete clumping results.")
setorder(clump_result, id, pval, rsid)
clump_keep <- clump_result[, .(SNP = rsid, phenotype_id = id)]
b_clumped <- merge(b, clump_keep, by = c("SNP", "phenotype_id"), all = FALSE)
if (nrow(b_clumped) != nrow(clump_result)) stop("Exposure merge did not preserve all index SNPs.")
setorder(b_clumped, phenotype_id, pval, SNP)
fwrite(clump_result, EXPOSURE_CLUMP_INDEX_FILE)
fwrite(b_clumped, EXPOSURE_CLUMPED_FILE)
message(sprintf("Finished: %d index SNP rows across %d exposures.",
                nrow(b_clumped), uniqueN(b_clumped$phenotype_id)))

