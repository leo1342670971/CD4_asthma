library(ieugwasr)
library(data.table)
library(dplyr)
input_file <- ""

# outcome file
outcome_dir <- ""

# outcome 
outcome_file_map <- c(
  "Asthma_GBMI" = "Asthma_GBMI.txt",
  "finngen_R12_J10_ASTHMA_EXMORE.gz" = "finngen_R12_J10_ASTHMA_EXMORE.gz.csv"
)

ld_bfile <- "E:/data/SMR/LD/bfile/EUR"
plink_bin <- "E:/data/SMR/LD/plink_win64/plink.exe"
output_csv <- "E:/data/SMR/LD_check_results.csv"

#change
window_size <- 1000000

input <- fread(input_file)

if(!all(c("exposure", "outcome", "SNP") %in% colnames(input))){
  stop("Input file must contain columns: exposure, outcome, SNP")
}

bim_file <- paste0(ld_bfile, ".bim")

if(!file.exists(bim_file)){
  stop("Cannot find bfile BIM file: ", bim_file)
}

eur_bim <- fread(bim_file, header = FALSE)
eur_snps <- unique(as.character(eur_bim$V2))

clean_chr <- function(x){
  gsub("^chr", "", as.character(x), ignore.case = TRUE)
}

input[, SNP := as.character(SNP)]
input <- unique(input, by = c("exposure", "outcome", "SNP"))

outcome_cache <- list()

load_outcome_data <- function(outcome_name){
  
  outcome_name <- as.character(outcome_name)
  
  if(!outcome_name %in% names(outcome_file_map)){
    stop(
      "Unknown outcome: ", outcome_name,
      ". Please add it to outcome_file_map."
    )
  }
  
  if(!is.null(outcome_cache[[outcome_name]])){
    return(outcome_cache[[outcome_name]])
  }
  
  outcome_path <- file.path(outcome_dir, outcome_file_map[[outcome_name]])
  
  if(!file.exists(outcome_path)){
    stop("Cannot find outcome file: ", outcome_path)
  }
  
  dat <- fread(outcome_path)
  
  if(!all(c("SNP", "chr", "pos", "pval") %in% colnames(dat))){
    stop(
      "Outcome file ", outcome_path,
      " must contain columns: SNP, chr, pos, pval"
    )
  }
  
  dat[, SNP := as.character(SNP)]
  dat[, chr_clean := clean_chr(chr)]
  dat[, pos := as.numeric(pos)]
  dat[, pval := as.numeric(pval)]
  
  dat <- dat[!is.na(SNP) & !is.na(chr_clean) & !is.na(pos)]
  
  outcome_cache[[outcome_name]] <<- dat
  
  return(dat)
}

results_list <- list()

for(i in 1:nrow(input)){
  
  row_i <- input[i, ]
  
  rsid <- as.character(row_i$SNP)
  exposure_id <- as.character(row_i$exposure)
  outcome_id <- as.character(row_i$outcome)
  
  cat(
    "Processing", i, "/", nrow(input),
    "| SNP:", rsid,
    "| outcome:", outcome_id, "\n"
  )
  
  outcome_dat <- tryCatch({
    load_outcome_data(outcome_id)
  }, error = function(e){
    
    results_list[[i]] <<- data.frame(
      SNP_query = rsid,
      exposure = exposure_id,
      outcome = outcome_id,
      query_outcome_chr = NA,
      query_outcome_pos = NA,
      region_start = NA,
      region_end = NA,
      SNP_ld_max = NA,
      R2 = NA,
      n_region_before_pfilter = NA,
      n_region_after_pfilter = NA,
      n_snps_input_ld = NA,
      n_snps_in_eur = NA,
      status = "outcome_file_error",
      note = e$message
    )
    
    return(NULL)
  })
  
  if(is.null(outcome_dat)){
    next
  }
  query_in_eur <- rsid %in% eur_snps
  
  if(!query_in_eur){
    results_list[[i]] <- data.frame(
      SNP_query = rsid,
      exposure = exposure_id,
      outcome = outcome_id,
      query_outcome_chr = NA,
      query_outcome_pos = NA,
      region_start = NA,
      region_end = NA,
      SNP_ld_max = NA,
      R2 = NA,
      n_region_before_pfilter = NA,
      n_region_after_pfilter = NA,
      n_snps_input_ld = NA,
      n_snps_in_eur = NA,
      status = "query_not_in_EUR_reference",
      note = "Query SNP is absent from EUR.bim, LD cannot be calculated"
    )
    next
  }
  
  query_outcome_pos <- outcome_dat[SNP == rsid]
  
  if(nrow(query_outcome_pos) == 0){
    results_list[[i]] <- data.frame(
      SNP_query = rsid,
      exposure = exposure_id,
      outcome = outcome_id,
      query_outcome_chr = NA,
      query_outcome_pos = NA,
      region_start = NA,
      region_end = NA,
      SNP_ld_max = NA,
      R2 = NA,
      n_region_before_pfilter = NA,
      n_region_after_pfilter = NA,
      n_snps_input_ld = NA,
      n_snps_in_eur = NA,
      status = "query_snp_not_found_in_outcome",
      note = "Query SNP was not found in outcome file, cannot define outcome-based 1000kb region"
    )
    next
  }
  query_outcome_pos <- query_outcome_pos[order(pval)][1]
  
  query_chr <- query_outcome_pos$chr_clean
  query_pos <- query_outcome_pos$pos
  
  region_start <- query_pos - window_size
  region_end <- query_pos + window_size

  region_snps_all <- outcome_dat[
    chr_clean == query_chr &
      pos >= region_start &
      pos <= region_end
  ]
  
  n_region_before_pfilter <- nrow(region_snps_all)
  
  if(n_region_before_pfilter == 0){
    results_list[[i]] <- data.frame(
      SNP_query = rsid,
      exposure = exposure_id,
      outcome = outcome_id,
      query_outcome_chr = query_chr,
      query_outcome_pos = query_pos,
      region_start = region_start,
      region_end = region_end,
      SNP_ld_max = NA,
      R2 = NA,
      n_region_before_pfilter = 0,
      n_region_after_pfilter = 0,
      n_snps_input_ld = NA,
      n_snps_in_eur = NA,
      status = "no_outcome_snp_in_1000kb_region",
      note = "No outcome SNP found in outcome-based 1000kb region"
    )
    next
  }

  region_snps <- region_snps_all[order(pval)]
  region_snps <- region_snps[pval < 1e-3]
  region_snps <- region_snps[!duplicated(SNP)]
  
  n_region_after_pfilter <- nrow(region_snps)
  
  if(n_region_after_pfilter == 0){
    results_list[[i]] <- data.frame(
      SNP_query = rsid,
      exposure = exposure_id,
      outcome = outcome_id,
      query_outcome_chr = query_chr,
      query_outcome_pos = query_pos,
      region_start = region_start,
      region_end = region_end,
      SNP_ld_max = NA,
      R2 = NA,
      n_region_before_pfilter = n_region_before_pfilter,
      n_region_after_pfilter = 0,
      n_snps_input_ld = NA,
      n_snps_in_eur = NA,
      status = "no_outcome_snp_after_pfilter",
      note = "Outcome SNPs exist in 1000kb region, but none passed pval < 1e-3"
    )
    next
  }
  
  if(nrow(region_snps) >= 500){
    region_snps <- region_snps[1:499]
  }
  
  snp_list <- unique(c(as.character(region_snps$SNP), rsid))
  
  n_snps_input_ld <- length(snp_list)
  n_snps_in_eur <- sum(snp_list %in% eur_snps)

  ld_mat <- NULL
  attempts <- 0
  
  while(attempts <= 5){
    
    ld_mat <- tryCatch({
      ld_matrix(
        variants = snp_list,
        with_alleles = TRUE,
        pop = "EUR",
        bfile = ld_bfile,
        plink_bin = plink_bin
      )
    }, error = function(e){
      warning(paste("LD calculation failed for SNP", rsid, ":", e$message))
      return(NULL)
    })
    
    if(!is.null(ld_mat)){
      break
    }
    
    attempts <- attempts + 1
    Sys.sleep(1)
  }
  
  if(is.null(ld_mat) || is.null(nrow(ld_mat))){
    results_list[[i]] <- data.frame(
      SNP_query = rsid,
      exposure = exposure_id,
      outcome = outcome_id,
      query_outcome_chr = query_chr,
      query_outcome_pos = query_pos,
      region_start = region_start,
      region_end = region_end,
      SNP_ld_max = rsid,
      R2 = 1,
      n_region_before_pfilter = n_region_before_pfilter,
      n_region_after_pfilter = n_region_after_pfilter,
      n_snps_input_ld = n_snps_input_ld,
      n_snps_in_eur = n_snps_in_eur,
      status = "ld_matrix_empty",
      note = "LD matrix is empty, output query SNP itself with R2=1 according to original logic"
    )
    next
  }

  query_pattern <- paste0("^", rsid, "($|_)")
  col.index <- grep(query_pattern, colnames(ld_mat))
  
  if(length(col.index) == 0){
    results_list[[i]] <- data.frame(
      SNP_query = rsid,
      exposure = exposure_id,
      outcome = outcome_id,
      query_outcome_chr = query_chr,
      query_outcome_pos = query_pos,
      region_start = region_start,
      region_end = region_end,
      SNP_ld_max = NA,
      R2 = NA,
      n_region_before_pfilter = n_region_before_pfilter,
      n_region_after_pfilter = n_region_after_pfilter,
      n_snps_input_ld = n_snps_input_ld,
      n_snps_in_eur = n_snps_in_eur,
      status = "query_not_in_ld_matrix",
      note = "Query SNP was in EUR.bim but not found in returned LD matrix"
    )
    next
  }
  
  col.index <- col.index[1]
  query_bfile_name <- colnames(ld_mat)[col.index]

  b <- (ld_mat[, col.index])^2
  names(b) <- rownames(ld_mat)
  b <- b[!is.na(b)]

  if(length(b) <= 1){
    results_list[[i]] <- data.frame(
      SNP_query = rsid,
      exposure = exposure_id,
      outcome = outcome_id,
      query_outcome_chr = query_chr,
      query_outcome_pos = query_pos,
      region_start = region_start,
      region_end = region_end,
      SNP_ld_max = query_bfile_name,
      R2 = 1,
      n_region_before_pfilter = n_region_before_pfilter,
      n_region_after_pfilter = n_region_after_pfilter,
      n_snps_input_ld = n_snps_input_ld,
      n_snps_in_eur = n_snps_in_eur,
      status = "only_query_snp_in_ld_matrix",
      note = "Only query SNP remained in EUR LD matrix"
    )
    next
  }

  b_no_self <- b[names(b) != query_bfile_name]
  
  if(length(b_no_self) == 0){
    results_list[[i]] <- data.frame(
      SNP_query = rsid,
      exposure = exposure_id,
      outcome = outcome_id,
      query_outcome_chr = query_chr,
      query_outcome_pos = query_pos,
      region_start = region_start,
      region_end = region_end,
      SNP_ld_max = query_bfile_name,
      R2 = 1,
      n_region_before_pfilter = n_region_before_pfilter,
      n_region_after_pfilter = n_region_after_pfilter,
      n_snps_input_ld = n_snps_input_ld,
      n_snps_in_eur = n_snps_in_eur,
      status = "no_other_snp_after_excluding_query",
      note = "No other SNP remained after excluding query SNP"
    )
    next
  }
  
  b_no_self <- b_no_self[order(b_no_self)]
  b_max <- b_no_self[length(b_no_self)]
  
  ld_snp_name <- names(b_max)
  ld_r2 <- as.numeric(b_max)
  
  results_list[[i]] <- data.frame(
    SNP_query = rsid,
    exposure = exposure_id,
    outcome = outcome_id,
    query_outcome_chr = query_chr,
    query_outcome_pos = query_pos,
    region_start = region_start,
    region_end = region_end,
    SNP_ld_max = ld_snp_name,
    R2 = ld_r2,
    n_region_before_pfilter = n_region_before_pfilter,
    n_region_after_pfilter = n_region_after_pfilter,
    n_snps_input_ld = n_snps_input_ld,
    n_snps_in_eur = n_snps_in_eur,
    status = "success",
    note = "LD calculated successfully using outcome-based 1000kb region"
  )
}
final_results <- bind_rows(results_list)

fwrite(final_results, output_csv)

cat("Local LD check finished! Results saved to:", output_csv, "\n")
