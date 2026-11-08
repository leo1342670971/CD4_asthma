

library(ieugwasr)
library(data.table)
library(dplyr)
# -----------------------------
# 参数设置
# -----------------------------
input_file <- ""

# outcome 文件所在路径
outcome_dir <- ""

# outcome 与文件名对应关系
outcome_file_map <- c(
  "Asthma_GBMI" = "Asthma_GBMI.txt",
  "finngen_R12_J10_ASTHMA_EXMORE.gz" = "finngen_R12_J10_ASTHMA_EXMORE.gz.csv"
)

ld_bfile <- "E:/data/SMR/LD/bfile/EUR"
plink_bin <- "E:/data/SMR/LD/plink_win64/plink.exe"
output_csv <- "E:/data/SMR/LD_check_results.csv"

# 上下游窗口大小：1000kb
window_size <- 1000000

# -----------------------------
# 读取 input
# -----------------------------
input <- fread(input_file)

# -----------------------------
# 检查 input 列名
# 现在不再依赖 chr / pos_start / pos_end，但保留检查可改为只检查必要列
# -----------------------------
if(!all(c("exposure", "outcome", "SNP") %in% colnames(input))){
  stop("Input file must contain columns: exposure, outcome, SNP")
}

# -----------------------------
# 读取 EUR.bim，检查 SNP 是否在本地 LD 参考面板中
# -----------------------------
bim_file <- paste0(ld_bfile, ".bim")

if(!file.exists(bim_file)){
  stop("Cannot find bfile BIM file: ", bim_file)
}

eur_bim <- fread(bim_file, header = FALSE)
eur_snps <- unique(as.character(eur_bim$V2))

# -----------------------------
# 统一 chr 格式
# -----------------------------
clean_chr <- function(x){
  gsub("^chr", "", as.character(x), ignore.case = TRUE)
}

# -----------------------------
# 不建议只按 SNP 去重
# 因为同一个 SNP 可能对应不同 exposure / outcome
# -----------------------------
input[, SNP := as.character(SNP)]
input <- unique(input, by = c("exposure", "outcome", "SNP"))

# -----------------------------
# outcome 文件缓存，避免重复 fread
# -----------------------------
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

# -----------------------------
# 结果列表
# -----------------------------
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
  
  # -----------------------------
  # 0. 根据 outcome 读取对应 outcome 文件
  # -----------------------------
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
  
  # -----------------------------
  # 1. 先检查 query SNP 是否在 EUR 参考面板中
  # -----------------------------
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
  
  # -----------------------------
  # 2. 在 outcome 中匹配 query SNP
  #    用 outcome 自己的 chr / pos 定义上下游 1000kb 区域
  # -----------------------------
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
  
  # 如果 outcome 中 query SNP 重复，保留 p 值最小的一条
  query_outcome_pos <- query_outcome_pos[order(pval)][1]
  
  query_chr <- query_outcome_pos$chr_clean
  query_pos <- query_outcome_pos$pos
  
  region_start <- query_pos - window_size
  region_end <- query_pos + window_size
  
  # -----------------------------
  # 3. 基于 outcome 坐标提取 query SNP 上下游 1000kb SNP
  # -----------------------------
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
  
  # -----------------------------
  # 4. 按 outcome pval 筛选 SNP
  #    先按 pval 升序，再保留 pval < 1e-3
  #    如果 >=500，取前499个
  # -----------------------------
  region_snps <- region_snps_all[order(pval)]
  region_snps <- region_snps[pval < 1e-3]
  
  # 去掉 outcome 中重复 SNP，保留 p 值最小的那条
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
  
  # -----------------------------
  # 5. 区域 SNP + query SNP
  # -----------------------------
  snp_list <- unique(c(as.character(region_snps$SNP), rsid))
  
  n_snps_input_ld <- length(snp_list)
  n_snps_in_eur <- sum(snp_list %in% eur_snps)
  
  # -----------------------------
  # 6. 本地 LD matrix 计算
  # -----------------------------
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
  
  # -----------------------------
  # 7. 如果 LD 矩阵为空
  # -----------------------------
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
  
  # -----------------------------
  # 8. 找 query SNP 在 LD 矩阵中的列
  #    with_alleles = TRUE 时列名可能是 rs71356699_T_C
  # -----------------------------
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
  
  # -----------------------------
  # 9. 取 query SNP 这一列，平方为 r2
  # -----------------------------
  b <- (ld_mat[, col.index])^2
  names(b) <- rownames(ld_mat)
  b <- b[!is.na(b)]
  
  # -----------------------------
  # 10. 如果 LD 矩阵只有 query SNP 一个 SNP
  # -----------------------------
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
  
  # -----------------------------
  # 11. 排除 query SNP 自身后，找最大 r2
  # -----------------------------
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

# -----------------------------
# 保存结果
# -----------------------------
final_results <- bind_rows(results_list)

fwrite(final_results, output_csv)

cat("Local LD check finished! Results saved to:", output_csv, "\n")
