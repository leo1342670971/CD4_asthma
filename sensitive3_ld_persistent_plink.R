library(data.table)
library(plinkbinr)
library(parallel)

# ============================================================
# 0. 参数设置
# ============================================================

WORK_DIR <- "~/lj74/CD4T"
setwd(WORK_DIR)

EXPOSURE_FILE <- "CD4all1.8.txt"

# 暴露筛选参数
P_THRESHOLD <- 0.05
F_THRESHOLD <- 5

# LD clumping参数
CLUMP_KB <- 10000
CLUMP_R2 <- 0.1
CLUMP_P  <- 1

# 1000G EUR参考面板共同前缀
EUR_BFILE <- "/data/h004/lj74/1000G.EUR/1000G_EUR"

# 使用稳定的持久化工作目录，避免/tmp临时文件被清理或无法读取
# 特意避免路径中出现空格
CLUMP_WORK_DIR <- "/data/h004/lj74/Rbao/plink_clump_sensitive3_parallel"

# 并行参数：账号目前可见8核时，建议先用4个并行PLINK任务
# 不建议一开始直接用满8核，因为多个PLINK进程会同时读取大型1000G参考面板，
# 速度可能受磁盘I/O限制，且内存占用会随并行数增加。
N_WORKERS <- 4L
PLINK_THREADS_PER_JOB <- 1L
PLINK_MEMORY_MB <- 4096L
PARALLEL_BATCH_SIZE <- N_WORKERS * 2L
STOP_ON_FIRST_FAILED_BATCH <- TRUE

# TRUE：删除旧进度，从头开始
# FALSE：读取completed_ids.txt并从上次完成位置继续
RESET_CLUMP_WORK <- TRUE

# 输出文件
EXPOSURE_FILTERED_FILE <- "CD4T_exposure_after_P_F_filter.csv"
EXPOSURE_CLUMP_INDEX_FILE <- "CD4T_exposure_LD_clump_index.csv"
EXPOSURE_CLUMPED_FILE <- "CD4T_exposure_after_LD_clumping.csv"

# 工作目录中的断点文件
CLUMP_PROGRESS_FILE <- file.path(CLUMP_WORK_DIR, "clump_progress.csv")
CLUMP_DONE_FILE <- file.path(CLUMP_WORK_DIR, "completed_ids.txt")
CLUMP_FAILED_FILE <- file.path(CLUMP_WORK_DIR, "failed_ids.csv")


# ============================================================
# 1. 辅助函数
# ============================================================

format_duration <- function(seconds) {
  if (is.na(seconds) || !is.finite(seconds) || seconds < 0) {
    return("NA")
  }
  
  seconds <- round(seconds)
  days <- seconds %/% 86400
  hours <- (seconds %% 86400) %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60
  
  sprintf(
    "%d天 %02d小时 %02d分钟 %02d秒",
    days, hours, minutes, secs
  )
}


# ============================================================
# 2. 初始化并检查PLINK工作目录
# ============================================================

if (RESET_CLUMP_WORK && dir.exists(CLUMP_WORK_DIR)) {
  message("RESET_CLUMP_WORK=TRUE：删除旧的clumping工作目录。")
  unlink(CLUMP_WORK_DIR, recursive = TRUE, force = TRUE)
}

dir.create(
  CLUMP_WORK_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(CLUMP_WORK_DIR)) {
  stop("无法创建PLINK工作目录：", CLUMP_WORK_DIR)
}

# 测试工作目录读写
write_test <- file.path(CLUMP_WORK_DIR, "write_test.txt")
writeLines("PLINK work directory write test", write_test)

if (!file.exists(write_test)) {
  stop("PLINK工作目录无法写入：", CLUMP_WORK_DIR)
}

unlink(write_test)

cat("PLINK持久化工作目录：", CLUMP_WORK_DIR, "\n")
cat("断点续跑模式：", ifelse(RESET_CLUMP_WORK, "从头开始", "继续上次任务"), "\n")


# ============================================================
# 3. 读取并筛选暴露数据
# ============================================================

b <- fread(EXPOSURE_FILE, header = TRUE)
setDT(b)

required_exp_cols <- c(
  "id", "SNP", "chr", "pos",
  "other allele", "effect allele",
  "maf", "beta", "se", "pval",
  "samplesize", "source"
)

missing_exp_cols <- setdiff(required_exp_cols, names(b))

if (length(missing_exp_cols) > 0) {
  stop(
    "暴露文件缺少以下列：",
    paste(missing_exp_cols, collapse = ", ")
  )
}

setnames(
  b,
  old = c("id", "other allele", "effect allele"),
  new = c("gene", "OA_exp", "EA_exp")
)

b[, `:=`(
  gene = trimws(as.character(gene)),
  SNP = trimws(as.character(SNP)),
  chr = as.character(chr),
  pos = as.integer(pos),
  OA_exp = toupper(trimws(as.character(OA_exp))),
  EA_exp = toupper(trimws(as.character(EA_exp))),
  maf = as.numeric(maf),
  beta = as.numeric(beta),
  se = as.numeric(se),
  pval = as.numeric(pval),
  samplesize = as.numeric(samplesize),
  source = trimws(as.character(source))
)]

b <- b[
  !is.na(gene) & gene != "" &
    !is.na(source) & source != "" &
    !is.na(SNP) & SNP != "" & SNP != "." &
    !is.na(beta) &
    !is.na(se) & se > 0 &
    !is.na(pval) &
    !is.na(chr) & !is.na(pos) &
    !is.na(EA_exp) & EA_exp != "" &
    !is.na(OA_exp) & OA_exp != ""
]

b[, F_statistic := (beta / se)^2]

b <- b[
  pval < P_THRESHOLD &
    is.finite(F_statistic) &
    F_statistic >= F_THRESHOLD
]

# gene × cell/state定义为一个暴露
b[, phenotype_id := paste(gene, source, sep = "__")]

# 固定暴露ID，避免后续format_data随机生成ID
b[, exposure_id := phenotype_id]

setorder(b, phenotype_id, SNP, pval)
b <- unique(b, by = c("phenotype_id", "SNP"))

cat("P/F筛选后暴露行数：", nrow(b), "\n")
cat("P/F筛选后暴露组合数：", uniqueN(b$phenotype_id), "\n")
cat("P/F筛选后唯一SNP数：", uniqueN(b$SNP), "\n")

if (nrow(b) == 0) {
  stop("P值和F统计量筛选后没有保留任何暴露SNP。")
}

fwrite(b, EXPOSURE_FILTERED_FILE)


# ============================================================
# 4. 检查PLINK及1000G EUR参考面板
# ============================================================

plink_path <- plinkbinr::get_plink_exe()

EUR_BED <- paste0(EUR_BFILE, ".bed")
EUR_BIM <- paste0(EUR_BFILE, ".bim")
EUR_FAM <- paste0(EUR_BFILE, ".fam")

cat("PLINK路径：", plink_path, "\n")
cat("BED文件：", EUR_BED, "\n")
cat("BIM文件：", EUR_BIM, "\n")
cat("FAM文件：", EUR_FAM, "\n")

if (!file.exists(plink_path)) {
  stop("找不到PLINK可执行文件：", plink_path)
}

if (!file.exists(EUR_BED)) {
  stop("找不到BED文件：", EUR_BED)
}

if (!file.exists(EUR_BIM)) {
  stop("找不到BIM文件：", EUR_BIM)
}

if (!file.exists(EUR_FAM)) {
  stop("找不到FAM文件：", EUR_FAM)
}

eur_snp <- fread(
  EUR_BIM,
  header = FALSE,
  select = 2,
  col.names = "rsid",
  showProgress = TRUE
)

eur_snp[, rsid := as.character(rsid)]
eur_snp <- unique(eur_snp[!is.na(rsid) & rsid != "" & rsid != "."])

cat("1000G EUR唯一rsID数：", nrow(eur_snp), "\n")


# ============================================================
# 5. 构建clumping输入并匹配参考面板
# ============================================================

clump_input <- unique(
  b[
    ,
    .(
      rsid = as.character(SNP),
      pval = as.numeric(pval),
      id = as.character(phenotype_id)
    )
  ]
)

clump_input <- clump_input[
  !is.na(rsid) &
    rsid != "" &
    rsid != "." &
    !is.na(pval) &
    is.finite(pval) &
    pval >= 0 &
    pval <= 1
]

cat("暴露LD clumping原始输入行数：", nrow(clump_input), "\n")
cat("暴露LD clumping原始暴露数：", uniqueN(clump_input$id), "\n")

n_before_panel <- nrow(clump_input)

# data.table连接，避免为超大向量反复执行%in%
setkey(eur_snp, rsid)
setkey(clump_input, rsid)

clump_input <- clump_input[
  eur_snp,
  on = "rsid",
  nomatch = 0L
]

# 恢复需要的列顺序并去重
clump_input <- unique(
  clump_input[
    ,
    .(
      rsid,
      pval,
      id
    )
  ]
)

setorder(clump_input, id, pval, rsid)

cat("匹配1000G EUR后clumping输入行数：", nrow(clump_input), "\n")
cat("未匹配参考面板的输入行数：", n_before_panel - nrow(clump_input), "\n")
cat("匹配1000G EUR后暴露数：", uniqueN(clump_input$id), "\n")

if (nrow(clump_input) == 0) {
  stop(
    "暴露SNP均未匹配到1000G EUR参考面板，",
    "请检查rsID或参考面板版本。"
  )
}

id_summary <- clump_input[
  ,
  .(
    n_SNP = uniqueN(rsid)
  ),
  by = id
]

setorder(id_summary, -n_SNP)

cat("每个暴露SNP数概况：\n")
print(summary(id_summary$n_SNP))
cat("SNP数最多的前10个暴露：\n")
print(head(id_summary, 10))

# 参考面板rsID已完成匹配，进入并行前释放超大对象，减少fork进程内存负担。
rm(eur_snp)
gc(verbose = FALSE)


# ============================================================
# 6. 初始化并行任务目录与断点信息
# ============================================================

CLUMP_INPUT_DIR <- file.path(CLUMP_WORK_DIR, "inputs")
CLUMP_OUTPUT_DIR <- file.path(CLUMP_WORK_DIR, "outputs")
CLUMP_RESULT_DIR <- file.path(CLUMP_WORK_DIR, "results")
CLUMP_DONE_DIR <- file.path(CLUMP_WORK_DIR, "done")
CLUMP_FAILURE_DIR <- file.path(CLUMP_WORK_DIR, "failed")
CLUMP_JOB_MAP_FILE <- file.path(CLUMP_WORK_DIR, "job_map.csv")

for (d in c(
  CLUMP_INPUT_DIR,
  CLUMP_OUTPUT_DIR,
  CLUMP_RESULT_DIR,
  CLUMP_DONE_DIR,
  CLUMP_FAILURE_DIR
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  
  if (!dir.exists(d)) {
    stop("无法创建任务目录：", d)
  }
}

available_cores <- parallel::detectCores(logical = TRUE)

if (is.na(available_cores) || available_cores < 1L) {
  available_cores <- 1L
}

N_WORKERS <- max(
  1L,
  min(
    as.integer(N_WORKERS),
    as.integer(available_cores)
  )
)

PARALLEL_BATCH_SIZE <- max(
  N_WORKERS,
  as.integer(PARALLEL_BATCH_SIZE)
)

cat("系统可见逻辑核心数：", available_cores, "\n")
cat("并行PLINK任务数：", N_WORKERS, "\n")
cat("每个PLINK任务线程数：", PLINK_THREADS_PER_JOB, "\n")
cat("每个PLINK任务内存上限：", PLINK_MEMORY_MB, "MB\n")
cat("每批提交任务数：", PARALLEL_BATCH_SIZE, "\n")

all_ids <- sort(unique(clump_input$id))

job_map <- data.table(
  global_index = seq_along(all_ids),
  file_tag = sprintf("job_%06d", seq_along(all_ids)),
  id = all_ids
)

fwrite(job_map, CLUMP_JOB_MAP_FILE)

# 每个成功任务都会产生独立.done文件；因此不会发生多个worker同时追加同一文件的问题。
done_files <- list.files(
  CLUMP_DONE_DIR,
  pattern = "^job_[0-9]+\\.done$",
  full.names = TRUE
)

# 只有.done标记和对应.result.csv同时存在时，才认定该任务已完成。
if (length(done_files) > 0L) {
  expected_result_files <- file.path(
    CLUMP_RESULT_DIR,
    sub(
      "\\.done$",
      ".result.csv",
      basename(done_files)
    )
  )
  
  done_files <- done_files[file.exists(expected_result_files)]
}

completed_ids <- character(0)

if (!RESET_CLUMP_WORK && length(done_files) > 0L) {
  completed_ids <- unique(
    trimws(
      unlist(
        lapply(
          done_files,
          readLines,
          warn = FALSE
        ),
        use.names = FALSE
      )
    )
  )
  
  completed_ids <- completed_ids[
    !is.na(completed_ids) & completed_ids != ""
  ]
}

todo_ids <- setdiff(all_ids, completed_ids)

cat("全部暴露组合数：", length(all_ids), "\n")
cat("断点中已经完成的暴露组合数：", length(completed_ids), "\n")
cat("本次需要处理的暴露组合数：", length(todo_ids), "\n")

if (length(todo_ids) == 0L) {
  message("所有暴露组合均已完成，直接汇总断点结果。")
}


# ============================================================
# 7. 定义单个暴露的PLINK clumping任务
#
# 并行策略：
# - 不是依赖PLINK --clump内部多线程；
# - 而是同时运行多个彼此独立的phenotype_id；
# - 每个worker只写自己的文件，避免并发写入冲突。
# ============================================================

run_one_clump <- function(current_id) {
  
  task_started <- Sys.time()
  
  map_row <- job_map[id == current_id]
  
  if (nrow(map_row) != 1L) {
    return(list(
      ok = FALSE,
      id = current_id,
      error = "无法在job_map中唯一定位当前暴露",
      n_input = NA_integer_,
      n_kept = NA_integer_,
      elapsed_seconds = 0
    ))
  }
  
  file_tag <- map_row$file_tag[[1]]
  
  current_dat <- clump_input[
    id == current_id,
    .(
      rsid,
      pval
    )
  ]
  
  current_dat <- unique(current_dat, by = "rsid")
  setorder(current_dat, pval, rsid)
  
  assoc_file <- file.path(
    CLUMP_INPUT_DIR,
    paste0(file_tag, ".assoc")
  )
  
  extract_file <- file.path(
    CLUMP_INPUT_DIR,
    paste0(file_tag, ".extract")
  )
  
  out_prefix <- file.path(
    CLUMP_OUTPUT_DIR,
    file_tag
  )
  
  plink_console_log <- paste0(out_prefix, ".console.log")
  plink_native_log <- paste0(out_prefix, ".log")
  clumped_file <- paste0(out_prefix, ".clumped")
  
  result_file <- file.path(
    CLUMP_RESULT_DIR,
    paste0(file_tag, ".result.csv")
  )
  
  done_file <- file.path(
    CLUMP_DONE_DIR,
    paste0(file_tag, ".done")
  )
  
  failure_file <- file.path(
    CLUMP_FAILURE_DIR,
    paste0(file_tag, ".failure.csv")
  )
  
  tryCatch({
    
    fwrite(
      current_dat[
        ,
        .(
          SNP = rsid,
          P = pval
        )
      ],
      assoc_file,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )
    
    fwrite(
      current_dat[
        ,
        .(
          rsid
        )
      ],
      extract_file,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    
    if (
      !file.exists(assoc_file) ||
      is.na(file.info(assoc_file)$size) ||
      file.info(assoc_file)$size <= 0
    ) {
      stop("无法创建PLINK clump输入文件：", assoc_file)
    }
    
    if (
      !file.exists(extract_file) ||
      is.na(file.info(extract_file)$size) ||
      file.info(extract_file)$size <= 0
    ) {
      stop("无法创建PLINK extract输入文件：", extract_file)
    }
    
    plink_args <- c(
      "--bfile", EUR_BFILE,
      "--extract", extract_file,
      "--clump", assoc_file,
      "--clump-snp-field", "SNP",
      "--clump-field", "P",
      "--clump-p1", format(CLUMP_P, scientific = FALSE),
      "--clump-p2", format(CLUMP_P, scientific = FALSE),
      "--clump-r2", format(CLUMP_R2, scientific = FALSE),
      "--clump-kb", as.character(CLUMP_KB),
      "--threads", as.character(PLINK_THREADS_PER_JOB),
      "--memory", as.character(PLINK_MEMORY_MB),
      "--out", out_prefix
    )
    
    status <- system2(
      command = plink_path,
      args = plink_args,
      stdout = plink_console_log,
      stderr = plink_console_log,
      wait = TRUE,
      env = c(
        paste0("OMP_NUM_THREADS=", PLINK_THREADS_PER_JOB),
        paste0("OPENBLAS_NUM_THREADS=", PLINK_THREADS_PER_JOB),
        paste0("MKL_NUM_THREADS=", PLINK_THREADS_PER_JOB)
      )
    )
    
    if (
      !identical(as.integer(status), 0L) ||
      !file.exists(clumped_file)
    ) {
      stop(
        "PLINK退出状态=", status,
        "；未生成预期的.clumped文件：", clumped_file
      )
    }
    
    clumped_dat <- fread(
      clumped_file,
      header = TRUE,
      fill = TRUE,
      data.table = TRUE
    )
    
    if (!"SNP" %in% names(clumped_dat)) {
      stop("PLINK .clumped结果中缺少SNP列：", clumped_file)
    }
    
    kept_rsids <- unique(
      as.character(
        clumped_dat[
          !is.na(SNP) & SNP != "",
          SNP
        ]
      )
    )
    
    if (length(kept_rsids) == 0L) {
      stop("PLINK运行成功，但没有读取到任何index SNP")
    }
    
    result_part <- current_dat[
      rsid %chin% kept_rsids,
      .(
        rsid,
        pval,
        id = current_id
      )
    ]
    
    setorder(result_part, pval, rsid)
    
    # 先写结果，再写.done标记；保证断点不会把未落盘结果视为已完成。
    fwrite(result_part, result_file)
    writeLines(current_id, done_file)
    
    elapsed_seconds <- as.numeric(
      difftime(
        Sys.time(),
        task_started,
        units = "secs"
      )
    )
    
    # 成功后清理大部分中间文件；result和done文件保留用于断点汇总。
    success_files <- c(
      assoc_file,
      extract_file,
      clumped_file,
      plink_console_log,
      plink_native_log,
      paste0(out_prefix, ".nosex")
    )
    
    unlink(
      success_files[file.exists(success_files)],
      force = TRUE
    )
    
    list(
      ok = TRUE,
      id = current_id,
      error = NA_character_,
      n_input = nrow(current_dat),
      n_kept = nrow(result_part),
      elapsed_seconds = elapsed_seconds,
      result_file = result_file
    )
    
  }, error = function(e) {
    
    elapsed_seconds <- as.numeric(
      difftime(
        Sys.time(),
        task_started,
        units = "secs"
      )
    )
    
    failure_row <- data.table(
      time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      id = current_id,
      n_SNP = nrow(current_dat),
      error = conditionMessage(e),
      assoc_file = assoc_file,
      extract_file = extract_file,
      plink_console_log = plink_console_log,
      plink_native_log = plink_native_log,
      expected_clumped_file = clumped_file
    )
    
    fwrite(failure_row, failure_file)
    
    list(
      ok = FALSE,
      id = current_id,
      error = conditionMessage(e),
      n_input = nrow(current_dat),
      n_kept = NA_integer_,
      elapsed_seconds = elapsed_seconds,
      failure_file = failure_file
    )
  })
}


# ============================================================
# 8. 分批并行运行、输出进度和ETA
# ============================================================

run_start <- Sys.time()
completed_this_run <- 0L
failed_this_run <- list()

if (length(todo_ids) > 0L) {
  
  batch_starts <- seq(
    from = 1L,
    to = length(todo_ids),
    by = PARALLEL_BATCH_SIZE
  )
  
  for (batch_number in seq_along(batch_starts)) {
    
    batch_start_index <- batch_starts[batch_number]
    batch_end_index <- min(
      batch_start_index + PARALLEL_BATCH_SIZE - 1L,
      length(todo_ids)
    )
    
    batch_ids <- todo_ids[
      batch_start_index:batch_end_index
    ]
    
    cat(
      sprintf(
        "\n========== 批次 %d/%d：提交%d个暴露，使用%d个并行worker ==========\n",
        batch_number,
        length(batch_starts),
        length(batch_ids),
        min(N_WORKERS, length(batch_ids))
      )
    )
    
    batch_started <- Sys.time()
    
    batch_results <- parallel::mclapply(
      X = batch_ids,
      FUN = run_one_clump,
      mc.cores = min(N_WORKERS, length(batch_ids)),
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
    
    # 防止mclapply中的worker异常直接返回try-error。
    batch_results <- lapply(
      seq_along(batch_results),
      function(j) {
        result_j <- batch_results[[j]]
        
        if (inherits(result_j, "try-error")) {
          return(list(
            ok = FALSE,
            id = batch_ids[j],
            error = as.character(result_j),
            n_input = NA_integer_,
            n_kept = NA_integer_,
            elapsed_seconds = NA_real_
          ))
        }
        
        result_j
      }
    )
    
    ok_vec <- vapply(
      batch_results,
      function(x) isTRUE(x$ok),
      logical(1)
    )
    
    n_ok <- sum(ok_vec)
    n_failed <- sum(!ok_vec)
    
    completed_this_run <- completed_this_run + n_ok
    
    if (n_failed > 0L) {
      failed_this_run <- c(
        failed_this_run,
        batch_results[!ok_vec]
      )
    }
    
    batch_elapsed <- as.numeric(
      difftime(
        Sys.time(),
        batch_started,
        units = "secs"
      )
    )
    
    total_elapsed <- as.numeric(
      difftime(
        Sys.time(),
        run_start,
        units = "secs"
      )
    )
    
    attempted_this_run <- min(
      batch_end_index,
      length(todo_ids)
    )
    
    avg_seconds_per_attempt <- total_elapsed / attempted_this_run
    remaining_attempts <- length(todo_ids) - attempted_this_run
    estimated_remaining <- avg_seconds_per_attempt * remaining_attempts
    
    cat(
      sprintf(
        paste0(
          "批次完成：成功=%d，失败=%d，批次耗时=%s\n",
          "本次已尝试=%d/%d；累计成功=%d；累计失败=%d\n",
          "总运行时间=%s；预计剩余=%s\n"
        ),
        n_ok,
        n_failed,
        format_duration(batch_elapsed),
        attempted_this_run,
        length(todo_ids),
        completed_this_run,
        length(failed_this_run),
        format_duration(total_elapsed),
        format_duration(estimated_remaining)
      )
    )
    
    if (n_failed > 0L) {
      cat("本批次失败的暴露：\n")
      
      for (failed_item in batch_results[!ok_vec]) {
        cat(
          "  - ",
          failed_item$id,
          "：",
          failed_item$error,
          "\n",
          sep = ""
        )
      }
      
      if (STOP_ON_FIRST_FAILED_BATCH) {
        stop(
          "并行clumping批次出现失败任务，已停止后续批次。",
          "失败详情保存在：",
          CLUMP_FAILURE_DIR
        )
      }
    }
    
    gc(verbose = FALSE)
  }
}


# ============================================================
# 9. 汇总所有独立结果文件
# ============================================================

result_files <- list.files(
  CLUMP_RESULT_DIR,
  pattern = "^job_[0-9]+\\.result\\.csv$",
  full.names = TRUE
)

if (length(result_files) == 0L) {
  stop(
    "没有找到任何clumping结果文件：",
    CLUMP_RESULT_DIR
  )
}

clump_result_list <- lapply(
  result_files,
  function(f) {
    tryCatch(
      fread(f),
      error = function(e) {
        stop(
          "无法读取clumping结果文件：",
          f,
          "\n错误：",
          conditionMessage(e)
        )
      }
    )
  }
)

clump_result <- rbindlist(
  clump_result_list,
  use.names = TRUE,
  fill = TRUE
)

required_clump_cols <- c("rsid", "pval", "id")
missing_clump_cols <- setdiff(required_clump_cols, names(clump_result))

if (length(missing_clump_cols) > 0L) {
  stop(
    "clumping汇总结果缺少以下列：",
    paste(missing_clump_cols, collapse = ", ")
  )
}

clump_result[, `:=`(
  rsid = as.character(rsid),
  pval = as.numeric(pval),
  id = as.character(id)
)]

clump_result <- unique(
  clump_result[
    !is.na(rsid) &
      rsid != "" &
      !is.na(id) &
      id != ""
  ],
  by = c("id", "rsid")
)

setorder(clump_result, id, pval, rsid)

# 更新兼容旧版的汇总进度文件和completed_ids文件。
fwrite(clump_result, CLUMP_PROGRESS_FILE)

completed_ids_final <- sort(unique(clump_result$id))
writeLines(completed_ids_final, CLUMP_DONE_FILE)

cat("\n全部可用LD clumping结果已汇总。\n")
cat("LD clumping后行数：", nrow(clump_result), "\n")
cat("LD clumping后暴露数：", uniqueN(clump_result$id), "\n")
cat("LD clumping后唯一SNP数：", uniqueN(clump_result$rsid), "\n")

if (nrow(clump_result) == 0L) {
  stop("暴露LD clumping后没有保留任何工具变量。")
}

if (uniqueN(clump_result$id) < length(all_ids)) {
  warning(
    "当前仅完成 ",
    uniqueN(clump_result$id),
    "/",
    length(all_ids),
    " 个暴露。若任务曾中断，请设置RESET_CLUMP_WORK=FALSE后继续。"
  )
}

fwrite(
  clump_result,
  EXPOSURE_CLUMP_INDEX_FILE
)


# ============================================================
# 10. 将clumping结果合并回完整暴露信息
# ============================================================

clump_keep <- unique(
  clump_result[
    ,
    .(
      SNP = as.character(rsid),
      phenotype_id = as.character(id)
    )
  ]
)

b_clumped <- merge(
  b,
  clump_keep,
  by = c("SNP", "phenotype_id"),
  all = FALSE
)

setorder(
  b_clumped,
  phenotype_id,
  pval,
  SNP
)

cat("合并回完整暴露信息后的行数：", nrow(b_clumped), "\n")
cat("合并回完整暴露信息后的暴露数：", uniqueN(b_clumped$phenotype_id), "\n")
cat("合并回完整暴露信息后的唯一SNP数：", uniqueN(b_clumped$SNP), "\n")

if (nrow(b_clumped) == 0) {
  stop("clumping结果与原始暴露数据合并后为空。")
}

fwrite(
  b_clumped,
  EXPOSURE_CLUMPED_FILE
)

cat("\n输出完成：\n")
cat("P/F筛选数据：", EXPOSURE_FILTERED_FILE, "\n")
cat("LD clumping索引：", EXPOSURE_CLUMP_INDEX_FILE, "\n")
cat("LD clumping完整暴露：", EXPOSURE_CLUMPED_FILE, "\n")
cat("工作及断点目录：", CLUMP_WORK_DIR, "\n")