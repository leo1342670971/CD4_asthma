
rm(list = ls())

setwd("C:/Users/Szq/Desktop/MR")

library(data.table)

## -----------------------------
## 1. 读取基因
## -----------------------------
gene_dt <- fread("基因药物查询.csv", header = TRUE)
stopifnot("gene" %in% names(gene_dt))

genes <- unique(na.omit(trimws(as.character(gene_dt$gene))))
genes <- genes[genes != ""]

cat("输入基因数：", length(genes), "\n")

## -----------------------------
## 2. 依赖包
## -----------------------------
pkgs <- c("httr2", "jsonlite", "stringr")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

scalar_chr <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  x <- as.character(unlist(x, use.names = FALSE))
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  x[1]
}

collapse_vec <- function(x, sep = "; ") {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  x <- as.character(unlist(x, use.names = FALSE))
  x <- unique(x[!is.na(x) & x != ""])
  if (!length(x)) return(NA_character_)
  paste(x, collapse = sep)
}

split_ids <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(character())
  x <- as.character(x)
  x <- unlist(strsplit(x, ";\\s*"), use.names = FALSE)
  x <- unique(x[!is.na(x) & x != ""])
  x
}

## -----------------------------
## 3. OpenTargets GraphQL 基础函数
## -----------------------------
base_url <- "https://api.platform.opentargets.org/api/v4/graphql"

ot_post <- function(query, variables = list()) {
  
  resp <- httr2::request(base_url) |>
    httr2::req_method("POST") |>
    httr2::req_headers(
      `Content-Type` = "application/json",
      `Accept` = "application/json"
    ) |>
    httr2::req_user_agent("R OpenTargets gene drug indication MOA clinical stage query") |>
    httr2::req_body_json(
      list(query = query, variables = variables),
      auto_unbox = TRUE
    ) |>
    httr2::req_timeout(90) |>
    httr2::req_retry(max_tries = 3) |>
    httr2::req_perform()
  
  out <- httr2::resp_body_json(resp, simplifyVector = FALSE)
  
  if (!is.null(out$errors)) {
    msg <- paste(
      vapply(out$errors, function(e) as.character(e$message %||% ""), character(1)),
      collapse = "\n"
    )
    stop("OpenTargets GraphQL error:\n", msg)
  }
  
  out$data
}

## -----------------------------
## 4. gene symbol -> Ensembl ID
## -----------------------------
q_search_target <- '
query($q: String!){
  search(
    queryString: $q,
    entityNames: ["target"],
    page: {index: 0, size: 10}
  ){
    hits{
      id
      object{
        ... on Target {
          id
          approvedSymbol
          approvedName
        }
      }
    }
  }
}
'

q_target_basic <- '
query($ensemblId: String!){
  target(ensemblId: $ensemblId){
    id
    approvedSymbol
    approvedName
  }
}
'

resolve_to_target <- function(gene) {
  
  gene_raw <- trimws(as.character(gene))
  gene_no_version <- sub("\\..*$", "", gene_raw)
  
  ## 输入已经是 ENSG
  if (grepl("^ENSG\\d+$", gene_no_version, ignore.case = TRUE)) {
    
    ensg <- toupper(gene_no_version)
    
    d <- tryCatch(
      ot_post(q_target_basic, list(ensemblId = ensg)),
      error = function(e) NULL
    )
    
    tinfo <- d$target %||% NULL
    
    if (is.null(tinfo)) {
      return(list(
        input_gene = gene_raw,
        ensemblId = ensg,
        approvedSymbol = gene_raw,
        match_type = "input_is_ensembl_no_info"
      ))
    }
    
    return(list(
      input_gene = gene_raw,
      ensemblId = scalar_chr(tinfo$id),
      approvedSymbol = scalar_chr(tinfo$approvedSymbol),
      match_type = "input_is_ensembl"
    ))
  }
  
  ## 输入是 gene symbol
  d <- tryCatch(
    ot_post(q_search_target, list(q = gene_raw)),
    error = function(e) NULL
  )
  
  hits <- d$search$hits %||% list()
  
  if (!length(hits)) {
    return(list(
      input_gene = gene_raw,
      ensemblId = NA_character_,
      approvedSymbol = gene_raw,
      match_type = "no_hit"
    ))
  }
  
  get_sym <- function(h) {
    obj <- h$object %||% NULL
    if (is.null(obj)) return(NA_character_)
    scalar_chr(obj$approvedSymbol)
  }
  
  get_id <- function(h) {
    obj <- h$object %||% NULL
    if (!is.null(obj) && !is.null(obj$id)) return(scalar_chr(obj$id))
    scalar_chr(h$id)
  }
  
  syms <- vapply(hits, get_sym, character(1))
  idx <- which(toupper(syms) == toupper(gene_raw))
  
  pick <- if (length(idx)) hits[[idx[1]]] else hits[[1]]
  
  list(
    input_gene = gene_raw,
    ensemblId = get_id(pick),
    approvedSymbol = get_sym(pick),
    match_type = if (length(idx)) "symbol_exact" else "first_hit"
  )
}

## -----------------------------
## 5. target -> drugs
## 这里用于判断某个 gene 是否有对应 drug
## 只保留 drugId、drugName、target 临床报告 ID 和 target-level clinical stage
## -----------------------------
q_target_drugs <- '
query($ensemblId: String!){
  target(ensemblId: $ensemblId){
    id
    approvedSymbol
    drugAndClinicalCandidates{
      rows{
        id
        maxClinicalStage
        drug{
          id
          name
        }
        diseases{
          diseaseFromSource
          disease{
            id
            name
          }
        }
        clinicalReports{
          id
          clinicalStage
        }
      }
    }
  }
}
'

fetch_target_drugs_simple <- function(ensemblId) {
  
  d <- ot_post(q_target_drugs, list(ensemblId = ensemblId))
  tinfo <- d$target %||% NULL
  
  if (is.null(tinfo)) return(data.table())
  
  rows <- tinfo$drugAndClinicalCandidates$rows %||% list()
  if (!length(rows)) return(data.table())
  
  res <- rbindlist(lapply(rows, function(r) {
    
    drug <- r$drug %||% list()
    reports <- r$clinicalReports %||% list()
    
    report_ids <- collapse_vec(lapply(reports, function(x) x$id %||% NA_character_))
    report_stages <- collapse_vec(lapply(reports, function(x) x$clinicalStage %||% NA_character_))
    
    data.table(
      drugId = scalar_chr(drug$id),
      drug = scalar_chr(drug$name),
      target_max_clinical_stage = scalar_chr(r$maxClinicalStage),
      target_report_ids = report_ids,
      target_report_clinical_stages = report_stages
    )
    
  }), fill = TRUE)
  
  res <- res[!is.na(drugId) & drugId != ""]
  unique(res)
}

## -----------------------------
## 6. drug -> indication + MOA
## indication 没有时后面填 NA
## MOA 没有时后面填 unknown
## -----------------------------
q_drug_full <- '
query($chemblId: String!){
  drug(chemblId: $chemblId){
    id
    name
    mechanismsOfAction {
      rows {
        mechanismOfAction
        targetName
        targets {
          id
          approvedSymbol
        }
        references {
          source
          urls
        }
      }
      uniqueActionTypes
      uniqueTargetTypes
    }
    parentMolecule {
      id
      name
    }
    childMolecules {
      id
      name
    }
    indications {
      rows {
        id
        maxClinicalStage
        disease {
          id
          name
        }
        clinicalReports {
          id
        }
      }
    }
  }
}
'

fetch_drug_full_simple <- function(drugId) {
  
  d <- ot_post(q_drug_full, list(chemblId = drugId))
  drug <- d$drug %||% NULL
  
  if (is.null(drug)) {
    return(list(
      indications = data.table(
        drugId = drugId,
        drug = NA_character_,
        clinical_indication = NA_character_,
        indication_max_clinical_stage = NA_character_,
        indication_report_ids = NA_character_
      ),
      moa = data.table(
        drugId = drugId,
        mechanisms_of_action = "unknown",
        moa_target_name = NA_character_,
        moa_target_symbol = NA_character_
      )
    ))
  }
  
  ## ---- 6.1 indications ----
  ind_rows <- drug$indications$rows %||% list()
  
  if (!length(ind_rows)) {
    
    ind_dt <- data.table(
      drugId = scalar_chr(drug$id %||% drugId),
      drug = scalar_chr(drug$name),
      clinical_indication = NA_character_,
      indication_max_clinical_stage = NA_character_,
      indication_report_ids = NA_character_
    )
    
  } else {
    
    ind_dt <- rbindlist(lapply(ind_rows, function(r) {
      
      dis <- r$disease %||% list()
      reports <- r$clinicalReports %||% list()
      
      data.table(
        drugId = scalar_chr(drug$id %||% drugId),
        drug = scalar_chr(drug$name),
        clinical_indication = scalar_chr(dis$name),
        indication_max_clinical_stage = scalar_chr(r$maxClinicalStage),
        indication_report_ids = collapse_vec(lapply(reports, function(x) x$id %||% NA_character_))
      )
      
    }), fill = TRUE)
    
    ind_dt <- unique(ind_dt)
  }
  
  ## ---- 6.2 mechanisms of action ----
  moa_obj <- drug$mechanismsOfAction %||% list()
  moa_rows <- moa_obj$rows %||% list()
  
  if (!length(moa_rows)) {
    
    moa_dt <- data.table(
      drugId = scalar_chr(drug$id %||% drugId),
      mechanisms_of_action = "unknown",
      moa_target_name = NA_character_,
      moa_target_symbol = NA_character_
    )
    
  } else {
    
    mechanism_vec <- lapply(moa_rows, function(r) r$mechanismOfAction %||% NA_character_)
    target_name_vec <- lapply(moa_rows, function(r) r$targetName %||% NA_character_)
    
    target_symbol_vec <- lapply(moa_rows, function(r) {
      targets <- r$targets %||% list()
      if (!length(targets)) return(NA_character_)
      vapply(targets, function(tg) scalar_chr(tg$approvedSymbol), character(1))
    })
    
    mechanisms_of_action <- collapse_vec(mechanism_vec)
    moa_target_name <- collapse_vec(target_name_vec)
    moa_target_symbol <- collapse_vec(target_symbol_vec)
    
    if (is.na(mechanisms_of_action) || mechanisms_of_action == "") {
      mechanisms_of_action <- "unknown"
    }
    
    moa_dt <- data.table(
      drugId = scalar_chr(drug$id %||% drugId),
      mechanisms_of_action = mechanisms_of_action,
      moa_target_name = moa_target_name,
      moa_target_symbol = moa_target_symbol
    )
  }
  
  list(
    indications = ind_dt,
    moa = moa_dt
  )
}

## -----------------------------
## 7. clinicalReports(ids) -> clinical stage
## 按你给的 ClinicalRecordsQuery 查询
## -----------------------------
q_clinical_records <- '
query($clinicalReportsIds: [String!]!) {
  clinicalReports(clinicalReportsIds: $clinicalReportsIds) {
    clinicalStage
    id
    source
    title
    trialOverallStatus
    trialStartDate
    type
  }
}
'

fetch_clinical_reports_simple <- function(report_ids) {
  
  report_ids <- unique(report_ids[!is.na(report_ids) & report_ids != ""])
  
  if (!length(report_ids)) {
    return(data.table(
      clinical_report_id = character(),
      clinical_stage = character()
    ))
  }
  
  d <- ot_post(q_clinical_records, list(clinicalReportsIds = report_ids))
  
  rows <- d$clinicalReports %||% list()
  
  if (!length(rows)) {
    return(data.table(
      clinical_report_id = character(),
      clinical_stage = character()
    ))
  }
  
  res <- rbindlist(lapply(rows, function(r) {
    data.table(
      clinical_report_id = scalar_chr(r$id),
      clinical_stage = scalar_chr(r$clinicalStage)
    )
  }), fill = TRUE)
  
  res <- res[!is.na(clinical_report_id) & clinical_report_id != ""]
  unique(res)
}

## -----------------------------
## 8. clinical report 缓存
## 避免同一个 clinicalReportId 反复查询
## -----------------------------
clinical_cache <- new.env(parent = emptyenv())

get_clinical_reports_cached <- function(report_ids) {
  
  report_ids <- unique(report_ids[!is.na(report_ids) & report_ids != ""])
  
  if (!length(report_ids)) {
    return(data.table(
      clinical_report_id = character(),
      clinical_stage = character()
    ))
  }
  
  missing_ids <- report_ids[
    !vapply(report_ids, exists, logical(1), envir = clinical_cache, inherits = FALSE)
  ]
  
  if (length(missing_ids)) {
    
    ## 分块查询，避免一次传太多 ID
    chunks <- split(missing_ids, ceiling(seq_along(missing_ids) / 100))
    
    for (chunk in chunks) {
      
      dt <- tryCatch(
        fetch_clinical_reports_simple(chunk),
        error = function(e) {
          data.table(
            clinical_report_id = character(),
            clinical_stage = character()
          )
        }
      )
      
      for (rid in chunk) {
        
        one <- dt[clinical_report_id == rid]
        
        if (!nrow(one)) {
          one <- data.table(
            clinical_report_id = rid,
            clinical_stage = NA_character_
          )
        }
        
        assign(rid, one, envir = clinical_cache)
      }
    }
  }
  
  out <- rbindlist(
    lapply(report_ids, function(rid) get(rid, envir = clinical_cache)),
    fill = TRUE
  )
  
  unique(out)
}

get_stage_from_reports <- function(report_ids) {
  
  report_ids <- unique(report_ids[!is.na(report_ids) & report_ids != ""])
  
  if (!length(report_ids)) return(NA_character_)
  
  dt <- get_clinical_reports_cached(report_ids)
  
  stages <- unique(dt$clinical_stage)
  stages <- stages[!is.na(stages) & stages != ""]
  
  if (!length(stages)) return(NA_character_)
  
  paste(stages, collapse = "; ")
}

## -----------------------------
## 9. 批量运行
## 不写任何大中间文件
## -----------------------------
result_list <- list()
failed_list <- list()

drug_cache <- new.env(parent = emptyenv())

for (i in seq_along(genes)) {
  
  g <- genes[i]
  cat(sprintf("[%d/%d] %s\n", i, length(genes), g))
  
  ## ---- 9.1 gene -> target ----
  mapped <- tryCatch(
    resolve_to_target(g),
    error = function(e) {
      failed_list[[length(failed_list) + 1]] <<- list(
        gene = g,
        step = "resolve_to_target",
        error = conditionMessage(e)
      )
      list(
        input_gene = g,
        ensemblId = NA_character_,
        approvedSymbol = g,
        match_type = "resolve_error"
      )
    }
  )
  
  ## 库里没有这个基因：不输出
  if (is.na(mapped$ensemblId) || mapped$ensemblId == "") {
    failed_list[[length(failed_list) + 1]] <- list(
      gene = g,
      step = "resolve_to_target",
      error = "No Ensembl target ID found"
    )
    next
  }
  
  Sys.sleep(0.1)
  
  ## ---- 9.2 target -> drugs ----
  target_drugs <- tryCatch(
    fetch_target_drugs_simple(mapped$ensemblId),
    error = function(e) {
      failed_list[[length(failed_list) + 1]] <<- list(
        gene = g,
        ensemblId = mapped$ensemblId,
        step = "fetch_target_drugs_simple",
        error = conditionMessage(e)
      )
      data.table()
    }
  )
  
  ## 这个基因没有对应药物：不输出
  if (!nrow(target_drugs)) next
  
  drug_ids <- unique(target_drugs$drugId)
  drug_ids <- drug_ids[!is.na(drug_ids) & drug_ids != ""]
  
  ## ---- 9.3 对每个 drugId 查 drug-level indication + MOA + clinical stage ----
  for (drug_id in drug_ids) {
    
    ## drug 缓存
    if (exists(drug_id, envir = drug_cache, inherits = FALSE)) {
      
      drug_full <- get(drug_id, envir = drug_cache)
      
    } else {
      
      Sys.sleep(0.1)
      
      drug_full <- tryCatch(
        fetch_drug_full_simple(drug_id),
        error = function(e) {
          failed_list[[length(failed_list) + 1]] <<- list(
            gene = g,
            drugId = drug_id,
            step = "fetch_drug_full_simple",
            error = conditionMessage(e)
          )
          list(
            indications = data.table(
              drugId = drug_id,
              drug = NA_character_,
              clinical_indication = NA_character_,
              indication_max_clinical_stage = NA_character_,
              indication_report_ids = NA_character_
            ),
            moa = data.table(
              drugId = drug_id,
              mechanisms_of_action = "unknown",
              moa_target_name = NA_character_,
              moa_target_symbol = NA_character_
            )
          )
        }
      )
      
      assign(drug_id, drug_full, envir = drug_cache)
    }
    
    ind_dt <- copy(drug_full$indications)
    moa_dt <- copy(drug_full$moa)
    
    ## 当前 gene 对应的 drug 名
    this_drug_name <- target_drugs[drugId == drug_id, unique(drug)]
    this_drug_name <- this_drug_name[!is.na(this_drug_name) & this_drug_name != ""]
    this_drug_name <- if (length(this_drug_name)) this_drug_name[1] else NA_character_
    
    ## target-level clinical stage/report ids，用于 indication 没有 report 时兜底
    this_target_stage <- target_drugs[drugId == drug_id, unique(target_max_clinical_stage)]
    this_target_stage <- this_target_stage[!is.na(this_target_stage) & this_target_stage != ""]
    this_target_stage <- if (length(this_target_stage)) paste(this_target_stage, collapse = "; ") else NA_character_
    
    this_target_report_ids <- target_drugs[drugId == drug_id, unique(target_report_ids)]
    this_target_report_ids <- split_ids(this_target_report_ids)
    
    this_target_report_stage <- get_stage_from_reports(this_target_report_ids)
    
    if (!nrow(ind_dt)) {
      ind_dt <- data.table(
        drugId = drug_id,
        drug = this_drug_name,
        clinical_indication = NA_character_,
        indication_max_clinical_stage = NA_character_,
        indication_report_ids = NA_character_
      )
    }
    
    ## 如果 drug name 缺失，用 target_drugs 的 drug name 补
    ind_dt[is.na(drug) | drug == "", drug := this_drug_name]
    
    ## 如果 MOA 缺失，按要求填 unknown
    if (!nrow(moa_dt)) {
      moa_dt <- data.table(
        drugId = drug_id,
        mechanisms_of_action = "unknown",
        moa_target_name = NA_character_,
        moa_target_symbol = NA_character_
      )
    }
    
    moa_dt[is.na(mechanisms_of_action) | mechanisms_of_action == "",
           mechanisms_of_action := "unknown"]
    
    ## 每个 drugId 的 MOA 已经 collapse 成一行
    moa_one <- moa_dt[1]
    
    ## ---- 9.4 生成最终小表 ----
    for (j in seq_len(nrow(ind_dt))) {
      
      report_ids_j <- split_ids(ind_dt$indication_report_ids[j])
      
      clinical_stage_j <- get_stage_from_reports(report_ids_j)
      
      ## clinicalReports 没查到时，依次用 indication max stage、target report stage、target max stage 兜底
      if (is.na(clinical_stage_j) || clinical_stage_j == "") {
        clinical_stage_j <- scalar_chr(ind_dt$indication_max_clinical_stage[j])
      }
      
      if (is.na(clinical_stage_j) || clinical_stage_j == "") {
        clinical_stage_j <- this_target_report_stage
      }
      
      if (is.na(clinical_stage_j) || clinical_stage_j == "") {
        clinical_stage_j <- this_target_stage
      }
      
      if (is.na(clinical_stage_j) || clinical_stage_j == "") {
        clinical_stage_j <- NA_character_
      }
      
      tmp <- data.table(
        gene = mapped$approvedSymbol %||% g,
        drugId = drug_id,
        drug = scalar_chr(ind_dt$drug[j]),
        clinical_indication = scalar_chr(ind_dt$clinical_indication[j]),
        mechanisms_of_action = scalar_chr(moa_one$mechanisms_of_action),
        moa_target_name = scalar_chr(moa_one$moa_target_name),
        moa_target_symbol = scalar_chr(moa_one$moa_target_symbol),
        clinical_stage = clinical_stage_j
      )
      
      result_list[[length(result_list) + 1]] <- tmp
    }
  }
}

## -----------------------------
## 10. 最终精简输出
## -----------------------------
if (length(result_list)) {
  
  final_dt <- rbindlist(result_list, fill = TRUE)
  
  ## 基因和药物必须存在；indication/MOA/clinical stage 可以为空或 unknown
  final_dt <- final_dt[
    !is.na(gene) & gene != "" &
      !is.na(drugId) & drugId != "" &
      !is.na(drug) & drug != ""
  ]
  
  ## 按你的要求：
  ## indication 没有则 NA
  ## MOA 没有则 unknown
  ## clinical stage 没有则 NA
  final_dt[clinical_indication == "", clinical_indication := NA_character_]
  final_dt[is.na(mechanisms_of_action) | mechanisms_of_action == "",
           mechanisms_of_action := "unknown"]
  final_dt[clinical_stage == "", clinical_stage := NA_character_]
  
  final_dt <- unique(
    final_dt,
    by = c(
      "gene",
      "drugId",
      "drug",
      "clinical_indication",
      "mechanisms_of_action",
      "moa_target_name",
      "moa_target_symbol",
      "clinical_stage"
    )
  )
  
  setorder(final_dt, gene, drug, clinical_indication)
  
} else {
  
  final_dt <- data.table(
    gene = character(),
    drugId = character(),
    drug = character(),
    clinical_indication = character(),
    mechanisms_of_action = character(),
    moa_target_name = character(),
    moa_target_symbol = character(),
    clinical_stage = character()
  )
}

fwrite(
  final_dt,
  "OpenTargets_gene_drug_indication_MOA_clinicalStage_simple.csv"
)

## -----------------------------
## 11. 保存很小的失败日志
## -----------------------------
if (length(failed_list)) {
  saveRDS(failed_list, "OpenTargets_failed_simple.rds")
  
  failed_dt <- rbindlist(failed_list, fill = TRUE)
  fwrite(failed_dt, "OpenTargets_failed_simple.csv")
}

cat("\n完成。\n")
cat("最终结果行数：", nrow(final_dt), "\n")
cat("输出文件：OpenTargets_gene_drug_indication_MOA_clinicalStage_simple.csv\n")

if (length(failed_list)) {
  cat("失败日志：OpenTargets_failed_simple.rds\n")
  cat("失败日志CSV：OpenTargets_failed_simple.csv\n")
}

print(head(final_dt, 30))
