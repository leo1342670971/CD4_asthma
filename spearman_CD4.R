
# 1. Settings -------------------------------------------------
work_dir <- if (.Platform$OS.type == "windows") {
 
}
background_name <- "01_all_eGene_cell_type_gene_deduplicated.csv"
significant_name <- "02_significant_eGene_cell_type_gene_deduplicated.csv"
figure_prefix <- "Figure5"  # Change to "Figure4" if needed.
run_gene_cluster_sensitivity <- TRUE

# The 01 file defines the background; the 02 file defines the significant set.
# MR filtering is not repeated here. Files 111.csv and FDR.csv are not inputs.
packages <- c("data.table", "stringr", "ggplot2", "ggrepel", "scales")
if (run_gene_cluster_sensitivity) packages <- c(packages, "geepack")
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install the missing packages first: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})
candidate_dirs <- c(work_dir, file.path(work_dir, "细胞类型富集"))
valid_dirs <- candidate_dirs[
  file.exists(file.path(candidate_dirs, background_name)) &
  file.exists(file.path(candidate_dirs, significant_name))
]
if (!length(valid_dirs)) {
  stop("Place BOTH 01 and 02 in work_dir or its 细胞类型富集 subfolder.")
}
data_dir <- valid_dirs[1]
background_file <- file.path(data_dir, background_name)
significant_file <- file.path(data_dir, significant_name)
out_dir <- tempfile(
  paste0("analysis_from_01_02_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
  tmpdir = data_dir
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Cannot create the output directory.")
save_csv <- function(x, filename) fwrite(x, file.path(out_dir, filename))

# 2. Read and identify gene-state records -----------------------
# Use Ensembl ID + FULL cell_type as the key.
# A missing gene symbol does not remove a record with a valid Ensembl ID.
clean_text <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}
read_gene_states <- function(path, exposure_column, fallback_symbol, role) {
  raw <- fread(path, encoding = "UTF-8", check.names = FALSE)
  required <- c(exposure_column, "cell_type")
  if (!all(required %in% names(raw))) {
    stop(role, " file is missing: ", paste(setdiff(required, names(raw)), collapse = ", "))
  }
  exposure <- clean_text(raw[[exposure_column]])
  id_matches <- str_extract_all(exposure, "ENSG[0-9]+")
  gene_id <- vapply(id_matches, function(ids) {
    ids <- unique(ids[!is.na(ids)])
    if (length(ids) == 1L) ids else NA_character_
  }, character(1))
  state <- clean_text(raw$cell_type)
  gene_symbol <- if ("gene" %in% names(raw)) {
    clean_text(raw$gene)
  } else if (fallback_symbol %in% names(raw)) {
    clean_text(raw[[fallback_symbol]])
  } else {
    rep(NA_character_, nrow(raw))
  }
  if (fallback_symbol %in% names(raw)) {
    fallback <- clean_text(raw[[fallback_symbol]])
    gene_symbol[is.na(gene_symbol)] <- fallback[is.na(gene_symbol)]
  }
  parsed <- data.table(
    source_csv_row = seq_len(nrow(raw)) + 1L, exposure_raw = exposure,
    gene_id = gene_id, cell_type = state, gene_symbol = gene_symbol
  )
  invalid <- parsed[is.na(gene_id) | is.na(cell_type)]
  save_csv(invalid, paste0("audit_", role, "_invalid_keys.csv"))
  if (nrow(invalid)) stop("Missing/ambiguous Ensembl ID or state; inspect the saved audit.")

  # Check that the state column agrees with the exposure when parseable.
  exposure_state <- str_match(
    exposure, "(?:^|:)([^:]+)_[0-9]+kb_combined[^:]*(?=:|$)"
  )[, 2]
  mismatch <- which(!is.na(exposure_state) & exposure_state != state)
  save_csv(parsed[mismatch], paste0("audit_", role, "_state_mismatches.csv"))
  if (length(mismatch)) stop("cell_type disagrees with exposure in the ", role, " file.")

  duplicate_keys <- parsed[, .N, by = .(gene_id, cell_type)][N > 1L]
  save_csv(duplicate_keys, paste0("audit_", role, "_duplicate_keys.csv"))
  save_csv(parsed[is.na(gene_symbol)], paste0("audit_", role, "_missing_symbols.csv"))
  pairs <- parsed[, .(
    gene_symbol = {
      symbols <- sort(unique(gene_symbol[!is.na(gene_symbol)]))
      if (length(symbols)) paste(symbols, collapse = ";") else NA_character_
    }
  ), by = .(gene_id, cell_type)]
  pairs[, gene_label := fifelse(is.na(gene_symbol), gene_id, gene_symbol)]
  list(pairs = pairs, raw_n = nrow(raw))
}
bg_input <- read_gene_states(background_file, "Exposure_all", "Gene_all", "background")
sig_input <- read_gene_states(significant_file, "Exposure", "Gene", "significant")
all_dt <- copy(bg_input$pairs)
significant_pairs <- sig_input$pairs
outside <- significant_pairs[!all_dt, on = .(gene_id, cell_type)]
save_csv(outside, "audit_significant_outside_background.csv")
if (nrow(outside)) {
  stop("Significant gene-state records are absent from 01. Resolve this before analysis.")
}
all_dt[, MR_sig := 0L]
all_dt[significant_pairs, on = .(gene_id, cell_type), MR_sig := 1L]
sig_dt <- all_dt[MR_sig == 1L]
if (uniqueN(all_dt$cell_type) < 2L) stop("At least two states are required.")
save_csv(all_dt, "01_background_gene_states_used.csv")
save_csv(sig_dt, "02_significant_gene_states_used.csv")

# 3. Rebuild state counts from the supplied gene-level inputs ---
# The denominator comes only from 01; the numerator comes from 02.
plot_df <- all_dt[, .(N_eGenes = .N, N_sig = sum(MR_sig)), by = cell_type]
plot_df[, prop_sig := N_sig / N_eGenes]
plot_df[, time := str_extract(cell_type, "(0h|16h|40h|5d|LA)$")]
plot_df[is.na(time), time := "Other"]
time_levels <- c("0h", "16h", "40h", "5d", "LA", "Other")
plot_df[, time := factor(time, levels = time_levels)]
# Simplify DISPLAY labels only, preserving full state identities.
plot_df[, display_state := str_replace_all(cell_type, "_(stim|unstim|uns)_", "_")]
ambiguous_labels <- plot_df[, .N, by = display_state][N > 1L, display_state]
plot_df[display_state %in% ambiguous_labels, display_state := cell_type]
plot_df[, cell_label := str_remove(display_state, "_(0h|16h|40h|5d|LA)$")]
save_csv(plot_df, "05_plot_data_eGenes_vs_significant_signals.csv")

# 4. Report unique genes separately from gene-state records -----
gene_reuse <- all_dt[, .(
  N_states_background = .N, N_states_significant = sum(MR_sig)
), by = gene_id]
setorder(gene_reuse, -N_states_significant, -N_states_background, gene_id)
save_csv(gene_reuse, "04_gene_reuse_across_states.csv")
audit_summary <- data.table(
  background_input_rows = bg_input$raw_n,
  significant_input_rows = sig_input$raw_n,
  background_gene_state_records = nrow(all_dt),
  significant_gene_state_records = nrow(sig_dt),
  unique_background_Ensembl_IDs = uniqueN(all_dt$gene_id),
  unique_significant_Ensembl_IDs = uniqueN(sig_dt$gene_id),
  background_duplicates_removed = bg_input$raw_n - nrow(all_dt),
  significant_duplicates_removed = sig_input$raw_n - nrow(sig_dt),
  significant_genes_in_multiple_states = sum(gene_reuse$N_states_significant > 1L),
  N_states = nrow(plot_df)
)
save_csv(audit_summary, "00_input_and_analysis_counts.csv")

# 5. Descriptive Spearman analysis and Figure B -----------------
fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                formatC(p, format = "f", digits = 3)))
}
time_colors <- c(
  "0h" = "#66C2A5", "16h" = "#FC8D62", "40h" = "#8DA0CB",
  "5d" = "#E78AC3", "LA" = "#A6D854", "Other" = "grey60"
)
sp <- cor.test(plot_df$N_sig, plot_df$N_eGenes,
               method = "spearman", exact = FALSE)
spearman_result <- data.table(
  method = "Spearman", N_states = nrow(plot_df),
  rho = as.numeric(sp$estimate), P_nominal = sp$p.value,
  interpretation = "Descriptive; nominal P does not account for dependence between states."
)
save_csv(spearman_result, "06_spearman_summary.csv")
# There is one global correlation test, not one per state.
# Report rho, not rho squared. N_sig is a subset of N_eGenes.
pB <- ggplot(plot_df, aes(N_sig, N_eGenes, color = time, label = cell_label)) +
  geom_point(size = 3.2, alpha = 0.9) +
  geom_text_repel(
    size = 4, max.overlaps = Inf, box.padding = 0.35,
    point.padding = 0.25, segment.color = NA, seed = 20261017,
    show.legend = FALSE
  ) +
  scale_color_manual(values = time_colors, drop = FALSE) +
  labs(
    x = "Number of MR-significant genes",
    y = "Number of eGenes",
    color = "Time\npoint",
    title = paste0("Spearman rho = ", sprintf("%.2f", as.numeric(sp$estimate))),
    subtitle = paste0("Nominal P = ", fmt_p(sp$p.value), " (exploratory)"),
    caption = paste(
      "Related states may share genes; the nominal P value does not account for this dependence.",
      "MR-significant genes are a subset of the background genes in each state.", sep = "\n"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(size = 9, hjust = 0))
ggsave(file.path(out_dir, paste0(figure_prefix, "B_scatter.pdf")),
       pB, width = 9, height = 8)
ggsave(file.path(out_dir, paste0(figure_prefix, "B_scatter.png")),
       pB, width = 9, height = 8, dpi = 300)

# 6. State-versus-other-states Fisher tests and pointwise CIs ---
# Original direction retained: alternative = "greater".
# The comparison pools gene-state records from all OTHER states.
# These records may share genes; see the separate GEE sensitivity.
z95 <- qnorm(0.975)
enrich_df <- copy(plot_df)
total_all <- sum(enrich_df$N_eGenes)
total_sig <- sum(enrich_df$N_sig)
enrich_df[, a := N_sig]
enrich_df[, b := N_eGenes - N_sig]
enrich_df[, c_other := total_sig - N_sig]
enrich_df[, d_other := total_all - total_sig - b]
# Wilson intervals for proportions: descriptive, not simultaneous.
enrich_df[, wilson_center :=
  (prop_sig + z95^2 / (2 * N_eGenes)) / (1 + z95^2 / N_eGenes)]
enrich_df[, wilson_halfwidth :=
  z95 * sqrt(prop_sig * (1 - prop_sig) / N_eGenes +
               z95^2 / (4 * N_eGenes^2)) / (1 + z95^2 / N_eGenes)]
enrich_df[, prop_CI_low := pmax(0, wilson_center - wilson_halfwidth)]
enrich_df[, prop_CI_high := pmin(1, wilson_center + wilson_halfwidth)]
enrich_df[, c("wilson_center", "wilson_halfwidth") := NULL]
fisher_results <- rbindlist(lapply(seq_len(nrow(enrich_df)), function(i) {
  row <- enrich_df[i]
  tab <- matrix(c(row$a, row$b, row$c_other, row$d_other), nrow = 2, byrow = TRUE)
  if (any(tab < 0)) stop("Negative 2x2 cell count.")
  if (any(colSums(tab) == 0)) {
    return(data.table(OR = NA_real_, OR_CI_low = NA_real_,
                      OR_CI_high = NA_real_, P = NA_real_))
  }
  ft_one <- fisher.test(tab, alternative = "greater")
  # Calculate a TWO-SIDED 95% OR interval separately from the one-sided P.
  ft_two <- fisher.test(tab, alternative = "two.sided", conf.level = 0.95)
  data.table(
    OR = unname(ft_one$estimate), OR_CI_low = ft_two$conf.int[1],
    OR_CI_high = ft_two$conf.int[2], P = ft_one$p.value
  )
}))
enrich_df <- cbind(enrich_df, fisher_results)
enrich_df[, FDR := p.adjust(P, method = "BH", n = nrow(enrich_df))]
enrich_df[, FDR_BY := p.adjust(P, method = "BY", n = nrow(enrich_df))]
enrich_df[, nominal_P_lt_0_05 := !is.na(P) & P < 0.05]
enrich_df[, BH_lt_0_05 := !is.na(FDR) & FDR < 0.05]
enrich_df[, sparse_2x2 := pmin(a, b, c_other, d_other) < 5]
enrich_df[, evidence := fcase(
  is.na(P), "Not estimable",
  BH_lt_0_05, "BH-adjusted P < 0.05",
  nominal_P_lt_0_05, "Nominal P < 0.05 only",
  default = "Nominal P >= 0.05"
)]
save_csv(enrich_df, "07_Fisher_counts_OR_CI_BH.csv")
save_csv(enrich_df[nominal_P_lt_0_05][order(P)], "08_nominal_enrichment.csv")
save_csv(enrich_df[BH_lt_0_05][order(FDR)], "09_BH_supported_enrichment.csv")

# 7. Figure C: original bar-chart layout, adding counts and CIs --
# The custom display order follows the supplied original script.
display_order <- c(
  "CD4_Naive_0h", "CD4_Naive_16h", "CD4_Naive_40h", "CD4_Naive_5d",
  "TN_0h", "TN_16h", "TN_40h", "TN_5d", "TN_LA", "TN_NFKB",
  "TN_IFN_16h", "TN_IFN_40h", "TN_IFN_5d", "TN_IFN_LA",
  "TN_cycling_40h", "TN_cycling_5d", "TN_HSP_5d", "HSP_16h",
  "CD4_Memory_0h", "CD4_Memory_16h", "CD4_Memory_40h", "CD4_Memory_5d",
  "TCM_0h", "TCM_16h", "TCM_40h", "TCM_5d", "TCM_LA",
  "TEM_0h", "TEM_16h", "TEM_40h", "TEM_5d", "TEM_LA",
  "TEMRA_0h", "TEMRA_16h", "TEMRA_40h", "TEMRA_5d", "TEMRA_LA",
  "TEM_HLApositive_40h", "TEM_HLApositive_5d", "TM_cycling_5d",
  "TM_ER-stress_40h", "T_ER-stress_5d", "TN2_40h",
  "nTreg_0h", "nTreg_16h", "nTreg_40h"
)
enrich_df[, plot_order := match(display_state, display_order)]
enrich_df[is.na(plot_order), plot_order := length(display_order) + 1L]
setorder(enrich_df, plot_order, cell_type)
enrich_df[, state_plot := factor(cell_type, levels = rev(cell_type))]
state_labels <- setNames(enrich_df$display_state, enrich_df$cell_type)
enrich_df[, annotation := paste0(
  N_sig, "/", N_eGenes, "   P=", fmt_p(P), "   BH=", fmt_p(FDR)
)]
pC <- ggplot(enrich_df, aes(state_plot, prop_sig, fill = evidence)) +
  geom_col(width = 0.70) +
  geom_errorbar(aes(ymin = prop_CI_low, ymax = prop_CI_high),
                width = 0.25, linewidth = 0.4, color = "black") +
  geom_text(aes(y = 1.03, label = annotation),
            hjust = 0, size = 2.7, color = "black") +
  scale_x_discrete(labels = state_labels) +
  scale_y_continuous(
    position = "right", limits = c(0, 1.70), breaks = seq(0, 1, 0.2),
    labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.01))
  ) +
  scale_fill_manual(
    values = c(
      "BH-adjusted P < 0.05" = "#7B3294",
      "Nominal P < 0.05 only" = "#E6C229",
      "Nominal P >= 0.05" = "grey75",
      "Not estimable" = "white"
    ),
    name = NULL
  ) +
  coord_flip(clip = "off") +
  labs(
    x = NULL, y = "Proportion of MR-significant genes among background eGenes",
    title = "Exploratory comparison across CD4+ T-cell states",
    subtitle = "Bars: proportions with pointwise 95% Wilson intervals | Labels: n/N, nominal P, BH-adjusted P",
    caption = paste0(
      "One-sided Fisher tests; each state versus all remaining states; ",
      nrow(enrich_df), " comparisons in the BH family.\n",
      "Fisher tests and Wilson intervals do not model gene reuse across states. ",
      "The gene-clustered sensitivity analysis is reported separately."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 10),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 9, hjust = 0),
    plot.margin = margin(8, 15, 8, 8)
  )
height_C <- max(8, 0.24 * nrow(enrich_df) + 2)
ggsave(file.path(out_dir, paste0(figure_prefix, "C_barplot_CI_BH.pdf")),
       pC, width = 14, height = height_C)
ggsave(file.path(out_dir, paste0(figure_prefix, "C_barplot_CI_BH.png")),
       pC, width = 14, height = height_C, dpi = 300)
save_csv(enrich_df, "10_FigureC_plot_data.csv")

# 8. Additional sensitivity accounting for within-gene reuse ----
# Binomial GEE: MR significance ~ current state versus other states.
# Cluster = Ensembl gene. Independence working correlation plus
# robust sandwich SE allows correlated records within each gene.
# The small-state estimates use asymptotic inference and can be unstable.
# This does not model LD between genes or all shared-donor errors.
sensitivity_status <- "Not requested"
if (run_gene_cluster_sensitivity) {
  gee_data <- all_dt[, .(gene_id, cell_type, MR_sig)]
  setorder(gee_data, gene_id, cell_type)  # Clusters must occupy contiguous rows.
  gee_data[, cluster_id := match(gene_id, unique(gene_id))]
  save_csv(gee_data, "11_gene_clustered_analysis_input.csv")

  gee_results <- rbindlist(lapply(plot_df$cell_type, function(current_state) {
    model_data <- as.data.frame(gee_data)
    model_data$in_state <- as.integer(model_data$cell_type == current_state)
    row <- enrich_df[cell_type == current_state]
    result <- data.table(
      cell_type = current_state, N_eGenes = row$N_eGenes, N_sig = row$N_sig,
      N_gene_clusters = uniqueN(gee_data$gene_id),
      OR_GEE = NA_real_, OR_GEE_CI_low = NA_real_, OR_GEE_CI_high = NA_real_,
      robust_SE_log_OR = NA_real_, P_GEE = NA_real_,
      sparse_2x2 = row$sparse_2x2, fit_status = "not fitted", fit_warnings = ""
    )
    if (any(unlist(row[, .(a, b, c_other, d_other)]) == 0)) {
      result[, fit_status := "not fitted: zero cell; separation possible"]
      return(result)
    }
    fit_warning_messages <- character()
    fit <- tryCatch(
      withCallingHandlers(
        geepack::geeglm(
          MR_sig ~ in_state, data = model_data, id = cluster_id,
          family = binomial(link = "logit"),
          corstr = "independence", std.err = "san.se",
          scale.fix = TRUE, scale.value = 1
        ),
        warning = function(w) {
          fit_warning_messages <<- c(fit_warning_messages, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    result[, fit_warnings := paste(unique(fit_warning_messages), collapse = "; ")]
    if (inherits(fit, "error")) {
      result[, fit_status := paste("fit failed:", conditionMessage(fit))]
      return(result)
    }
    not_converged <-
      (!is.null(fit$geese$error) && any(fit$geese$error != 0L)) ||
      (!is.null(fit$converged) && !isTRUE(fit$converged))
    if (not_converged) {
      result[, fit_status := "not converged"]
      return(result)
    }
    coefficients <- summary(fit)$coefficients
    beta <- coefficients["in_state", "Estimate"]
    robust_se <- coefficients["in_state", "Std.err"]
    if (!is.finite(beta) || !is.finite(robust_se) || robust_se <= 0) {
      result[, fit_status := "invalid coefficient or robust SE"]
      return(result)
    }
    result[, OR_GEE := exp(beta)]
    result[, OR_GEE_CI_low := exp(beta - z95 * robust_se)]
    result[, OR_GEE_CI_high := exp(beta + z95 * robust_se)]
    result[, robust_SE_log_OR := robust_se]
    # One-sided enrichment P, matching the original Fisher direction.
    result[, P_GEE := pnorm(beta / robust_se, lower.tail = FALSE)]
    result[, fit_status := "fitted"]
    result
  }))
  # Keep the full planned family size if some models cannot be estimated.
  gee_results[, FDR_BH_GEE := p.adjust(P_GEE, "BH", n = nrow(plot_df))]
  gee_results[, FDR_BY_GEE := p.adjust(P_GEE, "BY", n = nrow(plot_df))]
  save_csv(gee_results, "12_gene_clustered_GEE_sensitivity.csv")
  combined <- merge(
    enrich_df, gee_results,
    by = c("cell_type", "N_eGenes", "N_sig", "sparse_2x2"),
    all.x = TRUE, sort = FALSE
  )
  save_csv(combined, "13_Fisher_and_GEE_results.csv")
  sensitivity_status <- paste0(
    "GEE fitted for ", sum(gee_results$fit_status == "fitted"),
    "/", nrow(plot_df), " states."
  )
  if (any(gee_results$fit_status != "fitted")) {
    warning("Some GEE fits were unavailable; see fit_status. Do not count them as negative.")
  }
}


# 9. Summary and provenance -----------------------------------
provenance <- data.table(
  source_role = c("Background gene-state list", "Significant gene-state list"),
  path = c(normalizePath(background_file), normalizePath(significant_file)),
  md5 = unname(tools::md5sum(c(background_file, significant_file)))
)
save_csv(provenance, "source_files_and_MD5.csv")

writeLines(summary_lines, file.path(out_dir, "README_run_summary.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "R_sessionInfo.txt"))
print(audit_summary)
print(spearman_result)
print(enrich_df[order(P), .(
  cell_type, N_eGenes, N_sig, prop_sig, prop_CI_low, prop_CI_high,
  OR, OR_CI_low, OR_CI_high, P, FDR
)])
cat(paste(summary_lines, collapse = "\n"), "\n")
