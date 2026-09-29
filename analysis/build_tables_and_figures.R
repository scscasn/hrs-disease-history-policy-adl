#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
q3_dir <- file.path(run_dir, "q3")
table_dir <- file.path(q3_dir, "tables")
figure_dir <- file.path(q3_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_pair <- function(name) {
  read.csv(gzfile(file.path(run_dir, "derived", name)), stringsAsFactors = FALSE)
}

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_sd <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (sum(keep) < 2) return(NA_real_)
  m <- weighted_mean(x[keep], w[keep])
  sqrt(sum(w[keep] * (x[keep] - m)^2) / sum(w[keep]))
}

group_masks <- function(d) {
  target <- d$eligible_base == 1 & d$exposure_identified == 1 & d$follow_alive == 1
  list(
    "Survivor target" = target,
    "Outcome observed" = target & d$primary_outcome_observed == 1,
    "Alive outcome unknown" = target & d$primary_outcome_observed == 0
  )
}

continuous_row <- function(d, mask, variable, label, unit) {
  x <- d[[variable]][mask]
  w <- d$base_weight[mask]
  data.frame(
    characteristic = label,
    level = "Mean (weighted SD)",
    summary_type = "continuous",
    unweighted_n = sum(is.finite(x)),
    missing_n = sum(!is.finite(x)),
    estimate = weighted_mean(x, w),
    dispersion = weighted_sd(x, w),
    unit = unit,
    stringsAsFactors = FALSE
  )
}

categorical_rows <- function(d, mask, values, levels, label) {
  w <- d$base_weight[mask]
  x <- values[mask]
  denom <- sum(w[!is.na(x) & is.finite(w) & w > 0])
  do.call(rbind, lapply(levels, function(level_name) {
    hit <- !is.na(x) & x == level_name
    data.frame(
      characteristic = label,
      level = level_name,
      summary_type = "categorical",
      unweighted_n = sum(hit),
      missing_n = sum(is.na(x)),
      estimate = if (denom > 0) sum(w[hit & is.finite(w) & w > 0]) / denom else NA_real_,
      dispersion = NA_real_,
      unit = "proportion",
      stringsAsFactors = FALSE
    )
  }))
}

d <- read_pair("q2_hrs_w14_w15.csv.gz")
masks <- group_masks(d)
sex <- ifelse(d$sex_code == 1, "Male", ifelse(d$sex_code == 2, "Female", NA_character_))
race <- rep(NA_character_, nrow(d))
race[d$hispanic_code == 1 & !is.na(d$hispanic_code)] <- "Hispanic"
race[d$hispanic_code == 0 & d$race_code == 1 & !is.na(d$hispanic_code) & !is.na(d$race_code)] <- "Non-Hispanic White"
race[d$hispanic_code == 0 & d$race_code == 2 & !is.na(d$hispanic_code) & !is.na(d$race_code)] <- "Non-Hispanic Black"
race[d$hispanic_code == 0 & d$race_code == 3 & !is.na(d$hispanic_code) & !is.na(d$race_code)] <- "Other race/ethnicity"
proxy <- ifelse(d$base_proxy == 0, "Self interview", ifelse(d$base_proxy == 1, "Proxy interview", NA_character_))
changed <- ifelse(is.finite(d$final_count) & is.finite(d$cutoff_count), ifelse(d$final_count != d$cutoff_count, "Policy count changed", "No count change"), NA_character_)

table1 <- do.call(rbind, lapply(names(masks), function(group_name) {
  mask <- masks[[group_name]]
  pieces <- list(
    continuous_row(d, mask, "age", "Age", "years"),
    categorical_rows(d, mask, sex, c("Male", "Female"), "Sex"),
    categorical_rows(d, mask, race, c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other race/ethnicity"), "Race/ethnicity"),
    continuous_row(d, mask, "education_years", "Education", "years"),
    continuous_row(d, mask, "bmi", "Body mass index", "kg/m^2"),
    categorical_rows(d, mask, proxy, c("Self interview", "Proxy interview"), "Baseline interview type"),
    continuous_row(d, mask, "final_count", "Final-policy condition count", "conditions"),
    continuous_row(d, mask, "cutoff_count", "Cutoff-policy condition count", "conditions"),
    categorical_rows(d, mask, changed, c("No count change", "Policy count changed"), "Policy comparison")
  )
  out <- do.call(rbind, pieces)
  out$group <- group_name
  out$group_n <- sum(mask)
  out
}))
table1 <- table1[c("group", "group_n", "characteristic", "level", "summary_type", "unweighted_n", "missing_n", "estimate", "dispersion", "unit")]
write.csv(table1, file.path(table_dir, "Table1_baseline_characteristics_long.csv"), row.names = FALSE, na = "")

est <- read.csv(file.path(run_dir, "results", "Q2_ESTIMATES.csv"), stringsAsFactors = FALSE)
taylor <- read.csv(file.path(run_dir, "results", "Q2_TAYLOR_POLICY_SPECIFIC.csv"), stringsAsFactors = FALSE)

pull_est <- function(analysis, metric) {
  x <- est[est$analysis == analysis & est$metric == metric, , drop = FALSE]
  if (nrow(x) != 1) stop(paste("missing or duplicate estimate", analysis, metric))
  x
}

primary_rows <- list(
  transform(taylor[taylor$analysis == "W14_W15_primary" & taylor$policy == "final", ],
            result = "Per-condition prevalence ratio", policy = "Final retrospective", estimate = pr, conf_low = pr_low, conf_high = pr_high),
  transform(taylor[taylor$analysis == "W14_W15_primary" & taylor$policy == "cutoff", ],
            result = "Per-condition prevalence ratio", policy = "Index-wave cutoff", estimate = pr, conf_low = pr_low, conf_high = pr_high),
  transform(pull_est("W14_W15_primary", "delta_log_pr"), result = "Paired difference in log prevalence ratios", policy = "Cutoff minus final"),
  transform(pull_est("W14_W15_primary", "pr_ratio"), result = "Ratio of per-condition prevalence ratios", policy = "Cutoff divided by final"),
  transform(pull_est("W14_W15_primary", "stdprev_final"), result = "Standardized follow-up ADL-difficulty prevalence", policy = "Final retrospective"),
  transform(pull_est("W14_W15_primary", "stdprev_cutoff"), result = "Standardized follow-up ADL-difficulty prevalence", policy = "Index-wave cutoff"),
  transform(pull_est("W14_W15_categories", "delta_pd_1"), result = "Standardized prevalence-difference contrast versus 0 conditions", policy = "Paired policy difference: 1 condition"),
  transform(pull_est("W14_W15_categories", "delta_pd_2"), result = "Standardized prevalence-difference contrast versus 0 conditions", policy = "Paired policy difference: 2 conditions"),
  transform(pull_est("W14_W15_categories", "delta_pd_3plus"), result = "Standardized prevalence-difference contrast versus 0 conditions", policy = "Paired policy difference: 3+ conditions")
)
table2 <- do.call(rbind, lapply(primary_rows, function(x) {
  keep <- intersect(c("result", "policy", "estimate", "std_error", "conf_low", "conf_high", "df"), names(x))
  y <- x[1, keep, drop = FALSE]
  for (nm in setdiff(c("std_error", "conf_low", "conf_high", "df"), names(y))) y[[nm]] <- NA_real_
  y[c("result", "policy", "estimate", "std_error", "conf_low", "conf_high", "df")]
}))
table2$analysis_population <- "W14-to-W15 exposure-identified follow-up survivors"
table2$interpretation <- c(
  "Association estimate under final retrospectively revised history",
  "Association estimate under identifiable index-wave cutoff history",
  "Primary paired estimand; negative means a slightly shallower cutoff-policy slope",
  "Primary interpretation scale; prespecified practical-equivalence interval 0.9091 to 1.1000",
  "Secondary standardized prevalence under final-policy model",
  "Secondary standardized prevalence under cutoff-policy model",
  "Secondary categorical contrast; interval includes zero",
  "Secondary categorical contrast; interval includes zero",
  "Secondary categorical contrast; small negative policy difference"
)
write.csv(table2, file.path(table_dir, "Table2_primary_and_secondary_results.csv"), row.names = FALSE, na = "")

sensitivity_map <- data.frame(
  analysis = c(
    "W14_W15_primary", "W14_W15_weights_none", "W14_W15_weights_2.5_97.5",
    "W14_W15_complete_items", "W14_W15_self_only", "W14_W15_bmi_adjusted",
    "W14_W15_extreme_all0", "W14_W15_extreme_all1", "W14_W15_adverse_composite",
    "W14_W15_leave_out_hypertension", "W14_W15_leave_out_diabetes",
    "W14_W15_leave_out_cancer", "W14_W15_leave_out_lung_disease",
    "W14_W15_leave_out_heart_disease", "W14_W15_leave_out_stroke",
    "W14_W15_leave_out_psychiatric_problem", "W14_W15_leave_out_arthritis",
    "W14_W15_mnar_delta_0p5", "W14_W15_mnar_delta_1", "W14_W15_mnar_delta_2",
    "W14_W15_mnar_delta_3", "W15_W16_primary"
  ),
  display_order = 1:22,
  label = c(
    "Primary analysis", "No observation-weight truncation", "2.5th/97.5th truncation",
    "Complete five-item outcome", "Self interviews only", "BMI-adjusted outcome model",
    "Ambiguous F4/5 assigned 0", "Ambiguous F4/5 assigned 1", "Death/known-difficulty composite",
    "Leave out hypertension", "Leave out diabetes", "Leave out cancer", "Leave out lung disease",
    "Leave out heart disease", "Leave out stroke", "Leave out psychiatric problem", "Leave out arthritis",
    "MNAR odds multiplier 0.5", "MNAR odds multiplier 1", "MNAR odds multiplier 2",
    "MNAR odds multiplier 3", "Wave 15 to 16 time-horizon sensitivity"
  ),
  tier = c("Primary", rep("Sensitivity", 8), rep("Leave-one-condition-out", 8),
           rep("Sensitivity", 4), "Time-horizon sensitivity"),
  stringsAsFactors = FALSE
)

ratio <- est[est$metric == "pr_ratio" & est$analysis %in% sensitivity_map$analysis, ]
delta <- est[est$metric == "delta_log_pr" & est$analysis %in% sensitivity_map$analysis, c("analysis", "estimate")]
names(delta)[2] <- "delta_log_pr"
table3 <- merge(sensitivity_map, ratio, by = "analysis", all.x = TRUE, sort = FALSE)
table3 <- merge(table3, delta, by = "analysis", all.x = TRUE, sort = FALSE)
table3 <- table3[order(table3$display_order), ]
table3$equivalence_low <- 1 / 1.10
table3$equivalence_high <- 1.10
table3$equivalence_conclusion <- ifelse(table3$conf_low > table3$equivalence_low & table3$conf_high < table3$equivalence_high,
                                       "Inside prespecified practical-equivalence region", "Not fully inside region")
diagnostics <- fromJSON(file.path(run_dir, "results", "Q2_REPLICATE_DIAGNOSTICS.json"), simplifyVector = FALSE)
table3$target_n <- vapply(table3$analysis, function(x) as.numeric(diagnostics[[x]]$target_n), numeric(1))
table3$outcome_observed_n <- vapply(table3$analysis, function(x) as.numeric(diagnostics[[x]]$outcome_observed_n), numeric(1))
table3$target_change <- ifelse(
  grepl("extreme_all", table3$analysis), "Expanded F4/5 scenario-specific survivor domain",
  ifelse(table3$analysis == "W14_W15_adverse_composite", "Adverse-state target includes deaths and known ADL difficulty",
  ifelse(grepl("mnar_delta", table3$analysis), "Completed-outcome survivor target under the specified MNAR multiplier",
  ifelse(grepl("leave_out", table3$analysis), "Seven-condition count; primary survivor domain retained",
  ifelse(table3$analysis == "W14_W15_complete_items", "Outcome-observed set requires all five ADL items",
  ifelse(table3$analysis == "W14_W15_self_only", "Self-interview survivor domain",
  ifelse(table3$analysis == "W15_W16_primary", "Later non-equivalent time-horizon survivor domain",
         "Primary survivor target retained")))))))
primary_summary <- fromJSON(file.path(run_dir, "results", "Q2_PRIMARY_RESULTS.json"))
table3$result_version <- paste0("Q3-RECONCILED-v2 / Q2 rerun ", primary_summary$analysis_timestamp)
write.csv(table3, file.path(table_dir, "Table3_sensitivity_results.csv"), row.names = FALSE, na = "")
write.csv(table3, file.path(table_dir, "TableS3_sensitivity_results.csv"), row.names = FALSE, na = "")

qc14 <- fromJSON(file.path(run_dir, "derived", "q2_hrs_w14_w15_qc.json"))
qc15 <- fromJSON(file.path(run_dir, "derived", "q2_hrs_w15_w16_qc.json"))
flow <- data.frame(
  period = rep(c("Wave 14 to 15", "Wave 15 to 16"), each = 8),
  stage = rep(c("Positive-weight full frame", "Strict baseline eligible", "Exposure identified", "Follow-up survivor target", "Outcome observed", "Outcome positive", "Alive outcome unknown", "Deaths in identified baseline"), 2),
  count = c(
    qc14$full_weighted_frame, qc14$eligible_base, qc14$eligible_exposure_identified, qc14$identified_followup_survivors,
    qc14$identified_survivors_outcome_observed, qc14$identified_survivors_outcome_positive, qc14$alive_unknown_in_primary_target, qc14$deaths_in_identified_baseline,
    qc15$full_weighted_frame, qc15$eligible_base, qc15$eligible_exposure_identified, qc15$identified_followup_survivors,
    qc15$identified_survivors_outcome_observed, qc15$identified_survivors_outcome_positive, qc15$alive_unknown_in_primary_target, qc15$deaths_in_identified_baseline
  ),
  stringsAsFactors = FALSE
)
write.csv(flow, file.path(table_dir, "TableS1_cohort_flow_counts.csv"), row.names = FALSE)

condition_changes <- read.csv(file.path(run_dir, "results", "Q2_CONDITION_POLICY_CHANGES.csv"), stringsAsFactors = FALSE)
condition_changes$weighted_percent <- 100 * condition_changes$weighted_proportion
condition_changes$conf_low_percent <- 100 * condition_changes$conf_low
condition_changes$conf_high_percent <- 100 * condition_changes$conf_high
write.csv(condition_changes, file.path(table_dir, "TableS2_condition_policy_changes.csv"), row.names = FALSE)

plot_forest <- function() {
  dat <- table3[table3$tier != "Leave-one-condition-out", ]
  n <- nrow(dat)
  y <- rev(seq_len(n))
  par(mar = c(4.5, 12.5, 2.5, 1), family = "sans")
  plot(NA, xlim = c(0.90, 1.105), ylim = c(0.5, n + 0.8), xaxt = "n", yaxt = "n", xlab = "Ratio of per-condition prevalence ratios (cutoff / final)", ylab = "", bty = "n")
  axis(1, at = c(0.91, 0.95, 1.00, 1.05, 1.10), labels = c("0.91", "0.95", "1.00", "1.05", "1.10"))
  abline(v = 1, col = "#444444", lwd = 1.2)
  abline(v = c(1 / 1.10, 1.10), col = "#C44E52", lty = 2, lwd = 1.1)
  segments(dat$conf_low, y, dat$conf_high, y, col = ifelse(dat$tier == "Primary", "#005F73", "#4C78A8"), lwd = ifelse(dat$tier == "Primary", 2.5, 1.6))
  points(dat$estimate, y, pch = ifelse(dat$tier == "Primary", 18, 16), cex = ifelse(dat$tier == "Primary", 1.35, 1.0), col = ifelse(dat$tier == "Primary", "#005F73", "#4C78A8"))
  axis(2, at = y, labels = dat$label, las = 1, tick = FALSE, cex.axis = 0.72)
  mtext("Dashed lines: prespecified practical-equivalence bounds (1/1.10 and 1.10)", side = 3, line = 0.4, adj = 0, cex = 0.72)
}

plot_flow <- function() {
  par(mar = c(0.5, 0.5, 0.5, 0.5), family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  box_y <- c(0.90, 0.73, 0.56, 0.39, 0.22)
  labels <- c(
    sprintf("Positive-weight Wave 14 frame\nn = %s", format(qc14$full_weighted_frame, big.mark = ",")),
    sprintf("Strict baseline eligible\nAge >=65; no baseline five-task ADL difficulty\nn = %s", format(qc14$eligible_base, big.mark = ",")),
    sprintf("Eight-condition cutoff state identifiable\nn = %s", format(qc14$eligible_exposure_identified, big.mark = ",")),
    sprintf("Follow-up survivor target\nn = %s", format(qc14$identified_followup_survivors, big.mark = ",")),
    sprintf("Observed follow-up ADL state\nn = %s (%s positive; %s negative)", format(qc14$identified_survivors_outcome_observed, big.mark = ","), format(qc14$identified_survivors_outcome_positive, big.mark = ","), format(qc14$identified_survivors_outcome_observed - qc14$identified_survivors_outcome_positive, big.mark = ","))
  )
  for (i in seq_along(box_y)) {
    rect(0.23, box_y[i] - 0.055, 0.77, box_y[i] + 0.055, border = "#005F73", lwd = 1.5, col = "#E9F5F7")
    text(0.50, box_y[i], labels[i], cex = 0.78)
    if (i < length(box_y)) arrows(0.50, box_y[i] - 0.058, 0.50, box_y[i + 1] + 0.062, length = 0.08, col = "#555555")
  }
  text(0.81, 0.56, sprintf("Excluded: cutoff state not identified\nn = %s", format(qc14$eligible_base - qc14$eligible_exposure_identified, big.mark = ",")), adj = 0, cex = 0.72)
  segments(0.77, 0.56, 0.80, 0.56, col = "#555555")
  text(0.81, 0.39, sprintf("Deaths before follow-up interview\nn = %s", format(qc14$deaths_in_identified_baseline, big.mark = ",")), adj = 0, cex = 0.72)
  segments(0.77, 0.39, 0.80, 0.39, col = "#555555")
  text(0.81, 0.22, sprintf("Alive outcome unknown\nn = %s", format(qc14$alive_unknown_in_primary_target, big.mark = ",")), adj = 0, cex = 0.72)
  segments(0.77, 0.22, 0.80, 0.22, col = "#555555")
}

plot_changes <- function() {
  dat <- condition_changes[condition_changes$pair == "W14_W15", ]
  dat <- dat[order(dat$weighted_percent), ]
  labels <- gsub("_", " ", dat$condition)
  y <- seq_len(nrow(dat))
  par(mar = c(4.5, 8.5, 2.5, 1), family = "sans")
  plot(dat$weighted_percent, y, xlim = c(0, max(dat$conf_high_percent) * 1.12), ylim = c(0.5, nrow(dat) + 0.5), pch = 16, col = "#005F73", xaxt = "n", yaxt = "n", xlab = "Weighted percentage with policy-dependent condition status", ylab = "", bty = "n")
  segments(dat$conf_low_percent, y, dat$conf_high_percent, y, col = "#4C78A8", lwd = 1.5)
  points(dat$weighted_percent, y, pch = 16, col = "#005F73")
  axis(1)
  axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.78)
  mtext("Wave 14 baseline; descriptive construction frequencies only", side = 3, line = 0.4, adj = 0, cex = 0.76)
}

save_plot <- function(stem, width, height, plot_fun) {
  pdf(file.path(figure_dir, paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE)
  plot_fun()
  dev.off()
  png(file.path(figure_dir, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300, type = "cairo")
  plot_fun()
  dev.off()
}

save_plot("Figure1_cohort_flow", 9.0, 6.8, plot_flow)
save_plot("Figure2_policy_ratio_forest", 9.2, 7.4, plot_forest)
save_plot("FigureS1_condition_policy_changes", 7.5, 5.2, plot_changes)

manifest <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  source_files = c("derived/q2_hrs_w14_w15.csv.gz", "results/Q2_ESTIMATES.csv", "results/Q2_TAYLOR_POLICY_SPECIFIC.csv", "results/Q2_CONDITION_POLICY_CHANGES.csv"),
  tables = sort(list.files(table_dir, full.names = FALSE)),
  figures = sort(list.files(figure_dir, full.names = FALSE)),
  planned_main_elements = c("Table1_baseline_characteristics_long.csv", "Table2_primary_and_secondary_results.csv", "Figure2_policy_ratio_forest.pdf"),
  planned_supplementary_elements = c("TableS1_cohort_flow_counts.csv", "TableS2_condition_policy_changes.csv", "TableS3_sensitivity_results.csv", "Figure1_cohort_flow.pdf"),
  internal_not_for_submission = c("FigureS1_condition_policy_changes.pdf", "FigureS1_condition_policy_changes.png", "Table3_sensitivity_results.csv"),
  primary_claim = "The cutoff and final disease-history policies yield practically equivalent per-condition prevalence-ratio estimates in the selected survivor domain.",
  equivalence_ratio_bounds = c(1 / 1.10, 1.10),
  warning = "Wave 15-to-16 is a time-horizon sensitivity, not an independent replication; condition-specific change frequencies are descriptive and not outcome effects."
)
write_json(manifest, file.path(q3_dir, "Q3_TABLE_FIGURE_MANIFEST.json"), auto_unbox = TRUE, pretty = TRUE, digits = 10)

cat("Q3 tables and figures completed\n")
