#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "fail")

args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
conditions <- c("hypertension", "diabetes", "cancer", "lung_disease", "heart_disease", "stroke", "psychiatric_problem", "arthritis")

summarize_pair <- function(file, pair) {
  d <- read.csv(gzfile(file.path(run_dir, "derived", file)), stringsAsFactors = FALSE)
  frame <- d[d$full_weighted_frame == 1, , drop = FALSE]
  design <- svydesign(ids = ~raehsamp, strata = ~raestrat, weights = ~base_weight, data = frame, nest = TRUE)
  domain_mask <- frame$eligible_base == 1 & frame$exposure_identified == 1
  rows <- lapply(conditions, function(condition) {
    final <- frame[[paste0("final_", condition)]]
    cutoff <- frame[[paste0("cutoff_", condition)]]
    changed <- as.integer(final != cutoff)
    frame$.changed <- changed
    one_design <- update(design, .changed = changed)
    est <- svymean(~.changed, subset(one_design, domain_mask), na.rm = TRUE)
    ci <- confint(est, df = degf(subset(one_design, domain_mask)))
    data.frame(
      pair = pair,
      condition = condition,
      domain_n = sum(domain_mask),
      changed_n = sum(changed[domain_mask], na.rm = TRUE),
      weighted_proportion = as.numeric(coef(est)),
      conf_low = as.numeric(ci[1]),
      conf_high = as.numeric(ci[2]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

out <- rbind(
  summarize_pair("q2_hrs_w14_w15.csv.gz", "W14_W15"),
  summarize_pair("q2_hrs_w15_w16.csv.gz", "W15_W16")
)
write.csv(out, file.path(run_dir, "results", "Q2_CONDITION_POLICY_CHANGES.csv"), row.names = FALSE)
cat("CONDITION_CHANGE_SUMMARY_PASS\n")
print(out, row.names = FALSE)
