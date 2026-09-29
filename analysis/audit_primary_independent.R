#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})
options(survey.lonely.psu = "fail")

args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
d <- read.csv(gzfile(file.path(run_dir, "derived", "q2_hrs_w14_w15.csv.gz")), stringsAsFactors = FALSE)

stopifnot(!anyDuplicated(d$source_row))
stopifnot(sum(d$eligible_base == 1) == 7046)
stopifnot(sum(d$eligible_base == 1 & d$exposure_identified == 1) == 7008)
stopifnot(sum(d$primary_survivor_target == 1) == 6611)
stopifnot(sum(d$primary_outcome_observed == 1) == 5982)
stopifnot(sum(d$primary_outcome_observed == 1 & d$y_function == 1, na.rm = TRUE) == 638)
stopifnot(sum(d$eligible_base == 1 & d$exposure_identified == 1 & d$final_count != d$cutoff_count, na.rm = TRUE) == 252)

frame <- d[d$full_weighted_frame == 1, , drop = FALSE]
frame$a10 <- (frame$age - 75) / 10
frame$a10sq <- frame$a10^2
frame$sex_f <- factor(ifelse(frame$sex_code == 1, "M", ifelse(frame$sex_code == 2, "F", "Missing")), levels = c("M", "F", "Missing"))
race <- rep("Missing", nrow(frame))
race[frame$hispanic_code == 1 & !is.na(frame$hispanic_code)] <- "Hispanic"
race[frame$hispanic_code == 0 & frame$race_code == 1 & !is.na(frame$hispanic_code) & !is.na(frame$race_code)] <- "White"
race[frame$hispanic_code == 0 & frame$race_code == 2 & !is.na(frame$hispanic_code) & !is.na(frame$race_code)] <- "Black"
race[frame$hispanic_code == 0 & frame$race_code == 3 & !is.na(frame$hispanic_code) & !is.na(frame$race_code)] <- "Other"
frame$race_f <- factor(race, levels = c("White", "Black", "Hispanic", "Other", "Missing"))
frame$marital_f <- factor(ifelse(is.na(frame$marital_code), "Missing", paste0("Code", frame$marital_code)))
frame$proxy_f <- factor(ifelse(frame$base_proxy == 0, "Self", ifelse(frame$base_proxy == 1, "Proxy", "Missing")), levels = c("Self", "Proxy", "Missing"))
eligible_frame <- frame$eligible_base == 1
edu_median <- median(frame$education_years[eligible_frame], na.rm = TRUE)
bmi_median <- median(frame$bmi[eligible_frame], na.rm = TRUE)
frame$edu_missing <- as.integer(is.na(frame$education_years))
frame$bmi_missing <- as.integer(is.na(frame$bmi))
frame$edu_fill <- ifelse(is.na(frame$education_years), edu_median, frame$education_years)
frame$bmi_fill <- ifelse(is.na(frame$bmi), bmi_median, frame$bmi)

design <- svydesign(ids = ~raehsamp, strata = ~raestrat, weights = ~base_weight, data = frame, nest = TRUE)
rep_design <- as.svrepdesign(design, type = "Fay", fay.rho = 0.5, mse = TRUE)

primary_delta <- function(w, dat) {
  target <- dat$eligible_base == 1 & dat$exposure_identified == 1 & dat$follow_alive == 1
  observed <- target & !is.na(dat$y_function)
  td <- dat[target, , drop = FALSE]
  td$observed <- as.integer(!is.na(dat$y_function[target]))
  ww <- w[target]
  ww <- ww / mean(ww)
  obs_fit <- suppressWarnings(glm(
    observed ~ a10 + a10sq + sex_f + race_f + edu_fill + edu_missing + marital_f + proxy_f + bmi_fill + bmi_missing + cutoff_count,
    data = td,
    family = quasibinomial(link = "logit"),
    weights = ww,
    control = glm.control(maxit = 100, epsilon = 1e-8)
  ))
  stopifnot(obs_fit$converged)
  pr_obs <- as.numeric(predict(obs_fit, type = "response"))
  numerator <- weighted.mean(td$observed, w[target])
  sw <- numerator / pr_obs
  q <- quantile(sw[td$observed == 1], c(.01, .99), names = FALSE)
  sw <- pmax(q[[1]], pmin(q[[2]], sw))
  md <- td[td$observed == 1, , drop = FALSE]
  md$y <- md$y_function
  fit_w <- w[target][td$observed == 1] * sw[td$observed == 1]
  fit_w <- fit_w / mean(fit_w)
  final_fit <- suppressWarnings(glm(
    y ~ final_count + a10 + a10sq + sex_f + race_f + edu_fill + edu_missing + marital_f + proxy_f,
    data = md,
    family = quasipoisson(link = "log"),
    weights = fit_w,
    control = glm.control(maxit = 100, epsilon = 1e-8)
  ))
  cutoff_fit <- suppressWarnings(glm(
    y ~ cutoff_count + a10 + a10sq + sex_f + race_f + edu_fill + edu_missing + marital_f + proxy_f,
    data = md,
    family = quasipoisson(link = "log"),
    weights = fit_w,
    control = glm.control(maxit = 100, epsilon = 1e-8)
  ))
  stopifnot(final_fit$converged, cutoff_fit$converged)
  unname(coef(cutoff_fit)[["cutoff_count"]] - coef(final_fit)[["final_count"]])
}

audit <- withReplicates(rep_design, primary_delta, return.replicates = TRUE)
audit_estimate <- as.numeric(audit$theta)
audit_se <- sqrt(as.numeric(attr(audit$theta, "var")))
reported <- fromJSON(file.path(run_dir, "results", "Q2_PRIMARY_RESULTS.json"))
reported_estimate <- reported$primary_delta_log_pr$estimate
reported_se <- reported$primary_delta_log_pr$std_error

result <- list(
  audit_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  row_key_unique = TRUE,
  cohort_reconciliation = list(
    eligible_base = 7046,
    exposure_identified = 7008,
    survivor_target = 6611,
    outcome_observed = 5982,
    outcome_positive = 638,
    policy_count_changed = 252
  ),
  survey_design = list(
    frame_rows = nrow(frame),
    strata = length(unique(frame$raestrat)),
    psus = nrow(unique(frame[c("raestrat", "raehsamp")])),
    fay_replicates = ncol(weights(rep_design, type = "analysis")),
    degrees_of_freedom = degf(rep_design)
  ),
  independent_primary = list(estimate = audit_estimate, std_error = audit_se),
  reported_primary = list(estimate = reported_estimate, std_error = reported_se),
  absolute_difference = list(estimate = abs(audit_estimate - reported_estimate), std_error = abs(audit_se - reported_se)),
  tolerance = 1e-10,
  reproduced = abs(audit_estimate - reported_estimate) < 1e-10 && abs(audit_se - reported_se) < 1e-10,
  replicate_failure_count = sum(!is.finite(audit$replicates))
)

write_json(result, file.path(run_dir, "results", "Q2_INDEPENDENT_PRIMARY_REPRODUCTION.json"), auto_unbox = TRUE, pretty = TRUE, digits = 12)
if (!isTRUE(result$reproduced)) stop("independent primary reproduction exceeded tolerance")
cat("INDEPENDENT_REPRODUCTION_PASS\n")
print(result$independent_primary)
