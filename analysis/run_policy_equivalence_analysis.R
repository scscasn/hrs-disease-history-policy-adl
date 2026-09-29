#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survey)
  library(jsonlite)
})

options(survey.lonely.psu = "fail")

args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
dir.create(file.path(run_dir, "results"), showWarnings = FALSE, recursive = TRUE)

read_pair <- function(name) {
  read.csv(gzfile(file.path(run_dir, "derived", name)), stringsAsFactors = FALSE)
}

explicit_factor <- function(x, levels, labels) {
  x_chr <- ifelse(is.na(x), "MISSING", as.character(x))
  factor(x_chr, levels = c(levels, "MISSING"), labels = c(labels, "Missing"))
}

prepare <- function(d, pair_label) {
  d$pair_label <- pair_label
  d$age10 <- (d$age - 75) / 10
  d$age10sq <- d$age10^2
  d$sex <- explicit_factor(d$sex_code, c("1", "2"), c("Male", "Female"))
  race <- rep("Missing", nrow(d))
  race[!is.na(d$hispanic_code) & d$hispanic_code == 1] <- "Hispanic"
  race[!is.na(d$hispanic_code) & d$hispanic_code == 0 & d$race_code == 1] <- "White"
  race[!is.na(d$hispanic_code) & d$hispanic_code == 0 & d$race_code == 2] <- "Black"
  race[!is.na(d$hispanic_code) & d$hispanic_code == 0 & d$race_code == 3] <- "Other"
  d$race_eth <- factor(race, levels = c("White", "Black", "Hispanic", "Other", "Missing"))
  marital_levels <- c("1", "2", "3", "4", "5", "6", "7", "8")
  d$marital <- explicit_factor(d$marital_code, marital_levels, paste0("Code", marital_levels))
  d$base_proxy_f <- explicit_factor(d$base_proxy, c("0", "1"), c("Self", "Proxy"))

  eligible <- d$eligible_base == 1
  edu_med <- median(d$education_years[eligible], na.rm = TRUE)
  bmi_med <- median(d$bmi[eligible], na.rm = TRUE)
  d$edu_miss <- as.integer(is.na(d$education_years))
  d$bmi_miss <- as.integer(is.na(d$bmi))
  d$edu_imp <- ifelse(is.na(d$education_years), edu_med, d$education_years)
  d$bmi_imp <- ifelse(is.na(d$bmi), bmi_med, d$bmi)
  attr(d, "edu_median") <- edu_med
  attr(d, "bmi_median") <- bmi_med
  d
}

make_replicates <- function(d) {
  frame <- d[d$full_weighted_frame == 1, , drop = FALSE]
  design <- svydesign(
    ids = ~raehsamp,
    strata = ~raestrat,
    weights = ~base_weight,
    data = frame,
    nest = TRUE
  )
  rep_design <- as.svrepdesign(design, type = "Fay", fay.rho = 0.5, mse = TRUE)
  list(
    data = rep_design$variables,
    full = as.numeric(weights(rep_design, type = "sampling")),
    reps = as.matrix(weights(rep_design, type = "analysis")),
    scale = rep_design$scale,
    rscales = rep_design$rscales,
    df = degf(rep_design),
    design = design,
    rep_design = rep_design
  )
}

outcome_covars <- "age10 + age10sq + sex + race_eth + edu_imp + edu_miss + marital + base_proxy_f"
observation_covars <- paste0(outcome_covars, " + bmi_imp + bmi_miss")

covariates_for_spec <- function(spec, observation = FALSE) {
  covars <- if (isTRUE(spec$self_only)) {
    "age10 + age10sq + sex + race_eth + edu_imp + edu_miss + marital"
  } else {
    outcome_covars
  }
  if (observation) covars <- paste0(covars, " + bmi_imp + bmi_miss")
  covars
}

normalize_fit_weights <- function(w) {
  positive <- is.finite(w) & w > 0
  if (!any(positive)) stop("no positive fitting weights")
  out <- w
  out[positive] <- out[positive] / mean(out[positive])
  out
}

target_mask <- function(data, spec) {
  if (identical(spec$target_mode, "adverse")) {
    mask <- data$eligible_base == 1 & data$follow_iwstat %in% c(1, 4, 5)
  } else {
    mask <- data$eligible_base == 1 & data$follow_alive == 1
  }
  if (isTRUE(spec$extreme)) {
    mask <- mask & data$extreme_scenario_eligible == 1
  } else {
    mask <- mask & data$exposure_identified == 1
  }
  if (isTRUE(spec$self_only)) {
    mask <- mask & data$base_proxy == 0 & data$follow_proxy == 0
  }
  mask & !is.na(mask)
}

weight_components <- function(w, data, spec, outcome_name) {
  target <- target_mask(data, spec)
  y <- data[[outcome_name]]
  observed <- target & !is.na(y)
  idx <- which(target)
  if (length(idx) < 100 || sum(observed) < 80) stop("insufficient target/observed sample")

  obs_data <- data[idx, , drop = FALSE]
  obs_data$obs_ind <- as.integer(!is.na(y[idx]))
  obs_data$cutoff_current <- obs_data[[spec$cutoff_var]]
  obs_formula <- as.formula(paste("obs_ind ~", covariates_for_spec(spec, observation = TRUE), "+ cutoff_current"))
  obs_fit <- suppressWarnings(glm(
    obs_formula,
    data = obs_data,
    family = quasibinomial(link = "logit"),
    weights = normalize_fit_weights(w[idx]),
    control = glm.control(maxit = 100, epsilon = 1e-8)
  ))
  if (!isTRUE(obs_fit$converged)) stop("observation model did not converge")
  pred <- as.numeric(predict(obs_fit, newdata = obs_data, type = "response"))
  if (any(!is.finite(pred)) || any(pred <= 0) || any(pred >= 1)) stop("invalid observation probabilities")
  numerator <- weighted.mean(obs_data$obs_ind, w[idx])
  stabilized <- numerator / pred
  obs_in_target <- which(obs_data$obs_ind == 1)
  truncation <- spec$truncation
  if (identical(truncation, "none")) {
    lo <- -Inf
    hi <- Inf
  } else if (identical(truncation, "2.5_97.5")) {
    qs <- quantile(stabilized[obs_in_target], c(0.025, 0.975), names = FALSE, type = 7)
    lo <- qs[[1]]
    hi <- qs[[2]]
  } else {
    qs <- quantile(stabilized[obs_in_target], c(0.01, 0.99), names = FALSE, type = 7)
    lo <- qs[[1]]
    hi <- qs[[2]]
  }
  stabilized_trunc <- pmax(lo, pmin(hi, stabilized))
  analysis_w <- numeric(length(w))
  analysis_w[idx] <- w[idx] * stabilized_trunc
  analysis_w[!observed] <- 0
  list(
    target = target,
    observed = observed,
    y = y,
    analysis_w = analysis_w,
    numerator = numerator,
    stabilized = stabilized,
    stabilized_trunc = stabilized_trunc,
    limits = c(lo, hi),
    obs_fit = obs_fit
  )
}

fit_poisson <- function(data, rows, w, outcome, count_var, bmi_adjusted = FALSE, self_only = FALSE) {
  model_data <- data[rows, , drop = FALSE]
  model_data$.y <- outcome[rows]
  model_data$.count <- model_data[[count_var]]
  covars <- if (self_only) "age10 + age10sq + sex + race_eth + edu_imp + edu_miss + marital" else outcome_covars
  rhs <- paste(".count +", covars)
  if (bmi_adjusted) rhs <- paste(rhs, "+ bmi_imp + bmi_miss")
  fit <- suppressWarnings(glm(
    as.formula(paste(".y ~", rhs)),
    data = model_data,
    family = quasipoisson(link = "log"),
    weights = normalize_fit_weights(w[rows]),
    control = glm.control(maxit = 100, epsilon = 1e-8)
  ))
  if (!isTRUE(fit$converged)) stop(paste("outcome model failed:", count_var))
  beta <- unname(coef(fit)[[".count"]])
  if (!is.finite(beta)) stop(paste("non-finite exposure coefficient:", count_var))
  list(fit = fit, beta = beta)
}

core_estimator <- function(w, data, spec) {
  wc <- weight_components(w, data, spec, spec$outcome)
  rows <- wc$observed
  final <- fit_poisson(data, rows, wc$analysis_w, wc$y, spec$final_var, spec$bmi_adjusted, spec$self_only)
  cutoff <- fit_poisson(data, rows, wc$analysis_w, wc$y, spec$cutoff_var, spec$bmi_adjusted, spec$self_only)
  target_data <- data[wc$target, , drop = FALSE]
  target_data$.count <- target_data[[spec$final_var]]
  p_final <- as.numeric(predict(final$fit, newdata = target_data, type = "response"))
  target_data$.count <- target_data[[spec$cutoff_var]]
  p_cutoff <- as.numeric(predict(cutoff$fit, newdata = target_data, type = "response"))
  target_w <- w[wc$target]
  delta <- cutoff$beta - final$beta
  c(
    beta_final = final$beta,
    beta_cutoff = cutoff$beta,
    delta_log_pr = delta,
    pr_final = exp(final$beta),
    pr_cutoff = exp(cutoff$beta),
    pr_ratio = exp(delta),
    stdprev_final = weighted.mean(p_final, target_w),
    stdprev_cutoff = weighted.mean(p_cutoff, target_w)
  )
}

mnar_estimator <- function(w, data, spec, delta_odds) {
  wc <- weight_components(w, data, spec, spec$outcome)
  target <- wc$target
  observed <- wc$observed
  idx_obs <- which(observed)
  outcome_data <- data[idx_obs, , drop = FALSE]
  outcome_data$.y <- wc$y[idx_obs]
  outcome_data$.count <- outcome_data[[spec$cutoff_var]]
  impute_formula <- as.formula(paste(".y ~ .count +", covariates_for_spec(spec, observation = FALSE)))
  impute_fit <- suppressWarnings(glm(
    impute_formula,
    data = outcome_data,
    family = quasibinomial(link = "logit"),
    # The frozen MNAR specification calls for a survey-weighted shared
    # reference-outcome model.  Observation IP weights belong to the MAR
    # complete-case analysis, not to this completed-outcome reference model.
    weights = normalize_fit_weights(w[idx_obs]),
    control = glm.control(maxit = 100, epsilon = 1e-8)
  ))
  if (!isTRUE(impute_fit$converged)) stop("MNAR outcome-imputation model failed")
  target_data <- data[target, , drop = FALSE]
  target_data$.count <- target_data[[spec$cutoff_var]]
  p <- as.numeric(predict(impute_fit, newdata = target_data, type = "response"))
  if (any(!is.finite(p)) || any(p <= 0) || any(p >= 1)) stop("MNAR predictions outside (0,1)")
  odds <- p / (1 - p)
  p_delta <- delta_odds * odds / (1 + delta_odds * odds)
  y_complete <- rep(NA_real_, nrow(data))
  y_complete[target] <- ifelse(is.na(wc$y[target]), p_delta, wc$y[target])
  all_target_w <- numeric(length(w))
  all_target_w[target] <- w[target]
  final <- fit_poisson(data, target, all_target_w, y_complete, spec$final_var, FALSE, spec$self_only)
  cutoff <- fit_poisson(data, target, all_target_w, y_complete, spec$cutoff_var, FALSE, spec$self_only)
  td <- data[target, , drop = FALSE]
  td$.count <- td[[spec$final_var]]
  p_final <- as.numeric(predict(final$fit, newdata = td, type = "response"))
  td$.count <- td[[spec$cutoff_var]]
  p_cutoff <- as.numeric(predict(cutoff$fit, newdata = td, type = "response"))
  delta <- cutoff$beta - final$beta
  c(
    beta_final = final$beta,
    beta_cutoff = cutoff$beta,
    delta_log_pr = delta,
    pr_final = exp(final$beta),
    pr_cutoff = exp(cutoff$beta),
    pr_ratio = exp(delta),
    stdprev_final = weighted.mean(p_final, w[target]),
    stdprev_cutoff = weighted.mean(p_cutoff, w[target])
  )
}

category_estimator <- function(w, data, spec) {
  wc <- weight_components(w, data, spec, spec$outcome)
  rows <- wc$observed
  model_data <- data[rows, , drop = FALSE]
  model_data$.y <- wc$y[rows]
  make_cat <- function(x) factor(ifelse(x >= 3, "3+", as.character(x)), levels = c("0", "1", "2", "3+"))
  model_data$.cat_final <- make_cat(model_data$final_count)
  model_data$.cat_cutoff <- make_cat(model_data[[spec$cutoff_var]])
  f_final <- suppressWarnings(glm(
    as.formula(paste(".y ~ .cat_final +", covariates_for_spec(spec, observation = FALSE))),
    data = model_data,
    family = quasipoisson(link = "log"),
    weights = normalize_fit_weights(wc$analysis_w[rows]),
    control = glm.control(maxit = 100, epsilon = 1e-8)
  ))
  f_cutoff <- suppressWarnings(glm(
    as.formula(paste(".y ~ .cat_cutoff +", covariates_for_spec(spec, observation = FALSE))),
    data = model_data,
    family = quasipoisson(link = "log"),
    weights = normalize_fit_weights(wc$analysis_w[rows]),
    control = glm.control(maxit = 100, epsilon = 1e-8)
  ))
  if (!isTRUE(f_final$converged) || !isTRUE(f_cutoff$converged)) stop("category model failed")
  td <- data[wc$target, , drop = FALSE]
  tw <- w[wc$target]
  cats <- c("0", "1", "2", "3+")
  prev_final <- prev_cutoff <- numeric(length(cats))
  for (i in seq_along(cats)) {
    td$.cat_final <- factor(cats[[i]], levels = cats)
    td$.cat_cutoff <- factor(cats[[i]], levels = cats)
    prev_final[[i]] <- weighted.mean(as.numeric(predict(f_final, td, type = "response")), tw)
    prev_cutoff[[i]] <- weighted.mean(as.numeric(predict(f_cutoff, td, type = "response")), tw)
  }
  names(prev_final) <- cats
  names(prev_cutoff) <- cats
  ans <- numeric(0)
  for (cat in cats[-1]) {
    safe <- gsub("\\+", "plus", cat)
    pd_f <- prev_final[[cat]] - prev_final[["0"]]
    pd_c <- prev_cutoff[[cat]] - prev_cutoff[["0"]]
    ans[[paste0("pd_final_", safe)]] <- pd_f
    ans[[paste0("pd_cutoff_", safe)]] <- pd_c
    ans[[paste0("delta_pd_", safe)]] <- pd_c - pd_f
  }
  ans
}

run_replicated <- function(rep_obj, estimator, spec, label, extra = list()) {
  full <- do.call(estimator, c(list(w = rep_obj$full, data = rep_obj$data, spec = spec), extra))
  p <- length(full)
  reps <- matrix(NA_real_, nrow = ncol(rep_obj$reps), ncol = p)
  colnames(reps) <- names(full)
  errors <- rep(NA_character_, nrow(reps))
  for (j in seq_len(nrow(reps))) {
    val <- tryCatch(
      do.call(estimator, c(list(w = rep_obj$reps[, j], data = rep_obj$data, spec = spec), extra)),
      error = function(e) e
    )
    if (inherits(val, "error") || length(val) != p || any(!is.finite(val))) {
      errors[[j]] <- if (inherits(val, "error")) conditionMessage(val) else "non-finite or wrong-length estimate"
    } else {
      reps[j, ] <- val
    }
  }
  ok <- complete.cases(reps)
  failed <- which(!ok)
  failure_fraction <- length(failed) / nrow(reps)
  blocked <- failure_fraction > 0.05
  if (sum(ok) < 2) stop(paste("fewer than two successful replicates for", label))
  variance <- survey::svrVar(
    reps[ok, , drop = FALSE],
    scale = rep_obj$scale,
    rscales = rep_obj$rscales[ok],
    mse = TRUE,
    coef = full
  )
  se <- sqrt(diag(variance))
  inference_domain <- target_mask(rep_obj$data, spec)
  if (!identical(estimator, mnar_estimator)) {
    inference_domain <- inference_domain & !is.na(rep_obj$data[[spec$outcome]])
  }
  domain_df <- degf(subset(rep_obj$rep_design, inference_domain))
  crit <- qt(0.975, df = domain_df)
  estimates <- data.frame(
    analysis = label,
    metric = names(full),
    estimate = as.numeric(full),
    std_error = as.numeric(se),
    conf_low = as.numeric(full - crit * se),
    conf_high = as.numeric(full + crit * se),
    df = domain_df,
    scheduled_replicates = nrow(reps),
    successful_replicates = sum(ok),
    failed_replicates = length(failed),
    release_blocked = blocked,
    stringsAsFactors = FALSE
  )
  # The frozen primary estimand is delta_log_pr and pr_ratio is its
  # interpretation-scale transform.  Keep the ratio interval tied to the
  # same log-scale result object rather than constructing a second Wald
  # interval on the nonlinear ratio scale.
  delta_row <- which(estimates$metric == "delta_log_pr")
  ratio_row <- which(estimates$metric == "pr_ratio")
  if (length(delta_row) == 1L && length(ratio_row) == 1L) {
    estimates$estimate[ratio_row] <- exp(estimates$estimate[delta_row])
    estimates$conf_low[ratio_row] <- exp(estimates$conf_low[delta_row])
    estimates$conf_high[ratio_row] <- exp(estimates$conf_high[delta_row])
  }
  list(
    estimates = estimates,
    full = full,
    replicates = reps,
    diagnostics = list(
      analysis = label,
      target_n = sum(target_mask(rep_obj$data, spec)),
      outcome_observed_n = sum(target_mask(rep_obj$data, spec) & !is.na(rep_obj$data[[spec$outcome]])),
      target_mode = spec$target_mode,
      outcome = spec$outcome,
      final_count_variable = spec$final_var,
      cutoff_count_variable = spec$cutoff_var,
      scheduled = nrow(reps),
      successful = sum(ok),
      failed = length(failed),
      failure_fraction = failure_fraction,
      domain_degrees_of_freedom = domain_df,
      failed_indices = as.integer(failed),
      failed_messages = as.list(errors[failed]),
      release_blocked = blocked
    )
  )
}

run_point_only <- function(rep_obj, estimator, spec, label, extra = list()) {
  unit_weights <- rep(1, nrow(rep_obj$data))
  full <- do.call(estimator, c(list(w = unit_weights, data = rep_obj$data, spec = spec), extra))
  estimates <- data.frame(
    analysis = label,
    metric = names(full),
    estimate = as.numeric(full),
    std_error = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    df = NA_real_,
    scheduled_replicates = 0,
    successful_replicates = 0,
    failed_replicates = 0,
    release_blocked = FALSE,
    stringsAsFactors = FALSE
  )
  list(
    estimates = estimates,
    full = full,
    replicates = matrix(numeric(0), nrow = 0),
    diagnostics = list(
      analysis = label,
      scheduled = 0,
      successful = 0,
      failed = 0,
      failure_fraction = 0,
      failed_indices = list(),
      failed_messages = list(),
      release_blocked = FALSE,
      note = "Prespecified unweighted point-estimate sensitivity; no interval released."
    )
  )
}

taylor_policy_specific <- function(rep_obj, spec, label) {
  wc <- weight_components(rep_obj$full, rep_obj$data, spec, spec$outcome)
  dd <- rep_obj$data
  dd$taylor_weight <- wc$analysis_w
  dd$.y <- wc$y
  dd$.count_final <- dd[[spec$final_var]]
  dd$.count_cutoff <- dd[[spec$cutoff_var]]
  design <- svydesign(
    ids = ~raehsamp,
    strata = ~raestrat,
    weights = ~taylor_weight,
    data = dd,
    nest = TRUE
  )
  domain <- subset(design, wc$observed)
  reference_domain_df <- degf(subset(rep_obj$design, wc$observed))
  covars <- covariates_for_spec(spec, observation = FALSE)
  final_fit <- svyglm(as.formula(paste(".y ~ .count_final +", covars)), design = domain, family = quasipoisson(link = "log"))
  cutoff_fit <- svyglm(as.formula(paste(".y ~ .count_cutoff +", covars)), design = domain, family = quasipoisson(link = "log"))
  extract <- function(fit, term, policy) {
    b <- unname(coef(fit)[[term]])
    se <- sqrt(unname(vcov(fit)[term, term]))
    crit <- qt(0.975, df = reference_domain_df)
    data.frame(
      analysis = label,
      policy = policy,
      beta = b,
      beta_se = se,
      beta_low = b - crit * se,
      beta_high = b + crit * se,
      pr = exp(b),
      pr_low = exp(b - crit * se),
      pr_high = exp(b + crit * se),
      df = reference_domain_df,
      stringsAsFactors = FALSE
    )
  }
  rbind(extract(final_fit, ".count_final", "final"), extract(cutoff_fit, ".count_cutoff", "cutoff"))
}

full_diagnostics <- function(rep_obj, spec) {
  wc <- weight_components(rep_obj$full, rep_obj$data, spec, spec$outcome)
  sw_obs <- wc$stabilized[!is.na(wc$y[wc$target])]
  sw_trunc_obs <- wc$stabilized_trunc[!is.na(wc$y[wc$target])]
  list(
    target_n = sum(wc$target),
    observed_n = sum(wc$observed),
    positive_n = sum(wc$y[wc$observed] == 1),
    observation_rate_design_weighted = wc$numerator,
    stabilized_weight_summary_untruncated = as.list(setNames(as.numeric(quantile(sw_obs, c(0, .01, .025, .5, .975, .99, 1))), c("min", "p01", "p025", "median", "p975", "p99", "max"))),
    stabilized_weight_summary_truncated = as.list(setNames(as.numeric(quantile(sw_trunc_obs, c(0, .01, .025, .5, .975, .99, 1))), c("min", "p01", "p025", "median", "p975", "p99", "max"))),
    observation_model_converged = isTRUE(wc$obs_fit$converged)
  )
}

run_pair <- function(d, pair_label, include_full_set = FALSE) {
  rep_obj <- make_replicates(d)
  base_spec <- list(
    outcome = "y_function",
    final_var = "final_count",
    cutoff_var = "cutoff_count",
    truncation = "1_99",
    extreme = FALSE,
    self_only = FALSE,
    bmi_adjusted = FALSE,
    target_mode = "survivor"
  )
  analyses <- list()
  analyses[[paste0(pair_label, "_primary")]] <- run_replicated(rep_obj, core_estimator, base_spec, paste0(pair_label, "_primary"))
  if (include_full_set) {
    spec_none <- base_spec; spec_none$truncation <- "none"
    spec_mid <- base_spec; spec_mid$truncation <- "2.5_97.5"
    spec_complete <- base_spec; spec_complete$outcome <- "y_function_complete_items"
    spec_self <- base_spec; spec_self$self_only <- TRUE
    spec_bmi <- base_spec; spec_bmi$bmi_adjusted <- TRUE
    spec_zero <- base_spec; spec_zero$cutoff_var <- "cutoff_count_all0"; spec_zero$extreme <- TRUE
    spec_one <- base_spec; spec_one$cutoff_var <- "cutoff_count_all1"; spec_one$extreme <- TRUE
    analyses[[paste0(pair_label, "_weights_none")]] <- run_replicated(rep_obj, core_estimator, spec_none, paste0(pair_label, "_weights_none"))
    analyses[[paste0(pair_label, "_weights_2.5_97.5")]] <- run_replicated(rep_obj, core_estimator, spec_mid, paste0(pair_label, "_weights_2.5_97.5"))
    analyses[[paste0(pair_label, "_complete_items")]] <- run_replicated(rep_obj, core_estimator, spec_complete, paste0(pair_label, "_complete_items"))
    analyses[[paste0(pair_label, "_self_only")]] <- run_replicated(rep_obj, core_estimator, spec_self, paste0(pair_label, "_self_only"))
    analyses[[paste0(pair_label, "_bmi_adjusted")]] <- run_replicated(rep_obj, core_estimator, spec_bmi, paste0(pair_label, "_bmi_adjusted"))
    analyses[[paste0(pair_label, "_extreme_all0")]] <- run_replicated(rep_obj, core_estimator, spec_zero, paste0(pair_label, "_extreme_all0"))
    analyses[[paste0(pair_label, "_extreme_all1")]] <- run_replicated(rep_obj, core_estimator, spec_one, paste0(pair_label, "_extreme_all1"))
    analyses[[paste0(pair_label, "_categories")]] <- run_replicated(rep_obj, category_estimator, base_spec, paste0(pair_label, "_categories"))
    analyses[[paste0(pair_label, "_unweighted")]] <- run_point_only(rep_obj, core_estimator, base_spec, paste0(pair_label, "_unweighted"))
    spec_adverse <- base_spec; spec_adverse$outcome <- "adverse_state"; spec_adverse$target_mode <- "adverse"
    analyses[[paste0(pair_label, "_adverse_composite")]] <- run_replicated(rep_obj, core_estimator, spec_adverse, paste0(pair_label, "_adverse_composite"))
    condition_names <- c("hypertension", "diabetes", "cancer", "lung_disease", "heart_disease", "stroke", "psychiatric_problem", "arthritis")
    for (condition in condition_names) {
      final_var <- paste0("final_count_without_", condition)
      cutoff_var <- paste0("cutoff_count_without_", condition)
      rep_obj$data[[final_var]] <- rep_obj$data$final_count - rep_obj$data[[paste0("final_", condition)]]
      rep_obj$data[[cutoff_var]] <- rep_obj$data$cutoff_count - rep_obj$data[[paste0("cutoff_", condition)]]
      spec_loo <- base_spec; spec_loo$final_var <- final_var; spec_loo$cutoff_var <- cutoff_var
      label <- paste0(pair_label, "_leave_out_", condition)
      analyses[[label]] <- run_replicated(rep_obj, core_estimator, spec_loo, label)
    }
    for (delta_odds in c(0.5, 1, 2, 3)) {
      label <- paste0(pair_label, "_mnar_delta_", gsub("\\.", "p", as.character(delta_odds)))
      analyses[[label]] <- run_replicated(rep_obj, mnar_estimator, base_spec, label, list(delta_odds = delta_odds))
    }
  }
  list(
    analyses = analyses,
    design = list(
      rows = nrow(rep_obj$data),
      strata = length(unique(rep_obj$data$raestrat)),
      psus = nrow(unique(rep_obj$data[c("raestrat", "raehsamp")])),
      replicate_count = ncol(rep_obj$reps),
      df = rep_obj$df,
      scale = rep_obj$scale,
      rscales_unique = unique(rep_obj$rscales)
    ),
    primary_diagnostics = full_diagnostics(rep_obj, base_spec),
    replicate_weights = rep_obj$reps,
    taylor_primary = taylor_policy_specific(rep_obj, base_spec, paste0(pair_label, "_primary"))
  )
}

d1415 <- prepare(read_pair("q2_hrs_w14_w15.csv.gz"), "W14_W15")
d1516 <- prepare(read_pair("q2_hrs_w15_w16.csv.gz"), "W15_W16")

res1415 <- run_pair(d1415, "W14_W15", include_full_set = TRUE)
res1516 <- run_pair(d1516, "W15_W16", include_full_set = FALSE)

all_analyses <- c(res1415$analyses, res1516$analyses)
estimate_table <- do.call(rbind, lapply(all_analyses, function(x) x$estimates))
row.names(estimate_table) <- NULL
write.csv(estimate_table, file.path(run_dir, "results", "Q2_ESTIMATES.csv"), row.names = FALSE)
write.csv(rbind(res1415$taylor_primary, res1516$taylor_primary), file.path(run_dir, "results", "Q2_TAYLOR_POLICY_SPECIFIC.csv"), row.names = FALSE)

replicate_diagnostics <- lapply(all_analyses, function(x) x$diagnostics)
write_json(replicate_diagnostics, file.path(run_dir, "results", "Q2_REPLICATE_DIAGNOSTICS.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")

saveRDS(res1415$replicate_weights, file.path(run_dir, "derived", "q2_w14_w15_fay_analysis_weights.rds"), version = 3)
saveRDS(res1516$replicate_weights, file.path(run_dir, "derived", "q2_w15_w16_fay_analysis_weights.rds"), version = 3)

primary_row <- subset(estimate_table, analysis == "W14_W15_primary" & metric == "delta_log_pr")
ratio_row <- subset(estimate_table, analysis == "W14_W15_primary" & metric == "pr_ratio")
equiv <- log(1.10)
verdict <- if (isTRUE(primary_row$release_blocked)) {
  "INTERVAL_BLOCKED_REPLICATE_FAILURE"
} else if (primary_row$conf_low > -equiv && primary_row$conf_high < equiv) {
  "PRACTICAL_EQUIVALENCE_CRITERION_MET"
} else if (primary_row$estimate > equiv || primary_row$estimate < -equiv) {
  "POINT_ESTIMATE_EXCEEDS_10_PERCENT_THRESHOLD"
} else {
  "INCONCLUSIVE_RELATIVE_TO_EQUIVALENCE_REGION"
}

summary <- list(
  analysis_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  software = list(R = R.version.string, survey = as.character(packageVersion("survey")), jsonlite = as.character(packageVersion("jsonlite"))),
  design = list(W14_W15 = res1415$design, W15_W16 = res1516$design),
  primary_diagnostics = list(W14_W15 = res1415$primary_diagnostics, W15_W16 = res1516$primary_diagnostics),
  equivalence_bound_log = equiv,
  primary_delta_log_pr = as.list(primary_row[1, c("estimate", "std_error", "conf_low", "conf_high", "failed_replicates", "release_blocked")]),
  primary_pr_ratio = as.list(ratio_row[1, c("estimate", "std_error", "conf_low", "conf_high")]),
  prespecified_primary_interpretation = verdict,
  note = "Intervals are paired deterministic Fay-BRR intervals; every replicate refits the common observation model and both policy models."
)
write_json(summary, file.path(run_dir, "results", "Q2_PRIMARY_RESULTS.json"), auto_unbox = TRUE, pretty = TRUE, digits = 10, na = "null")
writeLines(capture.output(sessionInfo()), file.path(run_dir, "results", "Q2_SESSION_INFO.txt"))

cat("Completed Q2 locked analyses\n")
print(primary_row)
print(ratio_row)
cat("Primary interpretation:", verdict, "\n")
