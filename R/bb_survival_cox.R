########################################################
## bb_survival_cox.R
########################################################

## Ordinary Cox fit on BB resample

fit_cox_bb <- function(
    data,
    formula
) {
  fit <- survival::coxph(
    formula,
    data = data,
    ties = "breslow",
    x = TRUE,
    model = TRUE
  )
  bh <- survival::basehaz(
    fit,
    centered = FALSE
  )
  list(
    coef = stats::coef(fit),
    basehaz = bh
  )
}

## Run BB-Cox Gibbs posterior

run_bb_cox <- function(
    data,
    formula,
    M = 400,
    lambda = NULL,
    omega_fn = omega_ridge
) {
  raw_draws <- generate_bb_draws(
    fit_fn = fit_cox_bb,
    data = data,
    M = M,
    formula = formula
  )
  raw_draws <- Filter(
    Negate(is.null),
    raw_draws
  )
  M_eff <- length(raw_draws)
  if(M_eff == 0) {
    stop("All BB-Cox draws failed")
  }
  beta_mat <- do.call(
    rbind,
    lapply(raw_draws, `[[`, "coef")
  )
  beta_mat <- as.matrix(beta_mat)
  hr_mat <- exp(beta_mat)
  if (is.null(lambda)) {
    lambda <- 1 / sqrt(M)
  }
  gibbs_obj <- compute_gibbs_weights(
    draws = split(beta_mat, row(beta_mat)),
    omega_fn = omega_fn,
    lambda = lambda
  )
  w <- gibbs_obj$weights
  beta_summary_raw <- posterior_summary(
    beta_mat,
    w
  )
  beta_summary <- tibble::tibble(
    term = colnames(beta_mat),
    beta = beta_summary_raw$median,
    beta_lower = beta_summary_raw$lower,
    beta_upper = beta_summary_raw$upper
  )
  hr_summary_raw <- posterior_summary(
    hr_mat,
    w
  )
  hr_summary <- tibble::tibble(
    term = colnames(hr_mat),
    hr = hr_summary_raw$median,
    hr_lower = hr_summary_raw$lower,
    hr_upper = hr_summary_raw$upper
  )
  times <- sort(unique(
    unlist(
      lapply(
        raw_draws,
        function(x) x$basehaz$time
      )
    )
  ))
  S0_mat <- matrix(
    NA_real_,
    nrow = length(times),
    ncol = M_eff
  )
  for(m in seq_len(M_eff)) {
    bh <- raw_draws[[m]]$basehaz
    H0 <- approx(
      x = bh$time,
      y = bh$hazard,
      xout = times,
      method = "constant",
      rule = 2,
      f = 0
    )$y
    S0_mat[, m] <- exp(-H0)
  }
  S0_summary <- weighted_summary_km(
    S0_mat,
    w
  )
  baseline_survival <- tibble::tibble(
    time = times,
    median = S0_summary$median,
    lower = S0_summary$lower,
    upper = S0_summary$upper
  )
  results <- dplyr::left_join(
    beta_summary |>
      dplyr::transmute(
        term,
        beta = beta
      ),
    hr_summary |>
      dplyr::rename(
        HR = hr,
        HR_lower = hr_lower,
        HR_upper = hr_upper
      ),
    by = "term"
  )
  structure(
    list(
      beta_samples = beta_mat,
      hr_samples = hr_mat,
      beta_summary = beta_summary,
      hr_summary = hr_summary,
      summary = results,
      baseline_survival_draws = S0_mat,
      baseline_survival = baseline_survival,
      time = times,
      diagnostics = tibble::tibble(
        draw = seq_len(M_eff),
        risk = gibbs_obj$risk,
        omega = gibbs_obj$omega,
        weight = w
      ),
      lambda = lambda,
      prior_fn = omega_fn,
      formula = formula
    ),
    class = "BayesBoots"
  )
}

########################################################
########################################################
