########################################################
## bb_survival_crr.R
########################################################

## Subdistribution Analysis of Competing Risks on BB resample
## Fine JP and Gray RJ (1999) using cmprsk

fit_crr_bb <- function(
    data,
    formula,
    ...
) {
  ##
  var_list <- all.vars(formula)
  data <- data[complete.cases(data),]
  time <- data[,var_list[1]]
  status <- data[,var_list[2]]
  f <- as.formula(paste(paste(var_list[1])," ~ ", paste(var_list[-c(1,2)], collapse= "+")))
  X <- model.matrix(f, data = data)[, -1, drop = FALSE]
  fit <- cmprsk::crr(
    ftime = time,
    fstatus = status,
    cov1 = X,
    ...
  )
  list(
    coef = fit$coef,
    time = time
  )
}

## Run BB-Cox Gibbs posterior

run_bb_crr <- function(
    data,
    formula,
    M = 400,
    lambda = NULL,
    omega_fn = prior_normal,
    ...
) {
  raw_draws <- generate_bb_draws(
    fit_fn = fit_crr_bb,
    data = data,
    M = M,
    formula = formula,
    ...
  )
  raw_draws <- Filter(
    Negate(is.null),
    raw_draws
  )
  M_eff <- length(raw_draws)
  if(M_eff == 0) {
    stop("All BB-CRR draws failed")
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
        function(x) x$time
      )
    )
  ))
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
      time = times,
      diagnostics = tibble::tibble(
        draw = seq_len(M_eff),
        risk = gibbs_obj$risk,
        omega = gibbs_obj$omega,
        weight = w
      ),
      lambda = lambda,
      prior_fn = omega_fn,
      data = data,
      formula = formula
    ),
    class = "BayesBoots"
  )
}

########################################################
########################################################
