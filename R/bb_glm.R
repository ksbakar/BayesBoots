########################################################
## bb_glm.R
########################################################

## GLM - logistic model

fit_glm_bb <- function(
    data,
    formula,
    family = gaussian()
) {

  fit <- glm(
    formula,
    data = data,
    family = family
  )

  stats::coef(fit)
}

## Run BB-GLM

run_bb_glm <- function(
    data,
    formula,
    family = gaussian(),
    M = 400,
    lambda = 1,
    omega_fn = omega_ridge
) {
  raw_draws <- generate_bb_draws(
    fit_fn = fit_glm_bb,
    data = data,
    M = M,
    formula = formula,
    family = family
  )
  beta_mat <- do.call(rbind, raw_draws)
  beta_mat <- as.matrix(beta_mat)
  or_mat <- exp(beta_mat)
  gibbs_obj <- compute_gibbs_weights(
    draws = split(beta_mat, row(beta_mat)),
    omega_fn = omega_fn,
    lambda = lambda
  )
  w <- gibbs_obj$weights
  beta_sum_raw <- posterior_summary(
    beta_mat,
    w
  )
  beta_summary <- tibble::tibble(
    term = colnames(beta_mat),
    beta = beta_sum_raw$median,
    beta_lower = beta_sum_raw$lower,
    beta_upper = beta_sum_raw$upper
  )
  or_sum_raw <- posterior_summary(
    or_mat,
    w
  )
  or_summary <- tibble::tibble(
    term = colnames(or_mat),
    OR = or_sum_raw$median,
    OR_lower = or_sum_raw$lower,
    OR_upper = or_sum_raw$upper
  )
  results <- dplyr::left_join(
    beta_summary |>
      dplyr::transmute(
        term,
        beta = beta
      ),
    or_summary |>
      dplyr::rename(
        OR = OR,
        OR_lower = OR_lower,
        OR_upper = OR_upper
      ),
    by = "term"
  )
  structure(
    list(
      beta_samples = beta_mat,
      or_samples = or_mat,
      beta_summary = beta_summary,
      or_summary = or_summary,
      summary = results,
      weights = w,
      diagnostics = tibble::tibble(
        draw = seq_len(nrow(beta_mat)),
        risk = gibbs_obj$risk,
        omega = gibbs_obj$omega,
        weight = w
      ),
      prior_fn = omega_fn,
      formula = formula
    ),
    class = "BayesBoots"
  )
}

########################################################
########################################################
