########################################################
## bb_survival_km.R
########################################################

## Ordinary KM fit on BB resample

fit_km_bb <- function(
    data,
    formula
) {
  fit <- survival::survfit(
    formula,
    data = data
  )
  if(is.null(fit$strata)) {

    return(list(
      time = fit$time,
      surv = fit$surv
    ))
  }
  strata_names <- rep(
    names(fit$strata),
    fit$strata
  )
  list(
    time = fit$time,
    surv = fit$surv,
    strata = strata_names
  )
}

## Weighted posterior summary for KM curves

weighted_summary_km <- function(
    S_mat,
    weights
) {
  med <- numeric(nrow(S_mat))
  low <- numeric(nrow(S_mat))
  up  <- numeric(nrow(S_mat))
  for(i in seq_len(nrow(S_mat))) {
    qs <- weighted_quantile(
      x = S_mat[i, ],
      w = weights,
      probs = c(0.5, 0.025, 0.975)
    )
    med[i] <- qs[1]
    low[i] <- qs[2]
    up[i]  <- qs[3]
  }
  list(
    median = med,
    lower = low,
    upper = up
  )
}

## BB-KM runner

run_bb_km <- function(
    data,
    formula = survival::Surv(time, status) ~ 1,
    M = 400,
    lambda = NULL,
    omega_fn = prior_smoothness
) {
  curves <- generate_bb_draws(
    fit_fn = fit_km_bb,
    data = data,
    M = M,
    formula = formula
  )
  curves <- Filter(
    Negate(is.null),
    curves
  )
  M_eff <- length(curves)
  if(M_eff == 0) {
    stop("All KM bootstrap draws failed")
  }
  if (is.null(lambda)) {
    lambda <- 1 / sqrt(M)
  }
  gibbs <- compute_gibbs_weights(
    draws = curves,
    omega_fn = omega_fn,
    lambda = lambda,
    distance_fn = loss_km_L2
  )
  w <- gibbs$weights
  has_strata <- !is.null(curves[[1]]$strata)
  if(!has_strata) {
    times <- sort(unique(
      unlist(lapply(curves, `[[`, "time"))
    ))
    S_mat <- matrix(
      NA_real_,
      nrow = length(times),
      ncol = M_eff
    )
    for(m in seq_len(M_eff)) {
      S_mat[, m] <- approx(
        x = curves[[m]]$time,
        y = curves[[m]]$surv,
        xout = times,
        method = "constant",
        rule = 2,
        f = 0
      )$y
    }
    summary <- weighted_summary_km(
      S_mat,
      w
    )
    return(structure(
      list(
        curves = curves,
        summary = tibble::tibble(
          time = times,
          median = summary$median,
          lower = summary$lower,
          upper = summary$upper
        ),
        strata = NULL,
        diagnostics = tibble::tibble(
          draw = seq_len(M_eff),
          risk = gibbs$risk,
          omega = gibbs$omega,
          weight = w
        )
      ),
      class = "BayesBoots"
    ))
  }
  strata_levels <- unique(
    unlist(lapply(curves, function(x) x$strata))
  )
  stratified_summary <- purrr::map_dfr(
    strata_levels,
    function(s) {
      curves_s <- lapply(curves, function(x) {
        if(is.null(x$strata)) return(NULL)
        idx <- which(x$strata == s)
        if(length(idx) == 0) return(NULL)
        list(
          time = x$time[idx],
          surv = x$surv[idx]
        )
      })
      curves_s <- Filter(Negate(is.null), curves_s)
      curves_s <- Filter(
        function(x) {
          length(x$time) > 0 &&
            length(x$surv) > 0
        },
        curves_s
      )
      times_s <- sort(unique(
        unlist(lapply(curves_s, `[[`, "time"))
      ))
      S_mat_s <- matrix(
        NA_real_,
        nrow = length(times_s),
        ncol = length(curves_s)
      )
      for(m in seq_along(curves_s)) {
        S_mat_s[, m] <- approx(
          x = curves_s[[m]]$time,
          y = curves_s[[m]]$surv,
          xout = times_s,
          method = "constant",
          rule = 2,
          f = 0
        )$y
      }
      summary_s <- weighted_summary_km(
        S_mat_s,
        w[seq_along(curves_s)]
      )
      tibble::tibble(
        strata = as.character(s),
        time = times_s,
        median = summary_s$median,
        lower = summary_s$lower,
        upper = summary_s$upper
      )
    }
  )
  structure(
    list(
      summary = stratified_summary,
      curves_data = curves,
      diagnostics = tibble::tibble(
        draw = seq_len(M_eff),
        risk = gibbs$risk,
        omega = gibbs$omega,
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

