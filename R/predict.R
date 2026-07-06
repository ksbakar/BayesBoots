
########################################################
## prediction BB Wrapper
########################################################

predict.BayesBoots <- function(object, newdata, times = NULL, ...) {

  class_type <- intersect(c("BB-GLM", "BB-Cox", "BB-CRR", "BB-KM"), class(object))[1]
  switch(
    class_type,
    "BB-Cox" = {
      beta_mat <- object$beta_samples
      S0_mat   <- object$baseline_survival_draws
      base_t   <- object$time
      formula  <- object$formula
      M <- nrow(beta_mat)
      X <- model.matrix(formula, data = newdata)
      # remove intercept if present in beta
      if ("(Intercept)" %in% colnames(X)) {
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
      }
      # align columns safely
      common_vars <- intersect(colnames(X), colnames(beta_mat))
      X <- X[, common_vars, drop = FALSE]
      beta_mat <- beta_mat[, common_vars, drop = FALSE]
      # linear predictor
      lp <- X %*% t(beta_mat)
      hr <- exp(lp)
      # time grid
      if (is.null(times)) {
        times <- base_t
      }
      # baseline survival interpolation
      S0_interp <- apply(S0_mat, 2, function(s0) {
        approx(
          x = base_t,
          y = s0,
          xout = times,
          method = "constant",
          rule = 2
        )$y
      })
      n <- nrow(X)
      Tn <- length(times)
      S <- array(NA_real_, dim = c(Tn, n, M))
      for (m in seq_len(M)) {
        S[,,m] <- outer(
          S0_interp[, m],
          hr[, m],
          FUN = function(s0, rr) s0 ^ rr
        )
      }
      out <- list(
        model = "BB-Cox",
        time = times,
        survival_median = apply(S, c(1,2), median),
        survival_lower  = apply(S, c(1,2), quantile, 0.025),
        survival_upper  = apply(S, c(1,2), quantile, 0.975),
        linear_predictor = lp
      )
      class(out) <- c("predBB", "BayesBoots", out$model)
      out
    },
    "BB-GLM" = {
      beta_mat <- as.matrix(object$beta_samples)
      formula  <- object$formula
      X <- model.matrix(formula, data = newdata)
      if ("(Intercept)" %in% colnames(X)) {
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
      }
      common_vars <- intersect(colnames(X), colnames(beta_mat))
      if (length(common_vars) == 0) {
        stop("No matching variables between X and beta_samples")
      }
      X <- X[, common_vars, drop = FALSE]
      beta_mat <- beta_mat[, common_vars, drop = FALSE]
      # linear predictor
      lp <- X %*% t(beta_mat)
      # response scale depends on GLM family (assume logistic default)
      prob <- 1 / (1 + exp(-lp))
      # odds ratios relative to baseline
      or <- exp(lp)
      out <- list(
        model = "BB-GLM",
        summary = tibble::tibble(
        probability_median = apply(prob, 1, median),
        probability_lower  = apply(prob, 1, quantile, 0.025),
        probability_upper  = apply(prob, 1, quantile, 0.975),
        odds_ratio_median  = apply(or, 1, median),
        odds_ratio_lower   = apply(or, 1, quantile, 0.025),
        odds_ratio_upper   = apply(or, 1, quantile, 0.975)),
        linear_predictor   = lp
      )
      class(out) <- c("predBB", "BayesBoots", out$model)
      out
    },
    stop("Prediction not implemented for model type: ", class_type)
  )
}

## prediction plots

plot_bb_cox_multi <- function(pred, ids = 1:3) {
  df <- do.call(rbind, lapply(ids, function(p) {
    data.frame(
      time = pred$time,
      survival = pred$survival_median[, p],
      id = paste0("ID ", p)
    )
  }))
  ggplot(df, aes(x = time, y = survival, color = id)) +
    geom_line(linewidth = 1.1) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Time",
      y = "Survival probability",
      color = "Group",
      title = "Predicted Survival Curves (BB-Cox)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

plot_bb_glm_prob <- function(pred) {
  if (pred$model != "BB-GLM") {
    stop("This function only works for BB-GLM objects")
  }
  pred <- pred$summary
  df <- data.frame(
    id = seq_along(pred$probability_median),
    median = pred$probability_median,
    lower  = pred$probability_lower,
    upper  = pred$probability_upper
  )
  ggplot2::ggplot(df, ggplot2::aes(x = id, y = median)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2
    ) +
    ggplot2::labs(
      title = "BB-GLM Predicted Probabilities",
      x = "",
      y = "Probability"
    ) +
    ggplot2::theme_minimal()
}

plot.predBB <- function(x, ids=1:3, ...) {
  # ids = only for BB-Cox model
  model_type <- x$model
  switch(
    model_type,
    "BB-Cox" = {
      if (is.null(ids)) {
        ids <- 1:min(3, ncol(x$survival_median))
      }
      plot_bb_cox_multi(
        pred = x,
        ids = ids
      )
    },
    "BB-GLM" = {
      plot_bb_glm_prob(x)
    },
    stop("Unknown prediction model type: ", model_type)
  )
}

########################################################
########################################################

