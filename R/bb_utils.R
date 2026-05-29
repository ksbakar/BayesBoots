########################################################
## bb_utils.R
########################################################

## Stable softmax

stable_softmax <- function(x) {
  x <- x - max(x)
  ex <- exp(x)
  s <- sum(ex)
  if(!is.finite(s) || s <= 0) {
    return(rep(1 / length(x), length(x)))
  }
  ex / s
}

## L2 loss

loss_L2 <- function(x1, x2) {
  mean((x1 - x2)^2, na.rm = TRUE)
}

## L2 distance for Kaplan-Meier curves

loss_km_L2 <- function(curve1, curve2) {
  times <- sort(unique(c(
    curve1$time,
    curve2$time
  )))
  s1 <- approx(
    x = curve1$time,
    y = curve1$surv,
    xout = times,
    method = "constant",
    rule = 2,
    f = 0
  )$y
  s2 <- approx(
    x = curve2$time,
    y = curve2$surv,
    xout = times,
    method = "constant",
    rule = 2,
    f = 0
  )$y
  mean((s1 - s2)^2, na.rm = TRUE)
}

## Empirical risk

compute_empirical_risk <- function(
    draws,
    distance_fn
) {
  M <- length(draws)
  risk <- numeric(M)
  D <- matrix(0, M, M)
  for(i in seq_len(M)) {
    for(j in i:M) {
      d <- distance_fn(draws[[i]], draws[[j]])
      D[i, j] <- d
      D[j, i] <- d
    }
  }
  rowMeans(D)
}

compute_empirical_risk2 <- function(
    draws,
    distance_fn = loss_L2
) {
  M <- length(draws)
  risk <- numeric(M)
  for(m in seq_len(M)) {
    risk[m] <- mean(
      sapply(draws, function(d) {
        distance_fn(draws[[m]], d)
      })
    )
  }
  risk
}

## Gibbs posterior weights

compute_gibbs_weights <- function(
    draws,
    omega_fn,
    lambda = 1,
    distance_fn = loss_L2
) {
  risk <- compute_empirical_risk(
    draws,
    distance_fn
  )
  omega <- sapply(draws, omega_fn)
  score <- -lambda * risk - omega
  weights <- stable_softmax(score)
  list(
    weights = weights,
    risk = risk,
    omega = omega
  )
}

## Weighted quantile

weighted_quantile <- function(
    x,
    w,
    probs
) {
  stopifnot(length(x) == length(w))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  w <- w / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) {
    x[which(cw >= p)[1]]
  })
}

## Posterior summary

posterior_summary <- function(
    beta_mat,
    weights
) {
  p <- ncol(beta_mat)
  out <- vector("list", p)
  for(j in seq_len(p)) {
    qs <- weighted_quantile(
      beta_mat[, j],
      weights,
      probs = c(0.025, 0.5, 0.975)
    )
    out[[j]] <- tibble::tibble(
      parameter = colnames(beta_mat)[j],
      mean = sum(beta_mat[, j] * weights),
      lower = qs[1],
      median = qs[2],
      upper = qs[3]
    )
  }
  dplyr::bind_rows(out)
}

## Bayesian bootstrap weights

bb_weights <- function(n) {
  w <- rexp(n)
  w / sum(w)
}

## Generate BB posterior draws

generate_bb_draws <- function(
    fit_fn,
    data,
    M = 400,
    n_resample = nrow(data),
    ...
) {
  draws <- vector("list", M)
  n <- nrow(data)
  for(m in seq_len(M)) {
    w <- bb_weights(n)
    idx <- sample(
      seq_len(n),
      size = n_resample,
      replace = TRUE,
      prob = w
    )
    boot_data <- data[idx, , drop = FALSE]
    draws[[m]] <- fit_fn(
      data = boot_data,
      ...
    )
  }
  draws
}

## print fn

print.BayesBoots <- function(x, ...) {
  cat("\n############################################\n")
  cat(" Bayesian Bootstrap Induced Survival Models \n")
  cat("############################################\n\n")
  cat("Model class:", paste(class(x)[class(x) != "BayesBoots"], collapse = " | "), "\n\n")
  if (!is.null(x$formula)) {
    cat("Call:\n")
    print(x$formula)
    cat("\n")
  }
  cat("############################################\n")
  invisible(x)
}

## summary fn

summary.BayesBoots <- function(object, ..., digits = 4) {

  class_type <- class(object)[1]
  cat("\n############################################\n")
  cat(" Bayesian Bootstrap Induced Survival Models \n")
  cat("############################################\n\n")
  cat("Model type: ", class_type, "\n", sep = "")
  cat("Call:\n")
  cat(" ", deparse(object$formula), "\n")
  switch(
    class_type,
    "BB-GLM" = {
      cat("\n BB Posterior estimates:\n")
      print(data.frame(object$summary), digits = digits)
    },
    "BB-Cox" = {
      cat("\n BB Posterior estimates:\n")
      print(data.frame(object$summary), digits = digits)
    },
    "BB-KM" = {
      cat("\n Survival Curve:\n")
      print(data.frame(object$summary), digits = digits)
    },
    stop("Unknown BB model type: ", class_type)
  )
  invisible(object)
}

## plot fn

extract_time_var <- function(formula) {
  # survival::Surv(time, status) ~ x
  lhs <- as.character(formula)[2]
  if (grepl("Surv", lhs)) {
    # extract inside Surv(...)
    inside <- sub("Surv\\((.*)\\)", "\\1", lhs)
    # split time, status
    parts <- strsplit(inside, ",")[[1]]
    time_var <- trimws(parts[1])
    return(time_var)
  }
  stop("No Surv() structure found in formula; cannot infer time variable")
}

plot.BayesBoots <- function(object,
                            type = NULL,
                            ...) {
  # type can take: "hr", "or", "density", "baseline"
  class_type <- class(object)[1]
  plot_forest <- function(df, est, lower, upper, main) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 required for plotting")
    }
    ggplot2::ggplot(df,
                    ggplot2::aes(x = .data[[est]],
                                 y = reorder(term, .data[[est]]))) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = .data[[lower]],
                     xmax = .data[[upper]]),
        height = 0.2
      ) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2) +
      ggplot2::labs(title = main, x = "", y = "") +
      ggplot2::theme_minimal()
  }
  plot_survival <- function(df, main) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 required for plotting")
    }
    ggplot2::ggplot(df,
                    ggplot2::aes(x = time, y = median)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        alpha = 0.2
      ) +
      ggplot2::labs(title = main, x = "Time", y = "Survival") +
      ggplot2::theme_minimal()
  }
  switch(
    class_type,
    "BB-KM" = {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required for plotting")
      }
      df <- object$summary
      time_var <- extract_time_var(object$formula)
      # rename dynamically if needed
      if (!time_var %in% names(df)) {
        # fallback: assume it's already named "time"
        time_var <- "time"
      }
      ggplot2::ggplot(df,
                      ggplot2::aes(
                        x = .data[[time_var]],
                        y = median,
                        colour = strata,
                        fill = strata
                      )) +

        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower, ymax = upper),
          alpha = 0.2,
          colour = NA
        ) +
        ggplot2::labs(
          title = "BB-KM Survival Curves",
          x = time_var,
          y = "Survival"
        ) +
        ggplot2::theme_minimal()
    },
    "BB-Cox" = {
      if (is.null(type)) type <- "hr"
      switch(
        type,
        "hr" = plot_forest(
          object$hr_summary,
          est = "hr",
          lower = "hr_lower",
          upper = "hr_upper",
          main = "BB-Cox Hazard Ratios"
        ),
        "baseline" = plot_survival(
          object$baseline_survival,
          main = "BB-Cox Baseline Survival"
        ),
        "density" = {
          if (!requireNamespace("ggplot2", quietly = TRUE)) {
            stop("ggplot2 required")
          }
          df <- as.data.frame(object$hr_samples)
          df_long <- tidyr::pivot_longer(
            df,
            cols = everything(),
            names_to = "term",
            values_to = "value"
          )
          ggplot2::ggplot(df_long, ggplot2::aes(x = value)) +
            ggplot2::geom_density() +
            ggplot2::facet_wrap(~ term, scales = "free") +
            ggplot2::labs(title = "BB-Cox HR Posterior Densities") +
            ggplot2::theme_minimal()
        },
        stop("Unknown type for BB-Cox: ", type)
      )
    },
    "BB-GLM" = {
      if (is.null(type)) type <- "or"
      switch(
        type,
        "or" = plot_forest(
          dplyr::filter(object$or_summary, term != "(Intercept)"),
          est = "OR",
          lower = "OR_lower",
          upper = "OR_upper",
          main = "BB-GLM Odds Ratios"
        ),
        "density" = {
          if (!requireNamespace("ggplot2", quietly = TRUE)) {
            stop("ggplot2 required")
          }
          df <- as.data.frame(object$or_samples)
          df_long <- tidyr::pivot_longer(
            df,
            cols = everything(),
            names_to = "term",
            values_to = "value"
          )
          ggplot2::ggplot(df_long, ggplot2::aes(x = value)) +
            ggplot2::geom_density() +
            ggplot2::facet_wrap(~ term, scales = "free") +
            ggplot2::labs(title = "BB-GLM OR Posterior Densities") +
            ggplot2::theme_minimal()
        },
        stop("Unknown type for BB-GLM: ", type)
      )
    },
    stop("Unknown BB model class: ", class_type)
  )
}

########################################################
########################################################
