########################################################
## priorBB.R
########################################################

## prior generic fn

priorBB <- function(prior_fn){
  prior_fn
}

## ridge prior

prior_ridge <- function(beta,
                        sigma2 = 1) {
  sum(beta^2) / (2 * sigma2)
}

## lasso prior

prior_lasso <- function(beta,
                        tau = 1) {
  tau * sum(abs(beta))
}

## elasticnet prior

prior_elastic <- function(beta,
                          alpha = 0.5,
                          tau = 1) {
  alpha * sum(abs(beta)) +
    (1 - alpha) * sum(beta^2)
}

## Gaussian / Normal prior

prior_normal <- function(beta,
                         mu = 0,
                         sigma2 = 1) {
  sum((beta - mu)^2) / (2 * sigma2)
}

## Multivariate Normal prior

prior_mvn <- function(beta,
                      mu = NULL,
                      sigma2 = 1) {
  p <- length(beta)
  if(is.null(mu)) {
    mu <- rep(0, p)
  }
  if(length(sigma2) == 1) {
    sigma2 <- rep(sigma2, p)
  }
  sum((beta - mu)^2 / (2 * sigma2))
}

## Student-t prior

prior_student <- function(beta,
                          mu = 0,
                          sigma2 = 1,
                          df = 3) {
  sum(
    ((df + 1) / 2) *
      log(
        1 + ((beta - mu)^2) /
          (df * sigma2)
      )
  )
}

## Cauchy prior

prior_cauchy <- function(beta,
                         mu = 0,
                         scale = 1) {
  sum(
    log(
      1 + ((beta - mu)^2) / (scale^2)
    )
  )
}

## Horseshoe-style prior

prior_horseshoe <- function(beta,
                            tau = 1,
                            eps = 1e-8) {
  sum(
    log(
      1 + (tau^2) / (beta^2 + eps)
    )
  )
}

## Spike-and-slab style prior

prior_spike_slab <- function(beta,
                             spike_var = 0.01,
                             slab_var = 1,
                             pi = 0.5) {
  dens <- pi * dnorm(
    beta,
    mean = 0,
    sd = sqrt(spike_var)
  ) +
    (1 - pi) * dnorm(
      beta,
      mean = 0,
      sd = sqrt(slab_var)
    )
  -sum(log(dens + 1e-12))
}

## Smoothness prior for KM curves

prior_smoothness <- function(curve) {
  s <- curve$surv
  if(length(s) < 4) {
    return(0)
  }
  sum(diff(s, differences = 2)^2)
}

## roughness Prior for KM (opposite behavior)

prior_roughness <- function(curve) {
  # roughness Prior
  s <- curve$surv
  if(length(s) < 4) {
    return(0)
  }
  sum(abs(diff(s)))
}

## hazrd prior for KM

prior_hazard_stability <- function(curve) {
  # Hazard stability prior
  s <- curve$surv
  t <- curve$time
  if(length(s) < 3) {
    return(0)
  }
  ds <- diff(s)
  dt <- diff(t)
  h <- -ds / (dt * s[-length(s)] + 1e-8)
  h <- h[is.finite(h)]
  if(length(h) < 2) {
    return(0)
  }
  var(h)
}

## early failure for KM

prior_early_failure <- function(curve) {
  # Early failure prior
  t <- curve$time
  s <- curve$surv
  if(max(t) <= 0) {
    return(0)
  }
  w <- exp(-t / max(t))
  sum(w * (1 - s))
}

########################################################
########################################################
