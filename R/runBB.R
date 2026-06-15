########################################################
## Unified BB Wrapper
########################################################

runBB <- function(
    data,
    model = c("BB-KM", "BB-Cox", "BB-GLM"),
    formula = NULL,
    family = gaussian(),
    M = 400,
    lambda = 1,
    prior_fn = NULL,
    ...
) {
  options(warn=-1)
  model <- match.arg(model)
  if(is.null(prior_fn)){stop("Define argument: 'prior_fn'.")}
  dispatch <- list(
    "BB-KM" = function() {
      if(is.null(formula)) {
        stop("BB-KM requires a formula, e.g. Surv(time,status) ~ x1 + x2")
      }
      run_bb_km(
        data = data,
        formula = formula,
        M = M,
        lambda = lambda,
        omega_fn = prior_fn
      )
    },
    "BB-Cox" = function() {
      if(is.null(formula)) {
        stop("BB-Cox requires a formula, e.g. Surv(time,status) ~ x1 + x2")
      }
      run_bb_cox(
        data = data,
        formula = formula,
        M = M,
        lambda = lambda,
        omega_fn = prior_fn
      )
    },
    "BB-GLM" = function() {
      if(is.null(formula)) {
        stop("BB-GLM requires a formula, e.g. y ~ x1 + x2")
      }
      run_bb_glm(
        data = data,
        formula = formula,
        family = family,
        M = M,
        lambda = lambda,
        omega_fn = prior_fn
      )
    }
  )
  out <- dispatch[[model]]()
  class(out) <- c(class(out), model)
  out
}

########################################################
########################################################
