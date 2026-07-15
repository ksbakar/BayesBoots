\name{runBB}
\alias{runBB}

\title{Bayesian Bootstrap and Gibbs Posterior Framework}

\description{
This function provides an interface for Bayesian Bootstrap (BB) and Gibbs posterior
inference across multiple model classes, including Kaplan-Meier survival curves (BB-KM),
Cox proportional hazards models (BB-Cox), Fine-Gray i.e., competing risk regression (BB-CRR)
and generalized linear models (BB-GLM).

The function generates Bayesian bootstrap posterior draws using Dirichlet weights and
optionally refines them using a Gibbs posterior weighting scheme based on empirical risk
and user-specified prior penalty functions.
}

\usage{
runBB(
  data,
  model = c("BB-KM", "BB-Cox", "BB-CRR","BB-GLM"),
  formula = NULL,
  family = gaussian(),
  M = 400,
  lambda = NULL,
  prior_fn = NULL,
  ...
)
}

\arguments{

\item{data}{
A data frame containing the variables required for the chosen model.
For BB-KM, must include \code{time} and \code{status}.
For BB-Cox, BB-CRR and BB-GLM, must include variables referenced in \code{formula}.
}

\item{model}{
A character string specifying the model type. Options include:
\code{"BB-KM"}, \code{"BB-Cox"}, \code{"BB-CRR"}, and \code{"BB-GLM"}.
}

\item{formula}{
An object of class \code{formula}. Required for BB-Cox and BB-GLM models.
For BB-Cox, typically of the form \code{Surv(time, status) ~ x1 + x2}.
}

\item{family}{
A description of the error distribution and link function to be used in BB-GLM.
Defaults to \code{gaussian()}.
}

\item{M}{
Number of Bayesian bootstrap draws to generate.
}

\item{lambda}{
Gibbs posterior temperature parameter controlling the influence of empirical risk. If NULL then will use 1/sqrt(M).
}

\item{prior_fn}{
A user-defined prior function applied to each draw in the Gibbs posterior.
Common choices include ridge, lasso, or smoothness priors.
See details in \link{priorBB}.
}

\item{...}{
Additional arguments passed to model-specific fitting functions.
See details in \link{Surv}, \link{coxph}, \link{crr}, \link{glm}.
}

}

\value{

A list object of class \code{BayesBoots} containing model-specific output.

\item{summary}{
Posterior summary generated from Bayesian bootstrap draws (structure depends on model).
}

\item{curves_data / beta_samples / hr_samples / or_samples}{
Posterior samples generated from Bayesian bootstrap draws (structure depends on model).
}

\item{diagnostics}{
A data frame containing:
\itemize{
  \item empirical risk values
  \item prior penalty values (\code{omega})
  \item Gibbs weights
  \item ... (depends on model)
}
}

}

\details{

The function implements a two-stage Bayesian learning procedure:

\enumerate{
\item \strong{Bayesian Bootstrap Stage}:
Dirichlet(1,...,1) weights are used to generate posterior-like resamples.

\item \strong{Gibbs Posterior Reweighting Stage}:
Each bootstrap draw is reweighted using a Gibbs posterior of the form:

\deqn{
w_m \propto \exp(-\lambda R(S_m) - \Omega(S_m))
}

where \eqn{R(S_m)} is an empirical risk function and \eqn{\Omega(S_m)} is a user-defined prior penalty.
}

This framework allows likelihood-free inference and is applicable across survival,
regression, and generalized linear models.

}

\references{

Wasserman, L. & Bhattacharya, A. (2014). \emph{Statistical decision theory and Bayesian learning}.

Rubin, D. B. (1981). The Bayesian bootstrap. \emph{Annals of Statistics}.

(Extended Gibbs posterior learning framework applied in modern likelihood-free inference literature.)
}

\examples{
\dontrun{

######################################################
## Load datasets
######################################################

library(survival)
head(lung)
ibrary(MASS)
data(Melanoma)

######################################################
## BB-KM Example (lung cancer survival)
######################################################

km_data <- na.omit(lung)

fit_km <- runBB(
  data = km_data,
  model = "BB-KM",
  formula = survival::Surv(time, status) ~ ph.ecog,
  prior_fn = prior_smoothness,
  M = 300
)

fit_km
summary(fit_km)
plot(fit_km)

######################################################
## BB-Cox Example (lung cancer covariates)
######################################################

cox_data <- na.omit(lung)

fit_cox <- runBB(
  data = cox_data,
  model = "BB-Cox",
  formula = survival::Surv(time, status) ~ age + sex + ph.ecog + wt.loss,
  prior_fn = prior_ridge,
  M = 300
)

fit_cox
summary(fit_cox)
plot(fit_cox)
plot(fit_cox, type="hr")
plot(fit_cox, type = "surv", group_var = "sex")
plot(fit_cox, type="surv", group_var = "ph.ecog")
plot(fit_cox, type="density")
plot(fit_cox, type="baseline")

######################################################
## BB-CCR Example (Melanoma dataset)
######################################################

library(MASS)
data(Melanoma)

# Quick look at the event distribution
table(Melanoma$status)
# 1 = Died from melanoma (event of interest),
# 2 = Still alive (censored),
# 3 = Died other causes (competing risk)

# for 'status' input the BB-CRR argument requires
# 0 = censored observations
# 1, 2, ... = competing event causes

dat_crr <- Melanoma
dat_crr$status_crr <- ifelse(dat_crr$status == 1, 1,
                          ifelse(dat_crr$status == 3, 2, 0))
table(dat_crr$status)
table(dat_crr$status_crr)

fit_crr <- runBB(
  data = dat_crr,
  model = "BB-CRR",
  formula = survival::Surv(time, status) ~ age + sex + thickness + ulcer,
  failcode = 1, # default is 1, here the event of interest
  prior_fn = prior_ridge,
  M = 300
)

fit_crr
summary(fit_crr)
plot(fit_crr)

######################################################
## BB-GLM Example (lung cancer data)
######################################################

## we convert the survival data for logistic model
## this is also know as the pooled logistic model

glm_data <- na.omit(lung)
glm_data$year <- glm_data$time/365.25
max_data_time <- max(glm_data$year)
glm_data <- survSplit(Surv(year, status) ~ ., data = glm_data,
                       cut = seq(0, max_data_time, by = 1),
                       episode = "time_period")

fit_glm <- runBB(
  data = glm_data,
  model = "BB-GLM",
  formula = status ~ time_period + age + sex + ph.ecog + wt.loss,
  family = binomial(),
  prior_fn = prior_ridge,
  M = 300
)

fit_glm
summary(fit_glm)
plot(fit_glm)
plot(fit_glm, type="or")
plot(fit_glm, type="density")

######################################################

}
}


\seealso{
\code{\link{priorBB}},
\code{\link{predict}}
}
