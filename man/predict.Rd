\name{predict}
\alias{predict}

\title{Posterior Prediction for Bayesian Bootstrap Models}

\description{
Generates posterior predictions from fitted Bayesian Bootstrap models produced by
\code{\link{runBB}}. Predictions are available for Bayesian Bootstrap Cox
proportional hazards models (\code{"BB-Cox"}) and Bayesian Bootstrap generalized
linear models (\code{"BB-GLM"}).

For Cox models, posterior survival curves are computed by combining posterior
samples of regression coefficients with posterior baseline survival estimates.

For GLM models, posterior predictive probabilities and odds ratios are obtained
from the posterior distribution of the regression coefficients.
}

\usage{
\method{predict}{BayesBoots}(object, newdata, times = NULL, ...)
}

\arguments{

  \item{object}{
    A fitted Bayesian Bootstrap model returned by \code{\link{runBB}}.
  }

  \item{newdata}{
    A data frame containing covariate values for prediction.
    Variable names must match those used when fitting the model.
  }

  \item{times}{
    Optional vector of time points at which survival probabilities should be
    evaluated for \code{"BB-Cox"} models. If \code{NULL}, the original time grid
    used during model fitting is employed.
  }

  \item{\dots}{
    Additional arguments (currently ignored).
  }

}

\details{

The prediction output depends on the fitted model type.

\strong{BB-Cox}

For Bayesian Bootstrap Cox models, posterior survival curves are calculated as

\deqn{
S(t \mid x) = S_0(t)^{\exp(x^\top \beta)}
}

where:

\itemize{
  \item \eqn{S_0(t)} is the posterior baseline survival function.
  \item \eqn{\beta} is a posterior draw of the regression coefficients.
  \item \eqn{x} is the covariate vector for a new observation.
}

Posterior summaries returned include:

\itemize{
  \item \code{survival_median}
  \item \code{survival_lower} (2.5\% posterior quantile)
  \item \code{survival_upper} (97.5\% posterior quantile)
  \item \code{linear_predictor}
}

\strong{BB-GLM}

For Bayesian Bootstrap generalized linear models, posterior predictions are
computed from the posterior coefficient distribution.

For logistic models, posterior probabilities are obtained via

\deqn{
p = \frac{1}{1+\exp(-x^\top\beta)}
}

Posterior odds ratios are computed as

\deqn{
OR = \exp(x^\top\beta)
}

Posterior summaries returned include:

\itemize{
  \item \code{probability_median}
  \item \code{probability_lower} (2.5\% posterior quantile)
  \item \code{probability_upper} (97.5\% posterior quantile)
  \item \code{odds_ratio_median}
  \item \code{odds_ratio_lower}
  \item \code{odds_ratio_upper}
  \item \code{linear_predictor}
}

The returned prediction object has class

\preformatted{
c("predBB", "BayesBoots", model)
}

where \code{model} is either \code{"BB-Cox"} or \code{"BB-GLM"}.
}

\value{

An object of class \code{predBB}.

For \code{"BB-Cox"} models:

\itemize{
  \item \code{model} : model identifier.
  \item \code{time} : prediction time grid.
  \item \code{survival_median} : posterior median survival probabilities.
  \item \code{survival_lower} : lower posterior credible limit.
  \item \code{survival_upper} : upper posterior credible limit.
  \item \code{linear_predictor} : posterior linear predictors.
}

For \code{"BB-GLM"} models:

\itemize{
  \item \code{model} : model identifier.
  \item \code{summary} : tibble containing posterior summaries of probabilities and odds ratios.
  \item \code{linear_predictor} : posterior linear predictors.
}

}

\examples{
\dontrun{

## Cox model prediction

fit_cox <- runBB(
  data = lung,
  model = "BB-Cox",
  formula = Surv(time, status) ~ age + sex,
  prior_fn = prior_ridge
)

pred_cox <- predict(
  fit_cox,
  newdata = lung[1:5, ]
)

plot(pred_cox, ids = 1:5)

## GLM prediction

fit_glm <- runBB(
  data = mydata,
  model = "BB-GLM",
  formula = y ~ x1 + x2,
  family = binomial(),
  prior_fn = prior_ridge
)

pred_glm <- predict(
  fit_glm,
  newdata = mydata
)

plot(pred_glm)

}
}

\seealso{
\code{\link{runBB}},
\code{\link{priorBB}}
}
