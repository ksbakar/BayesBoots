\name{priorBB}

\title{Collection of Prior Functions for Bayesian Bootstrap (BB) Survival Models}

\description{
This set of functions defines a variety of prior distributions and regularization penalties used in Bayesian regression models and survival analysis. The priors include standard shrinkage priors for regression coefficients as well as smoothness and shape constraints for survival curves such as Kaplan--Meier estimates.
}

\usage{
priorBB(prior_fn)
}

\arguments{
  \item{prior_fn}{An R function to define prior distribution. See below for details.}
}

\details{

The possible prior  built-in functions are:

\itemize{
  \item \code{prior_ridge}: L2 penalty (Gaussian shrinkage).
  \item \code{prior_lasso}: L1 penalty (Laplace shrinkage).
  \item \code{prior_elastic}: Combination of L1 and L2 penalties.
  \item \code{prior_normal}: Gaussian prior penalty on coefficients.
  \item \code{prior_mvn}: Component-wise Gaussian prior (diagonal covariance).
  \item \code{prior_student}: Heavy-tailed Student-t prior penalty.
  \item \code{prior_cauchy}: Very heavy-tailed Cauchy prior penalty.
  \item \code{prior_horseshoe}: Horseshoe-style global shrinkage prior.
  \item \code{prior_spike_slab}: Mixture prior combining sparse and dense components.
  \item \code{prior_smoothness}: Penalizes curvature in survival curves (second differences).
  \item \code{prior_roughness}: Penalizes total variation in survival curves.
  \item \code{prior_hazard_stability}: Penalizes variability in implied hazard rates.
  \item \code{prior_early_failure}: Encourages higher weight on early failures in survival curves.
}

Corresponding functional forms:

\itemize{
\item \code{prior_ridge(beta, sigma2 = 1)}
\item \code{prior_lasso(beta, tau = 1)}
\item \code{prior_elastic(beta, alpha = 0.5, tau = 1)}
\item \code{prior_normal(beta, mu = 0, sigma2 = 1)}
\item \code{prior_mvn(beta, mu = NULL, sigma2 = 1)}
\item \code{prior_student(beta, mu = 0, sigma2 = 1, df = 3)}
\item \code{prior_cauchy(beta, mu = 0, scale = 1)}
\item \code{prior_horseshoe(beta, tau = 1, eps = 1e-8)}
\item \code{prior_spike_slab(beta, spike_var = 0.01, slab_var = 1, pi = 0.5)}
\item \code{prior_smoothness(curve)}
\item \code{prior_roughness(curve))}
\item \code{prior_hazard_stability(curve))}
\item \code{prior_early_failure(curve))}
}

Where the arguments:

\itemize{
  \item \code{beta}: Numeric vector of regression coefficients.
  \item \code{sigma2}: Variance parameter for ridge and normal-type priors.
  \item \code{tau}: Regularization parameter controlling shrinkage strength.
  \item \code{alpha}: Mixing parameter for elastic net prior, balancing L1 and L2 penalties.
  \item \code{mu}: Mean vector or scalar location parameter. If NULL in \code{prior_mvn}, defaults to a zero vector.
  \item \code{df}: Degrees of freedom for the Student-t prior. Controls tail heaviness.
  \item \code{scale}: Scale parameter for the Cauchy prior.
  \item \code{eps}: Small constant to avoid division by zero in horseshoe-style prior.
  \item \code{spike_var}: Variance of the spike component in spike-and-slab prior.
  \item \code{slab_var}: Variance of the slab component in spike-and-slab prior.
  \item \code{pi}: Mixing probability for the spike component.
  \item \code{curve}: A survival curve object containing at least \code{surv} (survival probabilities), and for some functions \code{time} (time points).
}

}

\examples{
\dontrun{

beta <- rnorm(10)

prior_ridge(beta, sigma2 = 2)

prior_lasso(beta, tau = 0.5)

prior_elastic(beta, alpha = 0.7, tau = 1)

prior_normal(beta, mu = 0, sigma2 = 1)

prior_student(beta, df = 5)

curve <- list(
  surv = c(1, 0.9, 0.7, 0.5, 0.3),
  time = c(0, 1, 2, 3, 4)
)

prior_smoothness(curve)
prior_roughness(curve)
prior_hazard_stability(curve)
prior_early_failure(curve)
}
}

\seealso{
\code{\link{runBB}}
}
