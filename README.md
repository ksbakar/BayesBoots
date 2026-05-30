
# BayesBoots

`BayesBoots` is an R package for Bayesian bootstrap methods and related statistical inference tools.

## Installation

The development version of `BayesBoots` can be installed directly from GitHub using the `devtools` package.

### Step 1: Install `devtools`

If you do not already have `devtools` installed, run:

```r
install.packages("devtools")
```

### Step 2: Install `BayesBoots` from GitHub

```r
devtools::install_github("ksbakar/BayesBoots")
```

### Step 3: Load the Package

```r
library(BayesBoots)
```

## Example

```r
library(BayesBoots)

# Example code

######################################################
## Load survival package datasets
######################################################

library(survival)
head(lung)

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
plot(fit_cox, type="density")
plot(fit_cox, type="baseline")


######################################################
## BB-GLM Example (Primary Biliary Cirrhosis - PBC)
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


```

## Reporting Issues

If you encounter bugs or have feature requests, please open an issue on the GitHub repository.

## License

Please see the `LICENSE` file for details.

