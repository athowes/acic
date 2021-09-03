#' setwd("adam/")

#' cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442","#0072B2","#D55E00","#CC79A7", "#999999")

# param <- 50
# seed <- 87

param <- commandArgs()$param
seed <- commandArgs()$seed

df_full <- dgp_2016(input_2016, param, seed) %>%
  as_tibble() %>%
  bind_cols(input_2016)

#' #' Save a dataset for users who are allergic to R
#' write_csv(df_full, "87-50.csv")

#' #' The columns of df are:
#' #' * z: the treatment indicator
#' #' * y: the response variable
#' #' * y.0, y.1: the potential outcomes
#' #' * mu.0, mu.1: the expected potential outcomes
#' #' * e: the true propensity score
#' names(df_full)

#' Want to predict the SATT
#' * SATT: sample average treatment effect on treated
#' * SATC: sample average treatment effect on control
satt_truth <- df_full %>%
  filter(z == 1) %>%
  summarise(sat = mean(y.1 - y.0)) %>%
  as.numeric()

#' Save the true propensity scores in-case want to compare to them later
e_truth <- select(df_full, e)

#' Not needed for now
#' #' A plot to look at confounding
#' #' e.g. is the distribution of potential outcomes different depending on assignment
#' pdf("87-50.pdf", h = 3, w = 6.25)
#'
#' df_full %>%
#'   pivot_longer(
#'     cols = c("y.0", "y.1"),
#'     names_prefix = "y.",
#'     names_to = "potential_outcome",
#'     values_to = "value"
#'   ) %>%
#'   mutate(
#'     observed = ifelse(potential_outcome == z, TRUE, FALSE)
#'   ) %>%
#'   ggplot(aes(x = potential_outcome, y = value, col = observed)) +
#'     geom_jitter(alpha = 0.5) +
#'     scale_colour_manual(values = c("#D3D3D3", cbpalette[3])) +
#'     labs(x = "Control or treatment", y = "Potential outcome", col = "Observed?") +
#'     theme_minimal()
#'
#' dev.off()

#' Remove the columns that I shouldn't have access to
df <- select(df_full, -y.0, -y.1, -mu.0, -mu.1, -e)

# Calculate the SATT
calculate_satt <- function(df, fit, name) {
  df_treated <- filter(df, z == 1)
  ey0 <- mean(predict(fit, mutate(df_treated, z = 0), type = "response"))
  ey1 <- mean(predict(fit, df_treated, type = "response"))
  out <- list(ey1 - ey0)
  names(out) <- c(paste(name))
  return(out)
}

#' 1. Regression adjustment
fit <- glm(y ~ ., family = "gaussian", data = df)

(result1 <- calculate_satt(df, fit, name = "regression"))

#' 2. Regression adjustment with propensity scores

fit_treatment <- glm(z ~ ., family = binomial(link = "logit"), data = select(df, -y))

#' e_fitted is the propensity
df_e <- mutate(df, e_fitted = predict(fit_treatment, df, type = "response"))

fit_e <- glm(y ~ 1 + z + e_fitted, family = "gaussian", data = df_e)

(result2 <- calculate_satt(df_e, fit_e, name = "regression_propensity"))

#' 2a. Regression adjustment with oracle propensity scores

df_oracle_e <- mutate(df, e_truth = unlist(e_truth))

fit_oracle_e <- glm(y ~ 1 + z + e_truth, family = "gaussian", data = df_oracle_e)

(result2a <- calculate_satt(df_oracle_e, fit_oracle_e, name = "regression_oracle_propensity"))

#' 3. Inverse probability weighting

#' 1) Model the treatment / control exposure to obtain propensity scores
#'  * Want to predict if someone is going to be treated or not. This helps us understand
#'    if there are differences between the control and treatment groups
#'  * Want to look at the distribution of the scores
#'    * Is there enough overlap?
#'    * Very close to 0 or 1?
#' 2) Convert propensity scores to inverse probability weights
#' 3) Assess covariate balance in weighted sample
#'   * Want the treatment and control groups to have similar distributions
#'     over the covariates
#'   * There is an R package for this "cobalt: Covariate Balance Tables and Plots"
#' 4) Model the outcome

#' Use the treatment model from 2.
df_ipw <- df_e %>%
  mutate(ipw = (z / e_fitted) + ((1 - z) / (1 - e_fitted)))

#' Don't include ipw or e_fitted in the regression equation
fit_ipw <- glm(y ~ ., family = "gaussian", data = select(df_ipw, -ipw, -e_fitted), weights = df_ipw$ipw)

(result3 <- calculate_satt(df_ipw, fit_ipw, name = "ipw"))

#' 3a. Inverse probability weighting with oracle propensity scores

df_oracle_ipw <- df_oracle_e %>%
  mutate(ipw = (z / e_truth) + ((1 - z) / (1 - e_truth)))

#' Don't include ipw or e_fitted in the regression equation
fit_oracle_ipw <- glm(
  y ~ .,
  family = "gaussian",
  data = select(df_oracle_ipw, -ipw, -e_truth),
  weights = df_oracle_ipw$ipw
)

(result3a <- calculate_satt(df_oracle_ipw, fit_oracle_ipw, name = "oracle_ipw"))

#' 4. Inverse probability weighting (with machine learning methods)

#' Note: could use ML for the treatment model, response model, or both
#' That said, perhaps it's difficult to find ML methods which work with
#' weighted likelihoods for the response model (doing IPW)

#' #' Prepare packages required for SuperLearner
#' #' eXtreme Gradient Boosting
#' install.packages(
#'   "xgboost",
#'   repos = c("http://dmlc.ml/drat/", getOption("repos")),
#'   type = "source"
#' )

fit_treatment_superlearner <- SuperLearner(
  Y = df$z,
  X = select(df, -z, -y),
  cvControl = list(V = 3),
  SL.library = c("SL.glm", "SL.glmnet", "SL.xgboost"),
  method = "method.NNLS",
  family = "binomial"
)

df_superlearner <- df %>%
  mutate(
    e_fitted = predict(fit_treatment_superlearner, type = "response")$pred,
    ipw = (z / e_fitted) + ((1 - z) / (1 - e_fitted))
  )

#' Don't include ipw or e_fitted in the regression equation
fit_superlearner <- glm(
  y ~ .,
  family = "gaussian",
  data = select(df_superlearner, -ipw, -e_fitted),
  weights = df_superlearner$ipw
)

(result4 <- calculate_satt(df_superlearner, fit_superlearner, name = "superlearner"))

#' Save to results

dir.create("results")

bind_cols(result1, result2, result2a, result3, result3a, result4) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("method") %>%
  rename("satt" = "V1") %>%
  mutate(
    truth = satt_truth,
    error = satt - truth,
    seed = seed,
    param = param
  ) %>%
  write_csv(paste0("results/result-", seed, "-", param, ".csv"))
