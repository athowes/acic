setwd("adam/")

cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442","#0072B2","#D55E00","#CC79A7", "#999999")

#' #' Start by installing the aciccomp2016 R package from Github
#' if (require("remotes", quietly = TRUE) == FALSE) {
#'   install.packages("remotes")
#'   require("remotes")
#' }
#'
#' remotes::install_github("vdorie/aciccomp/2016")

library(aciccomp2016)
library(tidyverse)
library(broom)
library(SuperLearner)

#' The 87th dataset simulated from the 50th setting
df_full <- dgp_2016(input_2016, 50, 87) %>%
  as_tibble() %>%
  bind_cols(input_2016)

#' Save a dataset for users who are allergic to R
write_csv(df_full, "87-50.csv")

#' The columns of df are:
#' * z: the treatment indicator
#' * y: the response variable
#' * y.0, y.1: the potential outcomes
#' * mu.0, mu.1: the expected potential outcomes
#' * e: the true propensity score
names(df_full)

#' Want to predict the SATT
#' * SATT: sample average treatment effect on treated
#' * SATC: sample average treatment effect on control
satt_truth <- df_full %>%
  filter(z == 0) %>%
  summarise(sat = mean(y.1 - y.0)) %>%
  as.numeric()

#' Save the true propensity scores in-case want to compare to them later
e_truth <- select(df_full, e)

#' A plot to look at confounding
#' e.g. is the distribution of potential outcomes different depending on assignment
pdf("87-50.pdf", h = 3, w = 6.25)

df_full %>%
  pivot_longer(
    cols = c("y.0", "y.1"),
    names_prefix = "y.",
    names_to = "potential_outcome",
    values_to = "value"
  ) %>%
  mutate(
    observed = ifelse(potential_outcome == z, TRUE, FALSE)
  ) %>%
  ggplot(aes(x = potential_outcome, y = value, col = observed)) +
    geom_jitter(alpha = 0.5) +
    scale_colour_manual(values = c("#D3D3D3", cbpalette[3])) +
    labs(x = "Control or treatment", y = "Potential outcome", col = "Observed?") +
    theme_minimal()

dev.off()

#' Remove the columns that I shouldn't have access to
df <- select(df_full, -y.0, -y.1, -mu.0, -mu.1, -e)

#' 1. Regression adjustment
fit <- glm(y ~ ., family = "gaussian", data = df)

compare_to_truth <- function(df, fit, name) {
  df_treated <- filter(df, z == 1)
  ey0 <- mean(predict(fit, mutate(df_treated, z = 0), type = "response"))
  ey1 <- mean(predict(fit, df_treated, type = "response"))
  out <- list(ey1 - ey0, satt_truth)
  names(out) <- c(paste(name), "truth")
  return(out)
}

compare_to_truth(df, fit, name = "regression")

#' 2. Regression adjustment with propensity scores

fit_treatment <- glm(z ~ 1 + ., family = binomial(link = "logit"), data = select(df, -y))

df_e <- augment_columns(fit_treatment, df, type.predict = "response") %>%
  #' e_fitted is the propensity
  rename(e_fitted = .fitted) %>%
  select("z", "y", starts_with("x"), "e_fitted")

fit_e <- glm(y ~ ., family = "gaussian", data = df_e)

compare_to_truth(df_e, fit_e, name = "regression_propensity")

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

df_ipw <- df_e %>%
  mutate(ipw = (z / e_fitted) + ((1 - z) / (1 - e_fitted)))

#' Don't include ipw or e_fitted in the regression equation
fit_ipw <- glm(y ~ ., family = "gaussian", data = select(df_ipw, -ipw, -e_fitted), weights = df_ipw$ipw)

compare_to_truth(df_ipw, fit_ipw, name = "ipw")

#' 3. Inverse probability weighting (with machine learning methods)
#' TODO
