setwd("adam/")

cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442","#0072B2","#D55E00","#CC79A7", "#999999")

#' Start by installing the aciccomp2016 R package from Github
if (require("remotes", quietly = TRUE) == FALSE) {
  install.packages("remotes")
  require("remotes")
}

remotes::install_github("vdorie/aciccomp/2016")

library(aciccomp2016)
library(tidyverse)

#' The 87th dataset simulated from the 50th setting
df <- dgp_2016(input_2016, 50, 87) %>%
  as_tibble()

#' The columns of df are:
#' * z: the treatment indicator
#' * y: the response variable
#' * y.0, y.1: the potential outcomes
#' * mu.0, mu.1: the expected potential outcomes
#' * e: the true propensity score
names(df)

#' Want to the SATT
#' * SATT: sample average treatment effect on treated
#' * SATC: sample average treatement effect on control
satt_truth <- df %>%
  group_by(z) %>%
  summarise(sat = mean(y.1 - y.0)) %>%
  filter(z == 0) %>%
  select(sat) %>%
  as.numeric()

#' A plot to look at confounding
#' e.g. is the distribution of potential outcomes different depending on assignment
pdf("87-50.pdf", h = 3, w = 6.25)

df %>%
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

#' 1. Regression adjustment
fit <- glm(y ~ 1 + z, family = "gaussian", data = df)

#' Only the treated
df_treated <- df %>% filter(z == 1)

ey0 <- mean(predict(fit, df_treated %>% mutate(z = 0), type = "response"))
ey1 <- mean(predict(fit, df_treated, type = "response"))

list(
  "regression_adjustment" = ey1 - ey0,
  "truth" = satt_truth
)
