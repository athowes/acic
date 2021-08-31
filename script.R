##' Start by installing the aciccomp2016 R package from Github
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

#' Want to compute the sample average treatment effect on the treated
#' Not aiming to compute the population ATE
df %>%
  group_by(z) %>%
  summarise(ate = mean(y.1 - y.0))
