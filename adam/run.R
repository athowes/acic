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

for(p in c(30, 40, 50)) {
  for(s in 83:87) {
    commandArgs <- function(...) list(param = p, seed = s)
    source("script.R")
  }
}
