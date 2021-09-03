#' 3 parameter settings, and 5 datasets from each
dir.create("data")

for(param in c(30, 40, 50)) {
  for(seed in 83:87) {
    dgp_2016(input_2016, param, seed) %>%
      as_tibble() %>%
      bind_cols(input_2016) %>%
      write_csv(paste0("data/", seed, "-", param, ".csv"))
  }
}

