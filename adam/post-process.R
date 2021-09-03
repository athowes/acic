library(ggplot2)

df <- list.files("results/", full.names = TRUE) %>%
  map_df(~read_csv(.)) %>%
  mutate(abs_error = abs(error))

cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442","#0072B2","#D55E00","#CC79A7", "#999999")

fct_reorg <- function(fac, ...) {
  fct_recode(fct_relevel(fac, ...), ...)
}

df <- df %>%
  mutate(
    method = fct_reorg(method,
      "Regression adjustment" = "regression",
      "Regression with propensity" = "regression_propensity",
      "Inverse probability weights with logistic regression" = "ipw",
      "Inverse probability weights with SuperLearner" = "superlearner",
      "Regression with propensity (oracle)" = "regression_oracle_propensity",
      "Inverse probability weights (oracle)" = "oracle_ipw"
    )
  )

pdf("compare-methods.pdf", h = 5, w = 9)

ggplot(df, aes(x = method, y = abs_error, fill = method)) +
  geom_boxplot() +
  facet_grid(~param) +
  scale_fill_manual(values = cbpalette) +
  labs(
    title = "ACIC results",
    subtitle = "3 settings, 5 seeds",
    y = "abs(error)",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "lines")
  )

dev.off()
