library(ggplot2)
library(dplyr)
library(ggpubr)

# Tukey Lambda quantile function
Q_tukey <- function(p, lambda) {
  if (lambda == 0) return(log(p / (1 - p)))
  if (lambda == 1) return(2 * p - 1)
  (p^lambda - (1 - p)^lambda) / lambda
}

# Setup grid
p_grid <- seq(0.001, 0.999, length.out = 1000)

# Define remaining cases (excluding λ = 0.5)
cases <- list(
  list(lambda = -1, label = "λ = -1 (Cauchy)", ref_cdf = pcauchy, range = c(-10, 10)),
  list(lambda = 0, label = "λ = 0 (Logistic)", ref_cdf = plogis, range = c(-6, 6)),
  list(lambda = 0.14, label = "λ = 0.14 (Normal)", ref_cdf = pnorm, range = c(-4, 4)),
  list(lambda = 1, label = "λ = 1 (Uniform)", ref_cdf = function(x) punif(x, min = -1, max = 1), range = c(-2, 2))
)

# Create one plot per case
plots <- lapply(cases, function(case) {
  lam <- case$lambda
  label <- case$label
  ref_cdf <- case$ref_cdf
  xlim <- case$range
  
  # Tukey CDF
  x_tukey <- Q_tukey(p_grid, lam)
  df_tukey <- data.frame(x = x_tukey, cdf = p_grid, dist = "Tukey Lambda")
  
  # Reference CDF
  x_ref <- seq(min(x_tukey), max(x_tukey), length.out = 1000)
  df_ref <- data.frame(x = x_ref, cdf = ref_cdf(x_ref), dist = "Reference")
  
  # Combine
  df_all <- bind_rows(df_tukey, df_ref)
  
  # Plot
  ggplot(df_all, aes(x = x, y = cdf, color = dist, linetype = dist)) +
    geom_line(size = 1.1) +
    scale_color_manual(values = c("Tukey Lambda" = "steelblue", "Reference" = "darkred")) +
    scale_linetype_manual(values = c("Tukey Lambda" = "solid", "Reference" = "dashed")) +
    labs(title = label, x = "x", y = "CDF") +
    coord_cartesian(xlim = xlim) +
    theme_minimal() +
    theme(legend.position = "none")
})

# Arrange plots in a 2×2 grid
ggarrange(plotlist = plots, ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "bottom")
