library(ggplot2)
library(ggpubr)
library(readr)
library(grid)

# Read data
setwd("~/Documents/Projects/GLD_inference/GLD_mu_sigma")

results_mu_df <- read_csv("~/Documents/Projects/GLD_inference/GLD_mu_sigma/gld_mu_ci.csv")
results_sigma_df <- read_csv("~/Documents/Projects/GLD_inference/GLD_mu_sigma/gld_sigma_ci.csv")

## --------------------------- MU PLOTS --------------------------- ##
plot_cov_mu <- ggplot(results_mu_df, aes(x = n, y = coverage, color = factor(mu))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", size = 1) +
  ylim(0.5, 1) +
  labs(
    title = expression("Coverage vs Sample Size (" * mu * ")"),
    x = "Sample Size (n)",
    y = "Coverage",
    color = expression(mu)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray85")
  ) 

plot_width_mu <- ggplot(results_mu_df, aes(x = n, y = width, color = factor(mu))) +
  geom_line(size = 1) +
  labs(
    title = expression("Width vs Sample Size (" * mu * ")"),
    x = "Sample Size (n)",
    y = "Interval Width",
    color = expression(mu)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray85")
  ) 

mu_combined <- ggarrange(plot_cov_mu, plot_width_mu, ncol = 2, common.legend = TRUE, legend = "bottom")

mu_combined
ggsave( "gld_mu_ci.pdf", width = 11, height = 5, dpi = 300, units = "in" )

## --------------------------- SIGMA PLOTS --------------------------- ##
plot_cov_sigma <- ggplot(results_sigma_df, aes(x = n, y = coverage, color = factor(sigma))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", size = 1) +
  ylim(0.5, 1) +
  labs(
    title = expression("Coverage vs Sample Size (" * sigma * ")"),
    x = "Sample Size (n)",
    y = "Coverage",
    color = expression(sigma)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray85")
  ) 

plot_width_sigma <- ggplot(results_sigma_df, aes(x = n, y = width, color = factor(sigma))) +
  geom_line(size = 1) +
  labs(
    title = expression("Width vs Sample Size (" * sigma * ")"),
    x = "Sample Size (n)",
    y = "Interval Width",
    color = expression(sigma)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray85")
  ) 

sigma_combined <- ggarrange(plot_cov_sigma, plot_width_sigma, ncol = 2, common.legend = TRUE, legend = "bottom")

ggsave( "gld_sigma_ci.pdf", width = 11, height = 5, dpi = 300, units = "in" )

