results_mu_df = read.csv("~/Documents/Projects/GLD_inference/GLD_mu_sigma/gld_mu_ci.csv")


ggplot(results_mu_df, aes(x = n, y = coverage, color = factor(mu))) +
  geom_line(size = 1) +
  labs(
    title = "Average Coverage vs Sample Size",
    x = "Sample Size (n)",
    y = "Coverage",
    color = "Mu"
  ) +
  ylim(0.5, 1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", size = 1) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title.x = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 16),  # Larger x-axis tick labels
    axis.text.y = element_text(size = 16),  # Larger y-axis tick labels
    axis.ticks.length = unit(0.3, "cm"),   # Longer axis ticks
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1.5, "lines")
  )


ggplot(results_mu_df, aes(x = n, y = width, color = factor(mu))) +
  geom_line(size = 1) +
  labs(
    title = "Average Width vs Sample Size",
    x = "Sample Size (n)",
    y = "Interval Width",
    color = "Mu"
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title.x = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 16),  # Larger x-axis tick labels
    axis.text.y = element_text(size = 16),  # Larger y-axis tick labels
    axis.ticks.length = unit(0.3, "cm"),   # Longer axis ticks
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1.5, "lines")
  )


##########################################################################################


results_sigma_df = read.csv("~/Documents/Projects/GLD_inference/GLD_mu_sigma/gld_sigma_ci.csv")
ggplot(results_sigma_df, aes(x = n, y = coverage, color = factor(sigma))) +
  geom_line(size = 1) +
  labs(
    title = "Average Coverage vs Sample Size",
    x = "Sample Size (n)",
    y = "Coverage",
    color = "Sigma"
  ) +
  ylim(0.5, 1) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", size = 1) +
  theme_bw()+
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title.x = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 16),  # Larger x-axis tick labels
    axis.text.y = element_text(size = 16),  # Larger y-axis tick labels
    axis.ticks.length = unit(0.3, "cm"),   # Longer axis ticks
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1.5, "lines")
  )

ggplot(results_sigma_df, aes(x = n, y = width, color = factor(sigma))) +
  geom_line(size = 1) +
  labs(
    title = "Average Width vs Sample Size",
    x = "Sample Size (n)",
    y = "Interval Width",
    color = "Sigma"
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(face = "bold", size = 24),
    axis.title.x = element_text(face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 22),
    axis.text.x = element_text(size = 16),  # Larger x-axis tick labels
    axis.text.y = element_text(size = 16),  # Larger y-axis tick labels
    axis.ticks.length = unit(0.3, "cm"),   # Longer axis ticks
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(1.5, "lines")
  )

