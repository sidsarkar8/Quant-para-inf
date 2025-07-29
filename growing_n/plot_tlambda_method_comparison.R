library(ggplot2)
library(ggforce)
library(dplyr)


all_results_long = read.csv("/Users/sidsarkar/Documents/Projects/GLD_inference/growing_n/tlambda_method_comparison.csv")

all_results_long$source[all_results_long$source == "DW"] = "CDF-CI"
all_results_long$source[all_results_long$source == "quant(NP boot)"] = "NP boot + quant"
all_results_long$source[all_results_long$source == "quant(P boot)"] = "P boot + quant"
all_results_long$source[all_results_long$source == "lmom(NP boot)"] = "NP boot + L-m"
all_results_long$source[all_results_long$source == "lmom(P boot)"] = "P boot + L-m"


ggplot(all_results_long, aes(x = lambda, y = value, color = source)) +
  geom_line() +
  geom_point() +
  geom_hline(
    data = data.frame(metric = "coverage", yint = 0.95),
    aes(yintercept = yint),
    linetype = "dotted",
    color = "black",
    linewidth = 0.8
  ) +
  facet_grid(rows = vars(metric), cols = vars(n), scales = "free_y") +
  labs(
    title = "Comparison of methods",
    x = expression(lambda),
    y = "Value",
    color = "Method"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "gray85")
  )+
  facetted_pos_scales(
    y = list(
      metric == "coverage" ~ scale_y_continuous(limits = c(0.5, 1)),
      metric != "coverage" ~ scale_y_continuous()  
    )
  )


#########################################
# Add a highlight variable
all_results_long_dw <- all_results_long %>%
  mutate(highlight = ifelse(source == "CDF-CI", "CDF-CI", "Other"))



ggplot(all_results_long_dw, aes(x = lambda, y = value)) + 
  geom_line(aes(color = highlight, group = source, linewidth = highlight), alpha = 0.9) + 
  geom_point(aes(color = highlight), alpha = 0.9) + 
  geom_hline(
    data = data.frame(metric = "coverage", yint = 0.95),
    aes(yintercept = yint),
    linetype = "dotted",
    color = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_grid(rows = vars(metric), cols = vars(n), scales = "free_y") + 
  labs(
    title = "Comparison of methods",
    x = expression(lambda),
    y = "Value",
    color = "Method"
  ) +
  scale_color_manual(
    values = c("CDF-CI" = "#F8766D", "Other" = "gray70")
  ) +
  scale_linewidth_manual(
    values = c("CDF-CI" = 1.2, "Other" = 0.6), guide = "none"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "gray85")
  ) +
  facetted_pos_scales(
    y = list(
      metric == "coverage" ~ scale_y_continuous(limits = c(0.5, 1)),
      metric != "coverage" ~ scale_y_continuous()
    )
  )

#############################


param_methods <- c("P boot + quant", "P boot + L-m")

res_para <- all_results_long %>%
  mutate(highlight = ifelse(source %in% param_methods, "Parametric", "Other"))

ggplot(res_para, aes(x = lambda, y = value)) + 
  geom_line(aes(color = highlight, group = source, linewidth = highlight), alpha = 0.9) + 
  geom_point(aes(color = highlight), alpha = 0.9) + 
  geom_hline(
    data = data.frame(metric = "coverage", yint = 0.95),
    aes(yintercept = yint),
    linetype = "dotted",
    color = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_grid(rows = vars(metric), cols = vars(n), scales = "free_y") + 
  labs(
    title = "Comparison of methods",
    x = expression(lambda),
    y = "Value",
    color = "Method"
  ) +
  scale_color_manual(
    values = c("Parametric" = "orange", "Other" = "gray70")
  ) +
  scale_linewidth_manual(
    values = c("Parametric" = 1.2, "Other" = 0.6),
    guide = "none"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "gray85")
  ) +
  facetted_pos_scales(
    y = list(
      metric == "coverage" ~ scale_y_continuous(limits = c(0.5, 1)),
      metric != "coverage" ~ scale_y_continuous()
    )
  )


