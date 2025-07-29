library(ggplot2)
library(ggforce)
library(dplyr)

setwd("/Users/sidsarkar/Documents/Projects/GLD_inference/growing_n")

all_results_long = read.csv("/Users/sidsarkar/Documents/Projects/GLD_inference/growing_n/tlambda_method_comparison.csv")

all_results_long$source[all_results_long$source == "DW"] = "CDF-CI"
all_results_long$source[all_results_long$source == "quant(NP boot)"] = "NP boot + quant"
all_results_long$source[all_results_long$source == "quant(P boot)"] = "P boot + quant"
all_results_long$source[all_results_long$source == "lmom(NP boot)"] = "NP boot + L-m"
all_results_long$source[all_results_long$source == "lmom(P boot)"] = "P boot + L-m"


desired_order <- c("n = 50", "n = 100", "n = 250", "n = 1000")
all_results_long$n <- factor(all_results_long$n, levels = desired_order)


all_results_long$metric[all_results_long$metric == "coverage"] = "Coverage"
all_results_long$metric[all_results_long$metric == "width"] = "Width"

res_np = all_results_long[all_results_long$source %in% c("CDF-CI","NP boot + quant","NP boot + L-m"),]
res_np$source[res_np$source == "NP boot + quant"] = "Quantile" 
res_np$source[res_np$source == "NP boot + L-m"] = "L-moment" 

res_p = all_results_long[!(all_results_long$source %in% c("NP boot + quant","NP boot + L-m")),]
res_p$source[res_p$source == "P boot + quant"] = "Quantile" 
res_p$source[res_p$source == "P boot + L-m"] = "L-moment" 



#################################################################

p_non_para = ggplot(res_np, aes(x = lambda, y = value, color = source)) + 
  geom_line(linewidth = 1) + 
  geom_point() + 
  geom_hline(
    data = data.frame(metric = "Coverage", yint = 0.95),
    aes(yintercept = yint),
    linetype = "dotted",
    color = "black",
    linewidth = 1
  ) +
  facet_grid(rows = vars(metric), cols = vars(n), scales = "free_y") + 
  labs(
    #title = "Comparison of methods: Non parametric bootstrap",
    x = expression(lambda),
    y = "Value",
    color = "Method"
  ) +
  scale_color_manual(
    values = c(
      "CDF-CI" = "orange",
      "Quantile" = "#1f78b4",  # Blue
      "L-moment" = "#33a02c"     # Green
    )
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
  ) +
  facetted_pos_scales(
    y = list(
      metric == "Coverage" ~ scale_y_continuous(limits = c(0.5, 1)),
      metric != "Coverage" ~ scale_y_continuous()
    )
  ) 



p_non_para

ggsave( "tl_non_para.pdf", width = 10, height = 6, dpi = 300, units = "in" )


#################################################################################

p_para = ggplot(res_p, aes(x = lambda, y = value, color = source)) + 
  geom_line(linewidth = 1) + 
  geom_point() + 
  geom_hline(
    data = data.frame(metric = "Coverage", yint = 0.95),
    aes(yintercept = yint),
    linetype = "dotted",
    color = "black",
    linewidth = 1
  ) +
  facet_grid(rows = vars(metric), cols = vars(n), scales = "free_y") + 
  labs(
    #title = "Comparison of methods: Parametric bootstrap",
    x = expression(lambda),
    y = "Value",
    color = "Method"
  ) +
  scale_color_manual(
    values = c(
      "CDF-CI" = "orange",
      "Quantile" = "#1f78b4",  # Blue
      "L-moment" = "#33a02c"     # Green
    )
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
  ) +
  facetted_pos_scales(
    y = list(
      metric == "Coverage" ~ scale_y_continuous(limits = c(0.5, 1)),
      metric != "Coverage" ~ scale_y_continuous()
    )
  ) 

p_para
ggsave( "tl_para.pdf", width = 10, height = 6, dpi = 300, units = "in" )
