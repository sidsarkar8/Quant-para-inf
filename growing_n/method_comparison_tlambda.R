## Importing all the required libraries 
library(extraDistr)
library(dplyr)
library(ggplot2)
library(Lmoments)
library(ggh4x)

#----------------------------------------------------
## Including Source Files
source("Documents/Academic_Stuffs/Projects/AK/R codes/DW_DKW_CI_F.R")
source("Documents/Academic_Stuffs/Projects/AK/R codes/DW_DKW_CI_comparison_lambda.R")

lambda_values =  c(-1,-0.5,0,0.5,1)
n_values = c(50,100,250,1000)
#---------------------------------------------------------------
## 1st competing method : estimating lambda by quantile comparison
## Then applying parametric and nonparametric bootstrap for CI, their width and coverage
lambda_est2 = function(samp){
  act_q = c(0.25, 0.75)
  target_q = quantile(samp, act_q)
  lambda_vals = numeric(length(act_q))
  
  for (i in seq_along(act_q)) {
    if (act_q[i] != 0.5) {
      f = function(lambda) q(act_q[i], lambda) - target_q[i]
      res = tryCatch(uniroot(f, lower = -100, upper = 100, tol = 1e-10)$root, error = function(e) NA)
      lambda_vals[i] = res
    }
  }
  
  return(mean(lambda_vals, na.rm = TRUE))
}

lambda_ci_2 = function(samp, B = 1000){
  n = length(samp)
  boot_estimates = replicate(B, {
    samp_b = sample(samp, n, replace = TRUE)
    res = tryCatch(lambda_est2(samp_b), error = function(e) NA)
    res
  })
  
  boot_estimates = na.omit(boot_estimates)
  if (length(boot_estimates) < 2) return(c(NA, NA))
  return(quantile(boot_estimates, c(0.025, 0.975)))
}

lambda_ci_2_para = function(samp, B = 1000){
  n = length(samp)
  lambda_est = lambda_est2(samp)
  boot_estimates = replicate(B, {
    samp_b = rtlambda(n,lambda_est)
    res = tryCatch(lambda_est2(samp_b), error = function(e) NA)
  })
  
  boot_estimates = na.omit(boot_estimates)
  if (length(boot_estimates) < 2) return(c(NA, NA))
  return(quantile(boot_estimates, c(0.025, 0.975)))
}

ci_2_coverage_width = function(n, lambda, R = 100, B = 100){
  res = replicate(R, {
    samp = rtlambda(n, lambda)
    ci = lambda_ci_2(samp, B = B)
    coverage = (ci[1] <= lambda) && (lambda <= ci[2])
    width = ci[2] - ci[1]
    c(coverage, width)
  })
  
  coverage_rate = mean(res[1, ], na.rm = TRUE) * 100
  avg_width = mean(res[2, ], na.rm = TRUE)
  return(c(coverage_rate, avg_width))
}

ci_2_para_coverage_width = function(n, lambda, R = 100, B = 100){
  res = replicate(R, {
    samp = rtlambda(n, lambda)
    ci = lambda_ci_2_para(samp, B = B)
    coverage = (ci[1] <= lambda) && (lambda <= ci[2])
    width = ci[2] - ci[1]
    c(coverage, width)
  })
  
  coverage_rate = mean(res[1, ], na.rm = TRUE) * 100
  avg_width = mean(res[2, ], na.rm = TRUE)
  return(c(coverage_rate, avg_width))
}

results_quant_boot <- list()

for (lambda in lambda_values) {
  res = t(sapply(n_values, function(n) ci_2_coverage_width(n, lambda, R = 100, B = 500)))
  df = data.frame(
    lambda = lambda,
    n = n_values,
    coverage = res[, 1],
    width = res[, 2]
  )
  results_quant_boot[[as.character(lambda)]] = df
}

results_quant_boot_df <- do.call(rbind, results_quant_boot)


results_quant_boot_para <- list()

for (lambda in lambda_values) {
  res = t(sapply(n_values, function(n) ci_2_para_coverage_width(n, lambda, R = 100, B = 500)))
  df = data.frame(
    lambda = lambda,
    n = n_values,
    coverage = res[, 1],
    width = res[, 2]
  )
  results_quant_boot_para[[as.character(lambda)]] = df
}

results_quant_boot_para_df <- do.call(rbind, results_quant_boot_para)

#---------------------------------------------------------------
## 2nd Competing Method : Lmoments comparison
## Then applying parametric and nonparametric bootstrap for CI, their width and coverage
## lmoments boot (nonpara boot and para boot)

lambda_est6 = function(samp) {
  f = function(lambda){
    if(lambda == 0 || lambda <= -1) return(NA)
    else return(2/lambda * (2/(2+lambda)-1/(1+lambda)))
  }
  v = Lmoments(samp)[2]
  g = function(lambda) { return(f(lambda) - v) }
  
  lambda_est = tryCatch({
    uniroot(g, lower = -1 + 1e-10, upper = 10, tol = 1e-10)$root
  }, error = function(e) NA)
  
  return(lambda_est)
}

lambda_ci_6 = function(samp, B = 500){
  n = length(samp)
  boot_estimates = replicate(B, {
    samp_b = sample(samp, n, replace = TRUE)
    lambda_est6(samp_b)
  })
  boot_estimates = na.omit(boot_estimates)
  if(length(boot_estimates) < 2) return(c(NA, NA))
  else return(quantile(boot_estimates, c(0.025, 0.975)))
}

lambda_ci_6_para = function(samp, B = 500){
  n = length(samp)
  lambda_est = lambda_est6(samp)
  boot_estimates = replicate(B, {
    samp_b = rtlambda(n,lambda_est)
    lambda_est6(samp_b)
  })
  boot_estimates = na.omit(boot_estimates)
  if(length(boot_estimates) < 2) return(c(NA, NA))
  else return(quantile(boot_estimates, c(0.025, 0.975)))
}

ci_6_coverage_width = function(n, lambda, R = 100, B = 500){
  res = replicate(R, {
    samp = rtlambda(n, lambda)
    ci = lambda_ci_6(samp, B)
    covered = (ci[1] <= lambda) && (lambda <= ci[2])
    width = ci[2] - ci[1]
    c(covered, width)
  })
  
  coverage_rate = mean(res[1, ], na.rm = TRUE) * 100
  avg_width = mean(res[2, ], na.rm = TRUE)
  return(c(coverage_rate, avg_width))
}

ci_6_para_coverage_width = function(n, lambda, R = 100, B = 500){
  res = replicate(R, {
    samp = rtlambda(n, lambda)
    ci = lambda_ci_6_para(samp, B)
    covered = (ci[1] <= lambda) && (lambda <= ci[2])
    width = ci[2] - ci[1]
    c(covered, width)
  })
  
  coverage_rate = mean(res[1, ], na.rm = TRUE) * 100
  avg_width = mean(res[2, ], na.rm = TRUE)
  return(c(coverage_rate, avg_width))
}


results_lmom_boot <- list()

for (lambda in lambda_values) {
  res = t(sapply(n_values, function(n) ci_6_coverage_width(n, lambda, R = 100, B = 500)))
  df = data.frame(
    lambda = lambda,
    n = n_values,
    coverage = res[, 1],
    width = res[, 2]
  )
  results_lmom_boot[[as.character(lambda)]] = df
}

results_lmom_boot_df <- do.call(rbind, results_lmom_boot)


results_lmom_boot_para <- list()

for (lambda in lambda_values) {
  res = t(sapply(n_values, function(n) ci_6_para_coverage_width(n, lambda, R = 100, B = 500)))
  df = data.frame(
    lambda = lambda,
    n = n_values,
    coverage = res[, 1],
    width = res[, 2]
  )
  results_lmom_boot_para[[as.character(lambda)]] = df
}

results_lmom_boot_para_df <- do.call(rbind, results_lmom_boot_para)


#----------------------------------------------------
## DW CI
results_dw <- data.frame(lambda = numeric(), n = integer(), coverage = double(), width = double())
for (lambda in lambda_values) {
  for (n in n_values) {
    res <- coverage_lambda_dw(lambda, n, R = 100)
    results_dw <- rbind(results_dw, data.frame(lambda = lambda, n = n, coverage = res[1], width = res[2]))
  }
}



#----------------------------------------------------
## Comparing the results from various methods
results_dw$source <- "DW"
results_quant_boot_df$source <- "quant(NP boot)"
results_quant_boot_para_df$source <- "quant(P boot)"
results_lmom_boot_df$source <- "lmom(NP boot)"
results_lmom_boot_para_df$source <- "lmom(P boot)"
all_results <- bind_rows(results_dw, results_quant_boot_df, results_quant_boot_para_df, results_lmom_boot_df,results_lmom_boot_para_df)
all_results$source <- factor(all_results$source, levels = c("DW", "quant(NP boot)","quant(P boot)","lmom(NP boot)","lmom(P boot)"))
# ggplot(all_results, aes(x = n, y = width, color = source, linetype = factor(lambda))) +
#   geom_line(size = 1) +
#   scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash","dashed")) +
#   labs(
#     title = "Confidence Interval Width vs Sample Size",
#     x = "Sample Size (n)",
#     y = "Interval Width",
#     color = "Method",
#     linetype = "Lambda"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(face = "bold", size = 16),  # Bold main title
#     axis.title.x = element_text(face = "bold", size = 14),  # Bold x-axis title
#     axis.title.y = element_text(face = "bold", size = 14),  # Bold y-axis title
#     legend.title = element_text(face = "bold", size = 12),  # Bold legend titles
#     legend.text = element_text(size = 10),  # Larger legend text (indices)
#     legend.key.size = unit(1.5, "lines")  # Larger legend keys (symbols)
#   )
# ggplot(all_results, aes(x = n, y = coverage / 100, color = source, linetype = factor(lambda))) +
#   geom_line(size = 1) +
#   geom_hline(yintercept = 0.95, color = "black", size = 1.2, linetype = "longdash") +  
#   scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash")) +
#   labs(
#     title = "Confidence Interval Coverage vs Sample Size",
#     x = "Sample Size (n)",
#     y = "Coverage",
#     color = "Method",
#     linetype = "Lambda"
#   ) +
#   ylim(0, 1) +
#   theme_minimal()+
#   theme(
#     plot.title = element_text(face = "bold", size = 16),  # Bold main title
#     axis.title.x = element_text(face = "bold", size = 14),  # Bold x-axis title
#     axis.title.y = element_text(face = "bold", size = 14),  # Bold y-axis title
#     legend.title = element_text(face = "bold", size = 12),  # Bold legend titles
#     legend.text = element_text(size = 10),  # Larger legend text (indices)
#     legend.key.size = unit(1.5, "lines")  # Larger legend keys (symbols)
#   )

all_results$coverage = all_results$coverage/100
all_results_long <- pivot_longer(all_results, cols = c(coverage, width), names_to = "metric", values_to = "value")
all_results_long <- all_results_long %>%
  mutate(n = factor(n, levels = sort(unique(n)), labels = paste0("n = ", sort(unique(n)))))

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


#write.csv(all_results_long,"Documents/Academic_Stuffs/Projects/AK/tlambda_method_comparison.csv")
#all_results_long = read.csv("Documents/Academic_Stuffs/Projects/AK/tlambda_method_comparison.csv")
  