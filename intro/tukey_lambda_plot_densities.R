library(ggplot2)
library(dplyr)

# Tukey Lambda quantile and derivative
Q_tukey <- function(p, lambda) {
  if (lambda == 0) return(log(p / (1 - p)))
  (p^lambda - (1 - p)^lambda) / lambda
}

Q_prime_tukey <- function(p, lambda) {
  if (lambda == 0) return(1 / (p * (1 - p)))
  p^(lambda - 1) + (1 - p)^(lambda - 1)
}

# Setup
p_grid <- seq(0.001, 0.999, length.out = 1000)
lambda_vals <- c(-1, 0, 0.14, 1, 2)  # Now includes λ = 1

# Compute densities
densities <- bind_rows(lapply(lambda_vals, function(lam) {
  q_vals <- Q_tukey(p_grid, lam)
  f_vals <- 1 / Q_prime_tukey(p_grid, lam)
  df <- data.frame(x = q_vals, density = f_vals, lambda = as.factor(lam))
  
  # Add zero-density points outside support if bounded (λ ≥ 1)
  if (lam >= 1) {
    x_left <- Q_tukey(0, lam)
    x_right <- Q_tukey(1, lam)
    df <- bind_rows(
      data.frame(x = x_left, density = 0, lambda = as.factor(lam)),
      df,
      data.frame(x = x_right, density = 0, lambda = as.factor(lam))
    )
  }
  return(df)
}))

# Plot
ggplot(densities, aes(x = x, y = density, color = lambda, group = lambda)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_cartesian(xlim = c(-4, 4), ylim = c(0, 1)) +
  labs(title = "TL densities: different choice of lambda",
       subtitle = "-1 (heavy), 0 (logistic), 0.14 (normal), 1,2 (uniform)",
       x = "x", y = "Density") +
  theme_minimal() +
  theme(legend.title = element_text(size = 10))
