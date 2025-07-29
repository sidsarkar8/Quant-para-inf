library(ggplot2)

# Tukey Lambda quantile function
Q_tukey <- function(p, lambda) {
  if (lambda == 0) return(log(p / (1 - p)))  # Logistic
  (p^lambda - (1 - p)^lambda) / lambda
}

# Derivative of quantile function
Q_prime_tukey <- function(p, lambda) {
  if (lambda == 0) return(1 / (p * (1 - p)))  # Logistic
  p^(lambda - 1) + (1 - p)^(lambda - 1)
}

# Set grid and parameter
lambda <- 0.14  # Approximately normal
p_grid <- seq(0.001, 0.999, length.out = 1000)

# Tukey Lambda: get x and density
x_tukey <- Q_tukey(p_grid, lambda)
f_tukey <- 1 / Q_prime_tukey(p_grid, lambda)

# Standard normal density
x_norm <- seq(min(x_tukey), max(x_tukey), length.out = 1000)
f_norm <- dnorm(x_norm)

# Combine into a data frame for plotting
df <- rbind(
  data.frame(x = x_tukey, density = f_tukey, dist = "Tukey Lambda (λ = 0.14)"),
  data.frame(x = x_norm, density = f_norm, dist = "Normal(0,1)")
)

# Plot
ggplot(df, aes(x = x, y = density, color = dist)) +
  geom_line(size = 1.1) +
  labs(title = "Tukey Lambda vs. Normal Distribution",
       subtitle = "λ = 0.14 approximately matches Normal(0,1)",
       x = "x", y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())
