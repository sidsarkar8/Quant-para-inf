param_grid_100_unif = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Unif/chixi_unif_100.csv")
param_grid_500_unif = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Unif/chixi_unif_500.csv")
param_grid_1000_unif = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Unif/chixi_unif_1000.csv")
param_grid_2000_unif = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Unif/chixi_unif_2000.csv")
param_grid_5000_unif = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Unif/chixi_unif_5000.csv")


param_grid_100_nor = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Normal/chixi_norm_100.csv")
param_grid_500_nor = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Normal/chixi_norm_500.csv")
param_grid_1000_nor = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Normal/chixi_norm_1000.csv")
param_grid_2000_nor = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Normal/chixi_norm_2000.csv")
param_grid_5000_nor = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Normal/chixi_norm_5000.csv")



param_grid_100_lnor = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Lognormal/chixi_lnorm_100.csv")
param_grid_500_lnor = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Lognormal/chixi_lnorm_500.csv")
param_grid_1000_lnor = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Lognormal/chixi_lnorm_1000.csv")
param_grid_2000_lnor = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Lognormal/chixi_lnorm_2000.csv")
param_grid_5000_lnor = read.csv("~/Documents/Projects/GLD_inference/Chi-Xi_CI/Lognormal/chixi_lnorm_5000.csv")



param_grid_100_unif$dist = "Uniform"
param_grid_100_unif$n = "n=100"
param_grid_500_unif$dist = "Uniform"
param_grid_500_unif$n = "n=500"
param_grid_1000_unif$dist = "Uniform"
param_grid_1000_unif$n = "n=1000"
param_grid_2000_unif$dist = "Uniform"
param_grid_2000_unif$n = "n=2000"
param_grid_5000_unif$dist = "Uniform"
param_grid_5000_unif$n = "n=5000"

param_grid_100_nor$dist = "Normal"
param_grid_100_nor$n = "n=100"
param_grid_500_nor$dist = "Normal"
param_grid_500_nor$n = "n=500"
param_grid_1000_nor$dist = "Normal"
param_grid_1000_nor$n = "n=1000"
param_grid_2000_nor$dist = "Normal"
param_grid_2000_nor$n = "n=2000"
param_grid_5000_nor$dist = "Normal"
param_grid_5000_nor$n = "n=5000"

param_grid_100_lnor$dist = "Lognormal"
param_grid_100_lnor$n = "n=100"
param_grid_500_lnor$dist = "Lognormal"
param_grid_500_lnor$n = "n=500"
param_grid_1000_lnor$dist = "Lognormal"
param_grid_1000_lnor$n = "n=1000"
param_grid_2000_lnor$dist = "Lognormal"
param_grid_2000_lnor$n = "n=2000"
param_grid_5000_lnor$dist = "Lognormal"
param_grid_5000_lnor$n = "n=5000"

chi_xi_grand_result = rbind(param_grid_100_unif,param_grid_500_unif,param_grid_1000_unif,param_grid_2000_unif,param_grid_5000_unif,
                            param_grid_100_nor,param_grid_500_nor,param_grid_1000_nor,param_grid_2000_nor,param_grid_5000_nor,
                            param_grid_100_lnor,param_grid_500_lnor,param_grid_1000_lnor,param_grid_2000_lnor,param_grid_5000_lnor)

# Convert admissible to numeric/factor if needed
chi_xi_grand_result <- chi_xi_grand_result %>%
  mutate(
    admissible = factor(admissible),
    n = factor(n, levels = c("n=100", "n=500", "n=1000", "n=2000", "n=5000")),  # ordered numerically
    dist = factor(dist, levels = c("Uniform", "Normal", "Lognormal"))  # optional: your desired row order
  )


highlight_point <- data.frame(
  chi = c(0, 0, 0.2844),           # example x-coordinates
  xi = c(1/2 - 1/sqrt(5), 0.3661, 0.3583),           # example y-coordinates
  dist = factor(c("Uniform", "Normal", "Lognormal"),
                levels = levels(chi_xi_grand_result$dist))
)

ggplot(chi_xi_grand_result, aes(x = chi, y = xi)) +
  geom_raster(aes(fill = admissible)) +
  geom_point(
    data = highlight_point,
    aes(x = chi, y = xi),
    shape = 21,
    fill = "gold",
    color = "black",
    size = 2,
    stroke = 1.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c("FALSE" = "blue", "TRUE" = "red"),
    name = "CI",
    labels = c("FALSE" = "Outside CI", "TRUE" = "Inside CI")
  ) +
  #coord_fixed(ratio = 1) +
  facet_grid(dist ~ n) +
  labs(
   # title = expression("CI for " * chi * " and " * xi * " for various distributions"),
    x = expression(chi),
    y = expression(xi)
  ) +
# theme(legend.position="none") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    #legend.title = element_text(size = 14, face = "bold"),
    #legend.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray85")
  ) 

ggsave("~/Documents/Projects/GLD_inference/Chi-Xi_CI/ppt_chi_xi_ci.pdf",  
       width = 10, height = 6, dpi = 300, units = "in" )

#write.csv(chi_xi_grand_result,"~/Documents/Projects/GLD_inference/Chi-Xi_CI/chi_xi_CI.csv")
