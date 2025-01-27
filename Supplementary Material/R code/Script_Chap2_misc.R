# R code - Chapter 2 Master Thesis

# Load Packages
library(ggplot2) 
library(spatstat) 
library(evd)

## Spatial Point Processes -----------------------------------------------------

### Illustration of poisson point process
set.seed(1234)
X <- rpoispp(5)
plot(X, pch=16, size=0.5, main="")
X <- rpoispp(15)
plot(X, pch=16, size=0.5, main="") 
X <- rpoispp(30)
plot(X, pch=16, size=0.5, main="") 

## Extreme Value Theory --------------------------------------------------------------------------

### Visualization of Generalized Pareto Distribution
set.seed(1234)
fixed_scale <- 1
fixed_shape <- 0 
shape_values <- c(-0.5, 0, 0.5, 1) # Different shape parameters for Plot 1
scale_values <- c(0.5, 1, 2)       # Different scale parameters for Plot 2
x_values <- seq(0, 30, length.out = 1000)  # Range of x values for density calculation

#### GPD for different shape parameters
gpd_shape_data <- data.frame()

for (shape in shape_values) {
  density_values <- dgpd(x_values, loc = 0, scale = fixed_scale, shape = shape)
  gpd_shape_data <- rbind(gpd_shape_data, data.frame(x = x_values, density = density_values, shape = shape))
}

gpd_shape_data$shape <- factor(gpd_shape_data$shape)

ggplot(gpd_shape_data, aes(x = x, y = density, color = shape)) +
  geom_line(size = 1) +
  labs(title = "GPD for different shapes and fixed scale (σ = 1)",
       x = "x (capped)",
       y = "Density",
       color = "Shape (ξ)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0, 10)

#### GPD for different scale parameters
gpd_scale_data <- data.frame()

for (scale in scale_values) {
  density_values <- dgpd(x_values, loc = 0, scale = scale, shape = fixed_shape)
  gpd_scale_data <- rbind(gpd_scale_data, data.frame(x = x_values, density = density_values, scale = factor(scale)))
}

gpd_scale_data$shape <- factor(gpd_scale_data$scale)

ggplot(gpd_scale_data, aes(x = x, y = density, color = scale)) +
  geom_line(size = 1) +
  labs(title = "GPD for different scales and fixed shape (ξ = 0)",
       x = "x (capped)",
       y = "Density",
       color = "Scale (σ)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(0, 10)

