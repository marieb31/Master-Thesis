# R code - Chapter 4 - Simulation Results - Master Thesis

## load Packages ---------------------------------------------------------------
library(ggplot2)
library(graphicalExtremes) # For Hüsler Reiss model fitting
library(spatstat) # For point process functions
library(dplyr)

## Set Working Directory -------------------------------------------------------
setwd("C:/Users/Marie BERNIERE/Desktop/Master Arbeit/MasterArbeitRepo/Code-Scripts")

## Load Functions --------------------------------------------------------------
source("utils.R") # import functions from utils

## Define functions for estimation routine -------------------------------------
estimate_HR_parameters <- function(X_data, threshold_quantile, method = "ns", criterion = "bic") {
  # For a given data frame, threshold quantile and algorithm choice, 
  # returns best Hüsler-Reiss graph and related matrices
  
  Y <- data2mpareto(X_data, p = threshold_quantile) # standardization to Pareto scale
  colnames(Y) <- colnames(X_data)
  
  rholist <- seq(0.01,0.21,0.02) # penalization range
  eglearn_fit <- eglearn(Y, rholist = rholist, reg_method = method, complete_Gamma = TRUE)
  
  criterion_eglearn_fit <- data.frame(rholist = rholist,
                                      values = NA,
                                      graph_size = NA)
  
  for (i in 1:length(rholist)) {
    Gamma <- eglearn_fit$Gamma[[i]]
    Graph <- eglearn_fit$graph[[i]]
    criterion_eglearn_fit[i, "values"] <- loglik_HR(data = Y, Gamma = Gamma, graph = Graph)[criterion]
    criterion_eglearn_fit[i, "graph_size"] <- gsize(Graph)
  }
  
  criterion_graph <- ggplot(data = criterion_eglearn_fit, aes(x = rholist, y = values)) +
    geom_line() +
    geom_point(shape = 16, size = 3) +
    labs(x = "rho", y = criterion) +
    scale_x_continuous(breaks = rholist,
                       labels = round(rholist, 2),
                       sec.axis = sec_axis(trans = ~.,
                                           breaks = rholist,
                                           labels = criterion_eglearn_fit$graph_size,
                                           name = "Number of edges")) +
    theme_minimal() +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
          legend.text = element_text(size = 14), legend.title = element_text(size = 16))
  
  
  best_index <- which.min(criterion_eglearn_fit[,"values"])
  best_Gamma <- eglearn_fit$Gamma[[best_index]] # Estimated Graph Variogram 
  
  k <- dim(Y)[2]
  theta_k <- Gamma2Theta(best_Gamma, k = k, full=TRUE) # Corresponding Precision Matrix
  theta_k[k,] <- colSums(theta_k) # fill last row
  theta_k[, k] <- rowSums(theta_k) # fill last column
  colnames(theta_k) <- colnames(X_data)
  rownames(theta_k) <- colnames(X_data) 
  
  output <- list(criterion_graph = criterion_graph, best_Gamma = best_Gamma, precision_matrix = theta_k)
  
  return(output)
}

estimate_max_dist <- function(data, ID_db, subset_window) {
  xy_coords <- data.frame()
  
  for (id in colnames(data)) {
    location <- as.integer(get_loc_from_id(id, ID_db))
    n_points <- sum(data[,id] > 0)
    for (n in 1:n_points) {
      xy_coords <- rbind(xy_coords, location)
    }
  }
  
  colnames(xy_coords) <- c("x", "y")
  ppp <- ppp(xy_coords$x, xy_coords$y, window = subset_window)
  ppp <- rjitter(ppp, 1) # add random jitter to simulate underlying point process
  kest <- Kest(ppp)
  max_R <- round(max(kest$r),1)
  
  return(max_R)
}

## Scenario 1 ------------------------------------------------------------------
### Load Simulated Data
grid_data_S1 <- read.csv("data_grid_S1_10000.csv")[,-1]
ID_db_S1<- read.csv("ID_db_S1_10000.csv")[,-1]

mean_plot <- plot_means_on_grid(grid_data =grid_data_S1, ID_db = ID_db_S1)

### Mean Plot
mean_plot

### Mean Plot for Subset
mean_plot +
  coord_cartesian(xlim = c(12, 20), ylim=c(13,21)) 

### Test_Data 
# Filter the data points using a grid subset and removing all locations with no event at all
ID_subset_loc <- ID_db_S1[ID_db_S1$gridX > 12, ]$ID # focusing on one area
ID_subset_1 <- colnames(grid_data_S1)[colSums(grid_data_S1 > 1) != 0] # only locations with at least one event
ID_subset <- intersect(ID_subset_loc, ID_subset_1) 

test_data_S1<- grid_data_S1[,ID_subset]
dim(test_data_S1) # measures at 29 sites over 10000 days

# aggregate the data weekly
aggregated_test_data_S1 <- apply(test_data_S1, 2, function(x) tapply(x, (seq_along(x) - 1) %/% 7, sum))
aggregated_test_data_S1 <- as.data.frame(aggregated_test_data_S1)
dim(aggregated_test_data_S1)
X_S1 <- aggregated_test_data_S1

## Scenario 2 -------------------------------------------------------------------
### Load Simulated Data 
grid_data_S2 <- read.csv("data_grid_S2_10000.csv")[,-1]
ID_db_S2 <- read.csv("ID_db_S2_10000.csv")[,-1]

mean_plot <- plot_means_on_grid(grid_data =grid_data_S2, ID_db = ID_db_S2)

### Mean Plot
mean_plot

### Mean Plot for Subset
mean_plot +
  coord_cartesian(xlim = c(12, 20), ylim=c(13,21)) 

## Test_Data 
# Filter the data points using a grid subset and removing all locations with no event at all
ID_subset_loc <- ID_db_S1[ID_db_S2$gridX > 12, ]$ID # focusing on one area
ID_subset_1 <- colnames(grid_data_S2)[colSums(grid_data_S2 > 1) != 0] # only locations with at least one event
ID_subset <- intersect(ID_subset_loc, ID_subset_1) 

test_data_S2 <- grid_data_S2[,ID_subset]
dim(test_data_S2) # measures at 29 sites over 10000 days

# aggregate the data weekly
aggregated_test_data_S2 <- apply(test_data_S2, 2, function(x) tapply(x, (seq_along(x) - 1) %/% 7, sum))
aggregated_test_data_S2 <- as.data.frame(aggregated_test_data_S2)
dim(aggregated_test_data_S2)
X_S2 <- aggregated_test_data_S2

## Scenario 3 ------------------------------------------------------------------
### Load Simulated Data 
grid_data_S3 <- read.csv("data_grid_S3_10000.csv")[,-1]
ID_db_S3 <- read.csv("ID_db_S3_10000.csv")[,-1]

mean_plot <- plot_means_on_grid(grid_data =grid_data_S3, ID_db = ID_db_S3)

### Mean Plot
mean_plot

### Filter Locations
# Filter the data points removing all locations with no event at all
ID_subset <- colnames(grid_data_S3)[colSums(grid_data_S3 > 1) != 0] # only locations with at least one event
test_data_S3 <- grid_data_S3[,ID_subset]
dim(test_data_S3) # measures at 68 sites over 10000 days

# aggregate the data weekly
aggregated_test_data_S3 <- apply(test_data_S3, 2, function(x) tapply(x, (seq_along(x) - 1) %/% 7, sum))
aggregated_test_data_S3 <- as.data.frame(aggregated_test_data_S3)
dim(aggregated_test_data_S3)
X_S3 <- aggregated_test_data_S3

## Graphical Modeling Results --------------------------------------------------

### results for Scenario 1:
Y_S1 <- data2mpareto(X_S1, p = 0.9) # check how many observations for threshold 90% quantile
print(paste(dim(Y_S1)[1], "observations in Scenario 1" ))

results_S1_ns_90 <- estimate_HR_parameters(X_S1, threshold_quantile = 0.9 , method = "ns", criterion = "bic")
results_S1_glasso_90 <- estimate_HR_parameters(X_S1, threshold_quantile = 0.9 , method = "glasso", criterion = "bic")
results_S1_ns_70 <- estimate_HR_parameters(X_S1, threshold_quantile = 0.7 , method = "ns", criterion = "bic")

### results for Scenario 2:
Y_S2 <- data2mpareto(X_S2, p = 0.9) # check how many observations for threshold 90% quantile
print(paste(dim(Y_S2)[1], "observations in Scenario 2" ))

results_S2_ns_90 <- estimate_HR_parameters(X_S2, threshold_quantile = 0.9 , method = "ns", criterion = "bic")
results_S2_glasso_90 <- estimate_HR_parameters(X_S2, threshold_quantile = 0.9 , method = "glasso", criterion = "bic")
results_S2_ns_70 <- estimate_HR_parameters(X_S2, threshold_quantile = 0.7 , method = "ns", criterion = "bic")

### results for Scenario 3:
Y_S3 <- data2mpareto(X_S3, p = 0.9) # check how many observations for threshold 90% quantile
print(paste(dim(Y_S3)[1], "observations in Scenario 3" ))

results_S3_ns_90 <- estimate_HR_parameters(X_S3, threshold_quantile = 0.9 , method = "ns", criterion = "bic")
results_S3_glasso_90 <- estimate_HR_parameters(X_S3, threshold_quantile = 0.9 , method = "glasso", criterion = "bic")
results_S3_ns_70 <- estimate_HR_parameters(X_S3, threshold_quantile = 0.7 , method = "ns", criterion = "bic")


## Visualize the results -------------------------------------------------------
# The results can be visualized by specfifying the following inputs:
### Inputs:
results <- results_S1_ns_90 
ID_db <- ID_db_S1
test_data <- test_data_S1
subset_window <- owin(c(12,20), c(12,20))  # For scenario 1 + 2: owin(c(12,20), c(12,20)), for scenario 1 owin(c(0,10), c(0,10))

# Run for results
results$criterion_graph
best_Gamma <- results$best_Gamma
precision_matrix <- results$precision_matrix
chi_matrix <- Gamma2chi(best_Gamma)

trigger_points <- ID_db %>% filter(Fixed_trigger == 1 & ID %in% colnames(test_data))

g_res <- graph_from_pm(precision_matrix = precision_matrix, ID_db = ID_db)
g_res + theme(aspect.ratio = 1, 
              axis.text = element_text(size = 16), axis.title = element_text(size = 16),
              legend.text = element_text(size = 16), legend.title = element_text(size = 16),
              title = element_text(size = 14)) +
        geom_point(data = trigger_points, aes(x = gridX, y = gridY), color ="red", size = 3)

g_res_chi <- graph_from_pm(precision_matrix = precision_matrix, Chi_matrix = chi_matrix, ID_db=ID_db)
g_res_chi + theme(aspect.ratio = 1, 
                  axis.text = element_text(size = 16), axis.title = element_text(size = 16),
                  legend.text = element_text(size = 16), legend.title = element_text(size = 16),
                  title = element_text(size = 14)) +
            geom_point(data = trigger_points, aes(x = gridX, y = gridY), color ="red", size = 3)

# Determine a max distance
max_dist <- estimate_max_dist(data = test_data, ID_db = ID_db, subset_window = subset_window)

precision_scaled <- scale_precision_matrix(precision_matrix, ID_db, dist_max = max_dist)
g_res_scaled <- graph_from_pm(precision_scaled, ID_db = ID_db)
g_res_scaled + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        title = element_text(size = 14))+
  labs(title = paste("Graph from Scaled Precision Matrix (max dist =", max_dist, ")")) +
  geom_point(data = trigger_points, aes(x = gridX, y = gridY), color ="red", size = 3)

g_res_scaled_chi <- graph_from_pm(precision_scaled, Chi_matrix = chi_matrix, ID_db = ID_db)
g_res_scaled_chi + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        title = element_text(size = 14))+
  labs(title = paste("Graph from Scaled Precision Matrix (max dist =", max_dist, ")")) +
  geom_point(data = trigger_points, aes(x = gridX, y = gridY), color ="red", size = 3)

### Graph Size for the scaled version:
graphB_size <- sum(precision_scaled[upper.tri(precision_scaled, diag = FALSE)] != 0)
print(paste("The scaled precision matrix has", graphB_size, "non zero entries."))


## Additional Plot -------------------------------------------------------------
### Number of locations where quantile is strictly positive for scenario 1
quantiles_test <- c(0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
portion_above_0 <- rep(NA, length(quantiles_test))

for (i in 1:length(quantiles_test)) {
  q <- quantiles_test[i]
  portion_above_0[i] <- sum(apply(X_S1, 2, quantile, q) > 0) / 29 
}

ggplot(data.frame(quantiles_test, portion_above_0), aes(x = quantiles_test, y = portion_above_0)) +
  geom_bar(stat = "identity", fill="steelblue") +
  labs(x = "Quantile levels", 
       y = "Portion of locations") +
  theme_minimal()
