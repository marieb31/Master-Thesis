# R code - Chapter 4 - Data Simulation - Master Thesis

## Load Packages ---------------------------------------------------------------
library(ggplot2)
library(evd) # for EVT distributions

## Set directory where to save the simulated data sets -------------------------
setwd("C:/Users/Marie BERNIERE/Desktop/Master Arbeit/MasterArbeitRepo/Code-Scripts")

## Load useful functions ------------------------------------------------------- 
source("utils.R")

## Visualizing sampling distribution of magnitude marks ------------------------
### Marks are set to 0 in case there is no event at the given location
### If an event occurs, the marks are following a Generalized Pareto Distribution

set.seed(1234)
x_GPD <- rgpd(n = 10000, loc = 1, scale = 0.5 , shape = 0.5)

### Histogram
ggplot(data.frame(x = x_GPD), aes(x)) +
  geom_histogram(binwidth = 0.1, fill = "#336699", color="#6699CC") +
  xlab("Values - capped at 10") +
  xlim(0,10)+
  theme_minimal()+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 10)))

## Grid Data Simulator ---------------------------------------------------------

simulate_spatial_events <- function(N, d, n_trigger_loc, max_trigger = 1, spread_directions) {
  
  # creation of a database of all locations with an associated ID
  D <- d*d
  ID_db <- cbind(paste0("ID_", 1:D), expand.grid(1:d, 1:d))
  colnames(ID_db) <- c("ID", "gridX", "gridY")
  
  # Sample a set of fixed trigger locations and add to data base
  fixed_trigger_IDs <- sample(ID_db$ID, n_trigger_loc, replace=FALSE)
  ID_db$Fixed_trigger <- ifelse(ID_db$ID %in% fixed_trigger_IDs, 1, 0)
  
  # Init Simulation Matrix (grid_data)
  grid_data <- matrix(NA, nrow = N, ncol = D)
  colnames(grid_data) <- paste0("ID_", 1:D)
  
  # iterate over N time periods
  for (t in 1:N){
    
    # Is there an event triggered ?
    is_event <- rbinom(1, 1, prob = 0.2)
    
    # if no event is triggered - zero mark
    if (is_event == 0) {
      for (j in 1:D) {
        if (is.na(grid_data[t,j])) {
          grid_data[t, j] = 0
        }
      }
      next # go to next iteration
    }
    
    # An event is triggered: Epidemic Spread
    
    # How many events ?
    n_events_triggered <- sample(1:max_trigger, 1) # How many start locations are triggered
    event_start_IDs <- sample(fixed_trigger_IDs, 
                              n_events_triggered, 
                              replace = FALSE) # not allowing several events on same location
    
    # Loop over each Event
    for (start_ID in event_start_IDs) {
      
      # trigger location gets a value for current iteration t
      grid_data[t, start_ID] <- rgpd(n = 1, loc = 1, scale = 0.5 , shape = 0.5)
      
      # determine epidemic spread - spread distance
      spread_distance <- sample(1:3, 1) # min distance is 1 and max distance is 3
      
      # get distance to start location for all locations
      start_location <- get_loc_from_id(start_ID, ID_db)
      
      distances <- sqrt((ID_db["gridX"] - as.vector(start_location["gridX"]))^2 + 
                        (ID_db["gridY"] - as.vector(start_location["gridY"]))^2) 
      
      # Loop over distances
      for (k in 1:spread_distance) {
        
        ID_potential_locations <- ID_db[which(distances <= k & distances > (k-1)), "ID"]
        
        # Loop over locations within distance
        for (potential_location in ID_potential_locations) {
          
          if (t+k <= N) { # making sure we are within the N time iterations
            
            # get direction of the potentially affected location
            direction_from_start <- get_direction(start_location, potential_location, ID_db) 
            
            if (direction_from_start  %in% spread_directions) { # If within specified direction
              affected_loc <- potential_location
              grid_data[t+k, affected_loc] = rgpd(n = 1, loc = 1, scale = 0.5 , shape = 0.5)
            } 
            
            else {
              # assign random events to the other potential locations randomly
              is_event <- rbinom(1, 1, p=0.3) # probability of an event 0.3 
              if (is_event == 1) {
                grid_data[t+k, potential_location] = rgpd(n = 1, loc = 1, scale = 0.5 , shape = 0.5)
              }
            }
          } 
        }
        
      }
    }
    
    # Assign zero values to all locations that were not affected even though an event took place
    for (j in 1:D) {
      if (is.na(grid_data[t,j])) {
        grid_data[t, j] = 0
      }
    }
  }
  
  # Return final data
  results <- list(grid_data, ID_db)
  return(results) 
}

## Scenario Simulation ---------------------------------------------------------

# Scenario 1
set.seed(1234)
scenario_1 <- simulate_spatial_events(N = 10000, 
                                      d = 20, 
                                      n_trigger_loc = 3, 
                                      max_trigger = 1,
                                      spread_direction = c("NE", "SE", "NW", "SW"))

data_grid_S1 <- scenario_1[1]
ID_db_S1 <- scenario_1[2]

write.csv(data_grid_S1, file = "data_grid_S1_10000.csv")
write.csv(ID_db_S1, file = "ID_db_S1_10000.csv")

# Scenario 2
set.seed(1234)
scenario_2 <- simulate_spatial_events(N = 10000, 
                                      d = 20, 
                                      n_trigger_loc = 3, 
                                      max_trigger = 1,
                                      spread_direction = c("NE", "SE"))

data_grid_S2 <- scenario_2[1]
ID_db_S2 <- scenario_2[2]

write.csv(data_grid_S2, file = "data_grid_S2_10000.csv")
write.csv(ID_db_S2, file = "ID_db_S2_10000.csv")

# Scenario 3
set.seed(1234)
scenario_3 <- simulate_spatial_events(N = 10000, 
                                      d = 10,
                                      n_trigger_loc = 6, 
                                      max_trigger = 2,
                                      spread_direction = c("NE", "SE"))

data_grid_S3 <- scenario_3[1]
ID_db_S3 <- scenario_3[2]

write.csv(data_grid_S3, file = "data_grid_S3_10000.csv")
write.csv(ID_db_S3, file = "ID_db_S3_10000.csv")

