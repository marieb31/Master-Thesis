# R code - Chapter 4 + 5 - Useful functions - Master Thesis

## Load Packages ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(igraph)

## Functions: ------------------------------------------------------------------

### Get the location for a given ID:
get_loc_from_id <- function(ID, ID_db) {
  # for a given ID, returns the location on the grid in format (x,y)
  loc <- ID_db[ID_db$ID == ID, c("gridX", "gridY")]
  return(loc)
}

### Get the direction of a location: 
get_direction <- function(start_location, location_ID, ID_db) {
  # Returns the direction of a location from a given start location
  ID_subset <- ID_db[ID_db$ID == location_ID,]
  loc_X <- ID_subset[ "gridX"]
  loc_Y <- ID_subset["gridY"]
  start_X <- start_location["gridX"]
  start_Y <- start_location["gridY"]
  
  # classify and return the direction (North-East, North-West, South-East, South-West...)
  if (loc_X >= start_X & loc_Y > start_Y) {
    return("NE")
  } else if (loc_X < start_X & loc_Y >= start_Y) {
    return("NW")
  } else if (loc_X > start_X & loc_Y <= start_Y) {
    return("SE")
  } else if (loc_X <= start_X & loc_Y < start_Y) {
    return("SW")
  } else {
    return("center") # location = start_location
  }
}

### Plot empirical mean values on the gris
plot_means_on_grid <- function(grid_data, ID_db) {
  
  # For a simulated datasets, plots mean values per locations on the grid
  
  mean_values <- apply(grid_data, 2, mean)
  data_plot <- cbind(mean_values, ID_db)
  
  trigger_locs <- ID_db[ID_db$Fixed_trigger == 1,]
  
  ggplot() +
    geom_tile(data = data_plot, aes(x = gridX, y = gridY, fill = mean_values), color = "grey") +  
    geom_point(data = trigger_locs, aes(x = gridX, y = gridY, shape = "Event Start"), 
               size = 4, color = "black") +
    scale_fill_gradient(low = "white", high = "red") +  
    scale_shape_manual(values = c("Event Start" = 4)) +
    theme_minimal() +
    theme(panel.grid = element_blank(), aspect.ratio = 1, 
          axis.text = element_text(size = 14), axis.title = element_text(size = 16),
          legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
    labs(fill = "Mean Value", shape = "", x = "Grid X", y = "Grid Y")+
    scale_x_continuous(breaks = seq(1, max(ID_db$gridX), 2), labels = seq(1, max(ID_db$gridX), 2))+
    scale_y_continuous(breaks = seq(1, max(ID_db$gridX), 2), labels = seq(1, max(ID_db$gridX), 2))
}


### Scale Precision matrix for a given maximum distance
scale_precision_matrix <- function(precision_matrix, dist_max, ID_db) {
  # for a given max distance and precision matrix, returns the precision matrix
  # with 0 entries for locations which have a distance above dist_max
  
  ids <- colnames(precision_matrix) # IDs of all locations
  d <- nrow(precision_matrix) # number of dimensions
  
  # double loop over all pair combinaitions of locations
  for (i in 1:d) {
    loc_i <- get_loc_from_id(ids[i], ID_db) # get location of ID i
    for (j in 1:d) {
      if (i != j) {
        loc_j <- get_loc_from_id(ids[j], ID_db) # get location of ID j
        
        # calculate distance between the two locations
        distance <- as.integer(sqrt((loc_i[1] - loc_j[1])^2 + (loc_i[2] - loc_j[2])^2))
        
        if (distance > dist_max) {
          precision_matrix[i, j] <- 0
        }
      }
    }
  }
  
  return(precision_matrix)
}

### Plot a Graph given a precision matrix
graph_from_pm <- function(precision_matrix, Chi_matrix = NULL, ID_db) {
  # Plots on the grid the graphical model for a given precision matrix. 
  # If Chi_matrix is specified, the edges are represented with a color-code corresponding to the estimated extremal correlations
  
  IDs <- colnames(precision_matrix) 
  locations <- ID_db %>% filter(ID %in% IDs)

  edges <- which(precision_matrix != 0 & upper.tri(precision_matrix), arr.ind = TRUE)
  
  if (is.null(Chi_matrix)) {
    
    edge_list <- data.frame(from = rownames(precision_matrix)[edges[,1]],
                            to = colnames(precision_matrix)[edges[,2]])
    
    edge_list <- merge(edge_list, locations, by.x = "from", by.y = "ID", all.x = TRUE)
    edge_list <- merge(edge_list, locations, by.x = "to", by.y = "ID", suffixes = c("_from", "_to"), all.x = TRUE)
    
    ggplot(locations, aes(x = gridX, y = gridY)) +
      geom_segment(data = edge_list, 
                   aes(x = gridX_from, y = gridY_from, 
                       xend = gridX_to, yend = gridY_to),
                   color = "blue",
                   alpha = 0.2,
                   size = 1) + # plot the edges with no color code (all blue) but with transparency
      geom_point(size = 3, alpha=0.8) + # plot the nodes on the grid
      labs(title = "Graph from Precision Matrix") +
      theme_minimal()
    
  } 
  
  else {
    # Same as before but color code according to Chi_matrix
    edge_list <- data.frame(from = rownames(precision_matrix)[edges[,1]],
                            to = colnames(precision_matrix)[edges[,2]],
                            chi_entry = Chi_matrix[edges])
    
    edge_list <- merge(edge_list, locations, by.x = "from", by.y = "ID", all.x = TRUE)
    edge_list <- merge(edge_list, locations, by.x = "to", by.y = "ID", suffixes = c("_from", "_to"), all.x = TRUE)
    
    ggplot(locations, aes(x = gridX, y = gridY)) +
      geom_segment(data = edge_list, 
                   aes(x = gridX_from, y = gridY_from, 
                       xend = gridX_to, yend = gridY_to, 
                       color = chi_entry, # edge color based on corresponding extremal correlation
                       alpha = chi_entry), # edge transparency based on corresponding extremal correlation
                   size = 1) +
      geom_point(size = 3, alpha=0.8) + # plot the nodes on the grid
      scale_alpha_continuous(range = c(0.6, 1), guide = "none") + 
      scale_color_gradientn(colours = c("white", "green", "darkgreen")) +
      labs(title = "Graph from Precision Matrix", color = expression(hat(chi)[ij]), x = "Grid X", y = "Grid Y") +
      theme_minimal()
  }
 
}


