# R script - Chap 5 - Wildfire data application - Master Thesis

# Load Packages ----------------------------------------------------------------
library(readxl) # to read the data tables
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate) # for manipulating dates 
library(extRemes) # for extreme value analysis
library(spatstat) # for point process functions
library(sf) # for handling long/lat and utm conversion
library(graphicalExtremes) # for extremal graphical modeling
library(igraph)

# Set working directory --------------------------------------------------------
setwd("C:/Users/Marie BERNIERE/Desktop/Master Arbeit/MasterArbeitRepo/Code-Scripts")

# Load an define useful functions ----------------------------------------------
source("utils.R")

get_poblacion_from_ID <- function(ID) {
  # for a given ID of the northern subset, returns the corresponding PoblaciÓn
  # In some cases, several Poblacion are represented by 1 ID.
  # the poblacion which occurs the most for this ID is returned.
  counts <- clmfires_north_df[clmfires_north_df$ID == ID,] %>% 
            group_by(Población) %>% summarise(count = n())
  
  Población <- counts$Población[which.max(counts$count)]
  
  return(Población)
}

to_utm <- function(data) {
  # converts lon lat columns present in data to utm coordinates in km
  # works for zone 30 (which is the one corresponding to Spain)
  sf_data <- st_as_sf(data, coords = c("Longitud", "Latitud"), crs = 4326)
  
  sf_data_utm <- st_transform(sf_data, crs = "+proj=utm +zone=30 +datum=WGS84")
  utm_coords <- st_coordinates(sf_data_utm)
  
  utm_coords_km <- utm_coords / 1000
  
  return(utm_coords_km)
}

# load data --------------------------------------------------------------------
data(clmfires)
clmfires_df <- as.data.frame(clmfires)
LaMancha_Outline <- as.data.frame(vertices(clmfires$window)[1])

lon_lat_spain <- read_excel("C:/Users/Marie BERNIERE/Desktop/Master Arbeit/MasterArbeitRepo/Data/listado-longitud-latitud-municipios-espana.xls",
                            sheet = 1, range = "A3:I8115")
lon_lat_LaMancha <- lon_lat_spain %>% filter(Comunidad == "Castilla La Mancha")
lon_lat_LaMancha[, c("utm_x", "utm_y")] <- to_utm(lon_lat_LaMancha)

# Data Description -------------------------------------------------------------
## Summaries
summary(clmfires)
summary(clmfires_df)

## Are there duplicates in points
duplicates <- dim(clmfires_df)[1] - dim(unique(clmfires_df[, c("x", "y")]))[1] 
print(paste("There are", duplicates, "duplicated points in the dataset"))

## Visualize the data frame:
clmfires_df_head <- head(clmfires_df)
clmfires_df_head

## Image plots of underlying geographic covariates 
images_200pxl <- clmfires.extra$clmcov200
par(mar = c(1, 1, 1, 1))

image(images_200pxl$elevation, main = "Elevation", las = 2)
polygon(LaMancha_Outline, border = "white", lwd = 4)

image(images_200pxl$slope, main = "Slope", las = 2)
polygon(LaMancha_Outline, border = "white", lwd = 4)

image(images_200pxl$landuse, main = "Landuse", las = 2)
polygon(LaMancha_Outline, border = "white", lwd = 4)

## Visualizations of the marked Point Pattern 
ggplot(clmfires_df, aes(x = x, y = y, color = cause))+
  geom_point()+
  geom_polygon(data = LaMancha_Outline, aes(x = x, y = y), fill = NA, color = "black" ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = c(0.1,0.8), legend.text = element_text(size = 16), legend.title = element_text(size = 18),
        aspect.ratio = 1)+
  guides(color = guide_legend(title = "Fire Cause"))

ggplot(clmfires_df, aes(x = x, y = y, color = burnt.area))+
  geom_polygon(data = LaMancha_Outline, aes(x = x, y = y), fill = NA, color = "black" ) +
  geom_point(alpha=0.5) +
  scale_color_gradient(
    low = "blue", high = "red", trans = "log",
    name = "Burned Land (ha)", breaks = c(1, 10, 100, 1000) ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = c(0.1,0.8), legend.text = element_text(size = 16), legend.title = element_text(size = 18),
        aspect.ratio = 1)

par(mar=c(5.1, 4.1, 4.1, 2.1)) #reset to default graphical parameters

## Match point process to administrative infos
### using min and max values of the utm coordinates of the two data set 

ref_coord_min_x <- lon_lat_LaMancha[which.min(lon_lat_LaMancha$utm_x),]
ref_fire_min_x <- clmfires_df[which.min(clmfires$x),]
ref_coord_max_x <- lon_lat_LaMancha[which.max(lon_lat_LaMancha$utm_x),]
ref_fire_max_x <- clmfires_df[which.max(clmfires$x),]

ref_coord_min_y <- lon_lat_LaMancha[which.min(lon_lat_LaMancha$utm_y),]
ref_fire_min_y <- clmfires_df[which.min(clmfires$y),]
ref_coord_max_y <- lon_lat_LaMancha[which.max(lon_lat_LaMancha$utm_y),]
ref_fire_max_y <- clmfires_df[which.max(clmfires$y),]

diff_x <- (ref_coord_min_x$utm_x - ref_fire_min_x$x +
          ref_coord_max_x$utm_x - ref_fire_max_x$x +
          ref_coord_min_y$utm_x - ref_fire_min_y$x +
          ref_coord_max_y$utm_x - ref_fire_max_y$x )/4

diff_y <- (ref_coord_min_y$utm_y - ref_fire_min_y$y +
             ref_coord_max_y$utm_y - ref_fire_max_y$y +
             ref_coord_min_x$utm_y - ref_fire_min_x$y +
             ref_coord_max_x$utm_y - ref_fire_max_x$y )/4

lon_lat_LaMancha[,"utm_x_scaled"] <- lon_lat_LaMancha$utm_x - diff_x
lon_lat_LaMancha[,"utm_y_scaled"] <- lon_lat_LaMancha$utm_y - diff_y

## visual check - looks good enough
ggplot(data=clmfires_df, aes(x = x, y = y)) +
  geom_point()+
  geom_point(data=lon_lat_LaMancha, aes(x = utm_x_scaled, y = utm_y_scaled),
             color = "red")+
  geom_polygon(data = LaMancha_Outline, aes(x = x, y = y), fill = NA, color = "black" )
  
## match each point to closest administrative infos
infos <- c("Provincia", "Población")
for (i in 1: nrow(clmfires_df)) {
  dist <- sqrt((clmfires_df$x[i] - lon_lat_LaMancha$utm_x_scaled)^2 + 
               (clmfires_df$y[i] - lon_lat_LaMancha$utm_y_scaled)^2)
  clmfires_df[i,infos] <- lon_lat_LaMancha[which.min(dist), infos] 
}

ggplot(clmfires_df, aes(x = x, y = y, color = Provincia))+
  geom_point()+
  geom_polygon(data = LaMancha_Outline, aes(x = x, y = y), fill = NA, color = "black" ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = c(0.1,0.8), legend.text = element_text(size = 16), legend.title = element_text(size = 18),
        aspect.ratio = 1)+
  guides(color = guide_legend(title = "Provincia"))


# Preliminary Analysis ---------------------------------------------------------

## Visualization for domination of low burnt area
clmfires_df$burnt.area.cat1ha <- ifelse(clmfires_df$burnt.area < 1, "Below 1", "Above 1")
clmfires_df$burnt.area.cat10ha <- ifelse(clmfires_df$burnt.area < 10, "Below 10", "Above 10")

sum(clmfires_df$burnt.area < 1)/length(clmfires_df$burnt.area)
sum(clmfires_df$burnt.area < 10)/length(clmfires_df$burnt.area)

ggplot(clmfires_df, aes(x=x, y=y, color=burnt.area.cat10ha))+
  geom_polygon(data = LaMancha_Outline, aes(x = x, y = y), fill = NA, color = "black" ) +
  geom_point(alpha=0.3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = c(0.2,0.85), legend.text=element_text(size = 16), legend.title=element_text(size = 18),
        aspect.ratio = 1)+
  guides(color = guide_legend(title = "Burned Land (ha)", override.aes = list(size = 4)))


ggplot(clmfires_df, aes(x=x, y=y, color=burnt.area.cat1ha))+
  geom_polygon(data = LaMancha_Outline, aes(x = x, y = y), fill = NA, color = "black" ) +
  geom_point(alpha=0.3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = c(0.2,0.85), legend.text = element_text(size = 16), legend.title=element_text(size = 18),
        aspect.ratio = 1)+
  guides(color = guide_legend(title = "Burned Land (ha)", override.aes = list(size = 4)))


## Burnt Area Histogram
ggplot(clmfires_df, aes(x = burnt.area))+
  geom_histogram(binwidth = 1, fill = "#660000")+
  xlim(0,100)+
  ylim(0,1500)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))+
  xlab("Burned Land (capped at 100)")

ggplot(clmfires_df, aes(x = burnt.area))+
  geom_histogram(binwidth = 1, fill = "#660000")+
  xlim(2,100)+
  ylim(0,450)+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))+
  xlab("Burned Land (capped at 2 and 100)")

## How many fires above 100ha ? above 1000ha ?
n_ab100 <- sum(clmfires_df$burnt.area > 100)
n_ab1000 <- sum(clmfires_df$burnt.area > 1000)
n_tot <- length(clmfires_df$burnt.area)

print(paste(n_ab100, "fires burned more than 100ha, (",round(n_ab100/n_tot*100,2) ,"% of all fires)"))
print(paste(n_ab1000, "fires burned more than 1000ha, (",round(n_ab1000/n_tot*100,2) ,"% of all fires)"))


## Univariate modeling or burned land

### Empirical Quantiles
quantiles <- c(0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99) 
quantiles_burnt_area <- round(quantile(clmfires_df$burnt.area, probs = quantiles), 2)
quantiles_df <- data.frame(quantile = quantiles, value = quantiles_burnt_area)

ggplot(quantiles_df, aes(x = quantile, y = value))+
  geom_point()+
  geom_line()+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))+
  xlab("Quantile")+
  ylab("Burned Land")+
  scale_x_continuous(breaks = c(0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99) )+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 14),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 14))

### GP fit using 60 % quantile as threshold
q_60 <- quantile(clmfires_df$burnt.area , p = 0.6)
print(q_60)
GP_fitted <- fevd(x = burnt.area, data = clmfires_df, type="GP", threshold = q_60)
plot(GP_fitted, type = "probprob", main = "Probability Plot")
plot(GP_fitted, type = "qq", main = " Quantile Plot")
plot(GP_fitted, type = "qq", main = " Quantile Plot", ylim = c(0,100), xlim = c(0,100))

### GP fit using 80 % quantile as threshold
q_80 <- quantile(clmfires_df$burnt.area , p = 0.8)
print(q_80)
GP_fitted <- fevd(x = burnt.area, data = clmfires_df, type="GP", threshold = q_80)
plot(GP_fitted, type = "probprob", main = "Probability Plot")
plot(GP_fitted, type = "qq", main = " Quantile Plot")
plot(GP_fitted, type = "qq", main = " Quantile Plot", ylim = c(0,100), xlim = c(0,100))

### GP fit using 95 % quantile as threshold
q_95 <- quantile(clmfires_df$burnt.area , p = 0.95)
print(q_95)
GP_fitted <- fevd(x = burnt.area, data = clmfires_df, type="GP", threshold = q_95)
plot(GP_fitted, type = "probprob", main = "Probability Plot")
plot(GP_fitted, type = "qq", main = " Quantile Plot")
plot(GP_fitted, type = "qq", main = " Quantile Plot", ylim = c(0,100), xlim = c(0,100))

## Estimating return levels 
daily_aggregate <- as.data.frame(clmfires_df %>% group_by(date)  %>% summarize(total_burnt_area = sum(burnt.area)))
full_dates <- data.frame(date = seq(min(daily_aggregate$date), max(daily_aggregate$date), by = "day"))

complete_data <- full_dates %>% left_join(daily_aggregate, by = "date") %>%
                 mutate(total_burnt_area = replace_na(total_burnt_area, 0))

q_95_daily <- quantile(daily_aggregate$total_burnt_area, p = 0.95)
daily_fit<- fevd(daily_aggregate$total_burnt_area, time.units = "days", type = "GP", threshold = q_95_daily)

plot(daily_fit, type = "probprob", main = "Probability Plot")
plot(daily_fit, type = "qq", main = " Quantile Plot")

return_levels <- return.level(daily_fit, return.period = seq(2,10,1))
return_levels_df <- data.frame(return_period = seq(2,10), return_levels = as.vector(return_levels))

ggplot(data = return_levels_df, aes(x = return_period, y = return_levels))+
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(axis.title.y = element_text(margin = margin(r = 10))) +
  labs(x = "Return Period (Years)", y = "Burned Land (ha)")

## Fire Activity Yearly Development
clmfires_df$year <- year(clmfires_df$date)

yearly_aggregate <- as.data.frame(clmfires_df %>% group_by(year)  %>% summarize(total_burnt_area = sum(burnt.area)))
yearly_aggregate$count <- as.data.frame(clmfires_df %>% group_by(year) %>% summarize(count = length(burnt.area)))$count

ggplot(yearly_aggregate, aes(x = year))+
  geom_line(aes(y = total_burnt_area, color = "Total Burned Land (ha)"))+
  geom_point(aes(y = total_burnt_area, color = "Total Burned Land (ha)"))+
  geom_line(aes(y = count*10, color = "Fire Count"))+
  geom_point(aes(y = count*10, color = "Fire Count"))+
  scale_y_continuous(name = "Total Burned Area", sec.axis = sec_axis(trans=~./10, name = "Fire Count")) +
  theme_minimal()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 10), size = 16),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 16),
        axis.title.y.right =  element_text(margin = margin(t = 0, r = 10, b = 0, l = 20)),
        axis.text = element_text(size = 14), 
        legend.position=c(0.2,0.9), legend.background = element_rect(fill = "white", color = "lightgrey"),
        legend.title = element_blank())+
  scale_x_continuous(breaks = seq(1998, 2007, 1))


## Fire Activity Montly Development
clmfires_df$month <- month(clmfires_df$date)
monthly_aggregate <- as.data.frame(clmfires_df %>% group_by(month) %>% summarize(count = length(burnt.area)))

ggplot(monthly_aggregate, aes(x = month, y = count))+
  geom_line()+
  geom_point()+
  theme_minimal()+
  ylab("Fire Count")+
  xlab("Month")+
  scale_x_continuous(breaks = seq(1, 12, 1))+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

# Graphical Modeling -----------------------------------------------------------

## Pre processing
### Focus on northern region
clmfires_north <- subset(clmfires, x >= 150 & y >= 250)
plot(clmfires_north, use.marks = FALSE, pch = 16, cex = 0.5, main = NULL)

### Get data subset and smooth spatially using 10x10km grid cell
clmfires_north_df <- clmfires_df %>% filter(clmfires_df$x >= 150 & clmfires_df$y >= 250)

clmfires_north_df$x_rounded <- round(clmfires_north_df$x/10)*10
clmfires_north_df$y_rounded <- round(clmfires_north_df$y/10)*10

## visualize the smoothed grid
ggplot(data = clmfires_north_df, aes(x = x_rounded, y = y_rounded))+
  geom_point()+
  theme_minimal()+
  labs(x = " Grid Y", y = "Grid Y")

relevant_columns <- c("x_rounded", "y_rounded", "date", "burnt.area", "Población")
clmfires_north_df <- clmfires_north_df[, relevant_columns]

### get an ID for each unique location
clmfires_north_df_grouped <- clmfires_north_df %>% group_by(x_rounded, y_rounded) %>%
                              mutate(ID = paste0("ID_", cur_group_id())) 

clmfires_north_df <- ungroup(clmfires_north_df_grouped)

### complete date range for each of the locations (inclusion of 0 events)
date_range <- seq(min(clmfires_north_df$date), max(clmfires_north_df$date), by = "day")

id_list <- unique(clmfires_north_df$ID)
print(paste("In this subset, there are", length(id_list), "possible locations"))
output_df <- data.frame(date = date_range)
output_df[ , id_list] <- 0 # init data matrix with zeros

for (i in 1:nrow(clmfires_north_df)) {
  id_i <- clmfires_north_df$ID[i]
  date_i <- clmfires_north_df$date[i]
  burnt.area_i <- clmfires_north_df$burnt.area[i]
  
  date_idx <- which(output_df$date == date_i)
  output_df[date_idx, id_i] = output_df[date_idx, id_i] + burnt.area_i
}

## aggregate the data monthly
output_df$date <- as.Date(output_df$date)
output_df <- output_df %>% mutate(month = floor_date(date, unit = "month"))
monthly_output_df <- output_df %>% group_by(month) %>% summarise(across(starts_with("ID"), sum, na.rm = TRUE), .groups = "drop")

## Use only fire seasons -> from may to october   
fire_season_monthly_output_df <- monthly_output_df %>% filter((month(month) <= 10 & month(month) >= 5))    

## create ID_db:
ID_db <- data.frame(ID = names(fire_season_monthly_output_df)[-1],
                    gridX = 0, gridY = 0)

ID_db <- unique(clmfires_north_df[,c("x_rounded", "y_rounded", "ID")])
colnames(ID_db) <- c("gridX", "gridY", "ID")

## Number of fires across the locations
fire_count <- apply(fire_season_monthly_output_df[,-1], 2, function(x) {sum( x != 0)})

ggplot(as.data.frame(fire_count), aes(x = fire_count)) +
  geom_bar(fill = "steelblue") +
  labs(x = "Number of months with fires in location (capped at 40)", 
       y = "Number of locations") +
  theme_minimal()+
  xlim(0,40)

## Focus on a subset of locations that have a quantile values > 1
ID_80 <- apply(fire_season_monthly_output_df[,-1], 2, quantile, 0.8)
ID_90 <- apply(fire_season_monthly_output_df[,-1], 2, quantile, 0.9)

ID_80_ab1 <- names(ID_80[ID_80 > 1]) # 4 locations
ID_90_ab1 <- names(ID_90[ID_90 > 1]) # 17 locations

## Graphical modeling 
ID_subset <- ID_90_ab1
q_th <- 0.90
X <- as.data.frame(fire_season_monthly_output_df)[,ID_subset]
Y <- data2mpareto(X, p = q_th)
colnames(Y) <- colnames(X)

dim(X)
dim(Y)

rholist <- seq(0.01,0.11,0.01)
criterion <- "bic"
eglearn_fit <- eglearn(Y, rholist = rholist, reg_method = "ns", complete_Gamma = TRUE)

criterion_eglearn_fit <- data.frame(rholist = rholist,
                                    values = NA,
                                    graph_size = NA)
for (i in 1:length(rholist)) {
  Gamma <- eglearn_fit$Gamma[[i]]
  Graph <- eglearn_fit$graph[[i]]
  if (is.null(Gamma)) {
    criterion_eglearn_fit[i, "values"] <- NA
    criterion_eglearn_fit[i, "graph_size"] <- gsize(Graph)
  }
  else {
    criterion_eglearn_fit[i, "values"] <- loglik_HR(data = Y, Gamma = Gamma, graph = Graph)[criterion]
    criterion_eglearn_fit[i, "graph_size"] <- gsize(Graph)
  }
}

ggplot(data = criterion_eglearn_fit, aes(x = rholist, y = values)) +
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
colnames(theta_k) <- colnames(X)
rownames(theta_k) <- colnames(X) 

precision_matrix <- theta_k
chi_matrix <- Gamma2chi(best_Gamma)

g_res_chi <- graph_from_pm(precision_matrix = precision_matrix, Chi_matrix = chi_matrix, ID_db=ID_db)
g_res_chi + theme(aspect.ratio = 1, 
                  axis.text = element_text(size = 16), axis.title = element_text(size = 16),
                  legend.text = element_text(size = 16), legend.title = element_text(size = 16),
                  title = element_text(size = 14)) +
            scale_alpha_continuous(range=c(0.1,1), guide = "none")

### Determine a max distance
#### Kest function
kest_north <- Kest(clmfires_north, correction = "best") # using best estimate according to documentation
plot(kest_north, main = NULL)
max_dist_north <- round(max(kest_north$r),1)

kest_all <- Kest(clmfires, correction = "best")
plot(kest_all, main = NULL)
max_dist_all <- round(max(kest_all$r),1)

### Scale Precision Matrix
max_dist <- max_dist_north
precision_scaled <- scale_precision_matrix(precision_matrix, ID_db, dist_max = max_dist)

g_res_scaled_chi <- graph_from_pm(precision_scaled, Chi_matrix = chi_matrix, ID_db = ID_db)
g_res_scaled_chi + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        title = element_text(size = 14))+
  labs(title = paste0("Graph from Scaled Precision Matrix (max dist = ", max_dist, ")")) +
  scale_alpha_continuous(range=c(0.1,1), guide = "none")


### Match location to Poblacion
ID_db[, "Población"] <- apply(ID_db, 1, function(x) { get_poblacion_from_ID(x[3]) })

ID_db_gm <- ID_db[ID_db$ID %in% colnames(X),]
dim(ID_db_gm)

# the following list of colors was made using the table in: https://www.sthda.com/english/wiki/colors-in-r
# This website provides with color codes taken originally from http://www.visibone.com
list_of_colors <- c("#FFFFFF", "#000000", "#FFFF00", "#CC0033", "#66FF33",
                    "#FF99FF", "#FF3300", "#9933FF", "#0033FF", "9933FF",
                    "#006600", "#999999", "#663333", "#FF9933", "#99CCCC",
                    "#FFCC99", "#333333")

g_res_scaled_chi +
  geom_point(data = ID_db, aes(x = gridX, y = gridY), color="grey")+
  geom_point(data = ID_db_gm, aes(x = gridX, y = gridY, fill = Población),
             shape = 21, size = 4)+
  scale_fill_manual(values = list_of_colors, name = "Población")+
  theme(legend.position = "right",
        axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        title = element_text(size = 14))+
  guides(fill = guide_legend(ncol = 2))+
  labs(title = paste("Graph from Scaled Precision Matrix (max dist =", max_dist, ")"))+
  xlim(xmin, xmax) +
  ylim(ymin, ymax)


## Plausibility of the Edges ?

### Huete - Barriopedro
ID_huete <- ID_db_cluster_gm[ID_db_cluster_gm$Población == "Huete", ]$ID
ID_barri <- ID_db_cluster_gm[ID_db_cluster_gm$Población == "Barriopedro", ]$ID

ggplot(X, aes(x = ID_76, y = ID_70))+
  geom_point(color = "blue", size = 2, alpha = 0.7)+
  labs(x = "Huete",
       y = "Barriopedro",
       title = "Monthly Burned Land (ha)")+
  theme_minimal()+
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 14))

### Casar (El) - Mohernando
ID_casar <- ID_db_cluster_gm[ID_db_cluster_gm$Población == "Casar (El)", ]$ID
ID_moher <- ID_db_cluster_gm[ID_db_cluster_gm$Población == "Mohernando", ]$ID

ggplot(X, aes(x = ID_8, y = ID_30))+
  geom_point(color = "blue", size = 2, alpha = 0.7)+
  labs(x = "Casar (El)",
       y = "Mohernando",
       title = "Monthly Burned Land (ha) - capped at 30")+
  theme_minimal()+
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 14))+
  xlim(0,30)+
  ylim(0,30)


