#plot site locations
library(ggmap)
#library(readxl)
library(tidyverse)
library(vegan)
library(rgdal)
library(lubridate)

#get satellite map function
local_map <- function(x,y, buffer) {
  e1 <- c(min(x, na.rm = TRUE) - buffer, min(y, na.rm = TRUE) - buffer, 
          max(x, na.rm = TRUE) + buffer, max(y, na.rm = TRUE) + buffer)
  get_map(location = c(e1), source = "google", maptype = "satellite")
}

#Bring in data and fix coordinates
my.df <- read_csv('Chandeleur_LA_Seagrass Locations.csv')
range(my.df$latitude, na.rm = T); range(my.df$longitude, na.rm = T) #check the ranges                        
my.df$latitude[my.df$latitude == 20.85865] <- 29.85865              #Spot fix
my.df$longitude[my.df$longitude == 99.93511] <- 88.93511
range(my.df$latitude, na.rm = T); range(my.df$longitude, na.rm = T)

my.df$longitude <- as.numeric(my.df$longitude); my.df$latitude <- as.numeric(my.df$latitude)       #make sure class=numeric            #make lats/longs numeric
my.df$longitude <- ifelse(my.df$longitude > 0, my.df$longitude * -1, my.df$longitude)              #correct hemisphere


#get satellite map of the area
satmap <- local_map(my.df[,"longitude"], my.df[,"latitude"], 0.1) #Set buffer range (degrees) with number value

# my.df$species.present <- NA
# my.df$species.present <- ifelse(my.df$`species present` %in% c('hal', 'hw', 'hw and he', 'hw and rm', 'hw and sf', 'hw or rm'), 'H. wrightii', NA)
# my.df$species.present <- ifelse(my.df$`species present` %in% c('hw and he', 'hw and rm', 'hw and sf', 'hw or rm',
#                                                                'sf hw'), 'mixed', my.df$species.present)

#plot with satellite map
#just points
ggmap(satmap,
      extent = "device", # "panel" keeps in axes, etc.
      ylab = "Latitude",
      xlab = "Longitude",
      legend = "right") +
  geom_point(data = my.df,
             aes(x = longitude, y = latitude), colour = 'orange', shape = 4) 


#get satellite map of the area
satmap <- local_map(my.df[,"longitude"], my.df[,"latitude"], 0.25) #Set buffer range (degrees) with number value

#plot with satellite map
#just points
ggmap(satmap,
      extent = "device", # "panel" keeps in axes, etc.
      ylab = "Latitude",
      xlab = "Longitude",
      legend = "right") +
  geom_point(data = my.df,
             aes(x = longitude, y = latitude), colour = 'orange', shape = 4) 

#color by species
ggmap(satmap,
      extent = "device", # "panel" keeps in axes, etc.
      ylab = "Latitude",
      xlab = "Longitude",
      legend = "right") +
  geom_point(data = my.df,
             aes(x = longitude, y = latitude, colour = `species present`)) 