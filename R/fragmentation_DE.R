# Marsh fragmentation in lower Terrebonne Bay
# Adam A. Kemberling
# 6/13/2018

#--------------------------------------------------------------------------------------
#---------------------------------------- load packages -------------------------------
#--------------------------------------------------------------------------------------
#library(conflicted)
library(readxl)
library(tidyverse)
library(vegan)
library(rgdal)
library(ggThemeAssist)
library(extrafont)
library(car)
library(lme4)
library(MuMIn)

#font setup
#font_import()
#loadfonts(device="pdf")       #Register fonts for pdf bitmap output
#loadfonts(device="win")       #Register fonts for Windows bitmap output
fonts() 

#test push
#test2


#--------------------------------------------------------------------------------------
#----------------------------------- Import datasets ----------------------------------
#--------------------------------------------------------------------------------------
#Move up one level to data sources
setwd('..')


boga <- read_xlsx("Bogaert calculation.xlsx")
frag50 <- read_xlsx("Frag at 50m and SAV likelihood.xlsx")
marsh <- read_xlsx("noaa_blu_objective1_marsh_survey.xlsx")
water <- read_xlsx("noaa_blu_objective1_openwater_survey.xlsx")
frag_metrics <- read_xlsx("Site level Frag_Metrics.xlsx")
full.t <- read_csv("all_transect_data.csv")

#throws <- read_xlsx("throw_trap_pts.dbf.xlsx") #contains midpoint coordinates used to create high res boxes
#sites <- readOGR(dsn = "~/Users/adamkemberling/Dropbox (The Craboratory)/Terrebonne Fragmentation paper/Data/noaa_blu_basecamp.kml", layer = "noaa_blu_basecamp.kml")


#--------------------------------------------------------------------------------------
#--------- Pairing crms data with marsh and open water data by sampling time ----------
#--------------------------------------------------------------------------------------
crms311 <- read.csv("HYDROGRAPHIC_HOURLY_311.csv")
crms345 <- read.csv("HYDROGRAPHIC_HOURLY_345.csv")
crms369 <- read.csv("HYDROGRAPHIC_HOURLY_369.csv")
names(crms311)

crms <- bind_rows(crms311, crms345)
crms <- bind_rows(crms, crms369)
rm(crms311); rm(crms345); rm(crms369)

table(crms$Station.ID)

str(crms)
crms <- crms %>% select(Station.ID:Adjusted.Battery..V., Organization.Name, Comments, Latitude, Longitude)

library(lubridate)

crms <- crms %>%  mutate(Station.ID = str_sub(Station.ID, 6, 8),
                         Date = mdy(Date),
                         Time = Time,
                         temp_raw = Raw.Water.Temperature,
                         temp_adj = Adjusted.Water.Temperature,
                         conductance_raw = Raw.Specific.Conductance..uS.cm.,
                         conductance_adj = Adjusted.Specific.Conductance..uS.cm.,
                         salinity_raw = Raw.Salinity..ppt.,
                         salinity_adj = Adjusted.Salinity..ppt.,
                         water_lvl_raw = Raw.Water.Level..ft.,
                         water_lvl_adj = Adjusted.Water.Level..ft.,
                         marsh_mat_elev_raw = Raw.Marsh.Mat.Elevation..ft.,
                         marsh_mat_elev_adj = Adjusted.Marsh.Mat.Elevation.to.Datum..ft.) %>% 
  dplyr::select(Station.ID, Date, Time,Latitude, Longitude, temp_raw, temp_adj, conductance_raw, conductance_adj, salinity_raw, salinity_adj,
                water_lvl_raw, water_lvl_adj, marsh_mat_elev_raw, marsh_mat_elev_adj, Organization.Name, Comments)

str(crms)
#write_csv(crms, "crms_three_stations.csv")

#where are they
crms %>%
  group_by(Station.ID) %>%
  slice(1) %>% 
  ggplot(aes(Longitude, Latitude, label = Station.ID)) +  
  geom_point() +
  geom_label()


#--------------------------------------------------------------------------------------
#-------------------------- Marsh Survey Dataset "marsh" ------------------------------
#--------------------------------------------------------------------------------------
# Objectives:
'So with this dataset:

Possible Response Variables:
Richness
Diversity Indices

Fixed Effects:
Depth
Salinity
Distance From Edge
*Accompanying CRMS Data, ?What time frame

Random Effects:
Month
Fragmentation
Site
Station'

#--------------------------------------------------------------------------------------
#---------------------------- Marsh Data Manipulation ---------------------------------
#--------------------------------------------------------------------------------------
marsh <- read_xlsx("noaa_blu_objective1_marsh_survey.xlsx")

#make any class/formatting changes
str(marsh)

marsh <- marsh %>% mutate(
  correct_elev_m = as.numeric(correct_elev_m),
  depth_cm = as.numeric(depth_cm),
  crms_site = factor(crms_site, levels = c("345", "311", "369")),
  month = ifelse(lubridate::month(date) == 4, "April", "September"),
  perc_bare = rep(0, nrow(marsh)),
  frag = factor(toupper(frag), levels = c("L","M","H"))) %>% 
  filter(month == "April") 

#for loop
for (i in 1:nrow(marsh)){
  marsh$perc_bare[i] <- 100 - sum(marsh[i,12:34], na.rm = TRUE)
}

#map function


#--------------------------------------------------------------------------------------
#----------------------  Match bogaert scale to the subsites --------------------------
#--------------------------------------------------------------------------------------
#low resolution bogaert values, subsite level
frag_metrics <- read_xlsx("Site level Frag_Metrics.xlsx")
frag_metrics %>% group_by(Site) %>% summarise(number = n())

#make id to join by site and frag level
frag_metrics <- frag_metrics %>% mutate(
  crms_site = Site,
  crms_site = ifelse(crms_site == "C", 311, crms_site),
  crms_site = ifelse(crms_site == "N", 369, crms_site),
  crms_site = ifelse(crms_site == "S", 345, crms_site),
  frag = `frag cat`,
  subsite_id = str_c(crms_site, frag, sep = "-"))

#join them to marsh dataset
marsh <- marsh %>% mutate(subsite_id = str_c(crms_site, frag, sep = "-")) %>% 
  left_join(., dplyr::select(frag_metrics, subsite_id, bogaert), by = "subsite_id") 




#Top Species Counts
diversity.df <- read_csv("all_transect_data.csv")
full.transects <- diversity.df %>% filter(dist_from_edge_m == "full_transect")
full.transects[,6:28] <- ifelse(full.transects[,6:28] > 0, 1, 0)
full.transects <- full.transects %>% summarise(Sa = (sum(Sa)/ 189) * 100,
                             Jr = (sum(Jr)/ 189) * 100,
                             Ag = (sum(Ag)/ 189) * 100,
                             Sp = (sum(Sp)/ 189) * 100,
                             Ds = (sum(Ds)/ 189) * 100,
                             Bm = (sum(Bm)/ 189) * 100,
                             Sb = (sum(Sb)/ 189) * 100,
                             Pv = (sum(Pv)/ 189) * 100,
                             Is = (sum(Is)/ 189) * 100,
                             St = (sum(St)/ 189) * 100,
                             Sr = (sum(Sr)/ 189) * 100,
                             Sc = (sum(Sc)/ 189) * 100,
                             Zm = (sum(Zm)/ 189) * 100,
                             Ep = (sum(Ep)/ 189) * 100,
                             Aa = (sum(Aa)/ 189) * 100,
                             Pa = (sum(Pa)/ 189) * 100,
                             Na = (sum(Na)/ 189) * 100,
                             If = (sum(If)/ 189) * 100,
                             Bf = (sum(Bf)/ 189) * 100,
                             Vl = (sum(Vl)/ 189) * 100,
                             's?' = (sum(`s?`)/ 189) * 100,
                             pavi = (sum(pavi)/ 189) * 100,
                             ScAm = (sum(ScAm)/ 189) * 100)

species_freq <- t(full.transects)
species_freq <- data.frame("Species" = rownames(species_freq), "Percent_of_Total_Transects" = as.numeric(species_freq[,1]))
species_freq <- species_freq %>% arrange(desc(Percent_of_Total_Transects)) %>% mutate(Number_of_Transects = (Percent_of_Total_Transects/100) * 189)
View(species_freq)



#--------------------------------------------------------------------------------------
#----------------------- vegan package diversity indices ------------------------------
#--------------------------------------------------------------------------------------
#percent cover and species encounter matrices
pcov.mat <- as.matrix(marsh[,12:34], dimnames = list(marsh$renamed_station,colnames(marsh)), rownames.force = TRUE)   #make a matrix of the percent cover values
bin.mat <- ifelse(pcov.mat > 0, 1, 0) #binary species encounter matrix


#remove empty sites with:
pcov.mat <- pcov.mat[rowSums(pcov.mat)!= 0,]
bin.mat <- bin.mat[rowSums(bin.mat)!= 0,]


#standardize them?
t1.mat <- decostand(pcov.mat, method = "max")

#Species richness
marsh$Richness <- specnumber(marsh[,12:34])
marsh$Shannon <- diversity(marsh[,12:34], index =  "shannon")
marsh$Simpson <- diversity(marsh[,12:34], index =  "simpson")
marsh$Evenness <- marsh$Shannon/log(marsh$Richness)         #Pielou's evenness




#--------------------------------------------------------------------------------------
#-------------------------- species encounter summaries -------------------------------
#--------------------------------------------------------------------------------------
'
Marsh plants

Most abundant marsh plant species overall
Most abundant marsh plant species at each site
Most abundant marsh plant species at each Bogaert value
Number of marsh plant species overall
# marsh plant species at each site
# marsh plant species at each Bogaert


Table listing marsh plant species found by site and by Bogaert index.
'
#use the matrices that include the 2 bare sites
pcov.mat.full <- as.matrix(marsh[,12:34], dimnames = list(marsh$renamed_station,colnames(marsh)), rownames.force = TRUE)   #make a matrix of the percent cover values
bin.mat.full <- ifelse(pcov.mat.full > 0, 1, 0) #binary species encounter matrix

#1. Most abundant marsh plant overall
colSums(bin.mat.full)


library(magrittr)
as.data.frame(bin.mat.full) %>% dplyr::group_by(crms_site) %>% summarise(
   Sa = sum(pcov_sa),
   Jr = sum(pcov_jr),
   Ag = sum(pcov_ag),
   Sp = sum(pcov_sp),
   Ds = sum(pcov_ds),
   Bm = sum(pcov_bm),
   Sb = sum(pcov_sb),
   Pv = sum(pcov_pv),
   Is = sum(pcov_is),
   St = sum(pcov_st),
   Sr = sum(pcov_sr),
   Sc = sum(pcov_sc),
   Zm = sum(pcov_zm),
   Ep = sum(pcov_ep),
   Aa = sum(pcov_aa),
   Pa = sum(pcov_pa),
   Na = sum(pcov_na),
   If = sum(pcov_if),
   Bf = sum(pcov_bf),
   Vl = sum(pcov_vl),
   pavi = sum(pcov_pavi),
   ScAm = (sum(pcov_ScAm) + sum(`pcov_s?`))
)



t1 <- sort(colSums(bin.mat.full), decreasing = TRUE)

t1.df <- data.frame(species = names(t1), encounters = as.numeric(t1))


#2. By Site
colSums(bin.mat.full[marsh$crms_site == "311",])
colSums(bin.mat.full[marsh$crms_site == "345",])
colSums(bin.mat.full[marsh$crms_site == "369",])

#3. By Bogaert Value
binary.df <- marsh %>% dplyr::select(crms_site, subsite_id,bogaert) %>% bind_cols(., as.data.frame(bin.mat.full))
binary.df %>%  group_by(bogaert) %>%  
  summarise(Sa = sum(pcov_sa),
            Jr = sum(pcov_jr),
            Ag = sum(pcov_ag),
            Sp = sum(pcov_sp),
            Ds = sum(pcov_ds),
            Bm = sum(pcov_bm),
            Sb = sum(pcov_sb),
            Pv = sum(pcov_pv),
            Is = sum(pcov_is),
            St = sum(pcov_st),
            Sr = sum(pcov_sr),
            Sc = sum(pcov_sc),
            Zm = sum(pcov_zm),
            Ep = sum(pcov_ep),
            Aa = sum(pcov_aa),
            Pa = sum(pcov_pa),
            Na = sum(pcov_na),
            If = sum(pcov_if),
            Bf = sum(pcov_bf),
            Vl = sum(pcov_vl),
            's?' = sum(`pcov_s?`),
            pavi = sum(pcov_pavi),
            ScAm = sum(pcov_ScAm))

#maybe use this table for ordination plotting by frag


#4. number of species overall
n <- 23

#5. number of species by site
binary.df %>%  group_by(crms_site) %>%  
  summarise(Sa = sum(pcov_sa),
            Jr = sum(pcov_jr),
            Ag = sum(pcov_ag),
            Sp = sum(pcov_sp),
            Ds = sum(pcov_ds),
            Bm = sum(pcov_bm),
            Sb = sum(pcov_sb),
            Pv = sum(pcov_pv),
            Is = sum(pcov_is),
            St = sum(pcov_st),
            Sr = sum(pcov_sr),
            Sc = sum(pcov_sc),
            Zm = sum(pcov_zm),
            Ep = sum(pcov_ep),
            Aa = sum(pcov_aa),
            Pa = sum(pcov_pa),
            Na = sum(pcov_na),
            If = sum(pcov_if),
            Bf = sum(pcov_bf),
            Vl = sum(pcov_vl),
            's?' = sum(`pcov_s?`),
            pavi = sum(pcov_pavi),
            ScAm = sum(pcov_ScAm))
data.frame(crms_site = c(345,311,369), spec_obs = c(4, 8, 20))

#6. number of species by bogaert, use table from number 3.
data.frame(bogaert = c(sort(unique(marsh$bogaert))), spec_obs = c(3, 3, 12, 2, 4, 2, 8, 10, 15))

table(marsh$subsite_id, marsh$bogaert)


#--------------------------------------------------------------------------------------
#-------------------------------- data exploration plots ------------------------------
#--------------------------------------------------------------------------------------


#Simpson's diversity
ggplot(marsh, aes(dist_from_edge_m, Simpson)) +
  geom_point() +
  facet_grid(~ frag) +
  ggtitle("Simpsons, DFE")

ggplot(marsh, aes(depth_cm, Simpson)) +
  geom_point() +
  facet_grid(~ frag) +
  ggtitle("Simpsons, depth")

ggplot(marsh, aes(correct_elev_m, Simpson)) +
  geom_point() +
  facet_grid(~ frag) +
  ggtitle("Simpsons, elevation")

#Shannon diversity
ggplot(marsh, aes(dist_from_edge_m, Shannon)) +
  geom_point() +
  facet_grid(~ frag) +
  geom_smooth(method='lm',formula= y ~ x)

ggplot(marsh, aes(depth_cm, Shannon)) +
  geom_point() 

ggplot(marsh, aes(correct_elev_m, Shannon)) +
  geom_point() 

#Species Richness
#fixed effects
ggplot(marsh, aes(dist_from_edge_m, Richness)) +
  geom_jitter(height = 0.1, width = 0.1) +
  facet_grid(~ crms_site) +
  xlab("Distance from Edge (m)") +
  ylab("# of Species")

ggplot(marsh, aes(depth_cm, Richness)) +
  geom_jitter(height = 0.1, width = 0.1) +
  facet_grid(~ crms_site) +
  xlab("Water Depth (cm)") +
  ylab("# of Species")

ggplot(marsh, aes(correct_elev_m, Richness)) +
  geom_jitter(height = 0.1, width = 0.1) +
  facet_grid(~ crms_site) +
  xlab("Elevation") +
  ylab("# of Species")


#Random Effects
ggplot(marsh, aes(Richness)) +
  geom_histogram() +
  facet_grid(~ crms_site + frag) 


#Sampling Period Variation
ggplot(data = marsh, aes(month, Richness)) +
  geom_boxplot() +
  facet_grid(~ crms_site)

#Within Site Variation
ggplot(data = marsh, aes(crms_site, Richness)) +
  geom_boxplot() +
  facet_grid(~ frag)


#three-factor boxplots
ggplot(data = marsh, aes(frag, Richness)) +
  geom_boxplot(aes(fill = frag)) +
  facet_wrap(~ month + crms_site) +
  coord_flip() +
  ggtitle("Species Richness")

ggplot(data = marsh, aes(frag, Shannon)) +
  geom_boxplot(aes(fill = frag)) +
  facet_wrap(~ month + crms_site) +
  coord_flip() +
  ggtitle("Shannon Diversity")

ggplot(data = marsh, aes(frag, Simpson)) +
  geom_boxplot(aes(fill = frag)) +
  facet_wrap(~ month + crms_site) +
  coord_flip() +
  ggtitle("Simpsons Diversity")

#plot site locations
library(ggmap)
marsh$longitude <- as.numeric(marsh$longitude); marsh$latitude <- as.numeric(marsh$latitude)                  #make lats/longs numeric
     
#get satellite map function
local_map <- function(x,y, buffer) {
  e1 <- c(min(x, na.rm = TRUE) - buffer, min(y, na.rm = TRUE) - buffer, 
          max(x, na.rm = TRUE) + buffer, max(y, na.rm = TRUE) + buffer)
  get_map(location = c(e1), source = "google", maptype = "satellite")
}

#maps
satmap <- local_map(marsh[,"longitude"], marsh[,"latitude"], 0.05)
marsh_311 <- local_map(marsh[marsh$crms_site == 311, "longitude"], marsh[marsh$crms_site == 311, "latitude"], 0.03)
marsh_345 <- local_map(marsh[marsh$crms_site == 345, "longitude"], marsh[marsh$crms_site == 345, "latitude"], 0.03)
marsh_369 <- local_map(marsh[marsh$crms_site == 369, "longitude"], marsh[marsh$crms_site == 369, "latitude"], 0.03)
site.names <- data.frame(site = c("0369", "0311", "0345" ), 
                         longitude = c(-90.71492, -90.76659, -90.73892),
                         latitude = c(29.29509, 29.21496, 29.18148))

#plot with satellite map
ggmap(satmap,
      extent = "device", # "panel" keeps in axes, etc.
      ylab = "Latitude",
      xlab = "Longitude",
      legend = "right") +
  geom_point(data = marsh,
             aes(x = longitude, y = latitude, colour = frag)) +
  geom_label(data = site.names, aes(x = longitude, y = latitude, label = site)) +
  ggtitle("All Sites")

ggmap(marsh_311,
      extent = "device", # "panel" keeps in axes, etc.
      ylab = "Latitude",
      xlab = "Longitude",
      legend = "right") +
  geom_point(data = marsh[marsh$crms_site == 311,],
             aes(x = longitude, y = latitude, colour = frag)) +
  ggtitle("CRMS311")

ggmap(marsh_345,
      extent = "device", # "panel" keeps in axes, etc.
      ylab = "Latitude",
      xlab = "Longitude",
      legend = "right") +
  geom_point(data = marsh[marsh$crms_site == 345,],
             aes(x = longitude, y = latitude, colour = frag)) +
  ggtitle("CRMS345")

ggmap(marsh_369,
      extent = "device", # "panel" keeps in axes, etc.
      ylab = "Latitude",
      xlab = "Longitude",
      legend = "right") +
  geom_point(data = marsh[marsh$crms_site == 369,],
             aes(x = longitude, y = latitude, shape = frag)) +
  ggtitle("CRMS369")
  

#--------------------------------------------------------------------------------------
#------------------------------  Site Level Ordinations  ------------------------------
#--------------------------------------------------------------------------------------
#Dataset with full transect totals

full.t <- read_csv("all_transect_data.csv")

#make full matrices that include the 2 bare sites, but not full transects
pcov.mat.full <- as.matrix(full.t[,6:28], dimnames = list(full.t$transect_id, colnames(full.t[,6:28])), rownames.force = TRUE) #make a matrix of the percent cover values
bin.mat.full <- ifelse(pcov.mat.full > 0, 1, 0) #binary species encounter matrix


#species accumulation curve
pcov.mat <- full.t[full.t$dist_from_edge_m == "full_transect", 6:28]
plot(specaccum(pcov.mat), xlab = "# of Transects", ylab = "# of species")

#new matrix done by subsite, removing the full transect totals
subsite.pc <- full.t[full.t$dist_from_edge_m != "full_transect",] %>%  group_by(subsite_id) %>%
  summarise(Sa = sum(Sa),
            Jr = sum(Jr),
            Ag = sum(Ag),
            Sp = sum(Sp),
            Ds = sum(Ds),
            Bm = sum(Bm),
            Sb = sum(Sb),
            Pv = sum(Pv),
            Is = sum(Is),
            St = sum(St),
            Sr = sum(Sr),
            Sc = sum(Sc),
            Zm = sum(Zm),
            Ep = sum(Ep),
            Aa = sum(Aa),
            Pa = sum(Pa),
            Na = sum(Na),
            If = sum(If),
            Bf = sum(Bf),
            Vl = sum(Vl),
            pavi = sum(pavi),
            ScAm = sum(ScAm, `s?`))

subsite.mat <- subsite.pc[,2:23]

metadata <- subsite.pc %>% select(subsite_id) %>% 
  mutate(crms_site = str_sub(subsite_id, 1 , 3))



#               compare species richness between subsites
boxplot(specnumber(pcov.mat.full) ~ full.t$subsite_id, ylab = "# of species")
boxplot(specnumber(pcov.mat.full) ~ full.t$crms_site, ylab = "# of species")
boxplot(full.t$Richness ~ full.t$subsite_id, ylab = "# of species") 


#             calculate Bray-Curtiscomm <- as.data.frame(subsite.pc)
rownames(subsite.mat) <- metadata$subsite_id
subsite.bc.dist <- vegdist(subsite.mat, method = "bray")

# cluster communities using average-linkage algorithm
subsite.bc.clust <- hclust(subsite.bc.dist, method = "average")
plot(subsite.bc.clust, ylab = "Bray-Curtis dissimilarity", las = 2)


# # plot cluster diagramby bogaert scale
# rownames(subsite.mat) <- metadata$bog_round
# bogaert.bc.clust <-  vegdist(subsite.mat, method = "bray")
# 
# png("dendrogram_bogaert.png", height = 4, width = 5, units = "in", res = 300)
# par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial", xaxt = "n", yaxt = "n")
# plot(bogaert.bc.clust, ylab = "Bray-Curtis dissimilarity", las = 2)
# dev.off()


#                            ordination
# The metaMDS function automatically transforms data and checks solution
# robustness
subsite.bc.mds <- metaMDS(subsite.mat, dist = "bray")
# Assess goodness of ordination fit (stress plot)
stressplot(subsite.bc.mds)
# plot site scores as text
ordiplot(subsite.bc.mds, display = "sites", type = "text")
# automated plotting of results - tries to eliminate overlapping labels
ordipointlabel(subsite.bc.mds)


#             ordination plots are highly customizable set up the plotting area but

# Set up plot dimensions
png("NMDS_noaxes.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
#, xaxt = "n", yaxt = "n"
mds.fig <- ordiplot(subsite.bc.mds, type = "none", xlim = c(-2,1.5), ylim = c(-1.5,1), las = 1, xaxt = "n", yaxt = "n")
# plot just the samples, colour by habitat, pch=19 means plot a circle
points(mds.fig, "sites", pch = 16, select = metadata$crms_site == "345")
points(mds.fig, "sites", pch = 17, select = metadata$crms_site == "311")
points(mds.fig, "sites", pch = 15, select = metadata$crms_site == "369")
# add confidence ellipses around habitat types
ordiellipse(subsite.bc.mds, metadata$crms_site, conf = 0.95, label = FALSE)
#tick marks
axis(side=1, at= seq(-2.0, 1.5, by = 0.5), tick = F, las = 1)
axis(side=2, at= seq(-1.5, 1.0, by = 0.5), tick = F, las = 1)

#Add text
text(c(-1.9, -0.15, 1.25), c(-0.4, 0.1, -0.9), labels = c("CRMS 0345", "CRMS 0311", "CRMS 0369"))
dev.off()


# Set up plot dimensions
png("NMDS_axes.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
#, xaxt = "n", yaxt = "n"
mds.fig <- ordiplot(subsite.bc.mds, type = "none", xlim = c(-1.7,1.4), ylim = c(-1.5,1), las = 1)
# plot just the samples, colour by habitat, pch=19 means plot a circle
points(mds.fig, "sites", pch = 16, select = metadata$crms_site == "345")
points(mds.fig, "sites", pch = 17, select = metadata$crms_site == "311")
points(mds.fig, "sites", pch = 15, select = metadata$crms_site == "369")
# add confidence ellipses around habitat types
ordiellipse(subsite.bc.mds, metadata$crms_site, conf = 0.95, label = FALSE)

#Add text
text(c(-1.9, -0.15, 1.25), c(-0.4, 0.1, -0.9), labels = c("CRMS 0345", "CRMS 0311", "CRMS 0369"))
dev.off()


#--------------------------------------------------------------------------------------
#----------------------   Transect level ordination plot  -----------------------------
#--------------------------------------------------------------------------------------

pcov.mat <- full.t[full.t$dist_from_edge_m == "full_transect", 6:28]


#--------------------------------------------------------------------------------------
#------------------  ANOSIM to verify differences quantitatively ----------------------
#--------------------------------------------------------------------------------------
pcov.mat <- full.t[full.t$dist_from_edge_m != "full_transect",]
pcov.mat <- pcov.mat[pcov.mat$Richness > 0, 6:28]
metadata <- full.t[full.t$dist_from_edge_m != "full_transect", ]
metadata <- metadata[metadata$Richness > 0, c(1:5,7:47)]


#ANOSIM testing for differences among sites using a bray-curtis distance matrix
metadata$crms_site <- factor(metadata$crms_site, levels = c('345', '311', '369'))

site_ano <-anosim(pcov.mat, metadata$crms_site, permutations=5000)
summary(site_ano)
dev.off()

#plot 
plot(site_ano)

png("ANOSIM_site.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
plot(site_ano, las = 1)
dev.off()

# values in boxplots are mean ranks within and between groups
tapply(site_ano$dis.rank, site_ano$class.vec, mean)

# Can calculate R with these values
mean_ranks<-tapply(site_ano$dis.rank, site_ano$class.vec,mean)
between <- mean_ranks[1]
within <- mean(mean_ranks[2:4])
(between-within)/((36*35)/4)


#plot a frequency distribution ofthe permuted R values
hist(site_ano$perm)


#--------------------------------------------------------------------------------------
#---------------------  ANOSIM bogaert fragmentation  ---------------------------------
#--------------------------------------------------------------------------------------
frag_metrics <- read_xlsx("Site level Frag_Metrics.xlsx")
frag_metrics <- frag_metrics %>% mutate(
  crms_site = Site,
  crms_site = ifelse(crms_site == "C", 311, crms_site),
  crms_site = ifelse(crms_site == "N", 369, crms_site),
  crms_site = ifelse(crms_site == "S", 345, crms_site),
  frag = `frag cat`,
  subsite_id = str_c(crms_site, frag, sep = "-"))


bog.df <- data.frame(subsite_id = unique(frag_metrics$subsite_id))
bog.df <- left_join(bog.df, frag_metrics[,c('subsite_id', 'bogaert')], by = "subsite_id")

metadata <- left_join(metadata, bog.df, by = "subsite_id")
metadata$bogaert <- factor(metadata$bogaert, levels = sort(unique(metadata$bogaert), decreasing =  T))
metadata$bog_round <- round(as.numeric(as.character(metadata$bogaert)))
metadata$bog_round <- factor(metadata$bog_round, levels = sort(unique(metadata$bog_round), decreasing =  T))

bog_ano <- anosim(pcov.mat, metadata$bogaert, permutations = 1000)
summary(bog_ano)

#plot
png("ANOSIM_bogaert.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
anosim_labels <- c("Between",c(as.character(sort(unique(round(as.numeric(as.character(metadata$bogaert)), digits = 2)), decreasing = T))))
plot(bog_ano, las = 1, xaxt = "n")
axis(1, at=1:10, labels=anosim_labels, tick = TRUE)
dev.off()

#plot a frequency distribution ofthe permuted R values
hist(bog_ano$perm)




#--------------------------------------------------------------------------------------
#----------------------- Water Survey Dataset "water" ---------------------------------
#--------------------------------------------------------------------------------------
water <- read_xlsx("noaa_blu_objective1_openwater_survey.xlsx")
water[water$frag == "i", "frag"] <- "m" 
water[water$frag == "m" & water$crms_site == 345,] #check that issue is fixed

# #high resolution bogaert values, pairs with throwtrap locations
# frag50 <- read_xlsx("Frag at 50m and SAV likelihood.xlsx")
# frag50 %>% group_by(site) %>% summarise(value.count = n())

#Check data structures
str(water)

#distance between sensors
sp <- .15

#make light attenuation measurement 
'Measuring Light Attenuation in Shallow Coastal Systems
Ana C. Brito1, Alice Newton, Teresa F. Fernandes, and Paul Tett'
water <- water %>% mutate(
  top_par1 = as.numeric(top_par1),
  top_par2 = as.numeric(top_par2),
  bot_par1 = as.numeric(bot_par1),
  bot_par2 = as.numeric(bot_par2),
  month = ifelse(lubridate::month(date) == 4, "April", "September"),
  top_avg = ((top_par1 + bot_par1)/2),                                      #bulb 1 is top bulb 
  bot_avg = ((bot_par2 + bot_par2)/2),                                      #bulb 2 is low bulb
  Kd_surface = -(log10(top_par2 / top_par1)) / sp,                          #Kd at surface, per meter 
  Kd_bottom = -(log10(bot_par2 / bot_par1)) / sp,                           #Kd at depth, per meter
  Kd_bot_test = Kd_bottom / (6 + (2/3)),                                    #attenuation per .15m
  Kd_surf_test = Kd_surface / (6 + (2/3)),
  frag = toupper(frag),
  subsite_id = str_c(crms_site, frag, sep = "-"),
  new_long = as.numeric(new_long),
  new_lat = as.numeric(new_lat),
  sal = as.numeric(sal),
  correct_elev_m = as.numeric(correct_elev_m)
)

#--------------------------------------------------------------------------------------
#----------------------  Match bogaert scale to the subsites --------------------------
#--------------------------------------------------------------------------------------
#low resolution bogaert values, subsite level
frag_metrics <- read_xlsx("Site level Frag_Metrics.xlsx")

#make id to join by site and frag level
frag_metrics <- frag_metrics %>% mutate(
  crms_site = Site,
  crms_site = ifelse(crms_site == "C", 311, crms_site),
  crms_site = ifelse(crms_site == "N", 369, crms_site),
  crms_site = ifelse(crms_site == "S", 345, crms_site),
  frag = `frag cat`,
  subsite_id = str_c(crms_site, frag, sep = "-"))

#join them to water dataset
water <- left_join(water, dplyr::select(frag_metrics, subsite_id, bogaert), by = "subsite_id") 


# #kd seems odd
# water %>% dplyr::select(renamed_station, top_par1, top_par2, bot_par1, bot_par2, Kd_surface, Kd_bottom) %>% 
#   dplyr::filter(top_par2 > top_par1) 
# water %>%  dplyr::select(renamed_station, top_par1, top_par2, bot_par1, bot_par2) %>% 
#   filter(bot_par2 > bot_par1)
#   
# #highlighing subsets
# p <- ggplot(water, aes(Kd_surface, Kd_bottom)) +
#   geom_point() +
#   geom_point(data = filter(water, top_par2 > top_par1), aes(Kd_surface, Kd_bottom), colour = "red") +
#   geom_point(data = filter(water, bot_par2 > bot_par1), aes(Kd_surface, Kd_bottom), colour = "blue")
# ggplotly(p)
# 
# 
# #attenuation per .15m, not 1m
# ggplot(water, aes(Kd_surf_test, Kd_bot_test)) 
#   geom_point()


#environmental variable sub-site summaries
water %>% filter(month == "April") %>% group_by(crms_site) %>% 
  summarise(sal_min = min(sal, na.rm = TRUE),
            sal_max = max(sal, na.rm = TRUE),
            sal_mean = mean(sal, na.rm = TRUE),
            elev_min = min(correct_elev_m, na.rm = TRUE),
            elev_max = max(correct_elev_m, na.rm = TRUE),
            elev_mean = mean(correct_elev_m, na.rm = TRUE))

water %>% filter(month == "April") %>% group_by(crms_site, frag) %>% 
  summarise(sal_min = min(sal, na.rm = TRUE),
            sal_max = max(sal, na.rm = TRUE),
            sal_mean = mean(sal, na.rm = TRUE),
            elev_min = min(correct_elev_m, na.rm = TRUE),
            elev_max = max(correct_elev_m, na.rm = TRUE),
            elev_mean = mean(correct_elev_m, na.rm = TRUE))


#Check locations of open water 311-H and 345-H to see if they were switched
#make lats/longs numeric

#get satellite map function
library(ggmap)
local_map <- function(x,y, buffer) {
  e1 <- c(min(x, na.rm = TRUE) - buffer, min(y, na.rm = TRUE) - buffer, 
          max(x, na.rm = TRUE) + buffer, max(y, na.rm = TRUE) + buffer)
  ggmap::get_map(location = c(e1), source = "google", maptype = "satellite")
}

water <- water[!is.na(water$new_long),]
water <- water[!is.na(water$new_lat),]

#maps
satmap <- local_map(water[,"new_long"], water[,"new_lat"], 0.05)
water_311H <- local_map(water[water$subsite_id == "311-H", "new_long"], water[water$subsite_id == "311-H", "new_lat"], 0.03)
water_345H <- local_map(water[water$subsite_id == "345-H", "new_long"], water[water$subsite_id == "345-H", "new_lat"], 0.03)
site.names <- data.frame(site = c("0369", "0311", "0345" ), 
                         longitude = c(-90.71492, -90.76659, -90.73892),
                         latitude = c(29.29509, 29.21496, 29.18148))
#plot with satellite map
ggmap(satmap,
      extent = "device", # "panel" keeps in axes, etc.
      ylab = "Latitude",
      xlab = "Longitude",
      legend = "right") +
  geom_point(data = water[water$subsite_id == "311-H",],
             aes(x = new_long, y = new_lat, colour = subsite_id)) +
  geom_point(data = water[water$subsite_id == "345-H",],
             aes(x = new_long, y = new_lat, colour = subsite_id)) +
  labs(fill = "Subsite_id from datsheet", title = str_c('Comparing Datasheet Labels',' to Actual Locations'),
       subtitle = str_c('White labels show general site location, \n',
          'colored points are what/where the datasheets say \n',
                        'plotted with the coordinates that go with them. \n',
                        'No clear sign of a mixup')) +
  geom_label(data = site.names, aes(x = longitude, y = latitude, label = site))
  

#--------------------------------------------------------------------------------------
#-----------------------  Table for open water survey models  -------------------------
#--------------------------------------------------------------------------------------
water %>% mutate(
  bogaert = as.numeric(bogaert),
  temp_c = as.numeric(temp_c),
  sal = as.numeric(sal),
  do_mgl = as.numeric(do_mgl),
  ph = as.numeric(ph),
  turb_ntu = as.numeric(turb_ntu),
  chla_ugl = as.numeric(chla_ugl)) %>% 
  select(subsite_id, bogaert, temp_c, sal, do_mgl, ph, turb_ntu, chla_ugl) %>% 
  group_by(subsite_id) %>%
  summarise(bogaert = mean(bogaert),
            temp_mean = mean(temp_c, na.rm = T),
            temp_se = sd(temp_c, na.rm = T)/ n(),
            sal_mwan = mean(sal, na.rm = T),
            sal_se = sd(sal, na.rm = T)/ n(),
            do_mean = mean(do_mgl, na.rm = T),
            do_se = sd(do_mgl, na.rm = T)/ n(),
            ph_mean = mean(ph, na.rm = T),
            ph_se = sd(ph, na.rm = T)/ n(),
            turbidity_mean = mean(turb_ntu, na.rm = T),
            turbidity_se = sd(turb_ntu, na.rm = T)/ n(),
            chla_mean = mean(chla_ugl, na.rm = T),
            chla_se = sd(chla_ugl, na.rm = T)/ n())









#--------------------------------------------------------------------------------------
#--------------------------- Marsh Diversity GLM's  -----------------------------------
#--------------------------------------------------------------------------------------

mod.df <- read_csv("all_transect_data.csv")

#variable setup
mod.df$crms_site <- factor(mod.df$crms_site, levels = c("345", "311", "369"))

#Histograms without zero inflation 
mod.df %>% filter(Richness > 1) %>% ggplot(aes(Shannon)) + geom_histogram()
mod.df %>% filter(Richness > 1) %>% ggplot(aes(Simpson)) + geom_histogram()
mod.df %>% filter(Richness > 1) %>% ggplot(aes(Evenness)) + geom_histogram()
mod.df %>% filter(Richness > 1) %>% ggplot(aes(Richness)) + geom_histogram()

#Richness and Simpson's summaries
mod.df <- mod.df %>% filter(dist_from_edge_m != 'full_transect', Richness != 0)

site_divsum <- mod.df %>% group_by(crms_site) %>% summarise(`min Richness` = min(Richness, na.rm = T),
                                             `max Richness` = max(Richness, na.rm = T),
                                             `med Richness` = median(Richness, na.rm = T),
                                             `mean Richness` = mean(Richness, na.rm = T),
                                             `stderr Richness` = sd(Richness, na.rm = T)/ sqrt(n()),
                                             `min Simpson` = min(Simpson, na.rm = T),
                                             `max Simpson` = max(Simpson, na.rm = T),
                                             `med Simpson` = median(Simpson, na.rm = T),
                                             `mean Simpson` = mean(Simpson, na.rm = T),
                                             `stderr Simpson` = sd(Simpson, na.rm = T)/ sqrt(n()))
View(site_divsum)

subsite_divsum <- mod.df %>% group_by(crms_site, bogaert) %>% summarise(`min Richness` = min(Richness, na.rm = T),
                                             `max Richness` = max(Richness, na.rm = T),
                                             `med Richness` = median(Richness, na.rm = T),
                                             `mean Richness` = mean(Richness, na.rm = T),
                                             `stderr Richness` = sd(Richness, na.rm = T)/ sqrt(n()),
                                             `min Simpson` = min(Simpson, na.rm = T),
                                             `max Simpson` = max(Simpson, na.rm = T),
                                             `med Simpson` = median(Simpson, na.rm = T),
                                             `mean Simpson` = mean(Simpson, na.rm = T),
                                             `stderr Simpson` = sd(Simpson, na.rm = T)/ sqrt(n()))


#--------------------------------------------------------------------------------------
#-------------------------  Data Exploration, Shannon  --------------------------------
#--------------------------------------------------------------------------------------
# glm(Shannon ~ bogaert + correct_elev_m + ss_sal_mean + (1 | crms_site), data = mod.df, family = "gaussian")

# Shannon and Covariates
mod.df %>% ggplot(aes(Shannon)) + geom_histogram()
plot(mod.df$correct_elev_m, mod.df$Shannon)         #FE
plot(mod.df$ss_sal_mean, mod.df$Shannon)            #FE
boxplot(mod.df$Shannon ~ mod.df$crms_site)          #RE
plot(mod.df$bogaert, mod.df$Shannon)                #FE

head(table(mod.df$Shannon)) #Shannon index is super zero inflated, Sa dominates

#### xyplots and box and whisker plots
## Data Exploration
X.labs <- c(mod.df$correct_elev_m, mod.df$ss_sal_mean, mod.df$ss_sal_min, mod.df$ss_sal_max, mod.df$crms_site,  mod.df$bogaert)
CPUE.reps <- rep(mod.df$Shannon, 6)
I12 <- rep(c("Elevation","Sal mean","Sal min","Sal max","CRMS site", "Bogaert"), each = nrow(mod.df))
ID12 <- rep(I12, 6)
xyplot(CPUE.reps ~ X.labs | ID12, col = 1,
       strip = function(bg = 'white', ...)
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE,
                     x = list(relation = 'free'),
                     y = list(relation = 'same')),
       xlab = 'Explanatory Variables',
       ylab = 'Shannon Index',
       panel = function(x,y){
         panel.grid(h = -1, v = 2)
         panel.points(x, y, col = 1)
       })




#--------------------------------------------------------------------------------------
#------------------------------- Simpson models ---------------------------------------
#--------------------------------------------------------------------------------------
#Dataset with full transect totals
full.t <- read_csv("all_transect_data.csv")

new.df <- filter(full.t, dist_from_edge_m != "full_transect", Simpson != 1)
new.df$crms_site <- factor(new.df$crms_site, levels = c('345', '311', '369'))

# Simpson and Covariates
theme_set(theme_classic())
new.df %>% ggplot(aes(Simpson)) + geom_density(fill = 'gray40') + facet_wrap(~crms_site)
new.df %>% ggplot(aes(correct_elev_m, Simpson)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~crms_site)
new.df %>% ggplot(aes(ss_sal_mean, Simpson)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~crms_site)
new.df %>% ggplot(aes(bogaert, Simpson)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~crms_site)
new.df %>% ggplot(aes(crms_site, Simpson)) + geom_boxplot()

#Simpson index is super zero inflated, Sa dominates
head(table(new.df$Simpson)) 

#change labels
new.df$crms_site <- str_c("0", new.df$crms_site)
new.df$crms_site <- factor(new.df$crms_site, levels = c('0345', '0311', '0369'))

#Response- Simpsons; Predictor- Site
mod1 <- lm(Simpson ~ crms_site, data = new.df)
summary(mod1)
anova(mod1)
mod1 <- aov(Simpson ~ crms_site, data = new.df)
plot(TukeyHSD(mod1))
AIC(mod1)
new.df %>%  ggplot(aes(crms_site, Simpson)) + geom_boxplot()
plot(Simpson ~ crms_site, data = new.df, 
     xlab = "CRMS Station", ylab = 'Simpson Diversity Index Scores',
     ylim = c(0,1), las = 1)
text(1, 0.9, 'a')
text(2, 0.9, 'b')
text(3, 0.9, 'b')

#Response- Simpsons; Predictor- Fragmentation
mod2 <- lm(Simpson ~ bogaert, data = new.df)
summary(mod2)
AIC(mod2)

x_range <- seq(min(new.df$bogaert), max(new.df$bogaert), 0.1)
preddat <- predict(mod2, list(bogaert = x_range), type = 'response', se.fit = T)
plot(Simpson ~ bogaert, data = new.df, las = 1,
     xlab = "Bogaert Fragmentation Metric", ylab = "Simpson Diversity Index Scores")
lines(x_range, preddat$fit, lwd = 2)
lines(x_range, preddat$fit + 1.96 * preddat$se.fit, lwd = 2, lty = 2)
lines(x_range, preddat$fit - 1.96 * preddat$se.fit, lwd = 2, lty = 2)




#Response- Simpsons; Predictor- Site, Fragmentation, SitexFrag
mod3 <- lmer(Simpson ~ bogaert + (1 | crms_site), data = new.df)
summary(mod3)
AIC(mod3)
r.squaredGLMM(mod3)

plot(Simpson ~ bogaert, data = new.df, las = 1,
     xlab = "Bogaert Fragmentation Metric", ylab = "Simpson Diversity Index Scores")
for (i in 1:length(unique(new.df$crms_site))) {
  abline(coef(mod3)$crms_site[i,1],coef(mod3)$crms_site[i,2], lty = i, lwd = 2)
}


#Anovas
new.df$bog.f <- parse_factor(as.character(round(as.numeric(as.character(new.df$bogaert)), digits = 2)), 
  levels = ,c(as.character(sort(unique(round(as.numeric(as.character(new.df$bogaert)), digits = 2)), decreasing = T))))

mod4 <- lm(Simpson ~ bog.f, data=new.df)
summary(mod4)
Anova(mod4)
new.df %>% ggplot(aes(bog.f, Simpson)) + geom_boxplot() + facet_wrap(~crms_site)
mod4 <- aov(Simpson ~ bog.f, data=new.df)
TukeyHSD(mod4, ordered = T)
plot(TukeyHSD(mod4, ordered = T)) 
#boxplot
boxplot(Simpson ~ bog.f, data = new.df, las = 1, ylim = c(0,1),
        xlab = 'Bogaert Fragmentation Metric', ylab = 'Simpson Diversity Index Scores')
# text(1, 0.9, 'a')
# text(2, 0.9, 'a')
# text(3, 0.9, 'ab')
# text(4, 0.9, 'c')
# text(5, 0.9, 'ab')
# text(6, 0.9, 'c')
# text(7, 0.9, 'ab')
# text(8, 0.9, 'c')
# text(9, 0.9, 'c')


#two-way
new.df$dist_from_edge_m <- parse_factor(new.df$dist_from_edge_m, levels = c(0,2,5))
mod5 <- lm(Simpson ~ crms_site + dist_from_edge_m + crms_site * dist_from_edge_m, data = new.df)
summary(mod5)
anova(mod5)
mod5 <- aov(Simpson ~ crms_site + dist_from_edge_m + crms_site:dist_from_edge_m, data = new.df)
plot(TukeyHSD(mod5))
boxplot(Simpson ~ dist_from_edge_m, data = new.df, las = 1, ylim = c(0,1),
        xlab = 'Distance from Marsh Edge (m)', ylab = 'Simpson Diversity Index Scores')


#GLMM elevation
interaction.plot(new.df$correct_elev_m, new.df$crms_site, new.df$Simpson)
mod6 <- lmer(Simpson ~ correct_elev_m + (1 | crms_site), data = new.df)
summary(mod6)
AIC(mod6)
#anova(mod6)
r.squaredGLMM(mod6)

plot(Simpson ~ correct_elev_m, data = new.df, las = 1,
     xlab = "Marsh Elevation (m)", ylab = "Simpson Diversity Index Scores")
for (i in 1:length(unique(new.df$crms_site))) {
  abline(coef(mod6)$crms_site[i,1],coef(mod6)$crms_site[i,2], lty = i, lwd = 2)
}

#--------------------------------------------------------------------------------------
#-------------------------------  Richness models ----------------------------------------
#--------------------------------------------------------------------------------------
#Data Exploration plots
new.df %>% ggplot(aes(Richness)) + geom_histogram() + facet_wrap(~crms_site)
mean(new.df$Richness); sd(new.df$Richness)
new.df %>% ggplot(aes(bogaert, Richness)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~crms_site)
new.df %>% ggplot(aes(correct_elev_m, Richness)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~crms_site)
new.df %>% ggplot(aes(ss_sal_mean, Richness)) + geom_point() + geom_smooth(method = 'lm') + facet_wrap(~crms_site)
new.df %>% ggplot(aes(ss_sal_min, Richness)) + geom_point()
new.df %>% ggplot(aes(ss_sal_max, Richness)) + geom_point()
new.df %>% ggplot(aes(crms_site, Richness)) + geom_boxplot()
new.df %>% ggplot(aes(subsite_id, Richness)) + geom_boxplot()




#Response- Richness; Predictor- Site
mod1 <- lm(Richness ~ crms_site, data = new.df)
summary(mod1)
anova(mod1)
#plot(TukeyHSD(mod1))
AIC(mod1)
boxplot(Richness ~ crms_site, data = new.df, las = 1, 
        xlab = "CRMS Station", ylab = 'Species Richness',
        ylim = c(0,6), las = 1)
text(1, 5.5, 'a')
text(2, 5.5, 'b')
text(3, 5.5, 'b')

#Response- Richness; Predictor- Fragmentation
mod2 <- lm(Richness ~ bogaert, data = new.df)
summary(mod2)
AIC(mod2)


x_range <- seq(min(new.df$bogaert), max(new.df$bogaert), 0.1)
preddat <- predict(mod2, list(bogaert = x_range), type = 'response', se.fit = T)
plot(Richness ~ bogaert, data = new.df, las = 1,
     xlab = "Bogaert Fragmentation Metric", ylab = "Species Richness")
lines(x_range, preddat$fit, lwd = 2)
lines(x_range, preddat$fit + 1.96 * preddat$se.fit, lwd = 2, lty = 2)
lines(x_range, preddat$fit - 1.96 * preddat$se.fit, lwd = 2, lty = 2)

#Response- Richness; Predictor- Site, Fragmentation, SitexFrag
mod3 <- lmer(Richness ~ bogaert + (1 | crms_site) + bogaert * (1 | crms_site), data = new.df)
summary(mod3)
AIC(mod3)
r.squaredGLMM(mod3)


plot(Richness ~ bogaert, data = new.df, las = 1,
     xlab = "Bogaert Fragmentation Metric", ylab = "Species Richness")
for (i in 1:length(unique(new.df$crms_site))) {
  abline(coef(mod3)$crms_site[i,1],coef(mod3)$crms_site[i,2], lty = i, lwd = 2)
}


#Anovas
mod4 <- lm(Richness ~ bog.f, data=new.df)
summary(mod4)
anova(mod4)
new.df %>% ggplot(aes(bog.f, Richness)) + geom_boxplot() + facet_wrap(~crms_site)
mod4 <- aov(Richness ~ bog.f, data=new.df)
TukeyHSD(mod4)
plot(TukeyHSD(mod4))
boxplot(Richness ~ bog.f, data = new.df, las = 1, ylim = c(0,6),
        xlab = 'Bogaert Fragmentation Score', ylab = 'Species Richness')


#two-way
new.df$dist_from_edge_m <- parse_factor(new.df$dist_from_edge_m, levels = c(0,2,5))
mod5 <- lm(Richness ~ crms_site + dist_from_edge_m + crms_site * dist_from_edge_m, data = new.df)
summary(mod5)
anova(mod5)
mod5 <- aov(Richness ~ crms_site + dist_from_edge_m + crms_site:dist_from_edge_m, data = new.df)
plot(TukeyHSD(mod5))

boxplot(Richness ~ dist_from_edge_m, data = new.df, las = 1, ylim = c(0,6),
        xlab = 'Bogaert Fragmentation Score', ylab = 'Species Richness')



#GLMM elevation
interaction.plot(new.df$correct_elev_m, new.df$crms_site, new.df$Richness)
mod6 <- lmer(Richness ~ correct_elev_m + (1 | crms_site), data = new.df)
summary(mod6)
AIC(mod6)
#anova(mod6)
r.squaredGLMM(mod6)


plot(Richness ~ correct_elev_m, data = new.df, las = 1,
     xlab = "Marsh Elevation (m)", ylab = "Species Richness")
for (i in 1:length(unique(new.df$crms_site))) {
  abline(coef(mod6)$crms_site[i,1],coef(mod6)$crms_site[i,2], lty = i, lwd = 2)
}





#------------- abandoned models ----
library(betareg)
# #Negative Binomial
# mod1 <- glmer.nb(Simpson ~ bogaert + correct_elev_m + ss_sal_mean + (1 | crms_site), data = new.df)
# library(zoib)
# mod2 <- zoib(Simpson ~ correct_elev_m + as.factor(crms_site)| 1 |correct_elev_m + as.factor(crms_site)| 1,
#      data = new.df, 
#      joint = FALSE, 
#      random = 12,
#      EUID = new.df$quadrat_id, 
#      zero.inflation = TRUE,
#      one.inflation = FALSE, 
#      n.iter = 1100, 
#      n.thin = 5,
#      n.burn=100)

# library(MuMIn)
# mod1 <- glmer(Richness ~ bogaert + correct_elev_m + ss_sal_mean + (1 | crms_site), data = new.df, family = poisson)
# summary(mod1)
# AIC(mod1)
# coef(mod1)
# r.squaredGLMM(mod1)
# 
# #best
# mod2 <- lmer(Richness ~ correct_elev_m  + ss_sal_mean + (1 | crms_site), data = new.df, family = poisson)
# AIC(mod2)
# r.squaredGLMM(mod2)
# 
# mod3 <- glmer(Richness ~ correct_elev_m  + ss_sal_mean + (1 | subsite_id), data = new.df, family = poisson)
# AIC(mod3)
# r.squaredGLMM(mod3)

#-----  3d plot of model 2, best model
# # Make data into grid
# #make a grid to predict the surface
# elev <- seq(min(new.df$correct_elev_m, na.rm = T)-0.5, max(new.df$correct_elev_m, na.rm = T)+0.5, .005)
# sal <- seq(min(new.df$ss_sal_mean, na.rm = T)-0.5, max(new.df$ss_sal_mean, na.rm = T)+0.5, .05 )
# 
# #try manually predicting each surface
# rich_grid <- expand.grid(elev, sal)
# #manually encode site level intercepts from coef(mod2) ranef(mod2)[1]
# CRMS345 <- -0.03061804
# CRMS311 <- -0.01619342
# CRMS369 <-  0.65217710 
# 
# rich_grid2CRMS345 <- predict(mod2, rich_grid, allow.new.levels=TRUE, re.form=NA)
# rich_grid2CRMS311 <- predict(mod2, rich_grid, allow.new.levels=TRUE, re.form=NA)
# rich_grid2CRMS369 <- predict(mod2, rich_grid, allow.new.levels=TRUE, re.form=NA)
# #fail
# 
# 
# #--------------  attempt 2 from stackoverflow
# # Predict elevation in two dimensions
# rich_grid2 <- expand.grid(elev, sal, group = unique(new.df$crms_site))
# rich_grid2$Richness <- predict(mod2, rich_grid, allow.new.levels=TRUE, re.form=NA)
# 
# 
# surf1 <- ( matrix(rich_grid2[rich_grid2$group == "345", ]$Richness, nrow = length(elev), ncol = length(sal)) )
# surf2 <- ( matrix(rich_grid2[rich_grid2$group == "311", ]$Richness, nrow = length(elev), ncol = length(sal)) )
# surf3 <- ( matrix(rich_grid2[rich_grid2$group == "369", ]$Richness, nrow = length(elev), ncol = length(sal)) )
# 
# group <- c(rep("345", nrow(surf1)), rep("311", nrow(surf2)), rep('369', nrow(surf3)))
# 
# 
# #plot them
# library(plotly)
# plot_ly(z = ~surf1, type="surface") %>%
#   add_surface(z = surf2, surfacecolor=surf2) %>% 
#   add_surface(z = surf3, surfacecolor=surf3) %>% 
#   layout(
#     title = "3d Surface of Model 2",
#     scene = list(
#       xaxis = list(title = "Elevation (m)"), 
#       yaxis = list(title = "Mean Salinity"),
#       zaxis = list(title = "Predicted Richness"))) 
# #awful


#--------------------------------------------------------------------------------------
#------------------------------ Marsh Elevation GLM -----------------------------------
#--------------------------------------------------------------------------------------
#Dataset with full transect totals
full.t <- read_csv("all_transect_data.csv")

new.df <- filter(full.t, dist_from_edge_m != "full_transect", !is.na(correct_elev_m)) %>% 
  mutate(Dist_from_edge = as.numeric(dist_from_edge_m),
         bog_round = round(as.numeric(as.character(bogaert))),
         bog_round = factor(bog_round, levels = sort(unique(bog_round), decreasing =  T))
  )

# GLM, no random effect
new.df %>% ggplot(aes(Dist_from_edge, correct_elev_m)) +
  geom_point() + geom_smooth(method = 'lm', colour = 'black') + theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.line = element_line(linetype = "solid"), 
    panel.background = element_rect(fill = NA))

mod1 <- lm(correct_elev_m ~ Dist_from_edge, data = new.df)

summary(mod1)
AIC(mod1)
coef(mod1)

# GLMM
new.df %>% ggplot(aes(Dist_from_edge, correct_elev_m)) +
  geom_point() + geom_smooth(method = 'lm', colour = 'black') +
  facet_wrap(~bog_round)


#Elevation is nice and normal, so basic model, gaussian distribution, maximum likelihood estimation
library(lme4)
mod2 <- lmer(correct_elev_m ~ Dist_from_edge + (1 | subsite_id), data = new.df, REML = F)
summary(mod2)


#to get p values
library(car)
Anova(mod2)
Anova(mod2, type = 2)

#to get r-squared
library(MuMIn)
r.squaredGLMM(mod2)

#Lower Aic scores from glmm

#-------------  plot GLM model

#predict over x range using model coefficients
x_range <- seq(min(new.df$Dist_from_edge, na.rm = T), max(new.df$Dist_from_edge, na.rm = T), 0.01)
preddat <- predict(mod1, list(Dist_from_edge = x_range), type = 'response', se.fit=TRUE)

dev.off()
png("Elevation_glm.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
plot(correct_elev_m ~ Dist_from_edge, data = new.df, las = 1,
     xlab = "Distance from Marsh Edge (m)", ylab = "Elevation (m)")
lines(x_range, preddat$fit, lwd = 2)
lines(x_range, preddat$fit + 1.96 * preddat$se.fit, lwd = 2, lty = 2)
lines(x_range, preddat$fit - 1.96 * preddat$se.fit, lwd = 2, lty = 2)
dev.off()


#----------- plot lme4 model
#extract coefficients from model
Vcov <- vcov(mod2, useScale = FALSE)                    #calculate variance covariance matrix
betas <- fixef(mod2)                                    #extract fixed-effects estimates
se <- sqrt(diag(Vcov))                                  #extract diagonal from matrix, take the sqrt
zval <- betas / se                                      #get z score
pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)        #get p value
cbind(betas, se, zval, pval)

#extract random intercepts
AIC(mod2)
coef(mod2)$subsite_id #column 1 is intercepts, column 2 is slope, which is the same for each in this mod
new.df  %>% group_by(subsite_id) %>% summarise(bogaert_value = mean(bogaert)) %>% arrange(desc(bogaert_value))


#Plot them all together
png("Elevation_glmm.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
plot(correct_elev_m ~ Dist_from_edge, data = new.df, las = 1,
     xlab = "Distance from Marsh Edge (m)", ylab = "Elevation (m)")
for (i in 1:length(unique(new.df$subsite_id))) {
  abline(coef(mod2)$subsite_id[i,1],coef(mod2)$subsite_id[i,2])
}
dev.off()





#---------------  Model 3 glm
mod3 <- glm(correct_elev_m ~ Dist_from_edge + bogaert, data = new.df)
summary(mod3)
AIC(mod3)
coef(mod3)
Anova(mod3)
RsquareAdj(mod3); r.squaredGLMM(mod3)

#make a grid to predict the surface
dist <- seq(min(new.df$Dist_from_edge, na.rm = T)-0.5, max(new.df$Dist_from_edge, na.rm = T)+0.5, .05 )
bog <- seq(min(new.df$bogaert, na.rm = T)-0.5, max(new.df$bogaert, na.rm = T)+0.5, .1 )
gom_grid <- expand.grid(dist, bog)
colnames(gom_grid) <- c("Dist_from_edge","bogaert")

# Predict elevation in two dimensions
gom_grid$correct_elev_m <- predict(mod3, gom_grid, type = "response")



# Make data into grid
library(reshape2)
Predicted_Elev <- acast(gom_grid, Dist_from_edge~bogaert, value.var = "correct_elev_m")

### PLOT flat surface with contours ###
# Get range of values
z.range <- c(min(gom_grid$correct_elev_m), max(gom_grid$correct_elev_m))
y.range <- c(min(gom_grid$bogaert), max(gom_grid$bogaert))
x.range <- c(min(gom_grid$Dist_from_edge), max(gom_grid$Dist_from_edge))

# Plot contour
library(oce)

png("Elevation_glm_3.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
plot( x = NA, y =NA,  xlab = "Distance from Marsh Edge (m)", ylab = 'Bogaert Fragmentation Score',xlim = x.range, ylim = y.range,  main = "Predicted Elevation")
image(dist, bog, Predicted_Elev, col = oceColorsTemperature(100), cex.lab = 1.5, cex.axis = 1.4, add=T, zlim = c(min(Predicted_Elev), max(Predicted_Elev)))
contour(dist, bog, Predicted_Elev, add=T, color = "black", axes = F, xlab=NA,  ylab=NA, main = NA)
dev.off()
#or use plotly for 3d
library(plotly)

#Set font attributes
f <- list(
  family = "Arial",
  size = 18,
  color = "#7f7f7f")

# volcano is a numeric matrix that ships with R
plot_ly(z = ~Predicted_Elev) %>% add_surface() %>% 
  layout(
    title = "3d Surface of Model 3",
    scene = list(
      xaxis = list(title = "Distance from Edge", titlefont = f), 
      yaxis = list(title = "Bogaert Score", titlefont = f),
      zaxis = list(title = "Predicted Elevation", titlefont = f))) 


#Elevation summaries
elev_summ <- full.t %>% 
  filter(dist_from_edge_m != 'full_transect') %>% 
  group_by(crms_site, dist_from_edge_m) %>% 
  summarise(mean_elev = mean(correct_elev_m, na.rm = T), 
            elev_se = sd(correct_elev_m, na.rm = T)/n())

View(elev_summ)

overal_elev <- full.t %>% 
  group_by(dist_from_edge_m) %>% 
  summarise(mean_elev = mean(correct_elev_m, na.rm = T), 
            elev_se = sd(correct_elev_m, na.rm = T)/n())

View(overal_elev)
#--------------------------------------------------------------------------------------
#------------------------ Under Water Elevation ~ Fragmentation -----------------------
#--------------------------------------------------------------------------------------
water <- read_xlsx("noaa_blu_objective1_openwater_survey.xlsx")

#make light attenuation measurement 
sp <- .15
water <- water %>% mutate(
  month = ifelse(lubridate::month(date) == 4, "April", "September"),
  sal = as.numeric(sal),
  correct_elev_m = as.numeric(correct_elev_m),
  crms_site = factor(crms_site, levels = c("369", "311", "345")),
  frag = toupper(frag),
  subsite_id = str_c(crms_site, frag, sep = '-'),
  depth_cm = as.numeric(depth_cm),
  top_par1 = as.numeric(top_par1),                                          #bulb 1 is top bulb
  top_par2 = as.numeric(top_par2),                                          #bulb 2 is low bulb
  bot_par1 = as.numeric(bot_par1),
  bot_par2 = as.numeric(bot_par2),                                                                       
  Kd_surface = -(log10(top_par1 / top_par2)) / sp,                          #attenuation per .15m for surface measurements
  Kd_bottom = -(log10(bot_par1 / bot_par2)) / sp,                           #Kd at bottom, per .15m                     #multiply light lost in 15cm by (total depth/15cm increments)
  SI_full = (bot_par2/top_par1) * 100) %>% 
  filter(month == "April") 

bog.df <- structure(list(subsite_id = c("369-L", "369-M", "369-H", "311-L", "311-M", "311-H", "345-L", "345-M", "345-H"), 
                         bogaert = c(157.409232427562,158.794835413732, 144.600763346679, 152.825951296738, 147.180047035677, 
                                     139.741186393264, 145.714196114844, 148.353319527527, 139.722070072008)), 
                    class = "data.frame", row.names = c(NA, -9L))


water <- left_join(water, bog.df, by = "subsite_id")


#summary stats on depth
water %>% summarise(
  min = min(correct_elev_m, na.rm = T), 
  max = max(correct_elev_m, na.rm = T), 
  mean = mean(correct_elev_m, na.rm = T), 
  se = (sd(correct_elev_m, na.rm = T)/ sqrt(length(water$correct_elev_m))))


water %>% ggplot(aes(1/bogaert, correct_elev_m)) + 
  geom_point() + 
  geom_smooth(method = 'lm', color = 'black')


wmod1 <- glm(correct_elev_m ~ bogaert, data = water)
summary(wmod1)
AIC(wmod1)
coef(wmod1)

#plot fit
plot(correct_elev_m ~ bogaert, data = water, las = 1,
     xlab = "Bogaert Fragmentation Metric", ylab = "Elevation (m)")
abline(coef(wmod1[1]), coef(wmod1[2]), lwd = 2)

#plot fitted vs. residuals
plot(fitted(wmod1), resid(wmod1), xlab = 'fitted vals', ylab = 'residuals')
abline(h = 0, lty = 2, lwd = 2, col = 'royalblue')




#-------------  plot GLM model
#predict over x range using model coefficients
x_range <- seq(min(water$bogaert, na.rm = T), max(water$bogaert, na.rm = T), 0.01)
preddat <- predict(wmod1, list(bogaert = x_range), type = 'response', se.fit=TRUE)

dev.off()

png("water_depth_glm.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
plot(correct_elev_m ~ bogaert, data = water, las = 1,
     xlab = "Bogaert Fragmentation Metric", ylab = "Elevation (m)")
lines(x_range, preddat$fit, lwd = 2)
lines(x_range, preddat$fit + 1.96 * preddat$se.fit, lwd = 2, lty = 2)
lines(x_range, preddat$fit - 1.96 * preddat$se.fit, lwd = 2, lty = 2)

dev.off()


#--------------------------------------------------------------------------------------
#------------------------------- %SI ~ bogaert ----------------------------------------
#--------------------------------------------------------------------------------------
hist(water$SI_full)
water$SI <- water$SI_full/100


wmod1 <- glm(SI ~ bogaert, data = water)
#plot fit
plot(SI ~ bogaert, data = water, las = 1,
     xlab = "Bogaert Fragmentation Metric", ylab = "Elevation (m)")
abline(coef(wmod1)[1], coef(wmod1)[2], lwd = 2)

#plot fitted vs. residuals
plot(fitted(wmod1), resid(wmod1), xlab = 'fitted vals', ylab = 'residuals')
abline(h = 0, lty = 2, lwd = 2, col = 'royalblue')


#or Complicated beta regressions
wmod2 <- betareg(SI ~ bogaert, data = water,link = "loglog")
wmod3 <- betareg(SI ~ bogaert, data = water)

library(ggplot2)
ggplot(water, aes(x = bogaert, y = SI)) +
  geom_point(size = 4, aes(fill = subsite_id), shape = 21) +
  scale_fill_grey() +
  geom_line(aes(y = predict(wmod2, water),
                colour = "log-log", linetype = "log-log")) +
  geom_line(aes(y = predict(wmod3, water), 
                colour = "logit", linetype = "logit")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_classic() + theme(text=element_text(size=18, family = "Arial"))



p <- ggplot(water, aes(x = bogaert, y = SI)) +
  geom_point(size = 4, aes(fill = subsite_id), shape = 21) +
  scale_fill_grey() +
  geom_line(aes(y = predict(wmod2, water),
                colour = "log-log", linetype = "log-log")) +
  geom_line(aes(y = predict(wmod3, water), 
                colour = "logit", linetype = "logit")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_classic() + theme(text=element_text(size=18, family = "Arial"))
ggsave("SI_betaregression.png", plot = p, height = 5, width = 5, units = "in", dpi = 300)




AIC(wmod1); AIC(wmod2); AIC(wmod3)


#Model 1 is better
summary(lm(SI ~ bogaert, data = water)) #lm gives you r2


#-------------  plot GLM model
wmod1 <- lm(SI ~ bogaert, data = water)

#predict over x range using model coefficients
x_range <- seq(min(water$bogaert, na.rm = T), max(water$bogaert, na.rm = T), by = 0.001)
preddat <- predict(wmod1, list(bogaert = x_range), type = 'response', se.fit=TRUE)

dev.off()

png("water_SI_glm.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
plot(SI ~ bogaert, data = water, las = 1,
     xlab = "Bogaert Fragmentation Metric", ylab = "% Surface Irradiance at Bottom")
lines(x_range, preddat$fit, lwd = 2)
lines(x_range, preddat$fit + 1.96 * preddat$se.fit, lwd = 2, lty = 2)
lines(x_range, preddat$fit - 1.96 * preddat$se.fit, lwd = 2, lty = 2)

dev.off()



#-------------  plot LM model, percentage SI
wmod1 <- lm(SI_full ~ bogaert, data = water)

x_range <- seq(min(water$bogaert, na.rm = T), max(water$bogaert, na.rm = T), by = 0.001)
preddat <- predict(wmod1, list(bogaert = x_range), type = 'response', se.fit=TRUE)


png("water_SI_percent_glm.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
plot(SI_full ~ bogaert, data = water, las = 1,
     xlab = "Bogaert Fragmentation Metric", ylab = "% Surface Irradiance at Bottom")
lines(x_range, preddat$fit, lwd = 2)
lines(x_range, preddat$fit + 1.96 * preddat$se.fit, lwd = 2, lty = 2)
lines(x_range, preddat$fit - 1.96 * preddat$se.fit, lwd = 2, lty = 2)
dev.off()






#--------------------------------------------------------------------------------------
#-----------------------------% Likelihood SAV ----------------------------------------
#--------------------------------------------------------------------------------------

frag50 <- read_xlsx("Frag at 50m and SAV likelihood.xlsx")

frag50 <- frag50 %>% 
  mutate(SAV = .$`Mean % likelihood SAV`) %>% 
  filter(!is.na(SAV))

frag50 %>% ggplot(aes(bogaert, SAV)) + geom_point()
frag50 %>% ggplot(aes(SAV)) + geom_density()


mod1 <- lm(SAV ~ bogaert, data = frag50)
summary(mod1)
AIC(mod1)

#Diagnose fit
x_range <- seq(min(frag50$bogaert), max(frag50$bogaert), by = 0.001)
preddat <- predict(mod1, list(bogaert = x_range), type = 'response', se.fit=TRUE)
ypred <- preddat$fit
# lines(x_range, preddat$fit, lwd = 2)
# lines(x_range, preddat$fit + 1.96 * preddat$se.fit, lwd = 2, lty = 2)
# lines(x_range, preddat$fit - 1.96 * preddat$se.fit, lwd = 2, lty = 2)


png("SAV_bog_glm.png", height = 4, width = 6, units = "in", res = 300)
par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
plot(SAV ~ bogaert, data = frag50, xlab = 'Bogaert Fragmentation Score', ylab = "% Likelihood SAV present")
lines(x_range, ypred, lwd = 2)
lines(x_range, preddat$fit + 1.96 * preddat$se.fit, lwd = 2, lty = 2)
lines(x_range, preddat$fit - 1.96 * preddat$se.fit, lwd = 2, lty = 2)
dev.off()

plot(fitted(mod1),resid(mod1), xlab = 'Fitted Values', ylab = 'Residuals')
abline(h = 0, lty = 2, lwd = 2, col = 'blue')

hist(resid(mod1))




#--------------------------------------------------------------------------------------
#----------------------------- species encounter heatmap ------------------------------
#--------------------------------------------------------------------------------------
full.t <- read_csv("all_transect_data.csv")

full.t <- filter(full.t, dist_from_edge_m != "full_transect") %>% 
  mutate(Dist_from_edge = as.numeric(dist_from_edge_m),
         bog_round = round(as.numeric(as.character(bogaert))),
         bog_round = factor(bog_round, levels = sort(unique(bog_round), decreasing =  T))
  )

bin.mat <- ifelse(full.t[,6:28] > 0, 1, 0) 
bin.mat <- as.tibble(bin.mat)
bin.mat <- bind_cols(full.t[,'crms_site'], bin.mat)
bin.mat$crms_site <- factor(bin.mat$crms_site, levels = c('345', '311', '369'))

library(magrittr)
bin.mat.site <- bin.mat %>% group_by(crms_site) %>% dplyr::summarise(
  Sa = sum(Sa),
  Jr = sum(Jr),
  Ag = sum(Ag),
  Sp = sum(Sp),
  Ds = sum(Ds),
  Bm = sum(Bm),
  Sb = sum(Sb),
  Pv = sum(Pv),
  Is = sum(Is),
  St = sum(St),
  Sr = sum(Sr),
  Sc = sum(Sc),
  Zm = sum(Zm),
  Ep = sum(Ep),
  Aa = sum(Aa),
  Pa = sum(Pa),
  Na = sum(Na),
  If = sum(If),
  Bf = sum(Bf),
  Vl = sum(Vl),
  PaVi = sum(pavi),
  ScAm = (sum(ScAm) + sum(`s?`))
)




bin.mat.site <- bin.mat %>% group_by(crms_site) %>% dplyr::summarise(
  `Spartina alterniflora` = sum(Sa),
  `Juncus roemarianus` = sum(Jr),
  `Avicennia germinans` = sum(Ag),
  `Spartina patens` = sum(Sp),
  `Distichlis spicata` = sum(Ds),
  `Batis maritima` = sum(Bm),
  `Salicornia bigelovii` = sum(Sb),
  `Paspalum vaginatum` = sum(Pv),
  `Ipomoea sagittata` = sum(Is),
  `Symphotrichum tenuifolium` = sum(St),
  `Schenoplectus robustus` = sum(Sr),
  `Spartina cynusuroides` = sum(Sc),
  `Zizaniopsis miliacea` = sum(Zm),
  `Eclipta prostrata` = sum(Ep),
  `Amaranthus australis` = sum(Aa),
  `Phragmites australis` = sum(Pa),
  `Nemophila aphylla` = sum(Na),
  `Iva frutescens` = sum(If),
  Bf = sum(Bf),
  Vl = sum(Vl),
  PaVi = sum(pavi),
  ScAm = (sum(ScAm) + sum(`s?`))
)


bin.mat.site[,2:23] <- ifelse(bin.mat.site[,2:23] > 1, 1, 0)


#remove unwanted empty columns
colSums(bin.mat.site[,2:23])
drop.cols <- c("Symphotrichum tenuifolium", 'Amaranthus australis', "Iva frutescens", "Bf", "Vl", "PaVi","ScAm")
bin.mat.site <- bin.mat.site %>% select(-one_of(drop.cols))




# png("pres_abs_tiles.png", height = 4, width = 6, units = "in", res = 300)
# par(mar = c(5,5,4,2), ps = 16, cex.lab = 1.125, cex.axis = 1, family = "Arial")
# image(t(as.matrix(bin.mat.site[c(3,1,2),23:2])), xlab = 'Species', ylab = 'CRMS Site', axes = F, col = c('gray60','white'))
# axis(1, at = seq(0,1, 1/21), labels = colnames(bin.mat.site)[2:23], tick = T, line = NA, las = 2)
# axis(2, at = seq(0,1, 1/2), labels = c('345','311','369'), tick = F, line = NA, las = 2)
# dev.off()

#or use ggplot
test.df <-bin.mat.site
test.df$crms_site <- factor(test.df$crms_site, levels = c('345', '311', '369'))
melt.t <- melt(test.df, id.vars = 'crms_site', value.name = "Presence")
melt.t$Ocurrence <- ifelse(melt.t$Presence == 1, 'Present', 'Absent')
melt.t$Ocurrence <- factor(melt.t$Ocurrence, levels = c("Present", "Absent"))
# library(plyr)
# melt.t <- ddply(melt.t, .(variable)) #rescale

#test plot
base_size <- 9
ggplot(melt.t, aes(variable, crms_site)) + 
  geom_tile(aes(fill = Ocurrence), colour = 'white') +
  scale_fill_manual(values = c('gray40', 'gray90')) +
  theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(face = 'italic', angle = 90, hjust = 1))


#Gray with white lines
p <- ggplot(melt.t, aes(variable, crms_site)) + 
  geom_tile(aes(fill = Ocurrence), colour = 'white') +
  scale_fill_manual(values = c('gray40', 'gray90')) +
  theme_classic(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text=element_text(size=12, family = "Arial"))
ggsave("pres_abs_heatmap_fixed.png", plot = p, height = 4, width = 6, units = "in", dpi = 300)


p <- ggplot(melt.t, aes(variable, crms_site)) +              #Change row indexes to exclude empty rows
  geom_tile(aes(fill = Ocurrence), colour = 'black') +
  scale_fill_manual(values = c('gray40', 'white')) +
  theme_classic(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(face = 'italic', angle = 90, hjust = 1, vjust = 0.5), text=element_text(size=12, family = "Arial"))
ggsave("pres_abs_heatmap_white.png", plot = p, height = 4, width = 6, units = "in", dpi = 300)


#Spin it and make it pretty
base_size <- 11
melt.t$crms_site <- recode(melt.t$crms_site, "345" = "0345", "311" = "0311", "369" = "0369")
p <- ggplot(melt.t, aes(crms_site, variable)) + 
  geom_tile(aes(fill = Ocurrence), colour = 'black') +
  scale_fill_manual(values = c('gray40', 'white')) +
  theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.y = element_text(face = 'italic', hjust = 1))
ggsave("pres_abs_heatmap_white_vert.png", plot = p, height = 4, width = 6, units = "in", dpi = 300)


names(full.t)
full.t %>% rename(Site = crms_site) %>% names(.)



