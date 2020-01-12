#### Final Plot Figures for Manuscript ####

####  Packages  ####
library(tidyverse)
library(patchwork)

#Data location
setwd('~/Dropbox/TBone_2019/data/')

#### 1. SAV ~ Bogaert  ####

#Data
frag50 <- read_xlsx("Frag at 50m and SAV likelihood.xlsx")

#Cleanup
frag50 <- frag50 %>% 
  mutate(SAV = .$`Mean % likelihood SAV`) %>% 
  filter(!is.na(SAV))

#Model
mod1 <- lm(SAV ~ bogaert, data = frag50)

#Predictions
x_range <- seq(min(frag50$bogaert), max(frag50$bogaert), by = 2)
preddat <- predict(mod1, list(bogaert = x_range), type = 'response', se.fit=TRUE)
mod1_preds <- data.frame(x = x_range,
                         predicted = preddat$fit,
                         ucl =  preddat$fit + 1.96 * preddat$se.fit,
                         lcl =  preddat$fit - 1.96 * preddat$se.fit
)

#Figure
(sav_bogaert <- ggplot(data = frag50, aes(bogaert, SAV)) +
    geom_point(shape = 1) + 
    geom_line(data = mod1_preds, aes(x, predicted), linetype = 1, size = 1) +
    geom_line(data = mod1_preds, aes(x, ucl), linetype = 2, size = 1) +
    geom_line(data = mod1_preds, aes(x, lcl), linetype = 2, size = 1) +
    #geom_smooth(formula = y ~ x, method = "glm", color = "black", se = TRUE) +
    scale_x_reverse() +
    theme(axis.text = element_text(size = 14, family = "Arial", color = "black"),
          axis.title = element_text(size = 14, family = "Arial", color = "black")) +
    labs(x = expression(paste("|",phi, "| `")),
         #x = "Bogaert Fragmentation Scale",
         y = "% Likelihood of SAV"))

ggsave(sav_bogaert, filename = "~/Dropbox/TBone_2019/figures/SAV_bog_lm.png", device = "png")


#### 2. Irradiance ~ Bogaert  ####

#Data
water <- read_xlsx("noaa_blu_objective1_openwater_survey.xlsx")


#Cleanup
#make light attenuation measurement 
sp <- .15
water <- water %>% mutate(
  month          = ifelse(lubridate::month(date) == 4, "April", "September"),
  sal            = as.numeric(sal),
  correct_elev_m = as.numeric(correct_elev_m),
  crms_site      = factor(crms_site, levels = c("369", "311", "345")),
  frag           = toupper(frag),
  subsite_id     = str_c(crms_site, frag, sep = '-'),
  depth_cm       = as.numeric(depth_cm),
  top_par1       = as.numeric(top_par1),                                          #bulb 1 is top bulb
  top_par2       = as.numeric(top_par2),                                          #bulb 2 is low bulb
  bot_par1       = as.numeric(bot_par1),
  bot_par2       = as.numeric(bot_par2),                                                                       
  Kd_surface     = -(log10(top_par1 / top_par2)) / sp,                          #attenuation per .15m for surface measurements
  Kd_bottom      = -(log10(bot_par1 / bot_par2)) / sp,                           #Kd at bottom, per .15m                     #multiply light lost in 15cm by (total depth/15cm increments)
  SI_full        = (bot_par2/top_par1) * 100) %>% 
  filter(month   == "April") 

bog.df <- structure(list(subsite_id = c("369-L", "369-M", "369-H", "311-L", "311-M", "311-H", "345-L", "345-M", "345-H"), 
                         bogaert = c(157.409232427562,158.794835413732, 144.600763346679, 152.825951296738, 147.180047035677, 
                                     139.741186393264, 145.714196114844, 148.353319527527, 139.722070072008)), 
                    class = "data.frame", row.names = c(NA, -9L))


water <- left_join(water, bog.df, by = "subsite_id")
water$SI <- water$SI_full/100

#Model
wmod1 <- lm(SI ~ bogaert, data = water)

#Predictions
x_range <- seq(min(water$bogaert, na.rm = T), 160, by = 2)
preddat <- predict(wmod1, list(bogaert = x_range), type = 'response', se.fit=TRUE)
mod1_preds <- data.frame(x = x_range,
                         predicted = preddat$fit,
                         ucl =  preddat$fit + 1.96 * preddat$se.fit,
                         lcl =  preddat$fit - 1.96 * preddat$se.fit
)

#Figure
(irrad_bogaert <- ggplot(data = water, aes(bogaert, SI)) +
    geom_point(shape = 1) + 
    geom_line(data = mod1_preds, aes(x, predicted), linetype = 1, size = 1) +
    geom_line(data = mod1_preds, aes(x, ucl), linetype = 2, size = 1) +
    geom_line(data = mod1_preds, aes(x, lcl), linetype = 2, size = 1) +
    #geom_smooth(formula = y ~ x, method = "glm", color = "black", se = TRUE) +
    scale_x_reverse() +
    theme(axis.text = element_text(size = 14, family = "Arial", color = "black"),
          axis.title = element_text(size = 14, family = "Arial", color = "black")) +
    labs(x = expression(paste("|",phi, "| `")),
         #x = "Bogaert Fragmentation Scale",
         y = "% Surface Irradiance at Bottom"))

ggsave(irrad_bogaert, filename = "~/Dropbox/TBone_2019/figures/irradiance_bog_lm.png", device = "png")


#### 3. Bathymmetry ~ Bogaert  ####

#Data
range(water$correct_elev_m, na.rm = T) #Bathymetry response

#Cleanup

#Model
bmod1 <- lm(correct_elev_m ~ bogaert, data = water)


#Predictions
x_range <- seq(min(water$bogaert, na.rm = T), 160, by = 2)
preddat <- predict(bmod1, list(bogaert = x_range), type = 'response', se.fit=TRUE)
mod1_preds <- data.frame(x = x_range,
                         predicted = preddat$fit,
                         ucl =  preddat$fit + 1.96 * preddat$se.fit,
                         lcl =  preddat$fit - 1.96 * preddat$se.fit
)

#Figure
(bathym_bogaert <- ggplot(data = water, aes(bogaert, correct_elev_m)) +
    geom_point(shape = 1) + 
    geom_line(data = mod1_preds, aes(x, predicted), linetype = 1, size = 1) +
    geom_line(data = mod1_preds, aes(x, ucl), linetype = 2, size = 1) +
    geom_line(data = mod1_preds, aes(x, lcl), linetype = 2, size = 1) +
    scale_x_reverse() +
    theme(axis.text = element_text(size = 14, family = "Arial", color = "black"),
          axis.title = element_text(size = 14, family = "Arial", color = "black")) +
    labs(
      x = expression(paste("|",phi, "| `")),
      #x = "Bogaert Fragmentation Scale",
      y = "Bathymetry (m)")
  )

#Figure
ggsave(bathym_bogaert, filename = "~/Dropbox/TBone_2019/figures/bathym_bog_lm.png", device = "png")


#### 4. Marsh Elevation ~ Bogaert * Distance from Marsh Edge + (1 | Location) ####

#Data

#Cleanup

#Model

#Predictions

#Figure
