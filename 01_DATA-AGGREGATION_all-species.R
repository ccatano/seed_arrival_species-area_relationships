##------------------------------------------------------------------------------
## Author: Christopher P. Catano                                               
##                                                                             
## Description: Aggregate species richness at each spatial scale for all species,
##              (arrived from the pool and those that recruited naturally)
##              
## Data:     CLE_data_clean.csv  
##           CLE_treatment_codes.csv
##------------------------------------------------------------------------------

rm(list = ls())

# Check for and install/load required packages
for (package in c('tidyverse')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}


#### 1. IMPORT DATA ####
# Raw species survey data
data.spe <- read.csv("raw_data/CLE_data_clean.csv")
data.spe <- data.spe[, -1]

col.nam <- c("Plot1", "Plot2", "Plot3", "Plot4", "Plot5")
data.spe[col.nam] <- sapply(sapply(data.spe[col.nam], as.character), as.numeric)
data.spe$Species <- as.character(data.spe$Species)
data.spe[ , col.nam][data.spe[ , col.nam] > 0] <- 1


#### 2. CALCULATE RICHNESS AT EACH SCALE ####

# 2.1) CALCULATE MEAN RICHNESS AT 1-m2
# Subset data to only include 1-m2 quadrat occurrences
m1_scale <- data.spe %>%
  filter(Scale == "1") 

# Sum total unique occurrences in each quadrat, then take mean occurrence across 
# all 5 quadrats on the transect
data_1m_rich <- m1_scale %>%
  select(-c(Species, Scale)) %>%
  group_by(Site, Half) %>% 
  summarise_each(funs(sum)) %>%
  gather(Plot, presence, Plot1:Plot5) %>%
  group_by(Site, Half) %>%
  summarize(rich1m = mean(presence, na.rm = T))


#-------------------------------------------------------------------------------
# 2.2) CALCULATE MEAN RICHNESS AT 25-m2

# Subset data to only include 1-m2 quadrat occurrences plus additional species 
# detected in the 25-m2 walkthrough around each quadrat
m5_scale <- data.spe %>%
  filter(Scale == "5" | Scale == "1") 

# Sum total unique occurrences in each 25-m2 sample, then take mean occurrence  
# across all 5 samples on the transect
data_5m_rich <- m5_scale %>%
  select(-c(Species, Scale)) %>%
  group_by(Site, Half) %>% 
  summarise_each(funs(sum)) %>%
  gather(Plot, presence, Plot1:Plot5) %>%
  group_by(Site, Half) %>%
  summarize(rich5m = mean(presence, na.rm = T))


#-------------------------------------------------------------------------------
# 3.3) CALCULATE MEAN RICHNESS AT 125-m2 (transect)

# Sum total occurrences of each species across all 5 25-m2 samples, then replace 
# all values > 1 (species occurred in multiple 25-m2 samples) with 1 to indicate 
# unique species occurrences at the transect scale. Sum unique species to get 
# transect richness
data_transect_rich <- m5_scale %>%
  gather(Plot, presence, Plot1:Plot5) %>%
  group_by(Site, Half, Species) %>%
  summarize(speocc = sum(presence, na.rm = T)) %>%
  mutate(speocc = replace(speocc, speocc > 0, 1)) %>%
  group_by(Site, Half) %>%
  summarize(rich25m = sum(speocc, na.rm = T))


#-------------------------------------------------------------------------------
# 2.4) CALCULATE MEAN RICHNESS AT SITE SCALE (actual area varies by site)

# Calculate new species observed in walkthrough of each treatment not already 
# observed within 1-m2 quadrats or 25-m2 walkthroughs along the transect. Then 
# merge with species richness calculated at 125-m2 (transect) scale
site_scale <- data.spe %>%
  filter(Scale == "Half Site") %>%
  select(-c(Plot2, Plot3, Plot4, Plot5)) %>%
  rename(presence = Plot1) %>%
  group_by(Site, Half) %>%
  summarize(Addrichness = sum(presence, na.rm = T)) %>%
  left_join(data_transect_rich, by = c("Site", "Half"))

# Add new additional species richness to transect richness to get total richness 
# at the site scale for each site x treatment combo
site_scale$siteRich <- site_scale$rich25m + site_scale$Addrichness

# Drop unnecessary columns
site_scale <- site_scale %>%
  select(-c(Addrichness)) %>%
  select(-c(rich25m))


#-------------------------------------------------------------------------------
#### 3. JOIN RICHNESS AT EACH SCALE INTO DATA SET FOR ANALYSIS  ####

SAR_allspecies <- data_1m_rich %>%
  left_join(data_5m_rich, by = c("Site", "Half")) %>%
  left_join(data_transect_rich, by = c("Site", "Half")) %>%
  left_join(site_scale, by = c("Site", "Half"))  %>%
  gather(scale, richness, rich1m:siteRich)

# Import treatment ID and site scale area
trt.data <- read.csv("raw_data/CLE_treatment_codes.csv")

SAR_allspecies <- SAR_allspecies %>%
  left_join(trt.data, by = c("Site", "Half", "scale")) %>%
  select(-c(Half, scale))

# Output data for analysis
write.csv(SAR_allspecies, "output/SAR_data_allspecies.csv")
