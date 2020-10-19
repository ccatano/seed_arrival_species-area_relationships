##------------------------------------------------------------------------------
## Author: Christopher P. Catano                                               
##                                                                             
## Description: Model effects of species pool size on species composition-soil
##              moisture relationships and spatial variation in species occupancy.
##              Produces FIGURE 2 from manuscript
##              
## Data:  CLE_data_clean.csv
##        CLE_treatment_codes.csv
##        seeded_species.csv
##        soil_moisture_CLE.csv
##------------------------------------------------------------------------------

rm(list = ls())

# Check for and install/load required packages
for (package in c('tidyverse','vegan', 'lme4', 'lmerTest', 'lattice', 'ggpubr')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}


#### 1. IMPORT & SUBSET DATA ####

# Import soil moisture data for each plot
soil <- read.csv("raw_data/soil_moisture_CLE.csv")

# Species composition data
data.spe <- read.csv("raw_data/CLE_data_clean.csv", row.names = 1)

# Make data at all sampling scales presence-absence
col.nam <- c("Plot1", "Plot2", "Plot3", "Plot4", "Plot5")
data.spe[col.nam] <- sapply(sapply(data.spe[col.nam], as.character), as.numeric)
data.spe$Species <- as.character(data.spe$Species)
data.spe[ , col.nam][data.spe[ , col.nam] > 0] <- 1

# Import treatment ID and site scale area
trt.data <- read.csv("raw_data/CLE_treatment_codes.csv")
trt.codes <- unique(trt.data[,c("Site", "Half", "species.richness.trt")])

# Join treatment codes with species composition data
data.spe <- data.spe %>%
  left_join(trt.codes, by = c("Site", "Half")) 

# Import and subset seeded species 
data.seed <- read.csv("raw_data/seeded_species.csv")

# Species from seed mix including codes for unknown sedges: CxGEN1 - CXGEN5)
data.seedHigh.nocx <- as.character(data.seed[data.seed$high.no.cxgen == "yes", "Species"])
data.seedLow <- as.character(data.seed[data.seed$low.rich.trt == "yes", "Species"])

# Subset to treat all unknown sedges as if not from seed mix
seeded.nocx <- data.spe[(data.spe$species.richness.trt == "High" & data.spe$Species %in% data.seedHigh.nocx) |
                         (data.spe$species.richness.trt == "Low" & data.spe$Species %in% data.seedLow) , ]
seeded.nocx <- seeded.nocx %>%
  select(-c(species.richness.trt))

# Subset data to only include 1-m2 quadrat occurrences plus additional species 
# detected in the 5-m2 walkthrough around each quadrat
m5_scale <- seeded.nocx %>%
  filter(Scale == "5" | Scale == "1") 

# sum counts for each species across 1 and 5m scales for each site half to get
# occurrence in each 5x5.
data_5m_comp <- m5_scale %>%
  select(-c(Scale)) %>%
  group_by(Species, Site, Half) %>% 
  summarise_each(funs(sum)) %>%
  gather(Plot, presence, Plot1:Plot5) %>%
  spread(key = Species, value = presence) %>%
  mutate_all(~replace(., is.na(.), 0))
  
# Calculate mean soil moisture over the multiple samples for each plot and merge
# with composition data
soil.mean <- soil %>%
  select(-c(date)) %>%
  group_by(Site, Half) %>% 
  summarise_each(funs(mean)) %>%
  gather(Plot, soil.moisture, Plot1:Plot5) %>%
  left_join(data_5m_comp, by = c("Site", "Half", "Plot"))
  
spe.env.data <- trt.data %>%
  select(-c("scale", "area")) %>%
  distinct() %>%
  left_join(soil.mean, by = c("Site", "Half")) %>%
  select(-c(Half))
data.spe.env <- spe.env.data[, -c(2)]


#### 2. MODEL SPECIES COMPOSITION-SOIL MOISTURE RELATIONSHIPS ####

# Combine site and treatment columns to loop over
data.2 <- data.spe.env %>%
  unite(Site_trt, Site, species.richness.trt, sep = "_")

# Calculate spatial variation in species occupancy for each replicate
replicates <- unique(data.2$Site_trt)
replicates <- as.data.frame(replicates)
LL <- lapply(1:nrow(replicates), 
             function(i) {
               site <- replicates[i, 1]
               data = data.2[data.2$Site_trt == site, ]	
               comp = data[,4:49]
               env = as.matrix(data[,3])
               dbRDA = capscale(comp ~ env, dist = "jaccard")
               RsquareAdj(dbRDA)[1]
                            }
)

result <- do.call(rbind, LL)
result <- as.data.frame(result)
result2 <- cbind(replicates, result)

RDA <- result2 %>%
  separate(col = "replicates", into = c("Site", "species.richness.trt"), sep = "_")

# Model effects of species pool size on composition-soil moisture relationships
RDA$species.richness.trt <- as.factor(RDA$species.richness.trt)
RDA$r.squared <- as.numeric((RDA$r.squared))
mod <- lmer(sqrt(r.squared) ~ species.richness.trt + (1 | Site), data = RDA)
summary(mod, ddf = "Kenward-Roger")
plot(mod, sqrt(abs(residuals(.))) ~ fitted(.))
qqmath(mod, id = 0.05)
dotplot(ranef(mod, condVar = TRUE))

# Save predicted values
RDA$eu.pred <- predict(mod, re.form = NA) 
RDA$species.richness.trt <- factor(RDA$species.richness.trt, levels = c("Low", "High"))

# Set plot theme
theme_set(theme_bw() +  
            theme(legend.position = "none", panel.grid.minor = element_blank(), 
                  panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5), 
                  text = element_text(size = 12), axis.text = element_text(size = 10)))

# plot result
(eu.pred <- ggplot(RDA, aes(x = species.richness.trt, y = sqrt(r.squared))) +
    geom_point(aes(color = species.richness.trt), size = 2, alpha = 0.5) + 
    geom_line(aes(group = Site), size = 0.25, color = "black", alpha = 0.2) +
    geom_line(colour = "black", aes(y = eu.pred, group = Site), size = 1, linetype = "dotted") +
    scale_x_discrete(labels = c('Small', 'Large')) +
    scale_color_manual(values=c("steelblue3","darkorange")) +
    xlab("Species pool size") +
    ylab(expression(paste("Composition-soil moisture (", italic(R^2),")"))) 
) 


#### 3. MODEL SPATIAL VARIATION IN SPECIES OCCUPANCY ####

# Calculate plot occupancy (proportion of plots on a transect occupied) by each
# core species (species common to both seed mixes) between treatments
core <- c("ANDGER", "BOUCUR", "ELYCAN", "KOEMAC", "SCHSCO", "CHAFAS", "CORLAN",
          "DALPUR", "ECHPUR", "LESCAP", "RATPIN", "RUDHIR")
m5_scale_core <- m5_scale[m5_scale$Species %in% core, ]

data_5m_core <- m5_scale_core %>%
  select(-c(Scale)) %>%
  group_by(Species, Site, Half) %>% 
  summarise_each(funs(sum)) %>%
  gather(Plot, presence, Plot1:Plot5) %>%
  spread(key = Species, value = presence) %>%
  mutate_all(~replace(., is.na(.), 0))

# Calculate differences in spatial aggregation of these species
data_5m_core2 <- trt.codes %>%
  left_join(data_5m_core, by = c("Site", "Half")) %>%
  select(-c(Half))

data_5m_core3 <- data_5m_core2 %>%
  unite(Site_trt, Site, species.richness.trt, sep = "_")

spe <- data_5m_core3[, c(3:12)]
spe$sum <- rowSums(spe)
trt <- data_5m_core3[, 1]
dis <- vegdist(spe, method = "jaccard")
beta <- betadisper(dis, group = trt)

df <- data.frame(spe.dispersion = beta$distances, 
                 group = beta$group)
df2 <- df %>%
group_by(group) %>%
  summarise(mean = mean(spe.dispersion))

df2 <- df2 %>%
  separate(col = "group", into = c("Site", "species.richness.trt"), sep = "_")
df2$species.richness.trt <- as.factor(df2$species.richness.trt)

# Model species pool effect on spatial variation in species occupancy
mod.disp <- lmer(mean ~ species.richness.trt + (1 | Site), data = df2)
summary(mod.disp, ddf = "Kenward-Roger")
plot(mod.disp, sqrt(abs(residuals(.))) ~ fitted(.))
qqmath(mod.disp, id = 0.05)
dotplot(ranef(mod.disp, condVar = TRUE))

# Save predicted values
df2$disp_pred <- predict(mod.disp, re.form = NA) 
df2$species.richness.trt <- factor(df2$species.richness.trt, levels = c("Low", "High"))

# plot result
(disp.pred <- ggplot(df2, aes(x = species.richness.trt, y = mean)) +
    geom_point(aes(color = species.richness.trt), size = 2, alpha = 0.5) + 
    geom_line(aes(group = Site), size = 0.25, color = "black", alpha = 0.2) +
    scale_color_manual(values=c("darkorange", "steelblue3")) +
    geom_line(colour = "black", aes(y = disp_pred, group = Site), size = 1) +
    scale_x_discrete(labels = c('Small', 'Large')) +
    xlab("Species pool size") +
    ylab("Spatial variation in composition \n of core seeded species")
)


#### 4. Make FIGURE 2 from manuscript ####
tiff("figures/newfigs_final.tiff", width = 5, height = 3, units = "in", res = 600)
ggarrange(eu.pred, disp.pred, ncol = 2, nrow = 1)
dev.off()
