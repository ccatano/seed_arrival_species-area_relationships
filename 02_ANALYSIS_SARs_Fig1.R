##------------------------------------------------------------------------------
## Author: Christopher P. Catano                                               
##                                                                             
## Description: Model effects of species pool size on species-area relationships
##              Produces results and Figure 1 from manuscript
##              
## Data required: SAR_data_allspecies.csv    
##                SAR_data_seededspecies_nocx.csv
##                SAR_data_nonseededspecies_nocx.csv
##------------------------------------------------------------------------------

rm(list = ls())

# Check for and install required packages
for (package in c('tidyverse', 'broom', 'lme4', 'lmerTest', 'emmeans', 'ggpubr',
                  'visreg', 'bootpredictlme4', 'emmeans', 'lattice')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}


#### 1. IMPORT AND JOIN DATA SETS ####

# All species
data.all <- read.csv("output/SAR_data_allspecies.csv", row.names = 1)

# Species present in the seed additions
data.seeded <- read.csv("output/SAR_data_seededspecies_nocx.csv", row.names = 1)

# Species not included in the seed additions
data.nonseeded <- read.csv("output/SAR_data_nonseededspecies_nocx.csv", row.names = 1)


#### 2. MODEL SPECIES POOL EFFECTS ON SPECIES-AREA RELATIONSHIPS ####

#---------------
# 2.1) All species
data.all$log.area <- log(data.all$area)
data.all$species.richness.trt <- as.factor(data.all$species.richness.trt)

z.mod.all <- lmer(richness ~ species.richness.trt*log.area + (1 | Site), data = data.all)
summary(z.mod.all, ddf = "Kenward-Roger")

# plots the Pearson residuals, residuals scaled by variance function, verses the 
# fitted values on the response scale.
plot(z.mod.all)
qqmath(z.mod.all, id = 0.05)
dotplot(ranef(z.mod.all, condVar = TRUE))

# Get slope of SARs for each treatment, and fixed effect
emtrends(z.mod.all, pairwise ~ species.richness.trt, var = "log.area")

# Get pairwise contrast of intercepts (differences in richness at 1 m2)
(emm <- emmeans(z.mod.all, "species.richness.trt", at = list(log.area = 0)))
pairs(emm)


#--------------------  
# 2.2) seeded species
data.seeded$log.area <- log(data.seeded$area)
data.seeded$species.richness.trt <- as.factor(data.seeded$species.richness.trt)

z.mod.seeded <- lmer(richness.seeded.nocx ~ species.richness.trt*log.area + (1 | Site), data = data.seeded)
summary(z.mod.seeded, ddf = "Kenward-Roger")

# plots the Pearson residuals, residuals scaled by variance function, verses the 
# fitted values on the response scale.
plot(z.mod.seeded)
qqmath(z.mod.seeded, id = 0.05)
dotplot(ranef(z.mod.seeded, condVar = TRUE))

# Get slope of SARs for each treatment, and fixed effect
emtrends(z.mod.seeded, pairwise ~ species.richness.trt, var = "log.area")

# Get pairwise contrast of intercepts (differences in richness at 1 m2)
(emm <- emmeans(z.mod.seeded, "species.richness.trt", at = list(log.area = 0)))
pairs(emm)


#-------------------
# 2.3) Non-seeded species
data.nonseeded$log.area <- log(data.nonseeded$area)
data.nonseeded$species.richness.trt <- as.factor(data.nonseeded$species.richness.trt)

z.mod.nonseeded <- lmer(richness.nonseeded.nocx ~ species.richness.trt*log.area + (1 | Site), data = data.nonseeded)
summary(z.mod.nonseeded, ddf = "Kenward-Roger")

# plots the Pearson residuals, residuals scaled by variance function, verses the 
# fitted values on the response scale.
plot(z.mod.nonseeded)
qqmath(z.mod.nonseeded, id = 0.05)
dotplot(ranef(z.mod.nonseeded, condVar = TRUE))

# Get slope of SARs for each treatment, and fixed effect
emtrends(z.mod.nonseeded, pairwise ~ species.richness.trt, var = "log.area")

# Get pairwise contrast of intercepts (differences in richness at 1 m2)
(emm <- emmeans(z.mod.nonseeded, "species.richness.trt", at = list(log.area = 0)))
pairs(emm)


#### 3. PLOT MODEL RESULTS ####

theme_set(theme_bw() +  
            theme(legend.position = "none", panel.grid.minor = element_blank(), 
                  panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5), 
                  text = element_text(size = 12), axis.text = element_text(size = 10)))

(all.spe <- visreg(z.mod.all, "log.area", by = "species.richness.trt", gg = T,
                   overlay = T, 
                   partial = F,
                   rug = F,
                   options(bootnsim = 100)) +
    xlab("ln(Area)") +
    ylab("Species richness") +
    scale_color_manual(values=c("darkorange", "steelblue3")) +
    ggtitle("All species") +
    theme(axis.title.y = element_blank())
)

(seeded.spe <- visreg(z.mod.seeded, "log.area", by = "species.richness.trt", gg = T,
                      overlay = T, 
                      partial = F,
                      rug = F,
                      options(bootnsim = 100)) +
    xlab("ln(Area)") +
    ylab("Species richness") +
    scale_color_manual(values=c("darkorange", "steelblue3")) +
    ggtitle("Seeded") +
    theme(axis.title.y = element_blank())
)

(nonseeded.spe <- visreg(z.mod.nonseeded, "log.area", by = "species.richness.trt", gg = T,
                      overlay = T, 
                      partial = F,
                      rug = F,
                      options(bootnsim = 100)) +
    xlab("ln(Area)") +
    ylab("Species richness") +
    scale_color_manual(values=c("darkorange", "steelblue3")) +
    ggtitle("Non-seeded") +
    theme(axis.title.y = element_blank())
)

# Combine plots to make Figure 1
plot <- ggarrange(all.spe, seeded.spe, nonseeded.spe, labels = c("(A)", "(B)", "(C)"), 
          ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")

#jpeg("figures/SAR_models_ancova_fig.jpg", width = 7, height = 3, units = "in", res = 600)
annotate_figure(plot,
                left = text_grob("Species richness", rot = 90),
                
)
#dev.off()