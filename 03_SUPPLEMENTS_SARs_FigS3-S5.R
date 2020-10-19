##------------------------------------------------------------------------------
## Author: Christopher P. Catano                                               
##                                                                             
## Description: Plot species-area relationships for each site x treatment combo.
##              Produces Supplemental figures S3-S5
##              
## Data:  SAR_data_allspecies_sampleplots.csv     
##        SAR_data_seeded_sampleplots.csv
##        SAR_data_nonseededspecies_sampleplots.csv
##------------------------------------------------------------------------------

rm(list = ls())

# Check for and install/load required packages
for (package in c('tidyverse')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}


#### 1. IMPORT DATA SETS ####

# All species
data.all <- read.csv("output/SAR_data_allspecies_sampleplots.csv", row.names = 1)

# Species present in the seed additions
data.seeded <- read.csv("output/SAR_data_seeded_sampleplots.csv", row.names = 1)

# Species not included in the seed additions
data.nonseeded <- read.csv("output/SAR_data_nonseededspecies_sampleplots.csv", row.names = 1)


#### 2. Visualize the species-area relationships for each site ####

theme_set(theme_bw() +  
            theme(legend.position = "bottom", legend.title.align = 0.5, legend.box.just = "center",
                  panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                  plot.title = element_text(hjust = 0.5), 
                  text = element_text(size = 12), axis.text = element_text(size = 10)))

# FIGURE S3
(SARs_all <- ggplot(data = data.all, aes(x = log(area), y = richness, color = species.richness.trt)) + 
   geom_smooth(method = "lm", se = F, size = 0.5, aes(group = species.richness.trt, 
                                                       colour = species.richness.trt)) +
   geom_point(alpha = 0.5) +
   facet_wrap(~Site) +
   labs(title = "All species") +
   ylab("Species richness") +
   xlab(expression(paste("Area (log"[e]^{}, " m"^2, ")"))) +
   scale_color_manual(name = "Species pool size",
                      labels = c("Large (70 species)", "Small (12 species)"),
                      values=c("darkorange", "steelblue3"))
)

# FIGURE S4
(SARs_seeded <- ggplot(data = data.seeded, aes(x = log(area), y = richness, color = species.richness.trt)) + 
    geom_smooth(method = "lm", se = F, size = 0.5, aes(group = species.richness.trt, 
                                                        colour = species.richness.trt)) +
    geom_point(alpha = 0.5) +
    facet_wrap(~Site) +
    labs(title = "Seeded species") +
    ylab("Species richness") +
    xlab(expression(paste("Area (log"[e]^{}, " m"^2, ")"))) +
    scale_color_manual(name = "Species pool size",
                       labels = c("Large (70 species)", "Small (12 species)"),
                       values=c("darkorange", "steelblue3")) 
)

# FIGURE S5
(SARs_resident <- ggplot(data = data.nonseeded, aes(x = log(area), y = richness, color = species.richness.trt)) + 
    geom_smooth(method = "lm", se = F, size = 0.5, aes(group = species.richness.trt, 
                                                        colour = species.richness.trt)) +
    geom_point(alpha = 0.5) +
    facet_wrap(~Site) +
    labs(title = "Non-seeded species") +
    ylab("Species richness") +
    xlab(expression(paste("Area (log"[e]^{}, " m"^2, ")"))) +
    scale_color_manual(name = "Species pool size",
                       labels = c("Large (70 species)", "Small (12 species)"),
                       values=c("darkorange", "steelblue3")) 
)