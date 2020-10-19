##------------------------------------------------------------------------------
## Author: Christopher P. Catano                                               
##                                                                             
## Description: Quantify and test for differences in abundance distributions
##              of species pool treatments (seed mixes)
##              
## Data:  seed_mix_SAD.csv                                                                      
##------------------------------------------------------------------------------

rm(list = ls())

# Check for and install required packages
for (package in c('tidyverse', 'lme4', 'lmerTest', 'lattice', 'ggpubr',
                  'codyn', 'vegan')) {
  if (!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}


#-------------------------------------------------------------------------------
# 1. IMPORT DATA
pools <- read.csv("raw_data/seed_mix_SAD.csv", stringsAsFactors = F)
str(pools)

pools <- pools %>%
  filter(!(Note == "unknown" | Note == "not seeded?")) %>%
  select(-c("Note", "oz.A", "lb.A"))

pools2 <- pools %>%
  group_by(seed.diversity) %>%
  arrange(seeds.m2) %>%
  mutate(total.seeds = sum(seeds.m2)) %>%
  mutate(proportion = seeds.m2 / total.seeds) %>%
  mutate(cumul = cumsum(proportion)) %>%
  mutate(log.per.abund = log10(proportion)) %>%
  mutate(log.abund = log10(seeds.m2+1)) %>%
  arrange(seed.diversity)

ranklow <- as.matrix(c(12:1))
rankhigh <- as.matrix(c(69:1))
rank <- rbind(rankhigh, ranklow) 
pools2$rank <- rank

theme_set(theme_bw() +  
            theme(legend.position = "bottom", legend.title.align = 0.5, legend.box.just = "center",
                  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
                  plot.title = element_text(hjust = 0.5), 
                  text = element_text(size = 12), axis.text = element_text(size = 10)))

# Empirical cumulative distribution function 
ggplot(pools2, aes(x = log.per.abund, y = cumul, colour = seed.diversity)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_line() +
  scale_color_manual(name = "Species pool size",
                     labels = c("High (70 species)", "Low (12 species)"),
                     values=c("darkorange", "steelblue3")) +
  xlab("Log(seeded abundance)") +
  ylab("Proportion") + 
  ggsave("figures/pool_ECDF.jpeg")


spe.site <- pools %>%
  select(-c(functional.group, focal)) %>%
  spread(key = species, value = seeds.m2) %>%
  mutate_all(~replace(., is.na(.), 0))

# Pielou's evenness: Shannon entropy / log(observed richness)
spe.site$S <- specnumber(spe.site[,2:70])
(spe.site$J <- (diversity(spe.site[,2:70])/log(spe.site$S)))
