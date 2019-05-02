#Figures and Tables


library(corrplot)

# Later can make an RDA file to save output from all of these
source("Analysis/thesis2019_Ch1_2-speciesrichness.R")
source("Analysis/thesis2019_Ch1_3-rarespecies.R")
source("Analysis/thesis2019_Ch1_4-PERMANOVA.R")
source("Analysis/thesis2019_Ch1_5-TimeSeries.R")




env.biwk.corr <- cor(pru.env.biwk %>% ungroup() %>% dplyr::select(-Year, -biweekly, -Station), use = "complete.obs")
env.ann.corr <- cor(pru.env.ann %>% ungroup() %>% dplyr::select(-Year,  -Station), use = "complete.obs")


corrplot.mixed(env.biwk.corr)
corrplot.mixed(env.ann.corr)

env.biwk.corr[upper.tri(env.biwk.corr)] <- 0
lower.tri(env.biwk.corr, diag = FALSE)



catchtable <- as_tibble(expand.grid(Year=2001:2018, 
                                    biweekly=1:4, 
                                    Species=unique(catchenviron$Species))) %>% 
  left_join(catchenviron, by = c("Species" = "Species", "Year" = "Year", "biweekly" = "biweekly")) %>% 
  # this first part makes sure to add in blank Year/biweek combos for species which are not frequently caught
  group_by(Species, Year, biweekly) %>% summarise(biwkCatch = sum(totcount)) %>% 
  mutate(biwkCatch = replace_na(biwkCatch, 0)) %>% #replace NAs with 0
  left_join(spplookup, by = c("Species" = "Species")) %>% ungroup() %>% dplyr::select(commonname, MeanCatch) %>%
  arrange(commonname) %>% rename(Species = commonname)



catchtable
