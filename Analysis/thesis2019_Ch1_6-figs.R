#Figures and Tables


library(corrplot)

# Later can make an RDA file to save output from all of these
source("Analysis/thesis2019_Ch1_2-speciesrichness.R")
source("Analysis/thesis2019_Ch1_3-rarespecies.R")
source("Analysis/thesis2019_Ch1_4-PERMANOVA.R")
source("Analysis/thesis2019_Ch1_5-TimeSeries.R")


#load(.RData)







env.biwk.corr <- cor(pru.env.biwk %>% ungroup() %>% dplyr::select(-Year, -biweekly, -Station), use = "complete.obs")
env.ann.corr <- cor(pru.env.ann %>% ungroup() %>% dplyr::select(-Year,  -Station), use = "complete.obs")


corrplot.mixed(env.biwk.corr)
corrplot.mixed(env.ann.corr)

env.biwk.corr[upper.tri(env.biwk.corr)] <- 0
lower.tri(env.biwk.corr, diag = FALSE)
