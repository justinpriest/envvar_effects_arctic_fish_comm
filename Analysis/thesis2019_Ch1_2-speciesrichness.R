###########
# SPECIES RICHNESS MODELING


# This is the first chunk of modeling after the script import. 



library(ggplot2)
library(viridis)

source("Analysis/thesis2019_Ch1_1-import&cleanup.R")


#############################
### Data Prep

spp_richness <- allcatch %>% filter(Species != "HYCS" & Species != "UNKN" & Station != 231) %>% group_by(Year) %>% 
  summarise(spp_yr = n_distinct(Species)) 

#summary(lm(spp_yr ~ Year, data = spp_richness))
spp_richness <- spp_richness %>% left_join(pru.env.ann %>% dplyr::select(-Station) %>% 
                                             mutate(Year = as.numeric(Year) + 2000) %>%
                                             group_by(Year) %>% summarise_all(mean), 
                                           by = c("Year" = "Year"))

summary(lm(spp_yr ~ Year + annwindspeed_kph + annwinddir + anndisch_cfs + annsal_ppt + anntemp_c + annwinddir_ew, data = spp_richness))

ggplot(data=spp_richness, aes(x=anntemp_c, y = spp_yr, color = Year)) + 
  geom_point(cex = 3) + geom_smooth(method = "lm", se=FALSE) +
  scale_color_viridis() + theme_bw() +
  geom_text(label = spp_richness$Year, nudge_y = 0.25, size =4)



summary(lm(spp_yr ~ Year + annsal_ppt - 1, data = spp_richness))


