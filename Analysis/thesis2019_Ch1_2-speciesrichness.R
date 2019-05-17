###########
# SPECIES RICHNESS MODELING


# This is the first chunk of modeling after the script import. 


library(ggplot2)
library(viridis)
library(MuMIn)
library(mgcv)
library(gratia) # devtools::install_github('gavinsimpson/gratia')
library(directlabels)

source(here::here("Analysis/thesis2019_Ch1_1-import&cleanup.R"))


#############################
### Data Prep

spp_richness <- allcatch %>% filter(Species != "HYCS" & Species != "UNKN" & Station != 231) %>% group_by(Year) %>% 
  summarise(spp_yr = n_distinct(Species)) 

spp_richness <- spp_richness %>% left_join(pru.env.ann %>% dplyr::select(-Station) %>% 
                                           mutate(Year = as.numeric(Year) + 2000) %>%
                                           group_by(Year) %>% summarise_all(mean), 
                                           by = c("Year" = "Year"))

#Set up standardized variables just in case to use later
spp_richness.std <- as.data.frame(scale(spp_richness)) %>% #scale to mu=0, sd=1
  dplyr::select(-annwinddir)


# Quick EDA
summary(lm(spp_yr ~ Year + annwindspeed_kph + annwinddir + anndisch_cfs + 
             annsal_ppt + anntemp_c + annwinddir_ew, data = spp_richness))

plotenv <- function(env_var){
  env_var1 <- rlang::sym(env_var)
  
  annotateloc <- min(spp_richness %>% dplyr::select(env_var)) +
          (max(spp_richness %>% dplyr::select(env_var)) - min(spp_richness %>% dplyr::select(env_var)))/5
  mod <- summary(glm(paste0("spp_yr ~ ", env_var), data = spp_richness))
  pval <- mod$coef[2,4]
  gof <- 1 - (mod$deviance / mod$null.deviance)

  p <- ggplot(data = spp_richness, aes(x = !! env_var1, y = spp_yr, color = Year)) +
    geom_point(cex = 3) +
    geom_smooth(method = "lm", se=FALSE) +
    scale_color_viridis() + theme_bw() +
    geom_text(label = spp_richness$Year, nudge_y = 0.25, size = 4) +
    annotate("text", x = annotateloc, y = 23, label = paste0("p-value is: ", round(pval,3))) +
    annotate("text", x = annotateloc, y = 23.5, label = paste0("GOF (1-dev/null.dev) is: ", round(gof,2))) 
  return(p)
}

plotenv("Year")
plotenv("annwindspeed_kph")
plotenv("anndisch_cfs")
plotenv("annsal_ppt") #Add variance of salinity?
plotenv("anntemp_c") 
plotenv("annwinddir_ew")


# SUMMARY: Individually, only Year is significant, with temp and salinity marginal. 
# So we conclude that species richness is significantly increasing over time
summary(glm(spp_yr ~ Year, data = spp_richness)) 




# Now let's do some model selection to best explain the results that we see 
options(na.action = "na.fail")
fullmodel <- glm(spp_yr ~ Year + annwindspeed_kph + anndisch_cfs + 
        annsal_ppt + anntemp_c + annwinddir_ew, data = spp_richness)
fullmodel1 <- glm(spp_yr ~ Year + annwindspeed_kph + anndisch_cfs + 
                  annsal_ppt + anntemp_c + Year:anntemp_c, data = spp_richness)
dredge(fullmodel1)

par(mar = c(3,5,6,4))
plot(dredge(fullmodel1), labAsExpr = TRUE)
subset(dredge(fullmodel1), delta < 4)

# SUMMARY: The top model is salinity + Year, though just year works well, and
# temp + wind speed is decent. 


#######################
## Does Species Richness increase over the season? 
## WEEKLY

# The next section will use the function addweeknum() will adds day of year and weeknumber columns
spp_richness.wk <- allcatch %>% filter(Species != "HYCS" & Species != "UNKN" & Station != 231) 
spp_richness.wk <- addweeknum(spp_richness.wk) %>% group_by(Year, week) %>% summarise(num_spp = n_distinct(Species)) 


effort.wk <- addweeknum(effort) %>% filter(Station != 231) %>%  
  group_by(Year, week) %>% summarise(weekeffort = sum(Effort_NetHrs, na.rm = TRUE)) %>%
  mutate(sample_proportion = weekeffort / (48*4*7)) #1344 is 4 nets with two cod ends, fishing 24/7 (7*24*4*2)
  # sample_proportion is the amount of sampling relative to 100% coverage (but can be slightly higher than 100%)
  

spp_richness.wk <- spp_richness.wk %>% left_join(effort.wk %>% dplyr::select(week, sample_proportion), 
                                                 by = c("Year" = "Year", "week" = "week")) %>%
  mutate(num_spp_adj = num_spp * sample_proportion)



ggplot(data = spp_richness.wk, aes(x=week, y = num_spp_adj, group=Year, color=Year)) + 
  geom_line() +
  geom_point() +
  scale_color_viridis() + theme_bw() 


summary(lm(num_spp ~ week + (Year), data=spp_richness.wk)) 



ggplot(addweeknum(pru.env.day) %>% group_by(Year, week) %>% summarise(meansal=mean(Salin_Mid, na.rm = TRUE)), 
       aes(x = week, y = meansal, group = Year, color = Year)) + 
  geom_smooth(method = "loess", se=FALSE, span=1.25)





library(gganimate) # devtools::install_github("thomasp85/gganimate")
library(gifski)
library(transformr)
# these each take about 30 seconds to create. Only use for presentations
richnessgif <- ggplot(spp_richness.wk, mapping = aes(x = week, y = num_spp_adj, group = Year, color = Year)) +
  #geom_line() + geom_point() +
  geom_smooth(method = "loess", se=FALSE, size = 2, span=1.25) +
  scale_x_continuous(breaks = 1:9) +
  scale_color_viridis() + theme_bw() +
  transition_states(Year, transition_length=1.5) +
  shadow_mark(size = 1, colour = 'grey') + ease_aes('linear') + 
  labs(title = 'Species Richness - Year: {closest_state}') + ylab("Weekly Species Richness")

weeklysalinity <- addweeknum(pru.env.day) %>%
  group_by(Year, week) %>% summarise(meansal = mean(Salin_Mid, na.rm = TRUE))

salinitygif <- ggplot(data = weeklysalinity, mapping = aes(x=week, y = meansal, color=Year)) + 
  geom_smooth(method = "loess", se=FALSE, size=2, span=1.25) +
  scale_x_continuous(breaks = 1:9) +
  scale_color_viridis() + theme_bw() +
  transition_states(Year, transition_length=1.5) +
  shadow_mark(size = 1, colour = 'grey') +
  ease_aes('linear') + 
  labs(title = 'Salinity - Year: {closest_state}') + ylab("Weekly Mean Salinity")

# You can visually see that the patterns between the two of them are very similar


# Create side by side. I turned this off because it takes a while to run 
# Guide: https://github.com/thomasp85/gganimate/wiki/Animation-Composition
# library(magick) # You need imagemagick installed
# richnessgif <- image_read(animate(richnessgif + theme(legend.position="none"), width = 300, height = 300))
# salinitygif <- image_read(animate(salinitygif, width = 300, height = 300))
# 
# combined_gif <- image_append(c(richnessgif[1], salinitygif[1]))
# for(i in 2:100){
#   combined <- image_append(c(richnessgif[i], salinitygif[i]))
#   combined_gif <- c(combined_gif, combined)
# }
# combined_gif
#save_animation(combined_gif, "combined.gif")



## SPECIES RICHNESS BIWEEKLY ##



spp_richness.biwk <- allcatch %>% filter(Species != "HYCS" & Species != "UNKN" & Station != 231) %>%
  addbiwknum() %>% group_by(Year, biweekly) %>% summarise(num_spp = n_distinct(Species)) 



effort.biwk <- addbiwknum(effort) %>% filter(Station != 231) %>%  
  group_by(Year, biweekly) %>% summarise(biweeklyeffort = sum(Effort_NetHrs, na.rm = TRUE)) %>%
  mutate(sample_proportion = biweeklyeffort / (48*4*15)) #2880 is 4 nets with two cod ends, fishing 24hrs for 15 days (15*24*4*2)
# sample_proportion is the amount of sampling relative to 100% coverage (but can be slightly higher than 100%)
# We're using 15 days as "full effort" sampling but e.g., July 16-31 is 16 days. This is why we scale for effort though


spp_richness.biwk <- spp_richness.biwk %>% left_join(effort.biwk %>% dplyr::select(biweekly, sample_proportion), 
                                                 by = c("Year" = "Year", "biweekly" = "biweekly")) %>%
  mutate(num_spp_adj = num_spp * sample_proportion)



ggplot(data = spp_richness.biwk, aes(x=biweekly, y = num_spp_adj, group=Year, color=Year)) + 
  geom_line() + geom_point() + scale_color_viridis() + theme_bw() 

ggplot(data = spp_richness.biwk, aes(x=Year, y = num_spp_adj, group=biweekly, color=biweekly)) + 
  geom_line() + geom_point() + scale_color_viridis() + theme_bw()


summary(glm(num_spp ~ biweekly + (Year), data=spp_richness.biwk)) 
summary(gam(num_spp ~ biweekly + (Year), data=spp_richness.biwk, gamma = 1.4))
summary(gam(num_spp ~ I(biweekly^2) + (Year), data=spp_richness.biwk, gamma = 1.4))
summary(gam(num_spp ~ biweekly + s(Year, k=-1), data=spp_richness.biwk, gamma = 1.4))
summary(gam(num_spp ~ biweekly^2 , data=spp_richness.biwk, gamma = 1.4))

summary(gam(num_spp ~ s(biweekly, k=3) + Year, data=spp_richness.biwk, gamma = 1.4))
testgam <- (gam(num_spp ~ s(Year, biweekly, k=6), data=spp_richness.biwk, gamma = 1.4))


dredge(gam(num_spp ~ s(biweekly, k=3) + Year, data=spp_richness.biwk, gamma = 1.4)) #TOP MODEL
dredge(gam(num_spp ~ biweekly*Year, data=spp_richness.biwk, gamma = 1.4))
dredge(gam(num_spp ~ biweekly + Year + I(biweekly^2), data=spp_richness.biwk, gamma = 1.4))
dredge(gam(num_spp ~ biweekly + s(Year), data=spp_richness.biwk, gamma = 1.4))


topmod.spprich <- gam(num_spp ~ s(biweekly, k=3) + Year, data=spp_richness.biwk, gamma = 1.4)


summary(topmod.spprich)
gam.check(topmod.spprich)
plot.gam(topmod.spprich, se=TRUE)
vis.gam(topmod.spprich, plot.type="contour", color="terrain")

vis.gam(testgam, plot.type="contour", color="terrain")


# Top model is species richness ~ s(biweekly, k=3) + Year
# GAM fits better than GLM
# We conclude that there is significant curvature in the model

gratia::draw(topmod.spprich)
gratia::appraise(topmod.spprich)


## Plotting top model
# Modified from pkg "gratia" https://github.com/gavinsimpson/gratia/blob/master/R/data-slice.R
modelpredict.gam <- function(object, var1, var2, n = 50) {
  #object is a gam model object, var1 is variable 1 (x axis), must be numeric
  #var2 is variable 2 (y axis), must be numeric, n is length out for pred intervals
  mf <- fix_offset(topmod.spprich, model.frame(topmod.spprich)) # offset not necessary for most models
  data1 <- seq(from = min(mf %>% dplyr::select(var1)), to = max(mf %>% dplyr::select(var1)), length.out = n)
  data2 <- seq(from = min(mf %>% dplyr::select(var2)), to = max(mf %>% dplyr::select(var2)), length.out = n)
  
  results <- crossing(data1, data2)
  names(results)[1:2] <- c(var1, var2)
  results$pred.gam <- predict.gam(topmod.spprich, newdata = results)
  print(results)
  
}


ggplot(modelpredict.gam(topmod.spprich, "Year", "biweekly"), aes(x=Year, y=biweekly, z=pred.gam)) + 
  geom_raster(aes(fill = pred.gam)) + geom_contour() +
  geom_dl(aes(label=..level..), method = list("top.pieces", cex=0.75), 
          stat="contour", breaks = seq(from=13, to=18.5, by=0.5)) +
  scale_fill_gradientn(colours= terrain.colors(6)) +
  scale_x_continuous(breaks = seq(from=2001, to=2018, by=2)) +
  labs(x="", y="Biweekly", fill = "Predicted\nSpecies\nRichness") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=12), 
        axis.text.x = element_text(angle = 35, hjust = 1)) 
#ggsave("plotexports/Fig_spprich.png", dpi = 300, width = 6.5, height = 4.33)




######## Station by Station effects, not used at this time but good summary to know
spp_richness.biwk.stn <- allcatch %>% filter(Species != "HYCS" & Species != "UNKN" & Station != 231) %>%
  addbiwknum() %>% group_by(Year, biweekly, Station) %>% summarise(num_spp = n_distinct(Species)) 


topmod.spprich220 <- gam(num_spp ~ s(biweekly, k=3) + Year, data=spp_richness.biwk.stn %>% filter(Station==220), gamma = 1.4)
topmod.spprich230 <- gam(num_spp ~ s(biweekly, k=3) + Year, data=spp_richness.biwk.stn %>% filter(Station==230), gamma = 1.4)
topmod.spprich214 <- gam(num_spp ~ s(biweekly, k=3), data=spp_richness.biwk.stn %>% filter(Station==214), gamma = 1.4)
topmod.spprich218 <- gam(num_spp ~ s(Year, k=4) + s(biweekly, k=3), data=spp_richness.biwk.stn %>% filter(Station==218), gamma = 1.4)

summary(topmod.spprich220)
summary(topmod.spprich230)
summary(topmod.spprich214)
summary(topmod.spprich218)

vis.gam(topmod.spprich220, plot.type="contour", color="terrain")
vis.gam(topmod.spprich218, plot.type="contour", color="terrain")
vis.gam(topmod.spprich230, plot.type="contour", color="terrain")
plot(topmod.spprich214, plot.type="contour", color="terrain")





