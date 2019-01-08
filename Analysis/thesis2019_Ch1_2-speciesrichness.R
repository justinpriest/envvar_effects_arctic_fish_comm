###########
# SPECIES RICHNESS MODELING


# This is the first chunk of modeling after the script import. 


library(ggplot2)
library(viridis)
library(MuMIn)

source("Analysis/thesis2019_Ch1_1-import&cleanup.R")


#############################
### Data Prep

spp_richness <- allcatch %>% filter(Species != "HYCS" & Species != "UNKN" & Station != 231) %>% group_by(Year) %>% 
  summarise(spp_yr = n_distinct(Species)) 

spp_richness <- spp_richness %>% left_join(pru.env.ann %>% dplyr::select(-Station) %>% 
                                             mutate(Year = as.numeric(Year) + 2000) %>%
                                             group_by(Year) %>% summarise_all(mean), 
                                           by = c("Year" = "Year"))

#Set up standardized variables just in case to use later
spp_richness.std <- as.data.frame(scale(spp_richness)) #scale to mu=0, sd=1
spp_richness.std <- spp_richness.std %>% dplyr::select(-annwinddir)


# Quick EDA
summary(lm(spp_yr ~ Year + annwindspeed_kph + annwinddir + anndisch_cfs + 
             annsal_ppt + anntemp_c + annwinddir_ew, data = spp_richness))

plotenv <- function(env_var){
  env_var1 <- rlang::sym(env_var)
  
  annotateloc <- min(spp_richness %>% dplyr::select(env_var)) +
          (max(spp_richness %>% dplyr::select(env_var)) - min(spp_richness %>% dplyr::select(env_var)))/5
  pval <- summary(lm(paste0("spp_yr ~ ", env_var), data = spp_richness))$coef[2,4]
  r2 <- summary(lm(paste0("spp_yr ~ ", env_var), data = spp_richness))$adj.r.squared

  p <- ggplot(data = spp_richness, aes(x = !! env_var1, y = spp_yr, color = Year)) +
    geom_point(cex = 3) +
    geom_smooth(method = "lm", se=FALSE) +
    scale_color_viridis() + theme_bw() +
    geom_text(label = spp_richness$Year, nudge_y = 0.25, size = 4) +
    annotate("text", x = annotateloc, y = 23, label = paste0("p-value is: ", round(pval,3))) + 
    annotate("text", x = annotateloc, y = 23.5, label = paste0("Adj R^2 is: ", round(r2,2))) 
  return(p)
}

plotenv("Year")
plotenv("annwindspeed_kph")
plotenv("anndisch_cfs")
plotenv("annsal_ppt")
plotenv("anntemp_c")
plotenv("annwinddir_ew")


# SUMMARY: Individually, only Year is significant, with temp and salinity marginal. 
# So we conclude that species richness is significantly increasing over time
summary(lm(spp_yr ~ Year, data = spp_richness)) 




# Now let's do some model selection to best explain the results that we see 
options(na.action = "na.fail")
fullmodel <- lm(spp_yr ~ Year + annwindspeed_kph + anndisch_cfs + 
        annsal_ppt + anntemp_c + annwinddir_ew, data = spp_richness)
fullmodel1 <- lm(spp_yr ~ Year + annwindspeed_kph + anndisch_cfs + 
                  annsal_ppt + anntemp_c + Year:anntemp_c, data = spp_richness)
dredge(fullmodel1)

par(mar = c(3,5,6,4))
plot(dredge(fullmodel1), labAsExpr = TRUE)
subset(dredge(fullmodel1), delta < 4)

# SUMMARY: The top model is salinity + Year, though just year works well, and
# temp + wind speed is decent. 


#######################
# Does Species Richness increase over the season? 

# This function adds the day of year column (if not already present), and then the weeknumber
addweeknum <- function(dataframename){
  dataframename <- dataframename %>% 
    mutate(day.of.year = yday(EndDate),
           week = ifelse(day.of.year <= 187, 1, 
                         ifelse(day.of.year > 187 & day.of.year <= 194, 2, 
                                ifelse(day.of.year > 194 & day.of.year <= 201, 3, 
                                       ifelse(day.of.year > 201 & day.of.year <= 208, 4, 
                                              ifelse(day.of.year > 208 & day.of.year <= 215, 5,
                                                     ifelse(day.of.year > 215 & day.of.year <= 222, 6,
                                                            ifelse(day.of.year > 222 & day.of.year <= 229, 7,
                                                                   ifelse(day.of.year > 229 & day.of.year <= 236, 8,
                                                                          ifelse(day.of.year > 236, 9, NA))))))))))
  
}


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


temp <- 
ggplot(data = addweeknum(pru.env.day) %>%
         group_by(Year, week) %>% summarise(meansal = mean(Salin_Mid, na.rm = TRUE)), 
       aes(x=week, y = meansal, group=Year, color=Year)) + 
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
