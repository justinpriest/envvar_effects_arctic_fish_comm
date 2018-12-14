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






