# COMMUNITY CHANGES IN TIME SERIES 

#################
### Assess nature of observed changes: linear trend, nonlinear trend, structural break

#################

library(mgcv)
library(strucchange)
library(ggsidekick) #devtools::install_github("seananderson/ggsidekick")
library(here)

source(here::here("code/arcticfishcomm-1_import&cleanup.R"))
source(here::here("code/arcticfishcomm-4_PERMANOVA.R"))


summary(lm(MDS1 ~ Year, data = nmdspoints.biwk)) # marginal
summary(lm(MDS2 ~ Year, data = nmdspoints.biwk)) # very significant
summary(lm(MDS3 ~ Year, data = nmdspoints.biwk)) # not significant


#now let's break it down by Station
nmdspoints.biwk %>% group_by(Station) %>% do(model = lm(MDS1 ~ Year, data = .)) %>% 
  tidy(model) # marginal at 218, not signif at any others
nmdspoints.biwk %>% group_by(Station) %>% do(model = lm(MDS2 ~ Year, data = .)) %>% 
  tidy(model) # significant at 220 and 214, marginal at 218
nmdspoints.biwk %>% group_by(Station) %>% do(model = lm(MDS3 ~ Year, data = .)) %>% 
  tidy(model) # significant at 220 and 214, marginal at 218


# visualize this for MDS1 & MDS2
ggplot(nmdspoints.biwk, aes(x=Year, y =MDS1, color = Station)) + 
  geom_point() + geom_smooth(method = "lm", se=FALSE)
ggplot(nmdspoints.biwk, aes(x=Year, y =MDS2, color = Station)) + 
  geom_point() + geom_smooth(method = "lm", se=FALSE)
ggplot(nmdspoints.biwk, aes(x=Year, y =MDS3, color = Station)) + 
  geom_point() + geom_smooth(method = "lm", se=FALSE)


#let's clean up the clutter and see if there is a non-linear trend
ggplot(nmdspoints.biwk %>% group_by(Year, Station) %>% summarise(MDS1 = mean(MDS1)), 
       aes(x=Year, y =MDS1, color = Station)) + 
  geom_line(cex=2.5) + geom_smooth(method = "lm", se=FALSE, cex=1.25) +
  scale_color_manual(values =  brewer.pal(4, "Set2")) +
  theme_bw() + theme(panel.grid.minor = element_blank()) 

ggplot(nmdspoints.biwk %>% group_by(Year, Station) %>% summarise(MDS2 = mean(MDS2)), 
       aes(x=Year, y =MDS2, color = Station)) + 
  geom_line(cex=2.5) + geom_smooth(method = "lm", se=FALSE, cex=1.25) +
  scale_color_manual(values =  brewer.pal(4, "Set2")) +
  theme_bw() + theme(panel.grid.minor = element_blank()) 
# looks mostly linear but we'll test that soon

ggplot(nmdspoints.biwk %>% group_by(Year, Station) %>% summarise(MDS3 = mean(MDS3)), 
       aes(x=Year, y =MDS3, color = Station)) + 
  geom_line(cex=2.5) + geom_smooth(method = "lm", se=FALSE, cex=1.25) +
  scale_color_manual(values =  brewer.pal(4, "Set2")) +
  theme_bw() + theme(panel.grid.minor = element_blank()) 



#Let's fit a nested effects linear model to account for Station effects by Year
model.lm1 <- lm(MDS1 ~ Station / Year, data = nmdspoints.biwk) # nested effect linear model
summary(model.lm1) # one marginal
model.lm2 <- lm(MDS2 ~ Station / Year, data = nmdspoints.biwk)
summary(model.lm2) # several signif
model.lm3 <- lm(MDS3 ~ Station / Year, data = nmdspoints.biwk)
summary(model.lm3) # nothing signif

plot(model.lm1) #diagnostics look decent
plot(model.lm2)
plot(model.lm3)
#summary: it's a better fit when we account for station differences

(nmdspoints.biwk %>% group_by(Station) %>% do(model = lm(MDS1 ~ Year, data = .)))$model


#Now we test for non-linear trends
summary(gam(MDS1 ~ Year + Station + biweekly, data = nmdspoints.biwk))  # library(mgcv)
summary(gam(MDS1 ~ s(Year) + Station + biweekly, data = nmdspoints.biwk))
summary(gam(MDS2 ~ Year + Station + biweekly, data = nmdspoints.biwk))
summary(gam(MDS2 ~ s(Year) + Station + biweekly, data = nmdspoints.biwk))
summary(gam(MDS3 ~ Year + Station + biweekly, data = nmdspoints.biwk))
summary(gam(MDS3 ~ s(Year) + Station + biweekly, data = nmdspoints.biwk)) #left these as gam() to compare deviance better
# summary: MDS1&2 have better fit with nonlinear, no diff MDS3 (measured by deviance)


#Plot predicted GAM output
gampredict.df <- expand.grid(Year=2001:2018, Station=c("214", "218", "220", "230"), biweekly=1:4)

MDS1.gam <- gam(MDS1 ~ s(Year) + Station + as.numeric(biweekly), data = nmdspoints.biwk)
gampredict.df$pred.MDS1 <- predict.gam(MDS1.gam, gampredict.df)
ggplot(gampredict.df, aes(x=Year, y=pred.MDS1, color=Station)) + 
  geom_point() + geom_smooth(se=FALSE, cex=2.5) + 
  scale_color_manual(values =  brewer.pal(4, "Set2")) 
summary(MDS1.gam)

MDS2.gam <- gam(MDS2 ~ s(Year) + Station + as.numeric(biweekly) -1, data = nmdspoints.biwk)
gampredict.df$pred.MDS2 <- predict.gam(MDS2.gam, gampredict.df)
ggplot(gampredict.df, aes(x=Year, y=pred.MDS2, color=Station)) +
  geom_point() + geom_smooth(se=FALSE, cex=2.5) + 
  scale_color_manual(values =  brewer.pal(4, "Set2")) 
summary(MDS2.gam)

MDS3.lm <- lm(MDS3 ~ Year + Station + as.numeric(biweekly) - 1, data = nmdspoints.biwk) # we can let this vary by stn
gampredict.df$pred.MDS3 <- predict(MDS3.lm, gampredict.df)
ggplot(gampredict.df, aes(x=Year, y=pred.MDS3, color=Station)) +
  geom_point() + geom_smooth(se=FALSE, cex=2.5) + 
  scale_color_manual(values =  brewer.pal(4, "Set2")) # Same output as earlier
summary(MDS3.lm)



#Finally, test for structural breaks
strucsummary <- function(nmdsdataframe=nmdspoints.biwk, station, mds = "MDS1"){
  #This function takes the nmds dataframe, filters it for a specific station,
  # turns it into a time series by year, then returns some summaries & plots
  # using the library "strucchange" 
  require(strucchange)
  require(dplyr)
  .datdf <- nmdsdataframe
  .timeseries <- .datdf %>% filter(Station == station) %>% dplyr::select(mds) %>% as.ts(mds)
  .fs.timeser <- Fstats(.timeseries ~1)
  print(plot(.fs.timeser))
  print(lines(breakpoints(.fs.timeser)))
  print(breakpoints(.fs.timeser))
  print("====================")
  print(summary(breakpoints(.timeseries ~ 1)))
}

strucsummary(nmdspoints.biwk, station=220, mds = "MDS1") 
strucsummary(nmdspoints.biwk, station=214, mds = "MDS1")
strucsummary(nmdspoints.biwk, station=218, mds = "MDS1")
strucsummary(nmdspoints.biwk, station=230, mds = "MDS1")
#Summary: For MDS1, all stations have 0 optimal TS breakpoints, though most have 1 breakpoint close behind

strucsummary(nmdspoints.biwk, station=220, mds = "MDS2") # one optimal breakpoint at 2016
strucsummary(nmdspoints.biwk, station=214, mds = "MDS2") # one optimal breakpoint at 2013
strucsummary(nmdspoints.biwk, station=218, mds = "MDS2") # one optimal breakpoint at 2016
strucsummary(nmdspoints.biwk, station=230, mds = "MDS2") # zero optimal breakpoints


strucsummary(nmdspoints.biwk, station=220, mds = "MDS3") # 
strucsummary(nmdspoints.biwk, station=214, mds = "MDS3") # 
strucsummary(nmdspoints.biwk, station=218, mds = "MDS3") # 
strucsummary(nmdspoints.biwk, station=230, mds = "MDS3") # 


############

# Focus on MDS3 if that's the axis that has the strongest time trend

temp5 <- nmdspoints.biwk %>% filter(Station == 220) %>% dplyr::select(MDS2) %>% as.ts(mds)
temp6 <- Fstats(temp5 ~1)
print(plot(temp6))
print(lines(breakpoints(temp6)))
summary(breakpoints(temp6))



# Fit all models (except structural change) within GAM framework:
# Now re-doing the same as above but for MDS1 & MDS3

### MDS1 ###
# Full model, smooth temporal trend by station:
mds1.gam1 <- gam(MDS1 ~ s(Year, by = Station) + Station + biweekly, data = nmdspoints.biwk)
summary(mds1.gam1)
visreg(mds1.gam1, "Year", by="Station")
visreg(mds1.gam1, "Station")
visreg(mds1.gam1, "biweekly")
# Trends differ but all decreasing

# Additive effects (same temporal trend across stations)
mds1.gam2 <- gam(MDS1 ~ s(Year) + Station + biweekly, data = nmdspoints.biwk)
summary(mds1.gam2)
visreg(mds1.gam2)

# No time trend:
mds1.gam3 <- gam(MDS1 ~ Station + biweekly, data = nmdspoints.biwk)
summary(mds1.gam3)
visreg(mds1.gam3)

# Compare the three models + null:
AIC(mds1.gam1, mds1.gam2, mds1.gam3, gam(MDS1 ~ 1, data = nmdspoints.biwk))
# second model also fits best


### MDS2 ###
# Full model, smooth temporal trend by station:
mds2.gam1 <- gam(MDS2 ~ s(Year, by = Station) + Station + biweekly, data = nmdspoints.biwk)
summary(mds2.gam1)
visreg(mds2.gam1, "Year", by="Station", gg = TRUE)
visreg(mds2.gam1, "Station")
visreg(mds2.gam1, "biweekly")
# Trends differ but all decreasing

# Additive effects (same temporal trend across stations)
mds2.gam2 <- gam(MDS2 ~ s(Year) + Station + biweekly, data = nmdspoints.biwk)
summary(mds2.gam2)
visreg(mds2.gam2)

# No time trend:
mds2.gam3 <- gam(MDS2 ~ Station + biweekly, data = nmdspoints.biwk)
summary(mds2.gam3)
visreg(mds2.gam3)

# Compare the three models + null:
AIC(mds2.gam1, mds2.gam2, mds2.gam3, gam(MDS2 ~ 1, data = nmdspoints.biwk))
# The second model fits best, suggesting a single, highly nonlinear trend over time.



### MDS3 ###
# Full model, smooth temporal trend by station:
mds3.gam1 <- gam(MDS3 ~ s(Year, by = Station) + Station + biweekly, data = nmdspoints.biwk)
summary(mds3.gam1)
visreg(mds3.gam1, "Year", by="Station")
visreg(mds3.gam1, "Station")
visreg(mds3.gam1, "biweekly")
# Trends differ but all decreasing

# Additive effects (same temporal trend across stations)
mds3.gam2 <- gam(MDS3 ~ s(Year) + Station + biweekly, data = nmdspoints.biwk)
summary(mds3.gam2)
visreg(mds3.gam2)

# No time trend:
mds3.gam3 <- gam(MDS3 ~ Station + biweekly, data = nmdspoints.biwk)
summary(mds3.gam3)
visreg(mds3.gam3)

# Compare the three models + null:
AIC(mds3.gam1, mds3.gam2, mds3.gam3, gam(MDS3 ~ 1, data = nmdspoints.biwk))
# all 3 models about the same. No time trend slightly best
