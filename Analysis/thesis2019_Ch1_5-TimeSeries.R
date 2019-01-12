# COMMUNITY CHANGES IN TIME SERIES 

#################
### Assess nature of observed changes: linear trend, nonlinear trend, structural break

#################

library(mgcv)
library(strucchange)
library(ggsidekick) #devtools::install_github("seananderson/ggsidekick")

source("Analysis/thesis2019_Ch1_4-PERMANOVA.R")


summary(lm(MDS1 ~ Year, data = nmdspoints.biwk)) # marginal
summary(lm(MDS2 ~ Year, data = nmdspoints.biwk)) # very significant
summary(lm(MDS3 ~ Year, data = nmdspoints.biwk)) # not significant


#now let's break it down by Station
nmdspoints.biwk %>% group_by(Station) %>% do(model = lm(MDS1 ~ Year, data = .)) %>% 
  tidy(model) # marginal at 218, not signif at any others
nmdspoints.biwk %>% group_by(Station) %>% do(model = lm(MDS2 ~ Year, data = .)) %>% 
  tidy(model) # significant at 220 and 214, marginal at 218

# visualize this for MDS1 & MDS2
ggplot(nmdspoints.biwk, aes(x=Year, y =MDS1, color = Station)) + 
  geom_point() + geom_smooth(method = "lm", se=FALSE)
ggplot(nmdspoints.biwk, aes(x=Year, y =MDS2, color = Station)) + 
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
summary(gam(MDS3 ~ s(Year) + Station + biweekly, data = nmdspoints.biwk))
# summary: MDS1&2 have better fit with nonlinear, no diff MDS3 (measured by deviance)



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

strucsummary(nmdspoints.biwk, station=220, mds = "MDS2") # one optimal breakpoint
strucsummary(nmdspoints.biwk, station=214, mds = "MDS2") # one optimal breakpoint
strucsummary(nmdspoints.biwk, station=218, mds = "MDS2") # one optimal breakpoint
strucsummary(nmdspoints.biwk, station=230, mds = "MDS2") # zero optimal breakpoints


strucsummary(nmdspoints.biwk, station=220, mds = "MDS3") # 
strucsummary(nmdspoints.biwk, station=214, mds = "MDS3") # 
strucsummary(nmdspoints.biwk, station=218, mds = "MDS3") # 
strucsummary(nmdspoints.biwk, station=230, mds = "MDS3") # 

