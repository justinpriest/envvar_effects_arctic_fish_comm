## I modified and added to your code starting at line 103
#  in: 'thesis2019_Ch1_2_-speciesrichness.R'

# This was inspired by your visualization of the salinity
# and species richness trends, but appropriately 
# accounts for effort. I used effort as an explanatory
# variable and it has the effect that you would expect
# (richness first increases with effort, then levels off)

# The model is perhaps a bit more complicated than what you did
# but includes smooth trends over time, with salinity and
# with sampling effort, as well as a linear trend over time
# and, to account for any 'unexplained' interannual variability,
# I also included a random year effect (using a smooth term
# with 'bs = "re"', which is not very well documented but 
# allows you to include a random effect in a 'gam' model)

# I started with the weekly data to be able to estimate 
# effects of effort more clearly:

spp_richness.wk <- spp_richness.wk %>% left_join(effort.wk %>% dplyr::select(Year, week, sample_proportion), 
                                                 by = c("Year" = "Year", "week" = "week")) 

weeklysalinity <- addweeknum(pru.env.day) %>%
  group_by(Year, week) %>% summarise(meansal = mean(Salin_Mid, na.rm = TRUE))

weekly <- left_join(spp_richness.wk, weeklysalinity)

weekly$Year.fac <- factor(weekly$Year)
weekly <- rename(weekly, effort = sample_proportion)

# GAM models with smooth seasonal term, smooth salinity term,
# random intercept for Year, and smooth trend with effort
wk.gam1 <- gam(num_spp ~ s(week, k=4) + s(meansal, k=4) 
               + s(Year.fac, bs="re") + s(effort, k=4)
               + Year, data=na.omit(weekly), 
               select=T, method="ML")
summary(wk.gam1)
par(mfrow=c(3,2), mar=c(4,4,1,1))
visreg(wk.gam1)

# It is somewhat difficult to untangle the effects of
# effort and salinity on species richness as they 
# are (moderately) correlated:
cor(weekly[,c("effort", "meansal")], use="pair")
# Both salinity (your animation) and effort are low
# in weeks 1 and 9, as is species richness:
plot(effort ~ week, data=weekly)
plot(num_spp ~ week, data=weekly)

# Can seasonal trend be 'explained' by salinity?
wk.gam2 <- gam(num_spp ~ s(meansal, k=4) 
               + s(Year.fac, bs="re") + s(effort, k=4)
               + Year, data=na.omit(weekly), 
               select=T, method="ML")
summary(wk.gam2)
par(mfrow=c(2,2))
visreg(wk.gam2)
# Salinity can account for a good part of the
# variability in species richness, but the R^2
# is considerably lower compared to the model with 
# a simple seasonal trend


# Moreover, salinity was essentially eliminated from
# the first model (edf=0), which the gam function can 
# do if 'select=T'. Hence let's remove salinity:
wk.gam3 <- gam(num_spp ~ s(week, k=4) + s(effort, k=4)
               + s(Year.fac, bs="re") + Year,
               data=na.omit(weekly), 
               select=T, method="ML")
summary(wk.gam3)
visreg(wk.gam3)
gam.vcomp(wk.gam3)  # Variances associated with terms
# The standard deviations suggest that the year-to-year
# variability is about half (0.69) of the residual 
# variability ('scale' = 1.22).

# The AIC also suggest that salinity does not result
# in an improvement and the model without the seasonal
# trend results in a much poorer fit:
AIC(wk.gam1, wk.gam3, wk.gam2)

# Because salinity is eliminated, we can fit the full data set:
# (Fit using REML, which gives better estimates but is not
# suitable for model comparisons)
wk.gam.all <- gam(num_spp ~ s(week, k=4) + s(effort, k=4)
               + s(Year.fac, bs="re") + Year,
               data=weekly, 
               select=T, method="REML")
summary(wk.gam.all)
visreg(wk.gam.all)
# The results are virtually identical, but the linear term
# is now significant at a 95% level (p = 0.0346)

# In summary, based on the weekly data, the results in your
# Fig. 1 are confirmed, which is great, and your rationale for
# using biweekly data is probably good. 

# However, rather than the adjusted species richness, you should
# include effort as an explanatory variable in your model,
# as well as a random effect for year (if you didn't?).


##########################################################
# Using biweekly data, continuing with your line 187

biweekly <- spp_richness.biwk %>% left_join(effort.biwk %>% dplyr::select(Year, biweekly, sample_proportion), 
                                                     by = c("Year" = "Year", "biweekly" = "biweekly"))
biweekly$Year.fac <- factor(biweekly$Year)
biweekly <- rename(biweekly, effort = sample_proportion)

# GAM models with smooth seasonal term, random intercept
# for Year, and effort
biwk.gam1 <- gam(num_spp ~ s(biweekly, k=4) + s(effort, k=4)
               + s(Year.fac, bs="re") 
               + Year, data=biweekly, 
               select=T, method="ML")
summary(biwk.gam1)
par(mfrow=c(2,2), mar=c(4,4,1,1))
visreg(biwk.gam1)

# For the biweekly model, the effort term is not significant,
# but richness does increase some with effort, as expected.
# This may suggest that - over two weeks of sampling - you 
# generally pick up all the species that are present, so you
# could safely remove effort from model:
biwk.gam2 <- gam(num_spp ~ s(biweekly, k=4)
                 + s(Year.fac, bs="re") 
                 + Year, data=biweekly, 
                 select=T, method="ML")
summary(biwk.gam2)
par(mfrow=c(2,2), mar=c(4,4,1,1))
visreg(biwk.gam2)

# One thing to note is that the uncertainty in the linear time
# trend is a lot larger than in your Fig. 1 and that is likely
# due to including the random effect:
biwk.gam3 <- gam(num_spp ~ s(biweekly, k=4)
                 + Year, data=biweekly, 
                 select=T, method="ML")
summary(biwk.gam3)
par(mfrow=c(2,2), mar=c(4,4,1,1))
visreg(biwk.gam3)
# Much tighter CI for both seasonal and long-term trend,
# but the model has strong residual patterns with some
# years having all negative and others all positive residuals
# --> that's why we need the random year effect!

AIC(biwk.gam1, biwk.gam2, biwk.gam3)
# The second model also the AIC-preferred model


