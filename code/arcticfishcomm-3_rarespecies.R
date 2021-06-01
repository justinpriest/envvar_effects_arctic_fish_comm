# Logistic binomial regression of rare species

library(here)
library(visreg)
library(scales)
source(here::here("code/arcticfishcomm-1_import&cleanup.R"))


# Data Import and Cleanup




#rarespp.biwk.pres %>% spread(Species, value = pres.abs) # turn wide to view it
rarespp.biwk.pres %>% group_by(Species) %>% summarise(instances = sum(pres.abs)) # summarize total catch instances




##############################
## DO RARE SPECIES INCREASE OVER TIME OR HAVE SEASONALITY


### All Rare Species Combined ###
summary(glm(pres.abs ~ Year + biweekly + Species + Station, family = "binomial", data = rarespp.biwk.pres))

biweeklogit <- glm(pres.abs ~ Year + biweekly, family = "binomial", data = rarespp.biwk.pres)

summary(biweeklogit) #biweekly is VERY significant


# We could restrict it to remove those REALLY rare species
# rareish <- rarespp.biwk.pres %>% group_by(Species) %>% summarise(instances = sum(pres.abs)) %>% filter(instances > 3) %>% 
#   dplyr::select(Species) %>% as.matrix() %>% as.vector()
# uncommonspp <- rarespp.biwk.pres %>% filter(Species %in% rareish)
# summary(glm(pres.abs ~ Year + biweekly, family = "binomial", data = uncommonspp))
# Doesn't really tell us much more than we already knew


plot(pres.abs ~ jitter(Year, amount = 0.4), data=rarespp.biwk.pres)
plot(pres.abs ~ jitter(biweekly, amount=0.4), data=rarespp.biwk.pres)


newdf <- expand.grid(Year=2001:2018, Station=c("220", "218", "214", "230"), biweekly = 1:4)
newdf <- newdf %>% mutate(pred.link = predict(biweeklogit, newdata = newdf, type = "link"),
                          pred.prob = plogis(pred.link))

plot(pres.abs ~ jitter(biweekly, amount=0.4), data=rarespp.biwk.pres)
lines(pred.prob ~ biweekly, data=newdf)

ggplot(data=newdf %>% filter(Station=="214", Year==2006), aes(x=biweekly, y = pred.prob)) + 
  geom_point() + geom_line()

ggplot(data=rarespp.biwk.pres, aes(x=biweekly, y=pres.abs)) + 
  geom_jitter(width = 0.45, height = 0, alpha = 0.4) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  ggtitle("All 15 Rare Species Combined") 



binom.analyzebyspp <- function(sppfilter){
  .raresppfilter <- rarespp.biwk.pres %>% filter(Species == sppfilter)
  .fit <- glm(pres.abs ~ Year + biweekly + Station, family = "binomial", 
      data = .raresppfilter)
  print(summary(.fit))
  # print(ggplot(data=.raresppfilter, aes(x=biweekly, y=pres.abs)) + 
  #         geom_jitter(width = 0.45, height = 0, alpha = 0.4) +
  #         geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  #         ggtitle(paste0(sppfilter, " by biweekly")) )
  # print(ggplot(data=.raresppfilter, aes(x=Year, y=pres.abs)) + 
  #       geom_jitter(width = 0.45, height = 0, alpha = 0.4) +
  #       geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  #       ggtitle(paste0(sppfilter, " by Year")) ) #keep code in case we need to reuse later
  # Compare the confidence bands from visreg. Much wider, and more accurate
  print(visreg(.fit, "Year", scale="response", jitter=T, gg=T, ylab="Probability of Occurrence") + ggtitle(sppfilter) )
  print(visreg(.fit, "biweekly", scale="response", jitter=T, gg=T, ylab="Probability of Occurrence") + ggtitle(sppfilter) )
}



binom.analyzebyspp("ARCH") #insig
binom.analyzebyspp("BRBT") #sig inc by year
binom.analyzebyspp("BRCS") #sig dec by year
binom.analyzebyspp("CHUM") #sig inc biwk 
binom.analyzebyspp("ELPT") #insig
binom.analyzebyspp("KPGL") #insig
binom.analyzebyspp("LIPA") #sig inc biwk
binom.analyzebyspp("LNSK") #insig
binom.analyzebyspp("RKGL") #insig
binom.analyzebyspp("SHSC") #insig
binom.analyzebyspp("SLEB") #insig
binom.analyzebyspp("SLSC") #sig inc by year
binom.analyzebyspp("SOCK") #insig
binom.analyzebyspp("WOLF") #insig
binom.analyzebyspp("WSGL") #insig



#include only species with significant changes in year (but also explore signif biweekly spp)

fit_rare <- glm(pres.abs ~ Year + biweekly + Species + Station - 1, family = binomial, 
           data = rarespp.biwk.pres %>% filter(Species %in% c("BRBT", "BRCS", "CHUM", "LIPA", "SLSC"))) 
summary(fit_rare)

fit_rare_BRBT <- glm(pres.abs ~ Year + biweekly + Station, family = binomial, 
                data = rarespp.biwk.pres %>% filter(Species == "BRBT")) 
summary(fit_rare_BRBT)

fit_rare_BRCS <- glm(pres.abs ~ Year + biweekly + Station, family = binomial, 
                     data = rarespp.biwk.pres %>% filter(Species == "BRCS")) 
summary(fit_rare_BRCS)

fit_rare_SLSC <- glm(pres.abs ~ Year + biweekly + Station, family = binomial, 
                     data = rarespp.biwk.pres %>% filter(Species == "SLSC")) 
summary(fit_rare_SLSC)


ggplot(data=rarespp.biwk.pres %>% filter(Species %in% c("BRBT", "BRCS", "CHUM", "LIPA", "SLSC")), 
       aes(x=biweekly, y=pres.abs, color=Species)) + 
  geom_jitter(width = 0.45, height = 0, shape=124, size =8) +
  ggtitle("new title") 





##########################

#Other exploratory plots




# ggplot() +
#   geom_jitter(data=rarespp.biwk.pres %>% filter(Species %in% c("BRBT", "BRCS", "CHUM", "LIPA", "SLSC")), 
#               aes(x=Year, y=pres.abs), width = 0.45, height = 0, shape=124, size =6) +
#   geom_line(data = predictedvals.year, aes(x=Year, y = predval, color = Species), cex=2)
# 
# 
# ggplot() +
#    # geom_jitter(data=rarespp.biwk.pres %>% filter(Species %in% c("BRBT", "BRCS", "CHUM", "LIPA", "SLSC")), 
#    #             aes(x=biweekly, y=pres.abs), width = 0.45, alpha = 0.8, height = 0, shape=124, size = 4) +
#   geom_line(data = predictedvals.biwk, aes(x=biweekly, y = predval, color = Species, linetype=Species), cex=2) +
#   scale_linetype_manual(values = c("dotdash", "solid","solid", "dotdash","solid")) +
#   scale_color_manual(values = c( "#000000", "#000000", "#bababa", "#4f4f4f", "#4f4f4f")) +
#   #scale_y_continuous(limits = c(0,0.25)) +
#   theme_bw()
# 
# 
# ggplot() +
#   # geom_jitter(data=rarespp.biwk.pres %>% filter(Species %in% c("BRBT", "BRCS", "CHUM", "LIPA", "SLSC")), 
#   #             aes(x=Year, y=pres.abs), width = 0.45, alpha = 0.8, height = 0, shape=124, size = 4) +
#   geom_line(data = predictedvals.year, aes(x=Year, y = predval, color = Species, linetype=Species), cex=2) +
#   scale_linetype_manual(values = c("dotdash", "solid","solid", "dotdash","solid")) +
#   scale_color_manual(values = c( "#000000", "#000000", "#bababa", "#4f4f4f", "#4f4f4f")) +
#   scale_y_continuous(limits = c(0,0.1)) +
#   theme_bw()
