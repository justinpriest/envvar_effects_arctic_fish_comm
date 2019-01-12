# Logistic binomial regression

source("Analysis/thesis2019_Ch1_1-import&cleanup.R")


# presence / absence


catchmatrix.biwk.pres <- catchenviron %>% group_by(Year, biweekly, Station, Species) %>% summarise(count = sum(totcount)) %>%
  spread(Species, value = count) %>% replace(., is.na(.), 0) %>% ungroup() %>% 
  arrange(Year, biweekly, Station) %>% # Just to double check order is correct
  gather(Species, pres.abs, -Year, -biweekly, -Station) %>% filter(!Species %in% keepspp) 
#mutate(Year = as.factor(Year), biweekly = as.factor(biweekly))

catchmatrix.biwk.pres$pres.abs[catchmatrix.biwk.pres$pres.abs > 0] <- 1 

catchmatrix.biwk.pres <- catchmatrix.biwk.pres %>% mutate(Year = as.double(Year) + 2000)


# analysis, all rare species combined
biweeklogit <- glm(pres.abs ~ Year + biweekly, family = "binomial", data = catchmatrix.biwk.pres)

summary(biweeklogit) #biweekly is VERY significant

plot(pres.abs ~ jitter(Year, amount = 0.4 ), data=catchmatrix.biwk.pres)
plot(pres.abs ~ jitter(biweekly, amount=0.4), data=catchmatrix.biwk.pres)


newdf <- expand.grid(Year=2001:2018, Station=c("220", "218", "214", "230"), biweekly = 1:4)

newdf <- newdf %>% mutate(pred.link = predict(biweeklogit, newdata = newdf, type = "link"),
                          pred.prob = plogis(pred.link))

plot(pres.abs ~ jitter(biweekly, amount=0.4), data=catchmatrix.biwk.pres)
lines(pred.prob ~ biweekly, data=newdf)

ggplot(data=newdf %>% filter(Station=="214", Year==2006), aes(x=biweekly, y = pred.prob)) + geom_line()

