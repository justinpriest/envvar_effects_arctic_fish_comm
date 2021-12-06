summary(lm(meandisch_cfs ~ Year, data = pru.env.biwk))
summary(lm(meandisch_cfs ~ Year, data = pru.env.day))
summary(lm(anndisch_cfs  ~ as.numeric(Year), data = pru.env.ann))


# RESULTS - Effects of Env Var on Sp Comp
pru.env.biwk %>%
  summarise(lotemp = min(Temp_Top, na.rm = TRUE),
            hitemp = max(Temp_Top, na.rm = TRUE),
            meantemp = mean(Temp_Top, na.rm = TRUE),
            sdtemp = sd(Temp_Top, na.rm = TRUE))

summary(lm(Temp_Top ~ Year, data = pru.env.biwk, na.action=na.omit))


pru.env.biwk %>%
  summarise(losal = min(Salin_Top , na.rm = TRUE),
            hisal = max(Salin_Top , na.rm = TRUE),
            meansal = mean(Salin_Top , na.rm = TRUE),
            sdsal = sd(Salin_Top , na.rm = TRUE))
summary(lm(Salin_Top ~ Year, data = pru.env.biwk, na.action=na.omit))


# requires running script 0
pru.env.ann_new <- pru.env.ann %>%
  mutate(disch_m3s = anndisch_cfs / 35.31466) # this converts it to m^3/s
summary(lm(disch_m3s  ~ as.numeric(Year), data = pru.env.ann_new))
summary(lm(anndisch_cfs  ~ as.numeric(Year), data = pru.env.ann_new))

min(pru.env.ann_new$disch_m3s)
max(pru.env.ann_new$disch_m3s)
mean(pru.env.ann_new$disch_m3s)

min(pru.env.biwk$meandisch_cfs) / 35.31
max(pru.env.biwk$meandisch_cfs) / 35.31
mean(pru.env.biwk$meandisch_cfs) / 35.31


summary(lm(annwinddir_ew ~ as.numeric(Year), data = pru.env.ann_new, na.action=na.omit))
summary(lm(annwinddir   ~ as.numeric(Year), data = pru.env.ann_new, na.action=na.omit))
summary(lm(annwindspeed_kph ~ as.numeric(Year), data = pru.env.ann_new, na.action=na.omit))



# RESULTS - Sp Richness and Rare Sp Presence

allcatch %>% filter(Species != "HYCS" & Species != "UNKN" & Station != 231) %>%
  addbiwknum() %>% 
  group_by(Year, Station, biweekly) %>% 
  summarise(num_spp = n_distinct(Species)) %>%
  ungroup() %>% 
  summarise(meanspp = mean(num_spp),
            lospp = min(num_spp),
            hispp = max(num_spp),
            sdspp = sd(num_spp))





#######################
# Creation of figures for Reviewer 1
pru.env.biwk %>%
  mutate(Station = as.factor(Station)) %>%
  ggplot(aes(x = Year, y = Temp_Top)) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_point(aes(shape = Station, fill = Station, color = Station), size = 2) +
  scale_shape_manual(values = c(1,21,21,21)) +
  scale_fill_manual(values = c("black", "lightgray", "gray40", "black")) +
  scale_color_manual(values = c("black", "gray40", "gray40", "black")) +
  theme_bw()

pru.env.biwk %>%
  mutate(Station = as.factor(Station)) %>%
  ggplot(aes(x = Year, y = Salin_Top)) +
  geom_point(aes(shape = Station, fill = Station, color = Station), size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_shape_manual(values = c(1,21,21,21)) +
  scale_fill_manual(values = c("black", "lightgray", "gray40", "black")) +
  scale_color_manual(values = c("black", "gray40", "gray40", "black")) +
  theme_bw()







topspcatch <- catchenviron %>%
  filter(Species %in% c("BDWF", "SFCD", "ARCD", "FHSC", "HBWF", "LSCS")) %>%
  group_by(Year, biweekly, Station, Species) %>%
  summarise(catch = sum(totcount)) 

lm(catch ~ Year, data = topspcatch %>% filter(Species == "BDWF")) %>% summary()


topspcatch %>%
  ggplot(aes(x = Year, y = catch)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~Species, scales = "free_y", ncol = 1)


