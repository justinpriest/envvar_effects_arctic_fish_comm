# CHANGES IN COMMUNITY STRUCTURE

# Script to run PERMANOVA models and analysis

# This was done at on an annual and a biweekly scale
# There are two sections, first Section 1 looks at 

library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(broom)
library(colorspace)
library(scales)


source("Analysis/thesis2019_Ch1_1-import&cleanup.R")

#this pulls in the following (relevant) dataframes:
head(catchenviron) # all catch data, left joined with environ data
head(pru.env.day)  # environmental data by day
head(pru.env.ann)  # annual summary of environmental data


head(catchmatrix)     # each row is a year/station, cols are species
head(catchmatrix.std) # standardize the above 0-1 (percent of max catch)
head(catchmatrix.day) # each row is a day/station, cols are species
head(catchmatrix.day.std)
head(catchmatrix.biwk.stdtrans)


#################
### Section 1: Test for changes in community assemblage structure
#################

################################
## ORDINATION & DISSIMILARITY ##



### BRAY-CURTIS DISTANCE
braydist <- vegdist(catchmatrix.std, method="bray") # This is a dissimilarity matrix

# library(gplots)
# heatmap.2(as.matrix(braydist))

totalNMDS <- metaMDS(braydist, k=3) # had to choose k=3 to get below 0.20

#Plot, put years by color, add environ vectors
plot(totalNMDS, display = "sites", type = "n", main = "Annual Results - colored by year grp")
text(totalNMDS, select=which(pru.env.ann$Year<2010), col="red")
text(totalNMDS, select=which(pru.env.ann$Year>2009), col="blue")
#unclear if any patterns

plot(totalNMDS, display = "sites", type = "n", main = "Annual Results - colored by station")
text(totalNMDS, select=which(pru.env.ann$Station==214), col="red")
text(totalNMDS, select=which(pru.env.ann$Station==218), col="blue")
text(totalNMDS, select=which(pru.env.ann$Station==220), col="black")
text(totalNMDS, select=which(pru.env.ann$Station==230), col="green")
#pretty distinct groupings by site

# From Franz: See how species are associated with sites
totalNMDS_raw <- metaMDS(catchmatrix.std, k=3)
plot(totalNMDS_raw$points, asp=1, pch=16, col=factor(pru.env.ann$Station))
text(totalNMDS_raw, display = "species", col=2)
# GRAY/RDWF assoc. w/ 230, DLVN w/ 214, PINK/PCHG/SFCD w/ 218 and CAPE/HBWF/LSCS w/ 220


nmdspoints <- as.data.frame(totalNMDS$points[1:71,]) # 71 is the num. year/stn combos (4 stns, 18 yrs (no 231 in 2001))
#row.names(nmdspoints) <- 2001:2017
nmdspoints$YearStn <- rownames(totalNMDS$points)
nmdspoints$Year <- as.numeric(substr(nmdspoints$YearStn, 1, 4))
nmdspoints$Station <- factor(substr(nmdspoints$YearStn, 5, 7))
nmdspoints$earlymidlate <- ifelse(nmdspoints$Year < 2007, "early", 
                                  ifelse(nmdspoints$Year >= 2013, "late", "mid"))
#this arbitrarily divides study in thirds, not the correct approach but good as EDA

ggplot(nmdspoints, aes(x=MDS1, y=MDS2)) + geom_point() + #aes(color=Station), cex=5
  geom_text(aes(label=YearStn, color=Station),hjust=.35, vjust=-.7, size=3)+
  theme_bw() + theme(panel.grid.minor = element_blank()) 

ggplot(nmdspoints, aes(x=MDS1, y=MDS2)) + geom_point() +
  geom_text(aes(label=YearStn, color=Year),hjust=.35, vjust=-.7, size=3)+
  theme_bw() + theme(panel.grid.minor = element_blank()) 
# I used earlymidlate as a color too, but pattern was too slight

ggplot(nmdspoints, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=Station), cex=5) + 
  scale_color_manual(values =  brewer.pal(4, "Set2")) +
  theme_bw() + theme(panel.grid.minor = element_blank()) 



mypal  <- colorRampPalette(brewer.pal(6, "Greens"))
mypal2 <- colorRampPalette(brewer.pal(6, "Greys"))
mypal3 <- colorRampPalette(brewer.pal(6, "Blues"))
mypal4 <- colorRampPalette(brewer.pal(6, "Reds"))

ggplot(nmdspoints, aes(x=MDS1, y=MDS2)) + geom_point() + 
  geom_text(aes(label=YearStn, color=interaction(Year, Station))) + 
  #scale_color_continuous(low = "#3fdeff", high = "#144f5b") +
  scale_colour_manual(values = c(mypal(18), mypal2(18), mypal3(18), mypal4(18))) +
  theme(legend.position = "none")



##########################
## nMDS axis regression ##

plot(MDS1 ~ Year, data = nmdspoints)
summary(lm(MDS1 ~ Year , data = nmdspoints)) # overall not significant
summary(lm(MDS2 ~ Year , data = nmdspoints)) # marginally significant

#now let's break it down by Station
#library(broom)
nmdspoints %>% group_by(Station) %>% do(model = lm(MDS1 ~ Year, data = .)) %>% 
  tidy(model) # significant trend at Station 230 & 214, marginal at 220 
nmdspoints %>% group_by(Station) %>% do(model = lm(MDS2 ~ Year, data = .)) %>% 
  tidy(model) # none significant

# visualize this for MDS1
ggplot(nmdspoints, aes(x=Year, y =MDS1, color = Station)) + 
  geom_point() + geom_smooth(method = "lm", se=FALSE)

# From Franz: "Axis 3 shows a very strong time trend and the trend is very consistent among the 4 sites"
ggplot(nmdspoints, aes(x=Year, y = MDS3, color = Station)) + 
  geom_point() + geom_smooth(se=FALSE)


######### BIWEEKLY
######### lets try with biweekly which is standardized and 4th root transformed
braydist.biwk <- vegdist(catchmatrix.biwk.stdtrans, method="bray")
totalNMDS.biwk <- metaMDS(braydist.biwk, k=3)  #not convergent with k=2

nmdspoints.biwk <- as.data.frame(totalNMDS.biwk$points[1:280,]) # 280 is the number of year/biwk/stn combos
nmdspoints.biwk$YearStn <- rownames(nmdspoints.biwk)
nmdspoints.biwk$Year <- as.numeric(substr(nmdspoints.biwk$YearStn, 1, 4))
nmdspoints.biwk$biweekly <- factor(substr(nmdspoints.biwk$YearStn, 5, 5))
nmdspoints.biwk$Station <- factor(substr(nmdspoints.biwk$YearStn, 6, 8))


ggplot(nmdspoints.biwk, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=Station), cex=5) +
  #geom_text(aes(label=YearStn, color=Station),hjust=.35, vjust=-.7, size=3)+
  scale_color_manual(values =  brewer.pal(4, "Set2")) +
  theme_bw() + theme(panel.grid.minor = element_blank()) 
# Looks like the sites ordinate out separately  



ggplot(nmdspoints.biwk, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=biweekly), cex=5) + 
  scale_color_manual(values = heat_hcl(6, power= c(1, 1.5))) + #library(colorspace)
  theme_bw() + theme(panel.grid.minor = element_blank()) 



# CORRELATIONS
Spp.cor <- data.frame(Species = as.character(""), MDS1.corr = as.numeric(0), 
                      MDS2.corr = as.numeric(0),  MDS3.corr = as.numeric(0),
                      stringsAsFactors = FALSE) # set up blank data frame
j = 1
for(i in colnames(catchmatrix.biwk.stdtrans)){
  .sppcor <- catchmatrix.biwk.stdtrans %>% gather(Species, abund) %>% filter(Species == i)
  # print(i)
  # print(cor(.sppcor$abund, nmdspoints.biwk$MDS1))
  # print(cor(.sppcor$abund, nmdspoints.biwk$MDS2))
  Spp.cor[j,1] <- i
  Spp.cor[j,2] <-cor(.sppcor$abund, nmdspoints.biwk$MDS1)
  Spp.cor[j,3] <-cor(.sppcor$abund, nmdspoints.biwk$MDS2)
  Spp.cor[j,4] <-cor(.sppcor$abund, nmdspoints.biwk$MDS3)
  j <- j+1
}
Spp.cor


ggplot(Spp.cor) + 
  geom_tile(aes(x="MDS1", y=Species, fill = MDS1.corr)) + 
  geom_tile(aes(x="MDS2", y=Species, fill = MDS2.corr)) +
  geom_tile(aes(x="MDS3", y=Species, fill = MDS3.corr)) + 
  scale_fill_gradient2(low = muted("darkred"), #muted requires "scales"
                       mid = "white", 
                       high = muted("cornflowerblue"), 
                       midpoint = 0) +
  labs(x = "", fill="Corr Coef")




#################
### Section 2: Quantify if / how changes are related to environmental variability
#################

## Bioenv ##
# Bioenv is a mantel type test: which combination of environmental var explain it best
bioenv(braydist ~ annwindspeed_kph + annwinddir + anndisch_cfs + 
         annsal_ppt + anntemp_c, pru.env.ann)
# salinity and temp are the best predictors

bioenv(braydist.biwk ~ biwkmeanspeed_kph + biwkmeandir + meandisch_cfs + 
         Salin_Top + Temp_Top + winddir_ew, pru.env.biwk)
# salinity, temp, and wind direction are best subset of env variables


## EnvFit ##

env.vectors.ann <- envfit(totalNMDS, pru.env.ann, permutations = 999)
env.vectors.ann # Year and Station centroids are significant
# also signif are temp, salin, and marginally wind direction

env.vectors.biwk <- envfit(totalNMDS.biwk, pru.env.biwk %>% 
                             dplyr::select(-winddir_ew, -Salin_Mid, -Temp_Mid), 
                           na.rm = TRUE, permutations = 999)
env.vectors.biwk # salin, wind, and year are signif





plot(totalNMDS) 
plot(env.vectors.ann, p.max = 0.02) # >pvals make plot too busy 
#way too busy to repeat for biweekly scale

ggplot(nmdspoints, aes(x=MDS1, y =MDS2)) + geom_point() +
  scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(-1, 1)) +
  geom_segment(data = data.frame(env.vectors.ann$vectors$arrows), aes(x=0, xend=NMDS1, y=0, yend=NMDS2))
#these segments don't match the arrows from before because of automatic scaling. I think. not sure


ggplot(nmdspoints, aes(x=MDS1, y =MDS2)) + geom_point(aes(color = Station)) +
  scale_x_continuous(limits = c(-0.43, 0.4)) + scale_y_continuous(limits = c(-0.4, 0.4)) +
  geom_segment(data = data.frame(env.vectors.ann$vectors$arrows) %>% 
                 cbind(r2=env.vectors.ann$vectors$r, pval =env.vectors.ann$vectors$pvals), 
               aes(x=0, xend=NMDS1 * r2, y=0, yend=NMDS2*r2))



ggplot(nmdspoints.biwk, aes(x=MDS1, y =MDS2)) + geom_point(aes(color = Station), cex =5) +
  scale_x_continuous(limits = c(-0.23, 0.25)) + scale_y_continuous(limits = c(-0.18, 0.25)) +
  geom_segment(data = data.frame(env.vectors.biwk$vectors$arrows) %>% 
                 cbind(r2=env.vectors.biwk$vectors$r, pval =env.vectors.biwk$vectors$pvals), 
               aes(x=0, xend=NMDS1 * (r2^0.4)/3, y=0, yend=NMDS2 * (r2^0.4)/3), cex =2) +#need to fix scale (mult by 5?)
  scale_color_manual(values =  brewer.pal(4, "Set2")) +
  theme_bw() + theme(panel.grid.minor = element_blank()) 


###############
## ANNUAL PERMANOVA ##

# following example from vegan tutorial, page 33
# http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf

betad <- betadiver(catchmatrix.std, "z")   # using Arrhenius z measure of beta diversity
boxplot(betadisper(betad, pru.env.ann$Year), main = "Annual")

adonis(betad ~ Year, pru.env.ann, perm=999) 
adonis(betad ~ Station, pru.env.ann, perm=999) 
adonis(betad ~ anntemp_c, pru.env.ann, perm=999) 
adonis(betad ~ annsal_ppt, pru.env.ann, perm=999) 
adonis(betad ~ anndisch_cfs, pru.env.ann, perm=999)     # not significant
adonis(betad ~ annwinddir_ew, pru.env.ann, perm=999)    # not significant
adonis(betad ~ annwindspeed_kph, pru.env.ann, perm=999) # not significant

adonis(betad ~ annsal_ppt + anntemp_c + Station + Year, pru.env.ann, perm=999) 
#note that terms are sequential so order matters!
adonis(betad ~ Station + annsal_ppt + anntemp_c, pru.env.ann, perm=999) 

boxplot(betadisper(betad, pru.env.ann$Station), main = "Annual")
#note that the western and eastern sites look similar


adonis(catchmatrix.std ~ annsal_ppt + annwinddir_ew + annwindspeed_kph + 
         anndisch_cfs + anntemp_c, data=pru.env.ann, perm=999)


adonis(braydist ~ Year + Station, data = pru.env.ann, perm = 9999)




######
## BIWEEKLY PERMANOVA ##

# because there are rows with NAs in the env data, we need to remove these
catchmatrix.biwk.stdtrans.sub <- catchmatrix.biwk.stdtrans[!rowSums(is.na(pru.env.biwk)) >0,]
pru.env.biwk.sub <- pru.env.biwk[!rowSums(is.na(pru.env.biwk)) >0,]

catchmatrix.biwk.cpue.stdtrans.sub <- catchmatrix.biwk.cpue.stdtrans[!rowSums(is.na(pru.env.biwk)) >0,]
pru.env.biwk.sub <- pru.env.biwk[!rowSums(is.na(pru.env.biwk)) >0,]


pru.env.biwk.std <- pru.env.biwk.sub
for (i in 4:ncol(pru.env.biwk.std)){ #starts at 4 to exclude Year, biweekly, and station cols
  pru.env.biwk.std[i] <- ((pru.env.biwk.sub[i]+1)^0.5)/max(((pru.env.biwk.sub[i]+1)^0.5))}
#using square root tranform



betad.biwk <- betadiver(catchmatrix.biwk.stdtrans.sub , "z") 
adonis(betad.biwk ~ Temp_Top + Salin_Top + winddir_ew + meandisch_cfs + Year + Station + biweekly, pru.env.biwk.std, perm=999) 
#winddir slightly better than speed. Temp signif if added first, but mostly captured by salin
# nonenviron explan var (Year, Stn, biweekly) are highly significant, esp seasonality (biweekly)
# But Franz recommends using the Bray-Curtis dissimilarity matrix
adonis2(catchmatrix.biwk.stdtrans.sub ~ Temp_Top + Salin_Top + meandisch_cfs + winddir_ew + Year + Station + biweekly, 
        pru.env.biwk.std, perm=999, by = "margin") 
# Seasonality, year, and station effects are the main contributors

#check whether CPUE changes anything
adonis2(catchmatrix.biwk.cpue.stdtrans.sub ~ Temp_Top + Salin_Top + meandisch_cfs + winddir_ew + Year + Station + biweekly, 
        pru.env.biwk.std, perm=999, by = "margin")
# same general trends, but the model fits better (resids dropped from 0.54 to 0.51), Stn R2 up


boxplot(betadisper(betad.biwk, pru.env.biwk.sub$Station), main = "Biweekly")
boxplot(betadisper(betad.biwk, pru.env.biwk.sub$Year), main = "Biweekly")
boxplot(betadisper(betad.biwk, pru.env.biwk.sub$biweekly), main = "Biweekly")



# Simper
summary(simper(catchmatrix, pru.env.ann$Station))
# ARCD, ARCS, BDWF, and LSCS account for most of the differences
summary(simper(catchmatrix.std, pru.env.ann$Station))
# HBWF, ARFL, DLVN, RDWF explain most of the standardized diffs

summary(simper(catchmatrix.biwk.stdtrans, pru.env.biwk.std$Year))
summary(simper(catchmatrix.biwk.stdtrans, pru.env.biwk.std$Station))
# Seems like THSB, RDWF, PINK, PCHG are most common? Hard to tell






