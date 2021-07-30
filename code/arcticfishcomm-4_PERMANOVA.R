# CHANGES IN COMMUNITY STRUCTURE

# Script to run PERMANOVA models and analysis

# Originally, this was done on both an annual and a biweekly scale
# Code here has been restricted to biweekly only results
# As such, hundreds of lines of code that looked at annual (even daily) results were cut.

# There are two sections: 
# Section 1 tests if species composition has changed
# Section 2 quantifies whether (and how) changes are related to environmental variability

library(vegan)
library(RColorBrewer)
library(broom)
library(colorspace)
library(scales)
library(here)

source(here::here("code/arcticfishcomm-1_import&cleanup.R"))

#this pulls in the following relevant dataframes:
head(pru.env.biwk) # biweekly summary of environmental data
head(catchmatrix.biwk.stdtrans) 
#each row is a biwk/station, cols are species
# This is then standardized to 0-1 (percent of max catch)


set.seed(7787) # Need to reproduce nMDS plots. Randomly chose this number
#################
### Section 1: Test for changes in community assemblage structure
#################

################################
## ORDINATION & DISSIMILARITY ##


# biweekly data is standardized and 4th root transformed

### BRAY-CURTIS DISTANCE
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

topcorrspp <- Spp.cor %>% 
  mutate(totcorr = (abs(MDS1.corr) + abs(MDS2.corr) + abs(MDS3.corr)),
         MDS12corr = (abs(MDS1.corr) + abs(MDS2.corr))) %>%
  top_n(9, MDS12corr) %>% pull(Species)

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
bioenv(braydist.biwk ~ as.numeric(Year) + as.numeric(biweekly) + biwkmeanspeed_kph + 
         biwkmeandir + meandisch_cfs + Salin_Top + Temp_Top + winddir_ew, 
       pru.env.biwk, metric = "euclidean") # old
# biweekly, salinity, temp are best subset of env variables
# If w/o biweekly, wind dir E-W is also important
bioenv(braydist.biwk ~ as.numeric(Year) + as.numeric(biweekly) + meandisch_cfs + 
         Salin_Top + Temp_Top + wind_vector, 
       pru.env.biwk, metric = "euclidean") # updated to include windvector
# does not change results

mantel(braydist.biwk, vegdist(pru.env.biwk %>% ungroup() %>%
         dplyr::select(biweekly, Salin_Top, Temp_Top), permutations = 999, na.rm = TRUE), method = "spearman")
# JTP June 2021: just had to add na.rm argument, double check this doesn't change results


## EnvFit ##
env.vectors.biwk <- envfit(totalNMDS.biwk, pru.env.biwk %>% 
                             mutate(Station = factor(Station)) %>%
                           #dplyr::select(-winddir_ew, -Salin_Mid, -Temp_Mid), 
                           dplyr::select(Year, biweekly, Station, Salin_Top, Temp_Top),   
                           na.rm = TRUE, permutations = 999)
env.vectors.biwk # salin, wind, and year are signif. temp marginal
# as a test, I included wind_vector (even though it wasn't one of the top params
#  from the Mantel test bioenv). It was signif p<0.001, R2 ~ 0.047




###############
## BIWEEKLY PERMANOVA ##
# Following example from vegan tutorial, page 33
# http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf


# because there are rows with NAs in the env data, we need to remove these
catchmatrix.biwk.stdtrans.sub <- catchmatrix.biwk.stdtrans[!rowSums(is.na(pru.env.biwk)) >0,]
pru.env.biwk.sub <- pru.env.biwk[!rowSums(is.na(pru.env.biwk)) >0,]

catchmatrix.biwk.cpue.stdtrans.sub <- catchmatrix.biwk.cpue.stdtrans[!rowSums(is.na(pru.env.biwk)) >0,]


# Original code. Not using once added wind_vector because negative values cause trans issues
# As per conversation with FJM (Sept 8), PERMANOVA does not require distributions
# pru.env.biwk.std <- pru.env.biwk.sub
# for (i in 4:ncol(pru.env.biwk.std)){ #starts at 4 to exclude Year, biweekly, and station cols
#   pru.env.biwk.std[i] <- ((pru.env.biwk.sub[i]+1)^0.5)/max(((pru.env.biwk.sub[i]+1)^0.5))}
#using square root tranform
# This standardizes each variable to itself


# Standardize variables to max obs
pru.env.biwk.std <- pru.env.biwk.sub
for (i in 4:ncol(pru.env.biwk.std)){ #starts at 4 to exclude Year, biweekly, and station cols
  pru.env.biwk.std[i] <- pru.env.biwk.sub[i]/max(((pru.env.biwk.sub[i])))}
# This standardizes each variable to itself



betad.biwk <- betadiver(catchmatrix.biwk.stdtrans.sub , "z") # using Arrhenius z measure of beta diversity
#adonis(betad.biwk ~ Temp_Top + Salin_Top + winddir_ew + meandisch_cfs + Year + Station + biweekly, pru.env.biwk.std, perm=999)
adonis(betad.biwk ~ Temp_Top + Salin_Top + wind_vector + meandisch_cfs + Year + Station + biweekly, pru.env.biwk.std, perm=999)
#winddir slightly better than speed. Temp signif if added first, but mostly captured by salin
# nonenviron explan var (Year, Stn, biweekly) are highly significant, esp seasonality (biweekly)
# But Franz recommends using the Bray-Curtis dissimilarity matrix
# adonis2(catchmatrix.biwk.stdtrans.sub ~ Temp_Top + Salin_Top + meandisch_cfs + winddir_ew + Year + Station + biweekly, 
#         pru.env.biwk.std, perm=999, by = "margin") # Turning this off for now for import speed
# Seasonality, year, and station effects are the main contributors

#check whether CPUE changes anything
adonis2(catchmatrix.biwk.cpue.stdtrans.sub ~ Temp_Top + Salin_Top + meandisch_cfs + winddir_ew + Year + Station + biweekly, 
        pru.env.biwk.std, perm=999, by = "margin")
# same general trends, but the model fits better (resids dropped from 0.54 to 0.51), Stn R2 up
# interaction btwn seasonal & station
adonis2(catchmatrix.biwk.cpue.stdtrans.sub ~ Temp_Top + Salin_Top + meandisch_cfs + wind_vector + Year + Station + biweekly, 
        pru.env.biwk.std, perm=999, by = "margin")


permanova.biwk <- adonis2(catchmatrix.biwk.cpue.stdtrans.sub ~ Year + Station + as.factor(biweekly) + #Station*biweekly + 
                          Temp_Top + Salin_Top, 
                          pru.env.biwk.std, perm=999, by = "margin")
permanova.biwk # Top PERMANOVA model. Explains 1-0.5254 = 47.5% of variation
# everything is very significant. 

# Explore using biweekly as a factor?

adonis(catchmatrix.biwk.cpue.stdtrans.sub ~ Temp_Top + Salin_Top + Station + as.factor(biweekly) + Year, 
       pru.env.biwk.std, perm=999, by = "terms")
# Sequential just to show how much effect Year has, once all other variables are accounted for

adonis(catchmatrix.biwk.cpue.stdtrans.sub ~ Station + as.factor(biweekly) + Year + Temp_Top + Salin_Top, 
       pru.env.biwk.std, perm=999, by = "terms")
# Sequential to see effect of env var, once spatio-temporal var accounted for: still some effects from temp/salin



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


# Final Analysis
EWsimper <- data.frame((summary(simper(catchmatrix.biwk.stdtrans, (pru.env.biwk.std %>% 
                        mutate(EW_stn = ifelse(Station==230 | Station == 214, "EastStn", "WestStn")))$EW_stn )))$EastStn_WestStn )
EWsimper






#finalcolors_BW <- c("#d6d6d6", "#575757", "#9b9b9b", "#7a7a7a") #ellipse colors
finalcolors_BW <- c("#9b9b9b", "#494949", "#9b9b9b", "#494949") #ellipse colors

# Final plot in B&W version
Fig3nMDS <- ggplot(nmdspoints.biwk %>%
         mutate(Station=factor(Station, levels=c("230", "214", "218", "220"))), # reorder factor by east-west
       aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(shape=Station, fill = Station), cex=4) + 
  scale_shape_manual(values=c(21, 21, 24, 24)) +
  scale_fill_manual(values =  c("#ffffff", "#494949", "#ffffff", "#494949")) +
  stat_ellipse(aes(group=Station, color=Station, lty=Station), size=2, show.legend = FALSE) + # turn on/off legend
  scale_linetype_manual(values=c(2,1,2,1)) + # can change linetypes here
  theme_bw() + theme(panel.grid.minor = element_blank(), 
                     text=element_text(family="Arial", size=12)) +
  coord_cartesian(xlim = c(-0.31, 0.31), ylim = c(-0.27, 0.21), expand = FALSE) + # added slight above and below
  geom_segment(data = data.frame(env.vectors.biwk$vectors$arrows) %>% 
                 cbind(r2=env.vectors.biwk$vectors$r, pval = env.vectors.biwk$vectors$pvals), 
               aes(x=0, xend=NMDS1 * (r2^0.5)/3, y=0, yend=NMDS2 * (r2^0.5)/3), cex = 2, arrow = arrow(length = unit(12, "points"))) +
  geom_label(data = data.frame(wascores(totalNMDS.biwk$points, w = catchmatrix.biwk.cpue)) %>% 
               mutate(species = rownames(.)) %>%
               filter(species %in% topcorrspp), 
             aes(x=MDS1, y=MDS2, label = species), family = "Arial") +
  scale_color_manual(values =  finalcolors_BW) + #ellipse color
  #annotate("text", x=0.21, y=0.15, label= "Salinity", family = "Arial") +
  #annotate("text", x=-0.05, y=0.01, label= "Year", family = "Arial") +
  #annotate("text", x=0.2, y=-0.13, label= "Biweekly", family = "Arial") +
  #annotate("text", x=0.02, y=-0.07, label= "Temp", family = "Arial") + 
  annotate("text", x=0.2, y=-0.25, label= "Stress = 0.156", family = "Arial") + 
  geom_label(aes(x=0.2, y=-0.14, label = "Biweekly"), fill = "black", color = "white", family = "Arial") +
  geom_label(aes(x=0.21, y=0.15, label = "Salinity"), fill = "black", color = "white", family = "Arial") +
  geom_label(aes(x=-0.05, y=0.01, label = "Year"), fill = "black", color = "white", family = "Arial") +
  geom_label(aes(x=0.02, y=-0.07, label = "Temp"), fill = "black", color = "white", family = "Arial") +
  annotate("text", x=0.05, y=-0.26, label= "230", family = "Arial") +
  annotate("text", x=0.08, y=-0.19, label= "214", family = "Arial") +
  annotate("text", x=-0.18, y=0.13, label= "218", family = "Arial") +
  annotate("text", x=-0.15, y=0.19, label= "220", family = "Arial") +
  geom_label(data=as.data.frame(env.vectors.biwk$factors$centroids) %>% mutate(Station = as.factor(substr(row.names(.), 8, 10) )),
             aes(x=NMDS1, y=NMDS2, label="X"), color="black", fill="light gray", cex=4, show.legend = FALSE)
Fig3nMDS

#ggsave(plot = Fig3nMDS, "plotexports/Fig3_biwknMDS_BW.png", dpi = 600, width = 7.5, height = 5)
#ggsave(plot = Fig3nMDS, "plotexports/Fig3_biwknMDS_BW.eps", dpi = 600, width = 7.5, height = 5, device = "eps")
# eps save won't work "family 'Arial' not included in postscript() device"






# #Color nMDS Plot
# finalcolors <- c("#b9a3c6", "#0063a0", "#8cc687", "#d86a6a")
# 
# ggplot(nmdspoints.biwk, aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color=Station), cex=5) + 
#   #stat_ellipse(aes(group=Station, color=Station), size=2, linetype=2, show.legend = FALSE) +
#   theme_bw() + theme(panel.grid.minor = element_blank(), 
#                      text=element_text(family="Times New Roman", size=12)) +
#   geom_segment(data = data.frame(env.vectors.biwk$vectors$arrows) %>% 
#                  cbind(r2=env.vectors.biwk$vectors$r, pval = env.vectors.biwk$vectors$pvals), 
#                aes(x=0, xend=NMDS1 * (r2^0.5)/3, y=0, yend=NMDS2 * (r2^0.5)/3), cex = 2, arrow = arrow(length = unit(12, "points"))) +
#   geom_label(data = data.frame(wascores(totalNMDS.biwk$points, w = catchmatrix.biwk.cpue)) %>% 
#                mutate(species = rownames(.)) %>%
#                filter(species %in% topcorrspp), 
#              aes(x=MDS1, y=MDS2, label = species), family = "Times New Roman") +
#   #scale_color_manual(values =  finalcolors) +
#   scale_color_manual(values =  brewer.pal(4, "Set2")) 
# # annotate("text", x=0.21, y=0.15, label= "Salinity", family = "Times New Roman") +
# # annotate("text", x=-0.05, y=0.01, label= "Year", family = "Times New Roman") +
# # annotate("text", x=0.2, y=-0.13, label= "Biweekly", family = "Times New Roman") +
# # annotate("text", x=0.02, y=-0.07, label= "Temp", family = "Times New Roman") 
# # geom_label(data=as.data.frame(env.vectors.biwk$factors$centroids) %>% mutate(Station = as.factor(substr(row.names(.), 8, 10) )), 
# #           aes(x=NMDS1, y=NMDS2, label="X", color=Station), fill="#e8e8e8", cex=5, show.legend = FALSE)
# 
# #ggsave("plotexports/Fig_biwknMDS.png", dpi = 300, width = 7.5, height = 5)


temp1 <- pru.env.biwk %>% 
  dplyr::select(Year, biweekly, Station, Salin_Top, Temp_Top) %>%
  mutate(Station = as.factor(Station))
temp1


temp <- envfit(totalNMDS.biwk, temp1,   
       na.rm = TRUE, permutations = 999)
temp
str(temp)

temp$factors$centroids




