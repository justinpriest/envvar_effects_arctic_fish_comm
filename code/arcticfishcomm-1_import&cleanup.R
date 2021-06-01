# Effects of Environmental Variables on a Nearshore Arctic Fish Community, 2001–2018

# Authors: Justin Priest, Franz Mueter, Scott Raborn, and Trent Sutton
# Last Major Update: June 2021
# Contact: justin.priest@alaska.gov




# DATA IMPORT AND CLEANUP


# This is the first script, importing data and cleaning it up for inclusion in this analysis. 
# For a full methodology of how the original LGL data and subsequent UAF data
#  were turned into the data present in this repository, please contact the author. 
# To use these data, contact Justin Priest
# Some of the data imported were created from the full database documented in file
# "arcticfishcomm-0_datacreation.R"


library(tidyverse)
#library(tidyr)
#library(tibble)
library(lubridate)
library(CircStats)
library(here)

source(here::here("code/functions.R"))






















###########################
#### EFFORT ####
###########################
# See end appendix for documenting how the effort file was created
effort <- read_csv("data/prudhoeeffort_2001-2018.csv") %>%
  mutate(EndDate = ymd(as.POSIXct(EndDate, format = "%m/%d/%Y")),
         StartDateTime = as.POSIXct(StartDateTime, format = "%m/%d/%Y %H:%M"),
         EndDateTime = as.POSIXct(EndDateTime, format = "%m/%d/%Y %H:%M"))

effort.biwk.stn <- addbiwknum(effort) %>% group_by(Year, biweekly, Station) %>% 
  summarise(Effort_NetHrs = sum(Effort_NetHrs, na.rm = TRUE)) 


effort.biwk <- addbiwknum(effort) %>% filter(Station != 231) %>%  
  group_by(Year, biweekly) %>% summarise(biweeklyeffort = sum(Effort_NetHrs, na.rm = TRUE)) %>%
  mutate(sample_proportion = biweeklyeffort / (48*4*15)) #2880 is 4 nets with two cod ends, fishing 24hrs for 15 days (15*24*4*2)
# sample_proportion is the amount of sampling relative to 100% coverage (but can be slightly higher than 100%)
# We're using 15 days as "full effort" sampling but e.g., July 16-31 is 16 days. This is why we scale for effort though


# Weekly effort is used to look at species richness increases over the season
effort.wk <- addweeknum(effort) %>% filter(Station != 231) %>%  
  group_by(Year, week) %>% summarise(weekeffort = sum(Effort_NetHrs, na.rm = TRUE)) %>%
  mutate(sample_proportion = weekeffort / (48*4*7)) #1344 is 4 nets with two cod ends, fishing 24/7 (7*24*4*2)
# sample_proportion is the amount of sampling relative to 100% coverage (but can be slightly higher than 100%)




##########################
#### SPECIES RICHNESS ####
##########################


spp_richness <- read.csv("data/prudhoespeciesrichness_annual_2001-2018.csv")



spp_richness.biwk <- read.csv("data/prudhoespeciesrichness_biwk_2001-2018.csv") %>% 
  left_join(effort.biwk %>% dplyr::select(biweekly, sample_proportion), 
            by = c("Year" = "Year", "biweekly" = "biweekly")) %>%
  mutate(num_spp_adj = num_spp * sample_proportion)





###########################
#### CATCH MATRICES ####
###########################

catchmatrix.biwk <- read.csv("data/prudhoecatchmatrixbiwk_2001-2018.csv") 
# %>%
#   rename("yrbiwkstn" = "X") %>% 
#   mutate(Year = as.numeric(substr(yrbiwkstn, 1, 4)),
#          biweekly = as.numeric(substr(yrbiwkstn, 5, 5)),
#          Station = as.numeric(substr(yrbiwkstn, 6, 8))) %>%
#   dplyr::select(Year, biweekly, Station, everything()) %>%
#   dplyr::select(-yrbiwkstn)


### Account for varying levels of effort ###
# Create CPUE dataframe (just for biweekly)
catchmatrix.biwk.cpue <- catchmatrix.biwk %>% gather(Species, catch, -Year, -biweekly, -Station) %>%
  left_join(effort.biwk.stn, by = c("Year" = "Year", "biweekly" = "biweekly", "Station" = "Station")) %>%
  mutate(CPUE_biwk = catch/(Effort_NetHrs/(24*2*14))) %>% #This is biweekly CPUE (2 nets fishing 24 hrs/day, for ~14 days)
  dplyr::select(-catch, -Effort_NetHrs) %>%
  spread(Species, value = CPUE_biwk)

# Need rownames for PERMANOVA
rownames(catchmatrix.biwk) <- paste0(catchmatrix.biwk$Year, catchmatrix.biwk$biweekly, catchmatrix.biwk$Station)
rownames(catchmatrix.biwk.cpue) <- paste0(catchmatrix.biwk$Year, catchmatrix.biwk$biweekly, catchmatrix.biwk$Station)

# PERMANOVA can't have any non-species columns, so drop those here
catchmatrix.biwk <- catchmatrix.biwk %>% dplyr::select(-Year, -biweekly, -Station)
catchmatrix.biwk.cpue <- catchmatrix.biwk.cpue %>% dplyr::select(-Year, -biweekly, -Station)




catchmatrix.biwk.stdtrans <- catchmatrix.biwk
for (i in 1:ncol(catchmatrix.biwk.stdtrans)){ #make sure year/stn cols already dropped
  catchmatrix.biwk.stdtrans[i] <- (catchmatrix.biwk[i]^0.25)/max((catchmatrix.biwk[i]^0.25))}
#using 4th root tranform

catchmatrix.biwk.cpue.stdtrans <- catchmatrix.biwk.cpue
for (i in 1:ncol(catchmatrix.biwk.cpue.stdtrans)){ #make sure year/stn cols already dropped
  catchmatrix.biwk.cpue.stdtrans[i] <- (catchmatrix.biwk.cpue[i]^0.25)/max((catchmatrix.biwk.cpue[i]^0.25))}

# make sure that 'catchmatrix' and catchmatrix.std are both set up in same order as pru.env.ann
#hist(catchmatrix.biwk.stdtrans$ARCS)




#########################
##### Rare Species ######
rarespp.biwk.pres <- read_csv("data/prudhoe_rarespeciesbiwk_2001-2018.csv") %>%
  mutate(Station = as.factor(Station))





####################################
##### Environmental Variables ######
pru.env.biwk <- read_csv("data/prudhoe_envdatabiweekly_2001-2018.csv")


# These are the "common" species. See script_0 for how this list was determined.
# Changes in this list shouldn't affect any analysis, but is used for plotting
keepspp <- c("ARCD", "ARCS", "ARFL", "BDWF", "CAPE", "DLVN", 
  "FHSC", "GRAY", "HBWF", "LSCS", "NNSB", "PCHG", 
  "PINK", "RBSM", "RDWF", "SFCD", "THSB")


