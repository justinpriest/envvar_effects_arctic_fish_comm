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




spp_richness
read.csv("data/speciesrichness_annual_2001-2018.csv")

read.csv("data/speciesrichness_biwk_2001-2018.csv")

















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



#############################
#### CATCH DATA ####
#############################

# Next, turn catch data into a 'wide' catch matrix (Species by Year/Station)
# drop enivron data too, just focusing on catch here
catchmatrix.all <- catchenviron %>% group_by(Year, Station, Species) %>% summarise(anncount = sum(totcount)) %>%
  spread(Species, value = anncount) %>% replace(., is.na(.), 0) %>% ungroup() #Make sure to replace NAs with 0
catchmatrix.all <- catchmatrix.all %>% arrange(Year, Station) # Just to double check order is correct

catchmatrix.day <- catchenviron %>% group_by(EndDate, Station, Species) %>% summarise(count = sum(totcount)) %>%
  spread(Species, value = count) %>% replace(., is.na(.), 0) %>% ungroup()
catchmatrix.day <- catchmatrix.day %>% arrange(EndDate, Station) # Just to double check order is correct

catchmatrix.biwk <- catchenviron %>% group_by(Year, biweekly, Station, Species) %>% summarise(count = sum(totcount)) %>%
  spread(Species, value = count) %>% replace(., is.na(.), 0) %>% ungroup()
catchmatrix.biwk <- catchmatrix.biwk %>% arrange(Year, biweekly, Station) # Just to double check order is correct


# Account for varying levels of effort


# Exclude rare species:
# Right now analysis is for only species >100 fish, all years combined. Change 100 to 0 to inc all spp
# the following code makes a list of the species to keep which have a threshold of 100 currently,
# then filters based on this list, then turns back into wide format
keepspp <- (catchmatrix.all %>% gather(Species, counts, -Year, -Station) %>%
              group_by(Species) %>% summarize(counts = sum(counts)) %>% filter(counts > 100))$Species
catchmatrix <- catchmatrix.all %>% gather(Species, counts, -Year, -Station) %>% filter(Species %in% keepspp) %>%
  spread(Species, value = counts)
catchmatrix.day <- catchmatrix.day %>% gather(Species, counts, -EndDate, -Station) %>% filter(Species %in% keepspp) %>%
  spread(Species, value = counts)
catchmatrix.biwk <- catchmatrix.biwk %>% gather(Species, counts, -Year, -biweekly, -Station) %>% filter(Species %in% keepspp) %>%
  spread(Species, value = counts)

# Add row names for ease of checking later
rownames(catchmatrix)      <- paste0(catchmatrix$Year, catchmatrix$Station)
rownames(catchmatrix.day)  <- paste0(year(catchmatrix.day$EndDate), yday(catchmatrix.day$EndDate), catchmatrix.day$Station)
rownames(catchmatrix.biwk) <- paste0(catchmatrix.biwk$Year, catchmatrix.biwk$biweekly, catchmatrix.biwk$Station)




# Create CPUE dataframe (just for biweekly right now)
catchmatrix.biwk.cpue <- catchmatrix.biwk %>% gather(Species, catch, -Year, -biweekly, -Station) %>%
  left_join(effort.biwk.stn, by = c("Year" = "Year", "biweekly" = "biweekly", "Station" = "Station")) %>%
  mutate(CPUE_biwk = catch/(Effort_NetHrs/(24*2*14))) %>% #This is biweekly CPUE (2 nets fishing 24 hrs/day, for ~14 days)
  dplyr::select(-catch, -Effort_NetHrs) %>%
  spread(Species, value = CPUE_biwk)

rownames(catchmatrix.biwk.cpue) <- paste0(catchmatrix.biwk$Year, catchmatrix.biwk$biweekly, catchmatrix.biwk$Station)



# Delete the Year & Stn cols (can't be present for PERMANOVA)
#pru.env.ann <- catchmatrix %>% dplyr::select(Year, Station)
catchmatrix <- catchmatrix %>% dplyr::select(-Year, -Station)
catchmatrix.all <- catchmatrix.all %>% dplyr::select(-Year, -Station)
catchmatrix.day <- catchmatrix.day %>% dplyr::select(-EndDate, -Station)
catchmatrix.biwk <- catchmatrix.biwk %>% dplyr::select(-Year, -biweekly, -Station)
catchmatrix.biwk.cpue <- catchmatrix.biwk.cpue %>% dplyr::select(-Year, -biweekly, -Station)

# standardize catches 0 to 1 (1 is max catch in a given year/station)
# note that the order corresponds to the now deleted Year/station combo. DON'T CHANGE ORDER
catchmatrix.std <- catchmatrix 
for (i in 1:ncol(catchmatrix.std)){ #make sure year/stn cols already dropped
  catchmatrix.std[i] <- catchmatrix[i]/max(catchmatrix[i])}

catchmatrix.day.std <- catchmatrix.day
for (i in 1:ncol(catchmatrix.day.std)){ #make sure year/stn cols already dropped
  catchmatrix.day.std[i] <- catchmatrix.day[i]/max(catchmatrix.day[i])}

catchmatrix.biwk.stdtrans <- catchmatrix.biwk
for (i in 1:ncol(catchmatrix.biwk.stdtrans)){ #make sure year/stn cols already dropped
  catchmatrix.biwk.stdtrans[i] <- (catchmatrix.biwk[i]^0.25)/max((catchmatrix.biwk[i]^0.25))}
#using 4th root tranform

catchmatrix.biwk.cpue.stdtrans <- catchmatrix.biwk.cpue
for (i in 1:ncol(catchmatrix.biwk.cpue.stdtrans)){ #make sure year/stn cols already dropped
  catchmatrix.biwk.cpue.stdtrans[i] <- (catchmatrix.biwk.cpue[i]^0.25)/max((catchmatrix.biwk.cpue[i]^0.25))}

# make sure that 'catchmatrix' and catchmatrix.std are both set up in same order as pru.env.ann
#hist(catchmatrix.biwk.stdtrans$ARCS)






##########################################
#APPENDIX

# Wind Data Creation
# The daily wind data was created using the following code:
# windhourly <- read.csv("../../Slope Project/Data/deadhorsewind_2001-2018_hourly.csv", header = TRUE, 
#          stringsAsFactors = FALSE) %>%
#   select(DATE, HOURLYWindSpeed, HOURLYWindDirection) %>% #remove extraneous data
#   rename(Date = DATE,
#          windhrly_mph = HOURLYWindSpeed,
#          windhrly_dir = HOURLYWindDirection) %>% 
#   mutate(Date = ymd(as.POSIXct(Date, format = "%m/%d/%Y")),
#          month = month(Date)) %>%
#   mutate_if(is.character, as.numeric) %>% select(Date, month, everything())
# 
# library(CircStats)
# winddaily <- windhourly %>% mutate(Year = year(Date)) %>% group_by(Date) %>% 
#   summarise(dailymeanspeed = mean(windhrly_mph, na.rm = TRUE),
#             dailymeandir = ((circ.mean(2*pi*na.omit(windhrly_dir)/360))*(360 / (2*pi))) %%360 )
# 
# write.csv(winddaily, file="deadhorsewind_2001-2018_summarized.csv")



# Effort Data Creation 
# effort <- full_join(all.len %>%
#                       distinct(EndDate, Net), 
#                     allcatch %>% distinct(EndDate, Net), 
#                     by = c("EndDate", "Net")) %>% 
#   left_join(all.len %>%
#               dplyr::select(EndDate, Net, StartDateTime, EndDateTime), 
#             by = c("EndDate", "Net")) %>% 
#   distinct(EndDate, Net, StartDateTime, EndDateTime) %>%
#   #this first mutate MANUALLY adds in skipped dates
#   mutate(StartDateTime=replace(StartDateTime, EndDate=="2018-07-10" & Net == "230N", 
#                                as.POSIXct("2018-07-09 09:25")), 
#          EndDateTime=replace(EndDateTime, EndDate=="2018-07-10" & Net == "230N", 
#                              as.POSIXct("2018-07-10 08:55")), 
#          
#          StartDateTime=replace(StartDateTime, EndDate=="2018-07-24" & Net == "220W", 
#                                as.POSIXct("2018-07-23 09:15")), 
#          EndDateTime=replace(EndDateTime, EndDate=="2018-07-24" & Net == "220W", 
#                              as.POSIXct("2018-07-24 14:45"))) %>% 
#   mutate(Year = year(EndDate),
#          Effort_NetHrs = as.numeric(EndDateTime - StartDateTime),
#          Station = substr(Net, 1,3)) %>%
#   arrange(EndDate, Net, Station)
# # A few NAs occur when we have catch but no lengths for that day 
# # (because the times are recorded in the length dataframe!)
# 
# write.csv(effort, "prudhoeeffort_2001-2018.csv")

