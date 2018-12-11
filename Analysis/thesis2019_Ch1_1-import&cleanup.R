###########
# Data Import and Cleanup

# This is the first script for Justin Priest's thesis analyses, Chapter 1
# Data is pulled from an .RData file, which was compiled using a script from
# "/Grad School/Slope Project/PB_DataImport2001-2016.R" which imports data from 
# the original excel file sent over by LGL Ecological Research. Data from 2017 
# and 2018 were entered into a Google Sheets doc, then exported to .csv files 
# which were imported in scripts "2017 import.R" and "2018 import.R"



library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)
library(CircStats)


# Read in the length, catch, and environmental data, plus species lookup table
load("Data/PrudhoeCatch&LengthDataset_2001-2018_Version11.Rdata")

# Clean up dataframe
allcatch <- allcatch %>% mutate(Net = factor(Net), 
                    Station = factor(substr(Net, 1, 3)))

env_allyears <- env_allyears %>% mutate(Station = factor(Station))


allcatch <- allcatch %>% add_column(Month = month(allcatch$EndDate), .after = 2)
env_allyears <- env_allyears %>% add_column(Month = month(env_allyears$Date), .after = 2)


# Separate out the temp and salinity. Drop the columns that are irrelevant
# Put data into 'long' format. NOTE: They used a slightly different definition 
# of the bottom in early years ('bottom 1.5'). I define this as analogous to 'bottom'
watertemps <- env_allyears %>% dplyr::select(-c(Salin_Top, Salin_Mid, Salin_Bot, Salin_Bot_1.5)) 
watersalin <- env_allyears %>% dplyr::select(-c(Temp_Top, Temp_Mid, Temp_Bot, Temp_Bot_1.5)) 

# if needed for later: watersalin %>% gather(depth, salin_ppt, -c(Year, Date, Month, Station))


##### WIND #####
# wind and air temp data from NOAA NCEI
# https://www.ncdc.noaa.gov/cdo-web/datasets/LCD/stations/WBAN:27406/detail
# https://www1.ncdc.noaa.gov/pub/data/cdo/documentation/LCD_documentation.pdf
# Station WBAN:27406
# These data were downloaded in 10-year increments (LCD CSV). 
# In appendix (at end of this doc, I summarized hourly data into a daily summary)


deadhorsewind <- read.csv("Data/deadhorsewind_2001-2018_daily_summarized.csv", header = TRUE, 
                          stringsAsFactors = FALSE) %>% 
  mutate(Date = ymd(as.POSIXct(Date, format = "%m/%d/%Y")),
         Year = year(Date),
         month = month(Date),
         dailymeanspeed = dailymeanspeed * 1.60934) %>%
  rename(dailymeanspeed_kph = dailymeanspeed) %>%
  filter(month == 7 | month == 8) %>%
  dplyr::select(Date, Year, month, everything()) # reorder month column



##### DISCHARGE #####
# Data from USGS https://waterdata.usgs.gov/nwis/uv/?site_no=15908000
# https://nwis.waterdata.usgs.gov/usa/nwis/uv/?cb_00060=on&format=rdb&site_no=15908000&period=&begin_date=2001-01-01&end_date=2018-09-01
sagdisch <- read.csv("Data/SagDischargeDaily_2001-2018.csv", header = TRUE) %>%
  mutate(datetime = as.POSIXct(paste0(date, " ", time), format = "%m/%d/%Y %H:%M"),
         Date = as_date(datetime),
         hour = hour(datetime) ) %>%
  dplyr::select(Date, hour, disch_cfs) %>%
  group_by(Date) %>% summarize(meandisch_cfs = mean(disch_cfs, na.rm = TRUE))



##### Catch #####
# Join the catch and environ data. 
catchenviron <- left_join(allcatch, watersalin %>% dplyr::select(-c(Month)), 
                          by = c("Year" = "Year", "EndDate" = "Date", "Station" = "Station")) %>% 
  left_join(watertemps %>% dplyr::select(-c( Month)), 
            by = c("Year" = "Year", "EndDate" = "Date", "Station" = "Station")) %>% 
  left_join(deadhorsewind %>% dplyr::select(-month), by = c("Year" = "Year", "EndDate" = "Date")) %>%
  left_join(sagdisch, by = c("EndDate" = "Date")) 

catchenviron <- catchenviron %>% 
  mutate(biweekly = ifelse(day.of.year <= 196, 1, 
                           ifelse(day.of.year > 196 & day.of.year <= 213, 2, #btwn July 15 and Aug 1
                                  ifelse(day.of.year > 213 & day.of.year <= 227, 3, #Aug 1 - 15
                                         ifelse(day.of.year > 227, 4, NA))))) # after Aug 15

###################################
# Join All Environmental Data

# This creates a dataframe summarizing annual environmental conditions at each site for each year
# and keep in same order as the corresponding catch dataframe. 
# The only tricky thing done is that I converted wind from degrees (0-360) to East-West using trig,
# before this I took the mean using circular averaging. Salin & temp are from midwater sampling
pru.env.ann <- catchenviron %>% dplyr::distinct(Year, Station) %>% arrange(Year, Station) %>%
  #the above section creates a dataframe of year/station combos. Need to decide what to do about Stn 231
  left_join(deadhorsewind %>% mutate(Year = year(Date)) %>% group_by(Year) %>% 
              summarise(annwindspeed_kph = mean(dailymeanspeed_kph, na.rm = TRUE),
                        annwinddir = ((circ.mean(2*pi*na.omit(dailymeandir)/360))*(360 / (2*pi))) %%360 ),
            by = c("Year" = "Year")) %>% 
  left_join(sagdisch %>% mutate(Year = year(Date), month = month(Date)) %>% 
              group_by(Year) %>% filter(month == 7 | month == 8) %>%
              summarise(anndisch_cfs = mean(meandisch_cfs, na.rm = TRUE)),by = c("Year" = "Year") ) %>% 
  left_join(watersalin %>% group_by(Year, Station) %>% summarise(annsal_ppt = mean(Salin_Mid, na.rm = TRUE)), 
            by = c("Year" = "Year", "Station" = "Station")) %>%
  left_join(watertemps %>% group_by(Year, Station) %>% summarise(anntemp_c = mean(Temp_Mid, na.rm = TRUE)), 
            by = c("Year" = "Year", "Station" = "Station")) %>%
  mutate(Year =  factor(Year, ordered = TRUE), # PERMANOVA will want it ordered
         annwinddir_ew = sin(annwinddir * pi / 180)) # this changes from polar coords to cartesian east-west (-1=W, 1=E)


#create similar setup on a daily scale
pru.env.day <- catchenviron %>% dplyr::distinct(EndDate, Station) %>% 
  mutate(Station=replace(Station, Station==231, 214)) %>% arrange(EndDate, Station) %>% 
  left_join(deadhorsewind %>% dplyr::select(-month), by = c("EndDate" = "Date")) %>%
  left_join(sagdisch, by = c("EndDate" = "Date")) %>%
  left_join(watersalin %>% dplyr::select(Date, Station, Salin_Top, Salin_Mid), 
            by = c("EndDate" = "Date", "Station" = "Station")) %>%
  left_join(watertemps %>% dplyr::select(Date, Station, Temp_Top, Temp_Mid), 
            by = c("EndDate" = "Date", "Station" = "Station")) %>%
  mutate(#EndDate =  factor(EndDate, ordered = TRUE), # PERMANOVA will want it ordered
    Year = year(EndDate),
    Station = factor(Station), 
    winddir_ew = sin(dailymeandir * pi / 180)) %>% # changes from degrees to cartesian east-west (-1=W, 1=E)
  arrange(EndDate, Station)


# finally, aggregate daily info to a biweekly scale
pru.env.biwk <- pru.env.day %>% 
  mutate(Year = year(EndDate), 
         day.of.year = yday(EndDate), 
         biweekly = factor(ifelse(day.of.year <= 196, 1, 
                                  ifelse(day.of.year > 196 & day.of.year <= 213, 2, #btwn July 15 and Aug 1
                                         ifelse(day.of.year > 213 & day.of.year <= 227, 3, #Aug 1 - 15
                                                ifelse(day.of.year > 227, 4, NA))))) ) %>% # after Aug 15
  group_by(Year, biweekly, Station) %>%
  summarise(biwkmeanspeed_kph = mean(dailymeanspeed_kph, na.rm = TRUE),
            biwkmeandir =  ((circ.mean(2*pi*na.omit(dailymeandir)/360))*(360 / (2*pi))) %%360,
            meandisch_cfs = mean(meandisch_cfs , na.rm=TRUE),
            Salin_Top = mean(Salin_Top, na.rm=TRUE),
            Salin_Mid = mean(Salin_Mid, na.rm=TRUE),
            Temp_Top = mean( Temp_Top, na.rm=TRUE),
            Temp_Mid = mean(Temp_Mid, na.rm=TRUE),
            winddir_ew = mean(winddir_ew, na.rm=TRUE))





#############################
### Catch Data

# now turn catch data into a 'wide' catch matrix (Species by Year/Station)
# drop enivron data too, just focusing on catch here
catchmatrix.all <- catchenviron %>% group_by(Year, Station, Species) %>% summarise(anncount = sum(totcount)) %>%
  spread(Species, value = anncount) %>% replace(., is.na(.), 0) %>% ungroup()
catchmatrix.all$Station[catchmatrix.all$Station == 231] <- 214 # Treat the 231 Station as the precursor to 214
catchmatrix.all <- catchmatrix.all %>% arrange(Year, Station) #was out of order in 2001 after station name change

catchmatrix.day <- catchenviron %>% group_by(EndDate, Station, Species) %>% summarise(count = sum(totcount)) %>%
  spread(Species, value = count) %>% replace(., is.na(.), 0) %>% ungroup()
catchmatrix.day$Station[catchmatrix.day$Station == 231] <- 214 # Treat the 231 Station as the precursor to 214
catchmatrix.day <- catchmatrix.day %>% arrange(EndDate, Station) #was out of order in 2001 after station name change

catchmatrix.biwk <- catchenviron %>% group_by(Year, biweekly, Station, Species) %>% summarise(count = sum(totcount)) %>%
  spread(Species, value = count) %>% replace(., is.na(.), 0) %>% ungroup()
catchmatrix.biwk$Station[catchmatrix.biwk$Station == 231] <- 214 # Treat the 231 Station as the precursor to 214
catchmatrix.biwk <- catchmatrix.biwk %>% arrange(Year, biweekly, Station) #was out of order in 2001 after station name change


# Exclude rare species:
# catchmatrix$Station <- as.numeric(catchmatrix$Station) #cheat to include it in shortcut
# catchmatrix <- catchmatrix[, which(colSums(catchmatrix) > 100)] 

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


rownames(catchmatrix) <- paste0(catchmatrix$Year, catchmatrix$Station)
rownames(catchmatrix.day) <- paste0(year(catchmatrix.day$EndDate), yday(catchmatrix.day$EndDate), catchmatrix.day$Station)
rownames(catchmatrix.biwk) <- paste0(catchmatrix.biwk$Year, catchmatrix.biwk$biweekly, catchmatrix.biwk$Station)



#pru.env.ann <- catchmatrix %>% dplyr::select(Year, Station)
catchmatrix <- catchmatrix %>% dplyr::select(-Year, -Station)
catchmatrix.all <- catchmatrix.all %>% dplyr::select(-Year, -Station)
catchmatrix.day <- catchmatrix.day %>% dplyr::select(-EndDate, -Station)
catchmatrix.biwk <- catchmatrix.biwk %>% dplyr::select(-Year, -biweekly, -Station)


# standardize catches 0 to 1 (1 is max catch in a given year/station)
#note that the order corresponds to the now deleted Year/station combo. DON'T CHANGE ORDER
catchmatrix.std <- catchmatrix 
for (i in 1:ncol(catchmatrix.std)){ #make sure year/stn cols already dropped
  catchmatrix.std[i] <- catchmatrix[i]/max(catchmatrix[i])}

catchmatrix.day.std <- catchmatrix.day
for (i in 1:ncol(catchmatrix.day.std)){ #starts at 3 to exclude Year and station cols
  catchmatrix.day.std[i] <- catchmatrix.day[i]/max(catchmatrix.day[i])}

catchmatrix.biwk.stdtrans <- catchmatrix.biwk
for (i in 1:ncol(catchmatrix.biwk.stdtrans)){ #starts at 3 to exclude Year and station cols
  catchmatrix.biwk.stdtrans[i] <- (catchmatrix.biwk[i]^0.25)/max((catchmatrix.biwk[i]^0.25))}
#using 4th root tranform

# make sure that 'catchmatrix' and catchmatrix.std are both set up in same order as pru.env.ann
#hist(catchmatrix.biwk.stdtrans$ARCS)


##########################################
#APPENDIX
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