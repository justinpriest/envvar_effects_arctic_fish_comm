# Effects of Environmental Variables on a Nearshore Arctic Fish Community, 2001–2018
# Authors: Justin Priest, Franz Mueter, Scott Raborn, and Trent Sutton
# Last Major Update: June 2021
# Contact: justin.priest@alaska.gov



# This script documents how data were created
# Because of data restrictions, not all of these data are public.
# As such, much of the code below won't run
# Begin with file "arcticfishcomm-1_import&cleanup.R"


# Read in the length, catch, and environmental data, plus species lookup table
load(here::here("data/PrudhoeCatch&LengthDataset_2001-2018_Version11.Rdata"))


#Will need this to calculate effort later 
all.len <- all.len %>% mutate(Net = paste0(Station, Side))

# Clean up dataframes
allcatch <- allcatch %>% mutate(Net = factor(Net), Station = factor(substr(Net, 1, 3)))
env_allyears <- env_allyears %>% mutate(Station = factor(Station))

allcatch <- allcatch %>% add_column(Month = month(allcatch$EndDate), .after = 2)
env_allyears <- env_allyears %>% add_column(Month = month(env_allyears$Date), .after = 2)

# Separate out the temp and salinity. Drop the columns that are irrelevant
# NOTE: They used a slightly different definition of the bottom in early years ('bottom 1.5'). 
# I define bottom 1.5 as analogous to 'bottom'
watertemps <- env_allyears %>% dplyr::select(-c(Salin_Top, Salin_Mid, Salin_Bot, Salin_Bot_1.5)) 
watersalin <- env_allyears %>% dplyr::select(-c(Temp_Top, Temp_Mid, Temp_Bot, Temp_Bot_1.5)) 



##### WIND #####
# wind and air temp data from NOAA NCEI
# https://www.ncdc.noaa.gov/cdo-web/datasets/LCD/stations/WBAN:27406/detail
# https://www1.ncdc.noaa.gov/pub/data/cdo/documentation/LCD_documentation.pdf
# Station WBAN:27406
# These data were downloaded in 10-year increments (LCD CSV). 
# In appendix (at end of this doc, I summarized hourly data into a daily summary)

deadhorsewind <- read.csv(here::here("data/deadhorsewind_2001-2018_daily_summarized.csv"), header = TRUE, 
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

sagdisch <- read.csv(here::here("data/SagDischargeDaily_2001-2018.csv"), header = TRUE) %>%
  mutate(datetime = as.POSIXct(paste0(date, " ", time), format = "%m/%d/%Y %H:%M"),
         Date = as_date(datetime),
         hour = hour(datetime) ) %>%
  dplyr::select(Date, hour, disch_cfs) %>%
  group_by(Date) %>% summarize(meandisch_cfs = mean(disch_cfs, na.rm = TRUE))



##### CATCH #####
# Join the catch and environ data. 
catchenviron <- left_join(allcatch, watersalin %>% dplyr::select(-c(Month)), 
                          by = c("Year" = "Year", "EndDate" = "Date", "Station" = "Station")) %>% 
  left_join(watertemps %>% dplyr::select(-c( Month)), 
            by = c("Year" = "Year", "EndDate" = "Date", "Station" = "Station")) %>% 
  left_join(deadhorsewind %>% dplyr::select(-month), by = c("Year" = "Year", "EndDate" = "Date")) %>%
  left_join(sagdisch, by = c("EndDate" = "Date")) %>%
  addbiwknum() %>% filter(Station != 231) 



###################################
# Join All Environmental Data

# This creates a dataframe summarizing annual environmental conditions at each site for each year
# and keep in same order as the corresponding catch dataframe. 
# The only tricky thing done is that I converted wind from degrees (0-360) to East-West using trig,
# before this I took the mean using circular averaging. Salin & temp are from midwater sampling
# Wind data were summarized before merging using code at the end of this script
pru.env.ann <- catchenviron %>% dplyr::distinct(Year, Station) %>% arrange(Year, Station) %>%
  #the above section creates a dataframe of year/station combos
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
  arrange(EndDate, Station) %>% 
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
  mutate(Year = year(EndDate)) %>% addbiwknum() %>% #add year and biweekly cols
  group_by(Year, biweekly, Station) %>%
  summarise(biwkmeanspeed_kph = mean(dailymeanspeed_kph, na.rm = TRUE),
            biwkmeandir = ((circ.mean(2*pi*na.omit(dailymeandir)/360))*(360 / (2*pi))) %%360,
            meandisch_cfs = mean(meandisch_cfs, na.rm=TRUE),
            Salin_Top = mean(Salin_Top, na.rm=TRUE),
            Salin_Mid = mean(Salin_Mid, na.rm=TRUE),
            Temp_Top = mean( Temp_Top, na.rm=TRUE),
            Temp_Mid = mean(Temp_Mid, na.rm=TRUE),
            winddir_ew = mean(winddir_ew, na.rm=TRUE),
            wind_vector = winddir_ew * biwkmeanspeed_kph) # Added Sept 2019, FJM advice








###################################
# Species Richness

spp_richness <- allcatch %>% filter(Species != "HYCS" & Species != "UNKN" & Station != 231) %>% group_by(Year) %>% 
  summarise(spp_yr = n_distinct(Species)) 

spp_richness <- spp_richness %>% left_join(pru.env.ann %>% dplyr::select(-Station) %>% 
                                             mutate(Year = as.numeric(Year) + 2000) %>%
                                             group_by(Year) %>% summarise_all(mean), 
                                           by = c("Year" = "Year"))


spp_richness.biwk <- allcatch %>% filter(Species != "HYCS" & Species != "UNKN" & Station != 231) %>%
  addbiwknum() %>% group_by(Year, biweekly) %>% summarise(num_spp = n_distinct(Species)) %>% 
  left_join(effort.biwk %>% dplyr::select(biweekly, sample_proportion), 
            by = c("Year" = "Year", "biweekly" = "biweekly")) %>%
  mutate(num_spp_adj = num_spp * sample_proportion)



##### WRITE CSVs #####


write.csv(spp_richness, "data/speciesrichness_annual_2001-2018.csv")
write.csv(spp_richness.biwk, "data/speciesrichness_biwk_2001-2018.csv")



