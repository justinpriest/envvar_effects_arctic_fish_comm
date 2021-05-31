##### FUNCTIONS #####

# This function adds the day of year column (if not already present), and then the weeknumber
addweeknum <- function(dataframename){
  dataframename <- dataframename %>% 
    mutate(day.of.year = yday(EndDate),
           week = ifelse(day.of.year <= 187, 1, 
                         ifelse(day.of.year > 187 & day.of.year <= 194, 2, 
                                ifelse(day.of.year > 194 & day.of.year <= 201, 3, 
                                       ifelse(day.of.year > 201 & day.of.year <= 208, 4, 
                                              ifelse(day.of.year > 208 & day.of.year <= 215, 5,
                                                     ifelse(day.of.year > 215 & day.of.year <= 222, 6,
                                                            ifelse(day.of.year > 222 & day.of.year <= 229, 7,
                                                                   ifelse(day.of.year > 229 & day.of.year <= 236, 8,
                                                                          ifelse(day.of.year > 236, 9, NA))))))))))
  
}

# Similarly, this function adds the day of year column, then the biweek number
addbiwknum <- function(dataframename){
  dataframename <- dataframename %>% 
    mutate(day.of.year = yday(EndDate),
           biweekly = ifelse(day.of.year <= 196, 1, # July 15 and before 
                             ifelse(day.of.year > 196 & day.of.year <= 212, 2, #btwn July 16 and 31
                                    ifelse(day.of.year > 212 & day.of.year <= 227, 3, #Aug 1 - 15
                                           ifelse(day.of.year > 227, 4, NA))))) # Aug 16 and after 
}

