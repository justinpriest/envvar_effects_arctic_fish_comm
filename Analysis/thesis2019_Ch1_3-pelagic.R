# PELAGIC VS BENTHIC 
# shelving this for now because it doesn't seem interesting


source("Analysis/thesis2019_Ch1_1-import&cleanup.R")



marinespp <- c("ARCD", "ARFL", "CAPE", "FHSC", "ELPT", "KPGL", "KPSF", "LIPA", "PCHG", "RBSM", 
               "RKGL", "SFCD", "SHSC", "SLEB", "WOLF", "WSGL")
freshspp <- c("ARCH", "BRBT", "GRAY", "LNSK", "SLSC")
anadspp <- c("ARCS", "BDWF", "BRCS", "CHUM", "DLVN", "HBWF", "HYCS", "LSCS", "NNSB", "PINK", "RDWF", "SOCK", "THSB")

spplookup <- spplookup %>% mutate(marine_fresh = if_else(Species %in% marinespp, "marine", 
                                                         if_else(Species %in% freshspp, "freshwater", 
                                                                 if_else(Species %in% anadspp, "anadromous", "NA"))))


spplookup <- spplookup %>% mutate(pelagic_benth = if_else(familygroup == "Cottidae" | familygroup == "Pleuronectidae" | 
                                                          familygroup == "Other", "benthic", "pelagic"), 
                     pelagic_benth = replace(pelagic_benth, Species == "NNSB" | Species == "THSB" | Species == "KPGL" |
                                             Species == "RKGL" | Species == "WSGL", "pelagic")) 
  

temp <- left_join(allcatch, spplookup %>% dplyr::select(-commonname, - familygroup), 
          by = c("Species" = "Species")) %>% 
  group_by(Year, pelagic_benth) %>% summarise(anncatch = sum(totcount))


ggplot(data = temp, aes(x=Year, y = anncatch, color = pelagic_benth)) + geom_line()
