# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com


############## Load libraries and dataset ##############

library(raster) # Geographic Data Analysis and Modeling, CRAN v3.5-2
library(terra) # Spatial Data Analysis, CRAN v1.4-20
library(sf) # Simple Features for R, CRAN v1.0-2
library(dplyr) # A Grammar of Data Manipulation, CRAN v1.0.5
library(tictoc) # Functions for timing R scripts, as well as implementations of Stack and List structures., CRAN v1.0

# Abundance ################
########### EXtract CEC data #############

# Read CEC 
cec <- raster("D:/Juan/Desktop/Universidad/PHD/Data/SoilGrids/CEC_mean_1k/CEC_mean_0-30_1k.tif")

# Read the point data
dat <- read.csv("D:/Juan/Desktop/Universidad/PHD/Chapters/Chapter4_InsectsNitrogen/Github/Nimpact_Invertebrates/Databases/AbundanceDataset_ArthNema_imputed.csv")

datsf <- st_as_sf(dat, coords = c("Longitude", "Latitude"), crs = "EPSG:4326", remove = FALSE)
datsf <- st_transform(datsf, crs(cec))
tic()
datsf$CEC <- terra::extract(cec, datsf)
toc()

# Check if any NA in CEC
any(is.na(datsf$CEC))

nas <- datsf %>% 
  dplyr::filter(is.na(CEC))

########### EXtract Nitrogen deposition data #############
# Read Ndep 
ndep <- raster("D:/Juan/Desktop/Universidad/PHD/Data/NitrogenDepositionGlobal/1984-2015/inorganic_tot_Ndep_kghayr_avrg_1984-2015.tif")

tic()
datsf$Ndep_1984_2015 <- terra::extract(ndep, datsf)
toc()

# Check if any NA in CEC
any(is.na(datsf$Ndep_1984_2015))

nas <- datsf %>% 
  dplyr::filter(is.na(Ndep_1984_2015))

########### EXtract Climate Research Unit data #############
# Read MATs per year 
matsfiles <- list.files("D:/Juan/Desktop/Universidad/PHD/Data/ClimateResearchUnit/MeanYears/TMP/", full.names = T)
mats <- stack(matsfiles)

# Read PREs (MAP) per year 
presfiles <- list.files("D:/Juan/Desktop/Universidad/PHD/Data/ClimateResearchUnit/MeanYears/PRE/", full.names = T)
pres <- stack(presfiles)

yearsInitial <- sort(unique(datsf$Study_Year_Initial))
yearsFinal <- sort(unique(datsf$Study_Year_Final))

datsf$CRU_MAT_Initial <- NA
datsf$CRU_MAT_Final <- NA
datsf$CRU_MAP_Initial <- NA
datsf$CRU_MAP_Final <- NA

datsf$UNIQID_Initial <- paste0(datsf$Source,datsf$Study,datsf$Study_Year_Initial)
datsf$UNIQID_Final <- paste0(datsf$Source,datsf$Study,datsf$Study_Year_Final)

tic("Initial year in:")
for(i in 1:length(yearsInitial)){
  dat_i <- datsf %>% filter(Study_Year_Initial %in% yearsInitial[i])
  
  year_i <- grep(yearsInitial[i], names(mats))
  
  mat_i <- mats[[year_i]]
  pre_i <- pres[[year_i]]
  
  dat_i$CRU_MAT_Initial <- terra::extract(mat_i, dat_i)
  dat_i$CRU_MAP_Initial <- terra::extract(pre_i, dat_i)
  
  datsf[datsf$UNIQID_Initial %in% dat_i$UNIQID_Initial, ]$CRU_MAT_Initial <- dat_i$CRU_MAT_Initial
  datsf[datsf$UNIQID_Initial %in% dat_i$UNIQID_Initial, ]$CRU_MAP_Initial <- dat_i$CRU_MAP_Initial
}
toc()


# Final year
tic("Final year in:")
for(i in 1:length(yearsFinal)){
  dat_i <- datsf %>% filter(Study_Year_Final %in% yearsFinal[i])
  
  year_i <- grep(yearsFinal[i], names(mats))
  
  mat_i <- mats[[year_i]]
  pre_i <- pres[[year_i]]
  
  dat_i$CRU_MAT_Final <- terra::extract(mat_i, dat_i)
  dat_i$CRU_MAP_Final <- terra::extract(pre_i, dat_i)
  
  datsf[datsf$UNIQID_Final %in% dat_i$UNIQID_Final, ]$CRU_MAT_Final <- dat_i$CRU_MAT_Final
  datsf[datsf$UNIQID_Final %in% dat_i$UNIQID_Final, ]$CRU_MAP_Final <- dat_i$CRU_MAP_Final
}
toc()

# Mean between years
datsf$CRU_MAT_mean <- (datsf$CRU_MAT_Initial+datsf$CRU_MAT_Final)/2
datsf$CRU_MAP_mean <- (datsf$CRU_MAP_Initial+datsf$CRU_MAP_Final)/2

# Check if any NA in CRU
any(is.na(datsf$CRU_MAT_Final))

nas <- datsf %>% 
  dplyr::filter(is.na(CRU_MAP_Final))

# Save csv
datsf <- as.data.frame(datsf)
datsf$geometry <- NULL
datsf$UNIQID_Initial <- NULL
datsf$UNIQID_Final <- NULL

write.csv(datsf, "D:/Juan/Desktop/Universidad/PHD/Chapters/Chapter4_InsectsNitrogen/Github/Nimpact_Invertebrates/Databases/AbundanceDataset_ArthNema_imputed_Env.csv", row.names = FALSE)


# Richness ################
########### EXtract CEC data #############

# Read the point data
datR <- read.csv("D:/Juan/Desktop/Universidad/PHD/Chapters/Chapter4_InsectsNitrogen/Github/Nimpact_Invertebrates/Databases/RichnessDataset_ArthNema_imputed.csv")

datsf <- st_as_sf(datR, coords = c("Longitude", "Latitude"), crs = "EPSG:4326", remove = FALSE)
datsf <- st_transform(datsf, crs(cec))
tic()
datsf$CEC <- terra::extract(cec, datsf)
toc()

# Check if any NA in CEC
any(is.na(datsf$CEC))

nas <- datsf %>% 
  dplyr::filter(is.na(CEC))

########### EXtract Nitrogen deposition data #############
tic()
datsf$Ndep_1984_2015 <- terra::extract(ndep, datsf)
toc()

# Check if any NA in CEC
any(is.na(datsf$Ndep_1984_2015))

nas <- datsf %>% 
  dplyr::filter(is.na(Ndep_1984_2015))

########### EXtract Climate Research Unit data #############
yearsInitial <- sort(unique(datsf$Study_Year_Initial))
yearsFinal <- sort(unique(datsf$Study_Year_Final))

datsf$CRU_MAT_Initial <- NA
datsf$CRU_MAT_Final <- NA
datsf$CRU_MAP_Initial <- NA
datsf$CRU_MAP_Final <- NA

datsf$UNIQID_Initial <- paste0(datsf$Source,datsf$Study,datsf$Study_Year_Initial)
datsf$UNIQID_Final <- paste0(datsf$Source,datsf$Study,datsf$Study_Year_Final)

tic("Initial year in:")
for(i in 1:length(yearsInitial)){
  dat_i <- datsf %>% filter(Study_Year_Initial %in% yearsInitial[i])
  
  year_i <- grep(yearsInitial[i], names(mats))
  
  mat_i <- mats[[year_i]]
  pre_i <- pres[[year_i]]
  
  dat_i$CRU_MAT_Initial <- terra::extract(mat_i, dat_i)
  dat_i$CRU_MAP_Initial <- terra::extract(pre_i, dat_i)
  
  datsf[datsf$UNIQID_Initial %in% dat_i$UNIQID_Initial, ]$CRU_MAT_Initial <- dat_i$CRU_MAT_Initial
  datsf[datsf$UNIQID_Initial %in% dat_i$UNIQID_Initial, ]$CRU_MAP_Initial <- dat_i$CRU_MAP_Initial
}
toc()


# Final year
tic("Final year in:")
for(i in 1:length(yearsFinal)){
  dat_i <- datsf %>% filter(Study_Year_Final %in% yearsFinal[i])
  
  year_i <- grep(yearsFinal[i], names(mats))
  
  mat_i <- mats[[year_i]]
  pre_i <- pres[[year_i]]
  
  dat_i$CRU_MAT_Final <- terra::extract(mat_i, dat_i)
  dat_i$CRU_MAP_Final <- terra::extract(pre_i, dat_i)
  
  datsf[datsf$UNIQID_Final %in% dat_i$UNIQID_Final, ]$CRU_MAT_Final <- dat_i$CRU_MAT_Final
  datsf[datsf$UNIQID_Final %in% dat_i$UNIQID_Final, ]$CRU_MAP_Final <- dat_i$CRU_MAP_Final
}
toc()

# Mean between years
datsf$CRU_MAT_mean <- (datsf$CRU_MAT_Initial+datsf$CRU_MAT_Final)/2
datsf$CRU_MAP_mean <- (datsf$CRU_MAP_Initial+datsf$CRU_MAP_Final)/2

# Check if any NA in CRU
any(is.na(datsf$CRU_MAT_Final))

nas <- datsf %>% 
  dplyr::filter(is.na(CRU_MAP_Final))

# Save csv
datsf <- as.data.frame(datsf)
datsf$geometry <- NULL
datsf$UNIQID_Initial <- NULL
datsf$UNIQID_Final <- NULL

write.csv(datsf, "D:/Juan/Desktop/Universidad/PHD/Chapters/Chapter4_InsectsNitrogen/Github/Nimpact_Invertebrates/Databases/RichnessDataset_ArthNema_imputed_Env.csv", row.names = FALSE)


