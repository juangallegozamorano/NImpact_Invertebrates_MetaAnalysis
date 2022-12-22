# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com


############## Load libraries and dataset ##############
library(metafor) # Meta-Analysis Package for R, CRAN v3.4-0
library(tidyverse) # Easily Install and Load the 'Tidyverse', CRAN v1.3.1
library(ggpubr) # 'ggplot2' Based Publication Ready Plots, CRAN v0.4.0
library(tictoc) # Functions for timing R scripts, as well as implementations of Stack and List structures., CRAN v1.0
library(plyr) # Tools for Splitting, Applying and Combining Data, CRAN v1.8.6
library(orchaRd) # Orchard Plots and Prediction Intervals for Meta-analysis, [github::itchyshin/orchard_plot] v0.0.0.9000

# Read data
dat <- read.csv("./Databases/RichnessDataset_ArthNema_imputed_Env.csv")

# Subset for the Arthropods
ArthO <- dat %>%
  dplyr::filter(Phylum %in% "Arthropoda")

# load necessary functions
source("Scripts/00_Functions.R")

# Arrange based on Control ID to calculate properly the var-covar matrix
ArthO <- ArthO %>%
  arrange(Control_ID)

# Create a unique identifier for each effect size, which corresponds to the observational id (i.e. residual variance in the model)
ArthO$RowID <-  1:length(ArthO$Source)

# Calculate var-covar matrix
calc.vRRDelta <- function(x) {
  v = matrix((x$Control_SDCorrected[1]^2 / (x$Control_N[1] * x$Control_MeanCorrected[1]^2) + 
                0.5*(x$Control_SDCorrected[1]^4/(x$Control_N[1]^2*x$Control_MeanCorrected[1]^4))), nrow = nrow(x), ncol = nrow(x))
  diag(v) = x$VAR
  return(v)
} 

ArthV <- bldiag(lapply(split(ArthO, ArthO$Control_ID), calc.vRRDelta))
is.positive.definite(ArthV) # TRUE!

# Check if variances are ok, if it's TRUE then is BAD
summary(ArthO$VAR)
vimaxmin <- max(ArthO$VAR) / min(ArthO$VAR) 
vimaxmin >= 1e7

# Transform continous variables
ArthO$logNadd <- log10(ArthO$Nitrogen_Added)
ArthO$logNadd2 <- ArthO$logNadd^2

ArthO$logNdep <- log10(ArthO$Ndep_1984_2015)
ArthO$logNdep2 <- ArthO$logNdep^2

ArthO$logDuration_Years <- log10(ArthO$Duration_Years)
ArthO$logCEC <- log10(ArthO$CEC)

# Merge Feeding guilds
table(ArthO$Feeding_Guild)
ArthO <- ArthO %>% 
  mutate(Feeding_GuildUsed = as.factor(case_when(Feeding_Guild %in% c("Parasite", "Parasitoid", "Plant parasite")~ "Parasite",
                                       Feeding_Guild %in% c("Nectar feeders", "Herbivore/ Fungivore", "Herbivore", "Fungivore")~ "Herbivore",
                                       Feeding_Guild %in% c("Detrivore")~ "Detritivore",
                                       Feeding_Guild %in% c("Omnivore", "Predator")~ "Predator-Omnivore",
                                       TRUE ~ Feeding_Guild)),
         Feeding_GuildCode = as.numeric(Feeding_GuildUsed))

table(ArthO$Feeding_GuildUsed)

# Create habitat
ArthO <- ArthO %>% 
  mutate(Habitat = case_when(
    str_detect(Vegetation_Type, "(?i)grass|prair|meadow|Lupinus|scrub|shrub|Bilberry|heath|marsh|Spartina|Milkweed|Peatland|Desertic") ~ "Non-forest",
    str_detect(Vegetation_Type, "(?i)forest|tree|pine|fern") ~ "Forest", 
    str_detect(Vegetation_Type, "(?i)plantat|orchard|tea|crop|agricul|culti|field|wheat|oat|Brussel|Cotton|Sugar Beet|Maize|Switchgrass|Miscanthus|Rosa|Soy|Sorghum bicolor|Cucumber|peas") ~ "Cultivated",
    TRUE ~ Vegetation_Type))

table(ArthO$Habitat)

# load Fixed effects
FixedEffects <- readRDS('Scripts/FixedEffects_Arthropods_Richness.rds')

length(FixedEffects)

############## Running models ##############
#Three levels of variation are:
# Level 1- effect size
# Level 2- within-study variation i.e RowID
# Level 3- between-study variation i.e. Source

RandomEffects <- list(~1|RowID, ~1|Source)
modelsArthR <- vector(mode = "list", length = length(FixedEffects))

for(i in 1:length(FixedEffects)){
  try({
  tic(paste("Model",i,"formula-",as.character(FixedEffects[[i]])[2]))

  modeli <- rma.mv(logRR,
                   mods = FixedEffects[[i]],
                   V = ArthV, 
                   data = ArthO, 
                   random = RandomEffects)

  toc()
  
  modelsArthR[[i]] <- modeli
  })
}

############## Extracting statictis from models ##############
# Creating empty data frame
I2df <- data.frame(ModelNr = rep(NA,length(FixedEffects)),
                   
                   FixedEffects = rep(NA,length(FixedEffects)),
                   RandomEffects = rep(NA,length(FixedEffects)),
                   
                   I2_total = rep(NA,length(FixedEffects)),
                   I2_RowID = rep(NA,length(FixedEffects)),
                   I2_Source = rep(NA,length(FixedEffects)),
                   
                   AIC = rep(NA,length(FixedEffects)),
                   BIC = rep(NA,length(FixedEffects)),
                   AICc = rep(NA,length(FixedEffects)),
                   mR2 = rep(NA,length(FixedEffects)),
                   cR2 = rep(NA,length(FixedEffects)))

for(i in 1:length(modelsArthR)){
  modeli <- modelsArthR[[i]]
  if(is.null(modeli)){
    I2df[i,]$ModelNr <- i
    I2df[i,]$FixedEffects <- paste0("logRR",as.character(FixedEffects[i]))
    next
  }
  I2_modeli <- i2_ml(modeli)
  
  I2df[i,]$ModelNr <- i
  I2df[i,]$FixedEffects <- paste0("logRR~",as.character(modeli$formula.mods)[2])
  I2df[i,]$RandomEffects <- "~1|RowID, ~1|Source"
  
  I2df[i,]$I2_total <- I2_modeli[1]
  I2df[i,]$I2_RowID <- I2_modeli[2]
  I2df[i,]$I2_Source <- I2_modeli[3]

  
  I2df[i,]$AIC <- modeli$fit.stats["AIC", "ML"]
  I2df[i,]$BIC <- modeli$fit.stats["BIC", "ML"]
  I2df[i,]$AICc <- modeli$fit.stats["AICc", "ML"]
  I2df[i,]$mR2 <- mR2.func(modeli)
  I2df[i,]$cR2 <- cR2.func(modeli)
}

I2df$I2_total <- round(I2df$I2_total, 3)
I2df$I2_RowID <- round(I2df$I2_RowID, 3)
I2df$I2_Source <- round(I2df$I2_Source, 3)
View(I2df) # View best model

write.csv(I2df, file = "./Results/Arthropods/Richness_Arthropods_model_selection_160.csv", row.names = F)

# Save models
save(modelsArthR, file = "./Results/Arthropods/Richness_modelsArth_160.RData")

##############  Checks of best model ############## 
# Check profile plot of best model
modeli <- modelsArthR[[1]] # Best model is the NULL model

profile(modeli)
