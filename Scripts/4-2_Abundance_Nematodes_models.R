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
dat <- read.csv("./Databases/AbundanceDataset_ArthNema_imputed_Env.csv")

#Subset for the Nematoda
Nema <- dat %>% 
  dplyr::filter(Phylum %in% "Nematoda")

# load necessary functions
source("Scripts/00_Functions.R")

# Arrange based on Control ID to calculate properly the var-covar matrix
Nema <- Nema %>%
  arrange(Control_ID)

# Create a unique identifier for each effect size, which corresponds to the observational id (i.e. residual variance in the model)
Nema$RowID <-  1:length(Nema$Source)

# Calculate var-covar matrix with Delta method
calc.vRRDelta <- function(x) {
  v = matrix((x$Control_SDCorrected[1]^2 / (x$Control_N[1] * x$Control_MeanCorrected[1]^2) + 
                0.5*(x$Control_SDCorrected[1]^4/(x$Control_N[1]^2*x$Control_MeanCorrected[1]^4))), nrow = nrow(x), ncol = nrow(x))
  diag(v) = x$VAR
  return(v)
} 


NemaV <- bldiag(lapply(split(Nema, Nema$Control_ID), calc.vRRDelta))
is.positive.definite(NemaV) # TRUE!

# Check if variances are ok, if it's TRUE then is BAD
summary(Nema$VAR)
vimaxmin <- max(Nema$VAR) / min(Nema$VAR) 
vimaxmin >= 1e7 # It has to be false, otherwise the models won't run

# Transform continuous variables that are skewed
Nema$logNadd <- log10(Nema$Nitrogen_Added)
Nema$logNadd2 <- Nema$logNadd^2

Nema$logNdep <- log10(Nema$Ndep_1984_2015)
Nema$logNdep2 <- Nema$logNdep^2

# Merge Feeding guilds
table(Nema$Feeding_Guild)
Nema <- Nema %>% 
  mutate(Feeding_GuildUsed = as.factor(case_when(str_detect(Feeding_Guild, "(?i)Preda|Omniv|Free|Microbi")~ "Predator-Omnivore",
                                                 str_detect(Feeding_Guild, "(?i)Parasit")~ "Parasite",
                                       TRUE ~ Feeding_Guild)),
         Feeding_GuildCode = as.numeric(Feeding_GuildUsed))

table(Nema$Feeding_GuildUsed)
table(Nema$Feeding_GuildCode)

# Create habitat
Nema <- Nema %>% 
  mutate(Habitat = case_when(
                             str_detect(Vegetation_Type, "(?i)grass|prair|meadow|Lupinus|scrub|shrub|Bilberry|heath|marsh|Spartina|Milkweed|Peatland|Desertic") ~ "Non-forest",
                             str_detect(Vegetation_Type, "(?i)forest|tree|pine|fern") ~ "Forest", 
                             str_detect(Vegetation_Type, "(?i)plantat|orchard|tea|crop|agricul|culti|field|wheat|oat|Brussel|Cotton|Sugar Beet|Maize|Switchgrass|Miscanthus|Rosa|Soy|Sorghum bicolor|Cucumber|peas") ~ "Cultivated",

                             TRUE ~ Vegetation_Type))
table(Nema$Habitat)

write.csv(Nema, "./Databases/AbundanceDataset_Nema_FeedingGuild_imputed.csv")

# load Fixed effects
FixedEffects <- readRDS('Scripts/FixedEffects_Nematodes_Abundance.rds')

length(FixedEffects)

############## Running models ##############
#Three levels of variation are:
# Level 1- effect size
# Level 2- within-study variation i.e RowID
# Level 3- between-study variation i.e. Source

RandomEffects <- list(~1|RowID, ~1|Source)
modelsNema <- vector(mode = "list", length = length(FixedEffects))

for(i in 1:length(FixedEffects)){
  
  tic(paste("Model",i,"formula-",as.character(FixedEffects[[i]])[2]))
  modeli <- rma.mv(logRR,
                   mods = FixedEffects[[i]],
                   V = NemaV, 
                   data = Nema, 
                   random = RandomEffects)
  toc()
  
  modelsNema[[i]] <- modeli
}


############## Extracting statictis from models ##############

# Create empty dataframe
I2dfNema <- data.frame(ModelNr = rep(NA,length(FixedEffects)),
                   
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

# Fill data frame with statictis of each model
for(i in 1:length(modelsNema)){
  modeli <- modelsNema[[i]]
  
  # Check if the model converged and have parameters
  if(is.null(modeli)){
    I2dfNema[i,]$ModelNr <- i
    I2dfNema[i,]$FixedEffects <- paste0("logRR",as.character(FixedEffects[i]))
    next
  }
  
  # Get all I2 from the model
  I2_modeli <- i2_ml(modeli)
  
  # Number of the model and fixed effects
  I2dfNema[i,]$ModelNr <- i
  I2dfNema[i,]$FixedEffects <- paste0("logRR~",as.character(modeli$formula.mods)[2])
  I2dfNema[i,]$RandomEffects <- "~1|RowID, ~1|Source"
  
  # Extract the different I2 for each random effect
  I2dfNema[i,]$I2_total <- I2_modeli[1]
  I2dfNema[i,]$I2_RowID <- I2_modeli[2]
  I2dfNema[i,]$I2_Source <- I2_modeli[3]
  
  # Extract ICs and R2
  I2dfNema[i,]$AIC <- modeli$fit.stats["AIC", "ML"]
  I2dfNema[i,]$BIC <- modeli$fit.stats["BIC", "ML"]
  I2dfNema[i,]$AICc <- modeli$fit.stats["AICc", "ML"]
  I2dfNema[i,]$mR2 <- mR2.func(modeli)
  I2dfNema[i,]$cR2 <- cR2.func(modeli)
}

# Round and save the dataframe
I2dfNema$I2_total <- round(I2dfNema$I2_total, 3)
I2dfNema$I2_RowID <- round(I2dfNema$I2_RowID, 3)
I2dfNema$I2_Source <- round(I2dfNema$I2_Source, 3)
View(I2dfNema)

write.csv(I2dfNema, file = "./Results/Nematode/Abundance_model_selection_256.csv", row.names = F)

# Save models
save(modelsNema, file = "./Results/Nematode/Abundance_modelsNema_256.RData")

# Save best model
modeliMAP <- modelsNema[[24]]
save(modeliMAP, file = "./Results/Nematode/Abundance_Bestmodel_MAP_FeedingGuild.RData")

##############  Checks of best model ############## 
# Check profile plot of best model
profile(modeliMAP)

# Check residuals of best model
resid <- residuals(modeliMAP) 
fitted <- fitted(modeliMAP) 

checkdf <- as.data.frame(cbind(resid,fitted))

ggplot(data = checkdf, aes(x= fitted, y = resid))+
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()

# Using Cook's distance (takes long time)
tic("Cooks distance in:")
cooks <- cooks.distance.rma.mv(modeliMAP, parallel = "multicore", ncpus = 8)
toc()

png(filename = "cooks_distance_plot.png",width = 800, height = 640, pointsize = 12, res = 120)
plot(cooks, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance")
dev.off()

