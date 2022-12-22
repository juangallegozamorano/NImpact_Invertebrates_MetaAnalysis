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

# loading phylogenetic matrix
load("./Results/Arthropods/Phylogeny/arth_phylo_Orders_cor.Rdata") #arth_phyloOrders_cor

# load necessary functions
source("Scripts/00_Functions.R")

############## Preparing dataset ##############
#Subset for the Arthropods with Order
ArthO <- dat %>% 
  dplyr::filter(Phylum %in% "Arthropoda") %>% 
  drop_na(Order)

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

# Check if variances are ok
summary(ArthO$VAR)
vimaxmin <- max(ArthO$VAR) / min(ArthO$VAR) 
vimaxmin >= 1e7 # It has to be false, otherwise the models won't run

# Transform continuous variables that are skewed
ArthO$logNadd <- log10(ArthO$Nitrogen_Added)
ArthO$logNadd2 <- ArthO$logNadd^2

ArthO$logNdep <- log10(ArthO$Ndep_1984_2015)
ArthO$logNdep2 <- ArthO$logNdep^2

# Merge Feeding guilds
table(ArthO$Feeding_Guild)
ArthO <- ArthO %>% 
  mutate(Feeding_GuildUsed = as.factor(case_when(Feeding_Guild %in% c("Parasite", "Parasitoid", "Plant-parasite")~ "Parasite",
                                       Feeding_Guild %in% c("Nectar feeders", "Herbivore/Fungivore","Fungi", "Herbivore", "Fungivore")~ "Herbivore-Fungivore",
                                       Feeding_Guild %in% c("Detrivore")~ "Detritivore",
                                       Feeding_Guild %in% c("Omnivore", "Predator")~ "Predator-Omnivore",
                                       TRUE ~ Feeding_Guild)),
         Feeding_GuildCode = as.numeric(Feeding_GuildUsed))

table(ArthO$Feeding_GuildUsed)

# Create habitat categories
ArthO <- ArthO %>% 
  mutate(Habitat = case_when(str_detect(Vegetation_Type, "(?i)grass|prair|meadow|Lupinus|Savan|scrub|shrub|Bilberry|heath|Milkweed|marsh|Spartina|Peatland|urban") ~ "Non-Forest",
                             str_detect(Vegetation_Type, "(?i)forest|tree|pine|fern|Old_second") ~ "Forest", 
                             str_detect(Vegetation_Type, "(?i)plantat|orchard|tea|crop|agricul|culti|field|wheat|Brussel|Cotton|Sugar Beet|Maize|Switchgrass|Miscanthus|Rosa|Soy") ~ "Cultivated",
                             TRUE ~ Vegetation_Type))
table(ArthO$Habitat)

# Create metamorphosis categories
ArthO <- ArthO %>% 
  mutate(Metamorphosis = case_when(str_detect(Order, "(?i)Coleoptera|Diptera|Hymenoptera|Lepidoptera|Mecoptera|Neuroptera|Raphidioptera") ~ "Complete",
                                   str_detect(Order, "(?i)Araneae|Blattodea|Dermaptera|Diplura|Geophilomorpha|Lithobiomorpha|Embioptera|Ephemeroptera|Hemiptera|Isopoda|Mesostigmata|Opiliones|Odonata|Orthoptera|Phasmatodea|Pseudoscorpiones|Psocoptera|Scolopendromorpha|Sarcoptiformes|Tetramerocerata|Thysanoptera|Trombidiformes|Polyxenida|Julida") ~ "Incomplete-Gradual",
                                   str_detect(Order, "(?i)Entomobryomorpha|Neelipleona|Poduromorpha|Symphypleona|Zygentoma|Protura") ~ "None",
                                   TRUE ~ Order))
table(ArthO$Metamorphosis)
ArthO$Metamorphosis <- as.factor(ArthO$Metamorphosis)

# Save database
write.csv(ArthO, "./Databases/AbundanceDataset_Arth_CompleteOrders_imputed.csv")

# load Fixed effects
FixedEffects <- readRDS('Scripts/FixedEffects_Arthropods_Abundance.rds')

length(FixedEffects)

############## Running models ##############
# Five levels of variation are:
# Level 1- effect size
# Level 2- Phylogenetic non-independece i.e. Order with phylogenetic matrix
# Level 3- within-study variation i.e RowID
# Level 4- between-study variation i.e. Source
# Level 5- Non-phylogenetic non-independece i.e. Order ID

# Create OrderID for Non-phylogenetic non-independece
ArthO$OrderID <- ArthO$Order

# Assign phylogenetic matrix to the Orders
phylocor <- list(Order= arth_phyloOrders_cor)

# Random effects
RandomEffects_Phylo <- list(~1|Order, ~1|RowID, ~ 1|Source, ~1|OrderID)

# Create a list to fill with all models
modelsArth_Phylo <- vector(mode = "list", length = length(FixedEffects))

for(i in 1:length(FixedEffects)){

  tic(paste("Phylo Model",i,"formula-",as.character(FixedEffects[[i]])[2]))
  modeli <- rma.mv(logRR,
                   mods = FixedEffects[[i]],
                   V = ArthV,
                   data = ArthO,
                   random = RandomEffects_Phylo,
                   R = phylocor)
  toc()

  modelsArth_Phylo[[i]] <- modeli
}


############## Extracting statictis from models ##############

# Create empty dataframe
I2df <- data.frame(ModelNr = rep(NA,length(FixedEffects)),
                   
                   FixedEffects = rep(NA,length(FixedEffects)),
                   RandomEffects = rep(NA,length(FixedEffects)),
                   
                   I2_total = rep(NA,length(FixedEffects)),
                   I2_RowID = rep(NA,length(FixedEffects)),
                   I2_Source = rep(NA,length(FixedEffects)),
                   I2_Order = rep(NA,length(FixedEffects)), 
                   I2_OrderID = rep(NA,length(FixedEffects)),
                   
                   PhyloYN = rep(NA,length(FixedEffects)),
                   AIC = rep(NA,length(FixedEffects)),
                   BIC = rep(NA,length(FixedEffects)),
                   AICc = rep(NA,length(FixedEffects)),
                   mR2 = rep(NA,length(FixedEffects)),
                   cR2 = rep(NA,length(FixedEffects)))

# Fill data frame with statictis of each model
for(i in 1:length(modelsArth_Phylo)){
  
  modeli <- modelsArth_Phylo[[i]]
  
  # Check if the model converged and have parameters
  if(is.null(modeli)){
    I2df[i,]$ModelNr <- i
    I2df[i,]$FixedEffects <- paste0("logRR",as.character(FixedEffects[i]))
    next
  }
  # Get all I2 from the model
  I2_modeli <- i2_ml(modeli)
  
  # Number of the model and fixed effects
  I2df[i,]$ModelNr <- i
  I2df[i,]$FixedEffects <- paste0("logRR~",as.character(modeli$formula.mods)[2])
  I2df[i,]$RandomEffects <- "~1|Order, ~1|RowID, ~1|Source, ~1|OrderID"
  
  # Extract the different I2 for each random effect
  I2df[i,]$I2_total <- I2_modeli[1]
  I2df[i,]$I2_Order <- I2_modeli[2]
  I2df[i,]$I2_RowID <- I2_modeli[3]
  I2df[i,]$I2_Source <- I2_modeli[4]
  I2df[i,]$I2_OrderID <- I2_modeli[5]
  

  I2df[i,]$PhyloYN <- "Yes"
  
  # Extract ICs and R2
  I2df[i,]$AIC <- modeli$fit.stats["AIC", "ML"]
  I2df[i,]$BIC <- modeli$fit.stats["BIC", "ML"]
  I2df[i,]$AICc <- modeli$fit.stats["AICc", "ML"]
  I2df[i,]$mR2 <- mR2.func(modeli)
  I2df[i,]$cR2 <- cR2.func_Phylo(modeli)
  
}

# Round and save the dataframe
I2df$I2_total <- round(I2df$I2_total, 2)
I2df$I2_Order <- round(I2df$I2_Order, 2)
I2df$I2_OrderID <- round(I2df$I2_OrderID, 2)
I2df$I2_RowID <- round(I2df$I2_RowID, 2)
I2df$I2_Source <- round(I2df$I2_Source, 2)
I2df$AICc <- round(I2df$AICc, 2)
I2df$BIC <- round(I2df$BIC, 2)
I2df$AIC <- round(I2df$AIC, 2)
View(I2df) # View to check the best model

write.csv(I2df, file = "./Results/Arthropods/Abundance_Arthropods_model_selection512_24102022.csv", row.names = F)

# Select best model and save
modeliMetamorMAT <- modelsArth_Phylo[[288]]
save(modeliMetamorMAT, file = "./Results/Arthropods/Abundance_Arthropods_MetamorphosisMAT_512phylo2.RData")

# Save all models
save(modelsArth_Phylo, file = "./Results/Arthropods/Abundance_modelsArth_OrderPhylo_OrderID512.RData")

write(modelsArth_Phylo, file = "./Results/Arthropods/Abundance_modelsArth_OrderPhylo_OrderID512.RData")
save(modelsArth_Phylo, file = "./Results/Arthropods/Abundance_modelsArth_OrderPhylo_OrderID512.RData")

##############  Checks of best model ############## 
# Check profile plot of best model
profile(modeliMetamorMAT)

# Check residuals of best model
resid <- residuals(modeliMetamorMAT) 
fitted <- fitted(modeliMetamorMAT) 

checkdf <- as.data.frame(cbind(resid,fitted))

ggplot(data = checkdf, aes(x= fitted, y = resid))+
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()

# Using Cook's distance (takes long time)
tic("Cooks distance in:")
cooks <- cooks.distance.rma.mv(modeliMetamorMAT, parallel = "multicore", ncpus = 8)
toc()

png(filename = "cooks_distance_plot.png",width = 800, height = 640, pointsize = 12, res = 120)
plot(cooks, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance")
dev.off()


