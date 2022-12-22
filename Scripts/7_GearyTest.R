# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com


############## Load libraries and dataset ##############
library(metafor) # Meta-Analysis Package for R, CRAN v3.4-0
library(tidyverse) # Easily Install and Load the 'Tidyverse', CRAN v1.3.1

# Arthropods Abundance ########
# Read data arthropods ####
ArthO <- read.csv("./Databases/AbundanceDataset_Arth_CompleteOrders_imputed.csv")
 
# Read phylogenetic matrix ####
load("./Results/Arthropods/Phylogeny/arth_phylo_Orders_cor.Rdata") #arth_phyloOrders_cor

# Load necessary functions ####
source("Scripts/00_Functions.R")

# Calculate varcovar matrix with Delta correction ####
ArthV <- bldiag(lapply(split(ArthO, ArthO$Control_ID), calc.vRRDelta))
is.positive.definite(ArthV) # TRUE!

# Assign phylogenetic matrix
phylocor <- list(Order= arth_phyloOrders_cor)

# Random effects
ArthO$OrderID <- ArthO$Order
RandomEffects_Phylo2 <- list(~1|Order, ~1|RowID, ~ 1|Source, ~1|OrderID)

#Fit intercept model and extract parameters for full model
mFull <- rma.mv(logRR, V=ArthV, random=RandomEffects_Phylo2, data=ArthO, R = phylocor)

# Check effect of uncorrected logRR
# Calculate varcovar matrix without Delta correction ####
ArthVuncorrected <- bldiag(lapply(split(ArthO, ArthO$Control_ID), calc.v))
is.positive.definite(ArthVuncorrected) # TRUE!

#Fit intercept model and extract parameters for full model without Delta correction
mFulluncorrected <- rma.mv(logRR_uncorrected, V=ArthVuncorrected, random=RandomEffects_Phylo2, data=ArthO, R = phylocor)

# Check for small sample sizes with Geary's test ####
# Get subset of data for which Geary's diagnostic > 3
ArthO$GearyTreat <- ArthO$Treatment_MeanCorrected / ArthO$Treatment_SDCorrected * ( (4*ArthO$Treatment_N^(3/2)) / (1+4*ArthO$Treatment_N) )
ArthO$GearyControl <- ArthO$Control_MeanCorrected / ArthO$Control_SDCorrected * ( (4*ArthO$Control_N^(3/2)) / (1+4*ArthO$Control_N) )
ArthO$GearyMin <- apply(data.frame(ArthO$GearyTreat,ArthO$GearyControl),1,min)
smallMeans <- ArthO$GearyMin <= 3 | is.na(ArthO$GearyMin)
largeMeans <- ArthO[!smallMeans, ]

# Create new variance covariance matrix for the selected subset
largeMeans$ControlIDsplit <- factor(largeMeans$Control_ID, levels<-unique(as.character(largeMeans$Control_ID)))
largeMeans_V <- metafor::bldiag(lapply(split(largeMeans,largeMeans$ControlIDsplit), calc.vRRDelta))
is.positive.definite(largeMeans_V) # TRUE!

#Fit intercept model and extract parameters for the subset
mGeary = rma.mv(logRR, V=largeMeans_V, random=RandomEffects_Phylo2, data=largeMeans)

# Check results ####
DataArth = data.frame(Type = c(paste('All-Delta (N=',mFull$k,')',sep=''), paste('All (N=',mFulluncorrected$k,')',sep=''), paste('Geary (N=',mGeary$k,')',sep='')),
                  Estimate = c(mFull$beta[1],mFulluncorrected$beta[1],mGeary$beta[1]),
                  CI_upper = c(mFull$ci.ub,mFulluncorrected$ci.ub,mGeary$ci.ub),
                  CI_lower = c(mFull$ci.lb,mFulluncorrected$ci.lb,mGeary$ci.lb))
DataCIArth = data.frame(CI = c(Data$CI_upper,Data$CI_lower),
                    Type = c(Data$Type, Data$Type))

ggplot(DataArth,aes(x=Type,y=Estimate)) +
  geom_point(size=2) +
  geom_line(data = DataCIArth, aes(x=Type,y=CI, group=Type)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", size=0.8) +    
  ylab(expression(LRR^Delta))+
  xlab("")+
  theme_bw()

# Nematodes Abundance ####
Nema <- read.csv("./Databases/AbundanceDataset_Nema_FeedingGuild_imputed.csv")

# Calculate varcovar matrix with Delta correction ####
NemaV <- bldiag(lapply(split(Nema, Nema$Control_ID), calc.vRRDelta))
is.positive.definite(NemaV) # TRUE!


#Fit intercept model and extract parameters for full model
RandomEffectsNema <- list(~1|RowID, ~1|Source)

mFullNema <- rma.mv(logRR, V = NemaV, data = Nema, random = RandomEffectsNema)


NemaVuncorrected <- bldiag(lapply(split(Nema, Nema$Control_ID), calc.v))
is.positive.definite(NemaVuncorrected) # TRUE!

#Fit intercept model and extract parameters for full model without Delta correction
mFullNemauncorrected <- rma.mv(logRR_uncorrected, V=NemaVuncorrected, random=RandomEffectsNema, data=Nema)

# Check for small sample sizes with Geary's test ####
# Get subset of data for which Geary's diagnostic > 3
Nema$GearyTreat <- Nema$Treatment_MeanCorrected / Nema$Treatment_SDCorrected * ( (4*Nema$Treatment_N^(3/2)) / (1+4*Nema$Treatment_N) )
Nema$GearyControl <- Nema$Control_MeanCorrected / Nema$Control_SDCorrected * ( (4*Nema$Control_N^(3/2)) / (1+4*Nema$Control_N) )
Nema$GearyMin <- apply(data.frame(Nema$GearyTreat,Nema$GearyControl),1,min)
smallMeansNema <- Nema$GearyMin <= 3 | is.na(Nema$GearyMin)
largeMeansNema <- Nema[!smallMeansNema, ]

# Create new variance covariance matrix for the selected subset
largeMeansNema$ControlIDsplit <- factor(largeMeansNema$Control_ID, levels<-unique(as.character(largeMeansNema$Control_ID)))
largeMeansNema_V <- metafor::bldiag(lapply(split(largeMeansNema,largeMeansNema$ControlIDsplit), calc.vRRDelta))
is.positive.definite(largeMeansNema_V) # TRUE!

#Fit intercept model and extract parameters for the subset
mGearyNema <- rma.mv(logRR, V=largeMeansNema_V, random=RandomEffectsNema, data=largeMeansNema)

# Check results
DataNema <- data.frame(Type = c(paste('All-Delta (N=',mFullNema$k,')',sep=''), paste('All (N=',mFullNemauncorrected$k,')',sep=''), paste('Geary (N=',mGearyNema$k,')',sep='')),
                  Estimate = c(mFullNema$beta[1],mFullNemauncorrected$beta[1],mGearyNema$beta[1]),
                  CI_upper = c(mFullNema$ci.ub,mFullNemauncorrected$ci.ub,mGearyNema$ci.ub),
                  CI_lower = c(mFullNema$ci.lb,mFullNemauncorrected$ci.lb,mGearyNema$ci.lb))
DataCINema <- data.frame(CI = c(DataNema$CI_upper,DataNema$CI_lower),
                    Type = c(DataNema$Type, DataNema$Type))

ggplot(DataNema,aes(x=Type,y=Estimate)) +
  geom_point(size=2) +
  geom_line(data = DataCINema, aes(x=Type,y=CI, group=Type)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", size=0.8) +    
  ylab(expression(LRR^Delta))+
  xlab("")+
  theme_bw()

# Arthropods Richness ####
# Read data
dat <- read.csv("./Databases/RichnessDataset_ArthNema_imputed_Env.csv")
# Subset for the Arthropods
ArthR <- dat %>%
  dplyr::filter(Phylum %in% "Arthropoda")
# Arrange based on Control ID to calculate properly the var-covar matrix
ArthR <- ArthR %>%
  arrange(Control_ID)

# Create a unique identifier for each effect size, which corresponds to the observational id (i.e. residual variance in the model)
ArthR$RowID <-  1:length(ArthR$Source)


ArthRV <- bldiag(lapply(split(ArthR, ArthR$Control_ID), calc.vRRDelta))
is.positive.definite(ArthV) # TRUE!

# Fit
RandomEffectsArthR <- list(~1|RowID, ~1|Source)

mFullArthR <- rma.mv(logRR, V = ArthRV, data = ArthR, random = RandomEffectsArthR)


# Check effect of uncorrected logRR
# Calculate varcovar matrix without Delta correction ####
ArthRVuncorrected <- bldiag(lapply(split(ArthR, ArthR$Control_ID), calc.v))
is.positive.definite(ArthVuncorrected) # TRUE!

#Fit intercept model and extract parameters for full model without Delta correction
mFullArthRuncorrected <- rma.mv(logRR_uncorrected, V=ArthRVuncorrected, data = ArthR, random=RandomEffectsArthR)

# Check for small sample sizes with Geary's test ####
# Get subset of data for which Geary's diagnostic > 3
ArthR$GearyTreat <- ArthR$Treatment_MeanCorrected / ArthR$Treatment_SDCorrected * ( (4*ArthR$Treatment_N^(3/2)) / (1+4*ArthR$Treatment_N) )
ArthR$GearyControl <- ArthR$Control_MeanCorrected / ArthR$Control_SDCorrected * ( (4*ArthR$Control_N^(3/2)) / (1+4*ArthR$Control_N) )
ArthR$GearyMin <- apply(data.frame(ArthR$GearyTreat,ArthR$GearyControl),1,min)
smallMeansArthR <- ArthR$GearyMin <= 3 | is.na(ArthR$GearyMin)
largeMeansArthR <- ArthR[!smallMeansArthR, ]

# Create new variance covariance matrix for the selected subset
largeMeansArthR$ControlIDsplit <- factor(largeMeansArthR$Control_ID, levels<-unique(as.character(largeMeansArthR$Control_ID)))
largeMeansArthR_V <- metafor::bldiag(lapply(split(largeMeansArthR,largeMeansArthR$ControlIDsplit), calc.vRRDelta))
is.positive.definite(largeMeansArthR_V) # TRUE!

#Fit intercept model and extract parameters for the subset
mGearyArthR = rma.mv(logRR, V=largeMeansArthR_V, random=RandomEffectsArthR, data=largeMeansArthR)

# Check results ####
Data = data.frame(Type = c(paste('All-Delta (N=',mFullArthR$k,')',sep=''), paste('All (N=',mFullArthRuncorrected$k,')',sep=''), paste('Geary (N=',mGearyArthR$k,')',sep='')),
                  Estimate = c(mFullArthR$beta[1],mFullArthRuncorrected$beta[1],mGearyArthR$beta[1]),
                  CI_upper = c(mFullArthR$ci.ub,mFullArthRuncorrected$ci.ub,mGearyArthR$ci.ub),
                  CI_lower = c(mFullArthR$ci.lb,mFullArthRuncorrected$ci.lb,mGearyArthR$ci.lb))
DataCI = data.frame(CI = c(Data$CI_upper,Data$CI_lower),
                    Type = c(Data$Type, Data$Type))

ggplot(Data,aes(x=Type,y=Estimate)) +
  geom_point(size=2) +
  geom_line(data = DataCI, aes(x=Type,y=CI, group=Type)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", size=0.8) +    
  ylab(expression(LRR^Delta))+
  xlab("")+
  theme_bw()

# Nematode Richness ####
# Subset for the Arthropods
NemaR <- dat %>%
  dplyr::filter(Phylum %in% "Nematoda")
# Arrange based on Control ID to calculate properly the var-covar matrix
NemaR <- NemaR %>%
  arrange(Control_ID)

# Create a unique identifier for each effect size, which corresponds to the observational id (i.e. residual variance in the model)
NemaR$RowID <-  1:length(NemaR$Source)


NemaRV <- bldiag(lapply(split(NemaR, NemaR$Control_ID), calc.vRRDelta))
is.positive.definite(ArthV) # TRUE!

# Fit
RandomEffectsNemaR <- list(~1|RowID, ~1|Source)

mFullNemaR <- rma.mv(logRR, V = NemaRV, data = NemaR, random = RandomEffectsNemaR)


# Check effect of uncorrected logRR
# Calculate varcovar matrix without Delta correction ####
NemaRVuncorrected <- bldiag(lapply(split(NemaR, NemaR$Control_ID), calc.v))
is.positive.definite(ArthVuncorrected) # TRUE!

#Fit intercept model and extract parameters for full model without Delta correction
mFullNemaRuncorrected <- rma.mv(logRR_uncorrected, V=NemaRVuncorrected, data = NemaR, random=RandomEffectsNemaR)

# Check for small sample sizes with Geary's test ####
# Get subset of data for which Geary's diagnostic > 3
NemaR$GearyTreat <- NemaR$Treatment_MeanCorrected / NemaR$Treatment_SDCorrected * ( (4*NemaR$Treatment_N^(3/2)) / (1+4*NemaR$Treatment_N) )
NemaR$GearyControl <- NemaR$Control_MeanCorrected / NemaR$Control_SDCorrected * ( (4*NemaR$Control_N^(3/2)) / (1+4*NemaR$Control_N) )
NemaR$GearyMin <- apply(data.frame(NemaR$GearyTreat,NemaR$GearyControl),1,min)
smallMeansNemaR <- NemaR$GearyMin <= 3 | is.na(NemaR$GearyMin)
largeMeansNemaR <- NemaR[!smallMeansNemaR, ]

# Create new variance covariance matrix for the selected subset
largeMeansNemaR$ControlIDsplit <- factor(largeMeansNemaR$Control_ID, levels<-unique(as.character(largeMeansNemaR$Control_ID)))
largeMeansNemaR_V <- metafor::bldiag(lapply(split(largeMeansNemaR,largeMeansNemaR$ControlIDsplit), calc.vRRDelta))
is.positive.definite(largeMeansNemaR_V) # TRUE!

#Fit intercept model and extract parameters for the subset
mGearyNemaR = rma.mv(logRR, V=largeMeansNemaR_V, random=RandomEffectsNemaR, data=largeMeansNemaR)

# Check results ####
Data = data.frame(Type = c(paste('All-Delta (N=',mFullNemaR$k,')',sep=''), paste('All (N=',mFullNemaRuncorrected$k,')',sep=''), paste('Geary (N=',mGearyNemaR$k,')',sep='')),
                  Estimate = c(mFullNemaR$beta[1],mFullNemaRuncorrected$beta[1],mGearyNemaR$beta[1]),
                  CI_upper = c(mFullNemaR$ci.ub,mFullNemaRuncorrected$ci.ub,mGearyNemaR$ci.ub),
                  CI_lower = c(mFullNemaR$ci.lb,mFullNemaRuncorrected$ci.lb,mGearyNemaR$ci.lb))
DataCI = data.frame(CI = c(Data$CI_upper,Data$CI_lower),
                    Type = c(Data$Type, Data$Type))

ggplot(Data,aes(x=Type,y=Estimate)) +
  geom_point(size=2) +
  geom_line(data = DataCI, aes(x=Type,y=CI, group=Type)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", size=0.8) +    
  ylab(expression(LRR^Delta))+
  xlab("")+
  theme_bw()

