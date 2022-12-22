# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com


############## Load libraries and dataset ##############

library(metafor)   # Meta-Analysis Package for R CRAN v3.4-0
library(metagear)  # Comprehensive Research Synthesis Tools for Systematic Reviews and Meta-Analysis CRAN v0.7
library(tidyverse) # Easily Install and Load the 'Tidyverse' CRAN v1.3.1

# Abundances ####
dat <- read.csv("D:/Juan/Desktop/Universidad/PHD/Chapters/Chapter2_InsectsNitrogen/Github/Nimpact_Invertebrates/Databases/AbundanceDataset_csv.csv")

# Remove studies with Deposition gradient
dat <- dat %>% 
  dplyr::filter(N_Added_Deposited == "Added")

# Filter to get only the treatments and not the controls
dat <- dat %>% 
  dplyr::group_by(Source, Study, Control_ID) %>% 
  dplyr::mutate(ControlYN = case_when(Nitrogen_Added == min(Nitrogen_Added) ~ "Control",
                                      TRUE ~ "Treatment"))
datT <- dat %>% 
  dplyr::filter(ControlYN == "Treatment") %>% 
  dplyr::ungroup()

# Remove observations where the Control AND the Treatment are 0
datT <- datT %>% 
  filter(!(Control_Mean == 0 & Treatment_Mean == 0))

# Remove studies where the control is not 0
ControlNoZero <- dat %>% 
  dplyr::filter(ControlYN == "Control") %>% 
  dplyr::filter(Nitrogen_Added != 0)

ControlNoZero <- unique(ControlNoZero$Study)

datT <- datT %>% 
  dplyr::filter(!(Study %in% ControlNoZero))

# Some summary numbers per phylum
table(datT$Phylum)

# Correct 0s in the Mean
datT <- datT %>% 
  dplyr::mutate(Control_MeanCorrected = case_when(Control_Mean == 0 ~ 1/(2 * Control_N * Control_ScalingConstant),
                                                  TRUE ~ Control_Mean),
                Treatment_MeanCorrected = case_when(Treatment_Mean == 0 ~ 1/(2 * Treatment_N * Treatment_ScalingConstant),
                                                    TRUE ~ Treatment_Mean))

# Impute SDs Check SDs that are 0 and convert to NA, better to impute all 0s because it's hard to believe that they are not 0
datT <- datT %>% 
  dplyr::mutate(Control_SDCorrected = case_when(Control_SD == 0 ~ NA_real_,
                                         TRUE ~ Control_SD),
         Treatment_SDCorrected = case_when(Treatment_SD == 0 ~ NA_real_,
                                           TRUE ~ Treatment_SD))

# Subset per phylum
Arth <- datT %>% 
  dplyr::filter(Phylum %in% "Arthropoda") %>% 
  as.data.frame()

Nema <- datT %>% 
  dplyr::filter(Phylum %in% "Nematoda")%>% 
  as.data.frame()

# Impute missing SDs separately
impute_missingness(Arth)
Arth_imp <- impute_SD(aDataFrame = Arth, columnSDnames = "Control_SDCorrected", columnXnames = "Control_MeanCorrected", method = "Bracken1992")

Arth_imp <- impute_SD(aDataFrame = Arth_imp, columnSDnames = "Treatment_SDCorrected", columnXnames = "Treatment_MeanCorrected", method = "Bracken1992")
impute_missingness(Arth_imp)

# Check if there any missing NAs
nas <- Arth_imp %>% 
  dplyr::filter(is.na(Control_SDCorrected)) %>% 
  dplyr::select(Source, contains("Treatment"), contains("Control"))

impute_missingness(Nema)
Nema_imp <- impute_SD(aDataFrame = Nema, columnSDnames = "Control_SDCorrected", columnXnames = "Control_MeanCorrected", method = "Bracken1992")

Nema_imp <- impute_SD(aDataFrame = Nema_imp, columnSDnames = "Treatment_SDCorrected", columnXnames = "Treatment_MeanCorrected", method = "Bracken1992")
impute_missingness(Nema_imp)

nas <- Nema_imp %>% 
  dplyr::filter(is.na(Control_SDCorrected) | is.na(Treatment_SDCorrected)) %>% 
  dplyr::select(Source, contains("Treatment"), contains("Control"))

# Join databases again
dat <- rbind(Arth_imp, Nema_imp)

# Calculate logRR with bias correction (Delta method according to Lajeunesse 2015)
dat$Ratio <- dat$Treatment_MeanCorrected/dat$Control_MeanCorrected
dat$logRR_uncorrected <- log(dat$Ratio)
dat$logRR <- dat$logRR_uncorrected + 
  0.5*(((dat$Treatment_SDCorrected^2)/(dat$Treatment_N*(dat$Treatment_MeanCorrected^2))) - ((dat$Control_SDCorrected^2)/(dat$Control_N*(dat$Control_MeanCorrected^2))))


dat$VAR_uncorrected <- (as.numeric(dat$Treatment_SDCorrected)^2)/(dat$Treatment_N*dat$Treatment_MeanCorrected^2) + (as.numeric(dat$Control_SDCorrected)^2)/(dat$Control_N*dat$Control_MeanCorrected^2)
dat$VAR <- dat$VAR_uncorrected + 0.5*(((dat$Treatment_SDCorrected^4)/((dat$Treatment_N^2)*(dat$Treatment_MeanCorrected^4))) + ((dat$Control_SDCorrected^4)/((dat$Control_N^2)*(dat$Control_MeanCorrected^4))))

write.csv(dat, "D:/Juan/Desktop/Universidad/PHD/Chapters/Chapter4_InsectsNitrogen/Github/Nimpact_Invertebrates/Databases/AbundanceDataset_ArthNema_imputed.csv", row.names = FALSE)

# Richness ####
datR <- read.csv("D:/Juan/Desktop/Universidad/PHD/Chapters/Chapter4_InsectsNitrogen/Github/Nimpact_Invertebrates/Databases/RichnessDataset_csv.csv")

# Remove studies with Deposition gradient
datR <- datR %>% 
  dplyr::filter(N_Added_Deposited == "Added")

# Remove papers with Raw data
raws <- datR %>% 
  dplyr::filter(str_detect(Data_Retrieved_From, "(?i)raw"))

datR <- datR %>% 
  dplyr::filter(!(Source %in% unique(raws$Source)))

# Filter to get only the treatments and not the controls
datR <- datR %>% 
  dplyr::group_by(Source, Study, Control_ID) %>% 
  dplyr::mutate(ControlYN = case_when(Nitrogen_Added == min(Nitrogen_Added) ~ "Control",
                                      TRUE ~ "Treatment"))
datRT <- datR %>% 
  dplyr::filter(ControlYN == "Treatment") %>% 
  dplyr::ungroup()


# Remove observations where the Control and the Treatment are 0
datRT <- datRT %>% 
  filter(!(Control_Mean == 0 & Treatment_Mean == 0))

# Remove studies where the control is not 0
ControlNoZero <- datR %>% 
  dplyr::filter(ControlYN == "Control") %>% 
  dplyr::filter(Nitrogen_Added != 0)

ControlNoZero <- unique(ControlNoZero$Source)

datRT <- datRT %>% 
  filter(!(Source %in% ControlNoZero))

# Summary numbers per phylum
table(datRT$Phylum)

# Correct 0s in the Mean
datRT <- datRT %>% 
  dplyr::mutate(Control_MeanCorrected = case_when(Control_Mean == 0 ~ 1/(2 * Control_N * Control_ScalingConstant),
                                                  TRUE ~ Control_Mean),
                Treatment_MeanCorrected = case_when(Treatment_Mean == 0 ~ 1/(2 * Treatment_N * Treatment_ScalingConstant),
                                                    TRUE ~ Treatment_Mean))

# Impute SDs Check SDs that are 0 and convert to NA, better to impute all 0s because it's hard to believe that they are not 0
datRT <- datRT %>% 
  mutate(Control_SDCorrected = case_when(Control_SD == 0 ~ NA_real_,
                                         TRUE ~ Control_SD),
         Treatment_SDCorrected = case_when(Treatment_SD == 0 ~ NA_real_,
                                           TRUE ~ Treatment_SD))


# Subset per phylum
ArthR <- datRT %>% 
  filter(Phylum %in% "Arthropoda") %>% 
  as.data.frame()

NemaR <- datRT %>% 
  filter(Phylum %in% "Nematoda") %>% 
  as.data.frame()


# Impute missing SDs separately
impute_missingness(ArthR)
ArthR_imp <- impute_SD(aDataFrame = ArthR, columnSDnames = "Control_SDCorrected", columnXnames = "Control_MeanCorrected", method = "Bracken1992")

ArthR_imp <- impute_SD(aDataFrame = ArthR_imp, columnSDnames = "Treatment_SDCorrected", columnXnames = "Treatment_MeanCorrected", method = "Bracken1992")
impute_missingness(ArthR_imp)

impute_missingness(NemaR)
NemaR_imp <- impute_SD(aDataFrame = NemaR, columnSDnames = "Control_SDCorrected", columnXnames = "Control_MeanCorrected", method = "Bracken1992")

NemaR_imp <- impute_SD(aDataFrame = NemaR_imp, columnSDnames = "Treatment_SDCorrected", columnXnames = "Treatment_MeanCorrected", method = "Bracken1992")
impute_missingness(NemaR_imp) 

# Join datRabases again
datR <- rbind(ArthR_imp, NemaR_imp)

# Calculate logRR with bias correction (Delta method according to Lajeunesse 2015)
datR$Ratio <- datR$Treatment_MeanCorrected/datR$Control_MeanCorrected
datR$logRR_uncorrected <- log(datR$Ratio)
datR$logRR <- datR$logRR_uncorrected + 
  0.5*(((datR$Treatment_SDCorrected^2)/(datR$Treatment_N*(datR$Treatment_MeanCorrected^2))) - ((datR$Control_SDCorrected^2)/(datR$Control_N*(datR$Control_MeanCorrected^2))))


datR$VAR_uncorrected <- (as.numeric(datR$Treatment_SDCorrected)^2)/(datR$Treatment_N*datR$Treatment_MeanCorrected^2) + (as.numeric(datR$Control_SDCorrected)^2)/(datR$Control_N*datR$Control_MeanCorrected^2)
datR$VAR <- datR$VAR_uncorrected + 0.5*(((datR$Treatment_SDCorrected^4)/((datR$Treatment_N^2)*(datR$Treatment_MeanCorrected^4))) + ((datR$Control_SDCorrected^4)/((datR$Control_N^2)*(datR$Control_MeanCorrected^4))))

write.csv(datR, "D:/Juan/Desktop/Universidad/PHD/Chapters/Chapter4_InsectsNitrogen/Github/Nimpact_Invertebrates/Databases/RichnessDataset_ArthNema_imputed.csv", row.names = FALSE)

