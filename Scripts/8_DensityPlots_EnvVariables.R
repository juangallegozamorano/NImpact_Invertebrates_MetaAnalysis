# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com


############## Load libraries and dataset ##############
library(tidyverse)
library(ggpubr)
library(tictoc)
library(plyr)

# Read data
dat <- read.csv("./Databases/AbundanceDataset_ArthNema_imputed_Env.csv")
datR <- read.csv("./Databases/RichnessDataset_ArthNema_imputed_Env.csv")

# Transform continuous variables
dat$logNadd <- log10(dat$Nitrogen_Added)
dat$logNadd2 <- dat$logNadd^2

dat$logNdep <- log10(dat$Ndep_1984_2015)
dat$logNdep2 <- dat$logNdep^2

datR$logNadd <- log10(datR$Nitrogen_Added)
datR$logNadd2 <- datR$logNadd^2

datR$logNdep <- log10(datR$Ndep_1984_2015)
datR$logNdep2 <- datR$logNdep^2

# Nitrogen added ####
# Abundance
ggplot(data=dat, aes(x=logNadd, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab(expression(paste(log[10],"(Nitrogen kg/ha/yr)")))+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")

ggsave(filename = "./Results/DensityVariables/Abundance_Nadd.pdf", width = 6, height =4, units = "in")

# Richness
ggplot(data=datR, aes(x=logNadd, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab(expression(paste(log[10],"(Nitrogen kg/ha/yr)")))+
  scale_x_continuous( breaks = seq(1, 3, by = 1))+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")
ggsave(filename = "./Results/DensityVariables/Richness_Nadd.pdf", width = 6, height =4, units = "in")

# Nitrogen deposition ####
# Abundance
ggplot(data=dat, aes(x=logNdep, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab(expression(paste(log[10],"(Nitrogen dep. kg/ha/yr)")))+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")

ggsave(filename = "./Results/DensityVariables/Abundance_Ndep.pdf", width = 6, height =4, units = "in")

# Richness
ggplot(data=datR, aes(x=logNdep, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab(expression(paste(log[10],"(Nitrogen dep. kg/ha/yr)")))+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")
ggsave(filename = "./Results/DensityVariables/Richness_Ndep.pdf", width = 6, height =4, units = "in")


# Duration ####
# Abundance
ggplot(data=dat, aes(x=Duration_Years, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab("Experiment duration (years)")+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")

ggsave(filename = "./Results/DensityVariables/Abundance_duration.pdf", width = 6, height =4, units = "in")
# Richness
ggplot(data=datR, aes(x=Duration_Years, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab("Experiment duration (years)")+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")

ggsave(filename = "./Results/DensityVariables/Richness_duration.pdf", width = 6, height =4, units = "in")

# MAT ####
# Abundance
ggplot(data=dat, aes(x=CRU_MAT_mean, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab("MAT (degrees Cº)")+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")

ggsave(filename = "./Results/DensityVariables/Abundance_MAT.pdf", width = 6, height =4, units = "in")

# Richness
ggplot(data=datR, aes(x=CRU_MAT_mean, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab("MAT (degrees Cº)")+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")
ggsave(filename = "./Results/DensityVariables/Richness_MAT.pdf", width = 6, height =4, units = "in")


# MAP ####
# Abundance
ggplot(data=dat, aes(x=CRU_MAP_mean, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab("MAP (mm/month)")+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")
ggsave(filename = "./Results/DensityVariables/Abundance_MAP.pdf", width = 6, height =4, units = "in")

# Richness
ggplot(data=datR, aes(x=CRU_MAP_mean, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab("MAP (mm/month)")+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")
ggsave(filename = "./Results/DensityVariables/Richness_MAP.pdf", width = 6, height =4, units = "in")



# CEC ####
# Abundance
ggplot(data=dat, aes(x=CEC, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab("CEC (cmol/kg)")+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")
ggsave(filename = "./Results/DensityVariables/Abundance_CEC.pdf", width = 6, height =4, units = "in")

# Richness
ggplot(data=datR, aes(x=CEC, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab("CEC (cmol/kg)")+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")
ggsave(filename = "./Results/DensityVariables/Richness_CEC.pdf", width = 6, height =4, units = "in")


# Feeding group ####
dat <- dat %>% 
  mutate(Feeding_GuildUsed = as.factor(case_when(Feeding_Guild %in% c("Parasite", "Parasitoid", "Plant parasite")~ "Parasite",
                                                 Feeding_Guild %in% c("Nectar feeders", "Herbivore/Fungivore","Fungi", "Herbivore", "Fungivore")~ "Herbivore-Fungivore",
                                                 Feeding_Guild %in% c("Detrivore")~ "Detritivore",
                                                 Feeding_Guild %in% c("Omnivore", "Predator")~ "Predator-Omnivore",
                                                 str_detect(Feeding_Guild, "(?i)Preda|Omniv|Free|Microbi")~ "Predator-Omnivore",
                                                 str_detect(Feeding_Guild, "(?i)Parasit")~ "Parasite",
                                                 TRUE ~ Feeding_Guild)),
         Feeding_GuildCode = as.numeric(Feeding_GuildUsed))

datR <- datR %>% 
  mutate(Feeding_GuildUsed = as.factor(case_when(Feeding_Guild %in% c("Nectar feeders", "Herbivore/Fungivore","Fungi", "Herbivore", "Fungivore")~ "Herbivore-Fungivore",
                                                 Feeding_Guild %in% c("Detrivore")~ "Detritivore",
                                                 str_detect(Feeding_Guild, "(?i)Preda|Omniv|Free|Microbi")~ "Predator-Omnivore",
                                                 str_detect(Feeding_Guild, "(?i)Parasit")~ "Parasite",
                                                 TRUE ~ Feeding_Guild)),
         Feeding_GuildCode = as.numeric(Feeding_GuildUsed))

# Abundance
DatSumFed <- dat %>% 
  group_by(Phylum,Feeding_GuildUsed ) %>% 
  tally() %>% 
  mutate(Label = case_when(Feeding_GuildUsed == "Detritivore" ~ "Det",
                           Feeding_GuildUsed == "Herbivore-Fungivore" ~ "Herb",
                           Feeding_GuildUsed == "Parasite" ~ "Par",
                           Feeding_GuildUsed == "Bacterivore" ~ "Bac",
                           Feeding_GuildUsed == "Predator-Omnivore" ~ "Pred",
                           Feeding_GuildUsed == "Unknown" ~ "Unk"))

ggplot(data=DatSumFed, aes(x=Label, y=n, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), alpha=.4)+
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab("Feeding guild")+
  ylab("Observations")+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")

ggsave(filename = "./Results/DensityVariables/Abundance_FeedingGuild.pdf", width = 6, height =4, units = "in")

# Richness
DatSumFedR <- datR %>% 
  group_by(Phylum,Feeding_GuildUsed ) %>% 
  tally() %>% 
  mutate(Label = case_when(Feeding_GuildUsed == "Detritivore" ~ "Det",
                           Feeding_GuildUsed == "Herbivore-Fungivore" ~ "Herb",
                           Feeding_GuildUsed == "Parasite" ~ "Par",
                           Feeding_GuildUsed == "Predator-Omnivore" ~ "Pred",
                           Feeding_GuildUsed == "Unknown" ~ "Unk"))

ggplot(data=DatSumFedR, aes(x=Label, y=n, group=Phylum, fill=Phylum, col = Phylum)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), alpha=.4)+
  scale_fill_manual(values = c("#6F3F97", "#FDBD63"))+
  scale_color_manual(values = c("#6F3F97", "#FDBD63"))+
  xlab("Feeding guild")+
  ylab("Observations")+
  theme_classic(base_size = 26)+
  theme(legend.position = "None")
ggsave(filename = "./Results/DensityVariables/Richness_FeedingGuild.pdf", width = 6, height =4, units = "in")

