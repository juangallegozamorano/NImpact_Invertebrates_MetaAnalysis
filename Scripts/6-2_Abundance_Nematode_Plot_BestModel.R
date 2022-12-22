# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com


############## Load libraries and dataset ##############

library(tidyverse) # Easily Install and Load the 'Tidyverse', CRAN v1.3.1
library(dummies) # Create dummy/indicator variables flexibly and efficiently, CRAN v1.5.6
library(metafor) # Meta-Analysis Package for R, CRAN v3.4-0
library(viridis) # Default Color Maps from 'matplotlib', CRAN v0.5.1
library(officer) # Manipulation of Microsoft Word and PowerPoint Documents, CRAN v0.3.18
library(flextable) # Functions for Tabular Reporting, CRAN v0.6.4

# Read data
Nema <- read.csv("./Databases/AbundanceDataset_Nema_FeedingGuild_imputed.csv")

Nema$logNadd <- log10(Nema$Nitrogen_Added)


I2dfNema <- read.csv("./Results/Nematode/Abundance_Nematodes_model_selection_256.csv")
# Best model
I2dfNema %>% 
   arrange(AICc) %>% 
   slice(1)

summary(Nema$CRU_MAP_mean)
summary(Nema$logNadd)

# Load the Best model: modeliMAP
load(file = "./Results/Nematode/Abundance_Bestmodel_MAP_FeedingGuild.RData")

# Export summary table
summary(modeliMAP)
sumtable <- data.frame(Moderator = rownames(modeliMAP$b),
                       Estimate = round(modeliMAP$b[,1],2),
                       Std.Error = round(modeliMAP$se,2),
                       CI.low = round(modeliMAP$ci.lb,2),
                       CI.up = round(modeliMAP$ci.ub,2),
                       p.value = round(modeliMAP$pval,2))

sumtable <- regulartable(sumtable) %>% 
  theme_booktabs() %>%  fontsize(part = "header", size= 12) %>% font(part = "all",fontname= "Times") %>% autofit()

read_docx() %>% body_add_flextable(sumtable) %>% print(target = "./Results/Nematode/Abundance_Nematodes_BestModel_MAPFeeding_Summary.docx")

# Check significance of moderators
anova(modeliMAP, btt = 2) # Nadd within refernce level, is significantly different from 0?
anova(modeliMAP, btt = 3:7) # FGuilds compared to reference level (Detritivores)
anova(modeliMAP, btt = 8) # MAP within refernce level, is significantly different from 0?
anova(modeliMAP, btt = 9:13) # FGuilds interacting with N
anova(modeliMAP, btt = 14) # MAP interacting with N

# Create grid for predictions
grid <- expand_grid(logNadd = seq(1.2,2.8, by = 0.1),
                    Feeding_GuildUsed = as.factor(c("Bacterivore","Fungivore","Herbivore","Parasite", "Predator-Omnivore", "Unknown")),
                    
                    MAP = seq(0,200, by = 5))

dummyFeeding <- dummy(grid$Feeding_GuildUsed, sep = ".")

newMods <- cbind(grid$logNadd, 
                 dummyFeeding,
                 grid$MAP, 
                 grid$logNadd*dummyFeeding,
                 grid$logNadd*grid$MAP)

colnames(newMods)

newMods <- newMods[, -c(2,9)] # Remove columns of the reference level

# Predict
predsMAP <- predict(modeliMAP, newmods = newMods, addx = T) # with addx=T you can see how the variables were used

# Convert the predictions into a dataframe and add the values of Nadd and Feeding Guild in a way that we can plot it later
predictionsMAP <- as.data.frame(predsMAP)
predictionsMAP$logNadd <- predictionsMAP$X.logNadd 
predictionsMAP$Nadd <- 10^(predictionsMAP$X.logNadd) 
predictionsMAP$MAP <- predictionsMAP$X.CRU_MAP_mean
predictionsMAP$PredPercent <- ((exp(predictionsMAP$pred))-1)*100
predictionsMAP$Feeding_GuildUsed <- grid$Feeding_GuildUsed


#### Plot ######
Nema %>% 
  dplyr::group_by(Feeding_GuildUsed) %>% 
  dplyr::summarize(MinNad = min(logNadd),
                   MeanNad = mean(logNadd),
                   MaxNad = max(logNadd))

Bac_predictions <- predictionsMAP %>% 
  filter(Feeding_GuildUsed == "Bacterivore", logNadd >= 1.2)
Fungi_predictions <- predictionsMAP %>% 
  filter(Feeding_GuildUsed == "Fungivore", logNadd >= 1.2)
Herb_predictions <- predictionsMAP %>% 
  filter(Feeding_GuildUsed == "Herbivore", logNadd >= 1.2)
Para_predictions <- predictionsMAP %>% 
  filter(Feeding_GuildUsed == "Parasite", logNadd >= 1.2)
PredOmn_predictions <- predictionsMAP %>% 
  filter(Feeding_GuildUsed == "Predator-Omnivore", logNadd >= 1.2)
Unk_predictions <- predictionsMAP %>% 
  filter(Feeding_GuildUsed == "Unknown", logNadd >= 1.4)

predictionsMAP <- rbind(Bac_predictions, Fungi_predictions, Herb_predictions, Para_predictions, PredOmn_predictions,Unk_predictions )
predictionsMAP <- predictionsMAP %>% 
  mutate(Feeding_GuildUsed = case_when(as.character(Feeding_GuildUsed) %in% "Parasite" ~ "Plant-parasites",
                                       TRUE ~ as.character(Feeding_GuildUsed))) %>% 
  filter(MAP >= 5)

Nema <- Nema %>% 
  mutate(Feeding_GuildUsed = case_when(as.character(Feeding_GuildUsed) %in% "Parasite" ~ "Plant-parasites",
                                       TRUE ~ as.character(Feeding_GuildUsed)))
# With CI
(MAP_plot_CI <- ggplot() +
    geom_hline(yintercept = 0, linetype="longdash", color = "grey")+
    geom_point(data=Nema, aes(x=logNadd, y=logRR, size = 1/VAR, col = CRU_MAP_mean), alpha = 0.2, shape=16)+

    geom_line(data = predictionsMAP, aes(logNadd, pred, col = MAP, group = MAP),size = 1, alpha = 0.7)+
    geom_line(data = predictionsMAP %>% filter(MAP == 80), aes(logNadd, pred), col = "black",size = 1)+
    geom_line(data = predictionsMAP %>% filter(MAP == 5), aes(logNadd, ci.ub, col = MAP, group = MAP),size = 1, linetype="dashed", alpha = 0.5)+
    geom_line(data = predictionsMAP %>% filter(MAP == 5), aes(logNadd, ci.lb, col = MAP, group = MAP),size = 1, linetype="dashed", alpha = 0.5)+
    geom_line(data = predictionsMAP %>% filter(MAP == 200), aes(logNadd, ci.ub, col = MAP, group = MAP),size = 1, linetype="dashed", alpha = 0.5)+
    geom_line(data = predictionsMAP %>% filter(MAP == 200), aes(logNadd, ci.lb, col = MAP, group = MAP),size = 1, linetype="dashed", alpha = 0.5)+

    
    # scale_linetype_manual(values=c("dotted", "dashed", "solid" ))+
    
    scale_color_gradient(low = "#F9C659", high = "#394998",limits = c(0,200), breaks = c(0, 50, 100, 150, 200))+
    # scale_color_viridis(option = "cividis", discrete = F, direction = -1, begin = 0.1, end = 0.9,
    #                     limits = c(0,200), breaks = c(0, 50, 100, 150, 200))+
    
    coord_cartesian(ylim=c(-3.5,3.5), xlim = c(0.5,3.5))+
    scale_x_continuous(breaks = seq(0,3, by = 1))+
    ylab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    xlab(expression(paste(log[10],"(Nitrogen kg/ha/year)")))+
    labs(color="MAP (mm/month)") +
    facet_wrap(~Feeding_GuildUsed, scales = "free")+
    theme_classic(base_size = 24)+
    guides(size = "none")+
    theme(strip.text = element_text(size = 20, face = "bold"),
          strip.background = element_blank(),
          axis.text = element_text(size=24), 
          legend.text = element_text(size=20), 
          legend.key.height = unit(1, "cm")))

pdf("Results/Nematode/Plots/Abundance_MAP_FeedingGuild_CI.pdf",  width = 15, height = 9)
MAP_plot_CI
dev.off()
