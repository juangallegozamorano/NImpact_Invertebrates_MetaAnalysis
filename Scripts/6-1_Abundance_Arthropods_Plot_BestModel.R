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
dat <- read.csv("./Databases/AbundanceDataset_Arth_CompleteOrders_imputed.csv")

#Subset for the Arthropods with Order
ArthO <- dat %>% 
  dplyr::filter(Phylum %in% "Arthropoda") %>% 
  drop_na(Order)

load(file = "./Results/Arthropods/Abundance_Arthropods_MetamorphosisMAT_512phylo2.RData")
modeliMetamorMAT <- modeliMetamorMATphylo2 # change variable name

# Export summary table
summary(modeliMetamorMAT)
sumtable <- data.frame(Moderator = rownames(modeliMetamorMAT$b),
                       Estimate = round(modeliMetamorMAT$b[,1],2),
                       Std.Error = round(modeliMetamorMAT$se,2),
                       CI.low = round(modeliMetamorMAT$ci.lb,2),
                       CI.up = round(modeliMetamorMAT$ci.ub,2),
                       p.value = round(modeliMetamorMAT$pval,2))

sumtable <- regulartable(sumtable) %>% 
  theme_booktabs() %>%  fontsize(part = "header", size= 12) %>% font(part = "all",fontname= "Times") %>% autofit()

read_docx() %>% body_add_flextable(sumtable) %>% print(target = "./Results/Arthropods/Abundance_Arthropods_BestModel_MetamorphosisMAT_Summary.docx")


# Wald test for the different moderators
anova(modeliMetamorMAT, btt=2:3) # If N is significant but only within the Reference level (Complete metamor)
anova(modeliMetamorMAT, btt=4) # If MAT is signifanctly different from 0 when N = 0
anova(modeliMetamorMAT, btt=5:6) # If Metamor groups are signifanctly different from 0 when N = 0
anova(modeliMetamorMAT, btt=7:8)# If different MAT differ from 0 when interacting with N (so if the interaction is significant)

anova(modeliMetamorMAT, btt=9:12)# If different Metamor levels differ from 0 when interacting with N (so if the interaction is significant)

# Create Grid
grid <- expand_grid(logNadd = seq(0.5,3.1, by = 0.01),
                    Metamorphosis = as.factor(c("Complete","Incomplete-Gradual","None")),
                    MAT =  seq(0,25, by = 0.5))

dummyMetamor <- dummy(grid$Metamorphosis, sep = ".")

newMods <- cbind(grid$logNadd, 
                 grid$logNadd^2,
                 grid$MAT,
                 dummyMetamor,
                 grid$logNadd*grid$MAT,
                 grid$logNadd^2*grid$MAT,
                 grid$logNadd*dummyMetamor,
                 grid$logNadd^2*dummyMetamor)


names(coef(modeliMetamorMAT)[-1])
colnames(newMods)

newMods <- newMods[, -c(4,9,12)] # Remove columns of the reference level
names(coef(modeliMetamorMAT)[-1])
colnames(newMods)

# Predict
preds <- predict(modeliMetamorMAT, newmods = newMods, addx = T) # with addx=T you can see how the variables were used

# Convert the predictions into a dataframe and add the values of Nadd and Feeding Guild in a way that we can plot it later
predictions <- as.data.frame(preds)
predictions$logNadd <- predictions$X.logNadd 
predictions$NitrogAdd <- 10^(predictions$logNadd)
predictions$MAT <- predictions$X.CRU_MAT_mean
predictions$PredPercent <- ((exp(predictions$pred))-1)*100
predictions <- predictions %>% 
   mutate(Metamorphosis = case_when(X.MetamorphosisIncomplete.Gradual == 1 ~ "Incomplete-Gradual",
                                    X.MetamorphosisNone == 1 ~ "None",
                                    TRUE ~ "Complete"))
ArthO %>% 
   dplyr::group_by(Metamorphosis) %>% 
   dplyr::summarize(MinT = min(CRU_MAT_mean),
             MeanT = mean(CRU_MAT_mean),
             MaxT = max(CRU_MAT_mean))

ArthO %>% 
  dplyr::group_by(Metamorphosis) %>% 
  dplyr::summarize(MinN = min(logNadd),
                   MeanN= mean(logNadd),
                   MaxN = max(logNadd))


completepredictions <- predictions %>% 
   dplyr::filter(Metamorphosis == "Complete", logNadd >= 0.85, logNadd <= 3.1)
incompletepredictions <- predictions %>% 
  dplyr::filter(Metamorphosis == "Incomplete-Gradual", logNadd >= 0.4)
nonepredictions <- predictions %>% 
  dplyr::filter(Metamorphosis == "None", logNadd >= 1.4, logNadd <= 2.42)

predictionsmeta <- rbind(completepredictions,incompletepredictions,nonepredictions)

#### Plot ######
# With CI
summary(ArthO$CRU_MAT_mean)

(MAT_plot <- ggplot() +
    geom_hline(yintercept = 0, linetype="longdash", color = "grey")+
    geom_point(data=ArthO, aes(x=logNadd, y=logRR, size = 1/VAR, color = CRU_MAT_mean), alpha = 0.1, shape=16)+
    
    geom_line(data = predictionsmeta, aes(logNadd, pred, col = MAT, group = MAT),size = 1, alpha = 0.7)+

    geom_line(data = predictionsmeta %>% filter(MAT == 10), aes(logNadd, pred), col = "black",size = 1)+
    geom_line(data = predictionsmeta %>% filter(MAT == 0), aes(logNadd, ci.ub, col = MAT, group = MAT),size = 1, linetype="dashed", alpha = 0.5)+
    geom_line(data = predictionsmeta %>% filter(MAT == 0), aes(logNadd, ci.lb, col = MAT, group = MAT),size = 1, linetype="dashed", alpha = 0.5)+
    geom_line(data = predictionsmeta %>% filter(MAT == 25), aes(logNadd, ci.ub, col = MAT, group = MAT),size = 1, linetype="dashed", alpha = 0.5)+
    geom_line(data = predictionsmeta %>% filter(MAT == 25), aes(logNadd, ci.lb, col = MAT, group = MAT),size = 1, linetype="dashed", alpha = 0.5)+
    
    scale_color_gradient(low = "#394998", high = "#F9C659",limits = c(0,25), breaks = c(0, 5, 10, 15, 20, 25))+
    
    scale_x_continuous(breaks = seq(0,3, by = 1))+
    coord_cartesian(ylim=c(-3.5,3.5), xlim = c(0.5,3.5))+
    ylab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    xlab(expression(paste(log[10],"(Nitrogen kg/ha/year)")))+
    labs(color="MAT (ºC)") +
    theme_classic(base_size = 40)+
    facet_wrap(~Metamorphosis, scales = "free")+
    guides(size = "none")+
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          axis.text = element_text(size=35), 
          axis.title.y=element_blank(),
          legend.title=element_blank(),
          axis.title.x=element_blank(),
          legend.text = element_text(size=25), 
          legend.key.height = unit(1, "cm")))

pdf("Results/Arthropods/Plots/Abundance_MetamormophosisMAT.pdf",  width = 20, height = 6)
MAT_plot
dev.off()
