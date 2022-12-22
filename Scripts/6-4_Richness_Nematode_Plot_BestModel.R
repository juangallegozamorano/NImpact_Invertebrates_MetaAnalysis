# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com


############## Load libraries and dataset ##############

library(tidyverse) # Easily Install and Load the 'Tidyverse', CRAN v1.3.1
library(dummies) # Create dummy/indicator variables flexibly and efficiently, CRAN v1.5.6
library(metafor) # Meta-Analysis Package for R, CRAN v3.4-0
library(officer) # Manipulation of Microsoft Word and PowerPoint Documents, CRAN v0.3.18
library(flextable) # Functions for Tabular Reporting, CRAN v0.6.4

# Read data
datR <- read.csv("./Databases/RichnessDataset_ArthNema_imputed_Env.csv")

#Subset for the Nematoda
NemaR <- datR %>% 
    dplyr::filter(Phylum %in% "Nematoda")

# Transform continous variables
NemaR$logNadd <- log10(NemaR$Nitrogen_Added)
NemaR$logNadd2 <- NemaR$logNadd^2

# load models
load(file = "./Results/Nematode/Richness_modelsNema_42.RData")

I2df_RichNema <- read.csv("./Results/Nematode/Richness_model_selection_42.csv")

I2df_RichNema %>% 
    arrange(AICc)%>% 
    slice(1)

# Best Model logNadd+logNadd2
modeli <- modelsNemaRich[[3]]
summary(NemaR$logNadd)

# Export summary table
summary(modeli)
anova(modeli)

sumtable <- data.frame(Moderator = rownames(modeli$b),
                       Estimate = round(modeli$b[,1],2),
                       Std.Error = round(modeli$se,2),
                       CI.low = round(modeli$ci.lb,2),
                       CI.up = round(modeli$ci.ub,2),
                       p.value = round(modeli$pval,2))



sumtable <- regulartable(sumtable) %>% 
    theme_booktabs() %>%  fontsize(part = "header", size= 12) %>% flextable::font(part = "all",fontname= "Times") %>% autofit()

read_docx() %>% body_add_flextable(sumtable) %>% print(target = "./Results/Nematode/Richness_Nematodes_BestModel_Nadd2_Summary.docx")


# Create grid for predictions
grid <- expand_grid(logNadd = seq(1.22,2.77, by = 0.01))

newMods <- cbind(grid$logNadd,
                 grid$logNadd^2)

#Predict
preds <- predict(modeli, newmods = newMods, addx = T) # with addx=T you can see how the variables were used

# Convert the predictions into a dataframe and add the values of Nadd and Feeding Guild in a way that we can plot it later
predictionsRichNema <- as.data.frame(preds)
predictionsRichNema$logNadd <- predictionsRichNema$X.logNadd 
predictionsRichNema$NitrogAdd <- 10^(predictionsRichNema$logNadd)
predictionsRichNema$PredPercent <- ((exp(predictionsRichNema$pred))-1)*100

#### Plots ######
(NplotRichNema <- ggplot() +
     geom_hline(yintercept = 0, linetype="longdash", color = "grey", size = 1)+
     geom_point(data=NemaR, aes(x=logNadd, y=logRR, size = 1/VAR), alpha = 0.2, shape=16, color="#FDBD63", show.legend=FALSE)+
     
     geom_line(data = predictionsRichNema, aes(logNadd, pred),size = 1.5, col = "#FDBD63")+
     geom_ribbon(data = predictionsRichNema, aes(logNadd, ymin = ci.lb, ymax = ci.ub), alpha = 0.2, fill = "#E9AC56")+
     
     coord_cartesian(ylim=c(-2.6,2.6), xlim = c(0.5,3.5))+
     scale_x_continuous(breaks = seq(0,3, by = 1))+
   
     ylab(expression(paste("Effect size (",lnRR^Delta, ")")))+
     xlab(expression(paste(log[10],"(Nitrogen kg/ha/yr)")))+
     theme_classic(base_size = 40)+
     theme(strip.text = element_text(size = 16, face = "bold"),
           axis.title.y=element_blank(),
           axis.title.x=element_blank(),
           strip.background = element_blank()))

pdf("Results/Nematode/Plots/Richness_Nema_Nadd_axis26.pdf",  width = 8, height = 7)
NplotRichNema
dev.off()

