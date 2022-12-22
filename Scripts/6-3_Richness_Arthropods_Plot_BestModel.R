# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com

############## Load libraries and dataset ##############

library(tidyverse)
library(dummies)
library(metafor)
library(viridis)
library(officer)
library(flextable)
library(magrittr)
# Read data
datR <- read.csv("./Databases/RichnessDataset_ArthNema_imputed_Env.csv")

#Subset for the Nematoda
ArthR <- datR %>%
        dplyr::filter(Phylum %in% "Arthropoda")

# Transform continous variables
ArthR$logNadd <- log10(ArthR$Nitrogen_Added)
ArthR$logNadd2 <- ArthR$logNadd^2

summary(ArthR$logNadd)

I2df_RichArth <- read.csv("./Results/Arthropods/Richness_Arthropods_model_selection_160.csv")

I2df_RichArth %>% 
   arrange(AICc) %>% 
   slice(1)

# load models
load(file = "./Results/Arthropods/Richness_modelsArth_160.RData")

# Best Model ~1 although the 3 of them are within 2 AIC units
modeliArthR <- modelsArthR[[1]]
anova(modeliArthR)
summary(ArthR$logNadd)

# Export summary table
summary(modeliArthR)

sumtable <- data.frame(Moderator = rownames(modeliArthR$b),
                       Estimate = round(modeliArthR$b[,1],2),
                       Std.Error = round(modeliArthR$se,2),
                       CI.low = round(modeliArthR$ci.lb,2),
                       CI.up = round(modeliArthR$ci.ub,2),
                       p.value = round(modeliArthR$pval,2))



sumtable <- regulartable(sumtable) %>% 
        theme_booktabs() %>%  
        fontsize(part = "header", size= 12) %>% 
        flextable::font(part = "all",fontname= "Times") %>% 
        autofit()

read_docx() %>% body_add_flextable(sumtable) %>% print(target = "./Results/Arthropods/Richness_Arthropods_BestModel_Null_Summary.docx")


# Create grid for predictions
grid <- expand_grid(logNadd = seq(1,3.1, by = 0.1))

newMods <- cbind(grid$logNadd)

#Predict
preds <- predict(modeliArthR) # with addx=T you can see how the variables were used

# Convert the predictions into a dataframe and add the values of Nadd and Feeding Guild in a way that we can plot it later
predictions <- as.data.frame(preds)
predictions <- rbind(predictions,predictions,predictions)
predictions$logNadd <- c(1,2,3.025)

#### Plot ######
# With CI
(NplotRichArth <- ggplot() +
         geom_hline(yintercept = 0, linetype="longdash", color = "grey", size = 1)+
         geom_point(data=ArthR, aes(x=logNadd, y=logRR, size = 1/VAR), alpha = 0.2, shape=16, color="#6F3F97", show.legend=FALSE)+
         
         geom_line(data = predictions, aes(logNadd, pred),size = 1.5, color="#6F3F97")+
         # geom_line(data = predictions, aes(logNadd, ci.ub),size = 1, linetype="dashed", alpha = 0.6, col = "#6F3F97")+
         # geom_line(data = predictions, aes(logNadd, ci.lb),size = 1, linetype="dashed", alpha = 0.6, col = "#6F3F97")+
         geom_ribbon(data = predictions, aes(logNadd, ymin = ci.lb, ymax = ci.ub), alpha = 0.2, fill="#6F3F97")+
         
         coord_cartesian(ylim=c(-2.6,2.6), xlim = c(0.5,3.5))+
         scale_x_continuous(breaks = seq(0,3, by = 1))+
         # scale_y_continuous(breaks = c(-1.5,0,1.5))+
         
         ylab(expression(paste("Effect size (",lnRR^Delta, ")")))+
         xlab(expression(paste(log[10],"(Nitrogen kg/ha/yr)")))+
         theme_classic(base_size = 40)+
         theme(strip.text = element_text(size = 16, face = "bold"),
               strip.background = element_blank(),
               axis.title.y=element_blank(),
               axis.title.x=element_blank(),
               legend.title = element_blank(),
               legend.text = element_text(size=12)))

pdf("Results/Arthropods/Plots/Richness_Arth_Nadd_axis26.pdf", width = 8, height =7)
NplotRichArth
dev.off()
