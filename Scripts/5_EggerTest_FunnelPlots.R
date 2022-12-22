library(metafor)
library(tidyverse)
library(ggpubr)
library(tictoc)
library(plyr)
# load necessary functions
source("Scripts/00_Functions.R")

# Abundance ####
# Read data
dat <- read.csv("./Databases/AbundanceDataset_ArthNema_imputed_Env.csv")

#Subset for the Arthropods with Order
Arth <- dat %>% 
  dplyr::filter(Phylum %in% "Arthropoda") 

#Subset for the Arthropods with Order
Nema <- dat %>% 
  dplyr::filter(Phylum %in% "Nematoda")


calc.vRRDelta <- function(x) {
  v = matrix((x$Control_SDCorrected[1]^2 / (x$Control_N[1] * x$Control_MeanCorrected[1]^2) + 
                0.5*(x$Control_SDCorrected[1]^4/(x$Control_N[1]^2*x$Control_MeanCorrected[1]^4))), nrow = nrow(x), ncol = nrow(x))
  diag(v) = x$VAR
  return(v)
} 

# Arrange based on Control ID to calculate properly the var-covar matrix
Arth <- Arth %>%
  arrange(Control_ID)

Nema <- Nema %>%
  arrange(Control_ID)

# For Arthropods
ArthV <- bldiag(lapply(split(Arth, Arth$Control_ID), calc.vRRDelta))
is.positive.definite(ArthV) # TRUE!

# For Nematodes
NemaV <- bldiag(lapply(split(Nema, Nema$Control_ID), calc.vRRDelta))
is.positive.definite(NemaV) # TRUE!

# Create a unique identifier for each effect size, which corresponds to the observational id (i.e. residual variance in the model)
Arth$RowID <-  1:length(Arth$Source)
Nema$RowID <-  1:length(Nema$Source)

# Set the Random Effects
RandomEffects <- list(~1|RowID, ~1|Source)

# Calculate seinv
Arth$se <- sqrt(Arth$VAR)
Arth$invse <- 1/Arth$se
Arth$e_n <- with(Arth, (4*(Control_N*Treatment_N)) / (Control_N + Treatment_N))
Arth$inv_n_tilda <-  with(Arth, (Control_N + Treatment_N)/(Control_N*Treatment_N))
Arth$sqrt_inv_n_tilda <-  with(Arth, sqrt(inv_n_tilda))

Nema$se <- sqrt(Nema$VAR)
Nema$invse <- 1/Nema$se
Nema$e_n <- with(Nema, (4*(Control_N*Treatment_N)) / (Control_N + Treatment_N))
Nema$inv_n_tilda <-  with(Nema, (Control_N + Treatment_N)/(Control_N*Treatment_N))
Nema$sqrt_inv_n_tilda <-  with(Nema, sqrt(inv_n_tilda))

# Egger test Abundance ####
# Abundance arthropods
Arth_egger <- rma.mv(logRR,
                   V = ArthV, 
                   mod = ~ 1 + sqrt_inv_n_tilda,
                   data = Arth, 
                   random = RandomEffects,
                   test = "t",
                   method = "REML", 
                   level=95, digits=4) 


# when p-val of mods is significant < 0.05 asimmetry (caused either by true heterogeneity or publication bias) is detected
summary(Arth_egger)# pval > 0.05 for mods, good!

# Abundance nematodes
Nema_egger <- rma.mv(logRR,
                     V = NemaV, 
                     mod = ~ 1 + sqrt_inv_n_tilda,
                     data = Nema, 
                     random = RandomEffects,
                     test = "t",
                     method = "REML", 
                     level=95, digits=4) 

# when p-val of mods is significant < 0.05 asimmetry (caused either by true heterogeneity or publication bias) is detected
summary(Nema_egger)# pval > 0.05 for mods, good!


# plots
# (Arth_plotegger <- ggplot(data = Arth) +
#     geom_point(aes(x=logRR, y=invse), 
#                shape=16, size=3,color= "#0034A0", alpha= 0.15)+
#     theme_classic()+
#     # scale_y_reverse()+
#     xlab(expression(paste("Effect size (",lnRR^Delta, ")")))+
#     ylab("1/SE")+
#     geom_vline(xintercept=0, size=1, linetype="solid", color="darkgray")+
#     geom_vline(xintercept=Arth_egger[[1]][1], size=1, linetype="dashed", color="black")+
#     theme(axis.text.x = element_text(size=20,angle = 0),
#           axis.text.y = element_text(size=20,angle = 0),
#           plot.title = element_blank(),
#           axis.title.x = element_text(size=20),
#           axis.title.y = element_text(size=20),
#           legend.position="none"))
# 
# (Nema_plotegger <- ggplot(data = Nema) +
#     geom_point(aes(x=logRR, y=invse), 
#                shape=16, size=3,color= "#0034A0", alpha= 0.15)+
#     theme_classic()+
#     # scale_y_reverse()+
#     xlab(expression(paste("Effect size (",lnRR^Delta, ")")))+
#     ylab("1/SE")+
#     geom_vline(xintercept=0, size=1, linetype="solid", color="darkgray")+
#     geom_vline(xintercept=Nema_egger[[1]][1], size=1, linetype="dashed", color="black")+
#     theme(axis.text.x = element_text(size=20,angle = 0),
#           axis.text.y = element_text(size=20,angle = 0),
#           plot.title = element_blank(),
#           axis.title.x = element_text(size=20),
#           axis.title.y = element_text(size=20),
#           legend.position="none"))
# 
# 
# pdf("Results/ArthAbundance_Eggerplot.pdf", width = 7, height = 6)
# funnel(x = Arth$logRR, 
#        sei = Arth$se,
#        yaxis = "seinv",
#        level = c(90, 95, 99),
#        ylab = "1/SE",
#        pch = 16,
#        col = alpha('#0034A0', 0.35),
#        shade = c("gray45", "gray65", "gray90") ,
#        # refline = 0,
#        steps = 5,
#        ylim = c(1, 60),
#        back =  "white",
#        lty = c("blank","dashed"),
#        xlab = expression(paste("Effect size (",lnRR^Delta, ")")), legend = "topright",
#        cex.lab = 1.5,
#        cex.axis = 1.5,
#        cex = 1.5)
# dev.off()

pdf("Results/ArthAbundance_Eggerplot_EffectiveNsize.pdf", width = 7, height = 6)
funnel(Arth$logRR, Arth$VAR, ni = Arth$e_n, yaxis="ni",
       #xlim = c(-3, 3),
       ylab = "Effective sample size",
       pch = 16,
       col = alpha('#0034A0', 0.35),
       shade = c("gray45", "gray65", "gray90") ,
       back =  "white",
       lty = c("blank","dashed"),
       xlab = "Effect size (lnRR)",
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex = 1.5)
dev.off()


pdf("Results/NemaAbundance_Eggerplot_EffectiveNsize.pdf", width = 7, height = 6)
funnel(Nema$logRR, Nema$VAR, ni = Nema$e_n, yaxis="ni",
       #xlim = c(-3, 3),
       ylab = "Effective sample size",
       pch = 16,
       col = alpha('#0034A0', 0.35),
       shade = c("gray45", "gray65", "gray90") ,
       back =  "white",
       lty = c("blank","dashed"),
       xlab = "Effect size (lnRR)",
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex = 1.5)
dev.off()

# Time-lag effects Abundance ####
# Arthropods Abundance
Arth$Source_Year.c <- as.vector(scale(Arth$Source_Year, scale = F))

# Application of Equation 23 from the main manuscript
Arth_timelag <- rma.mv(logRR, 
                       V = ArthV, 
                       mods=~1 + se + Source_Year.c,
                       random=RandomEffects, 
                       method="REML",
                       test="t", # using t dist rather than z 
                       data=Arth)

# Check summary
summary(Arth_timelag)

# Check time-lag effects 
predict.Arth_timelag <- predict(Arth_timelag,
                                newmods= cbind(mean(Arth$se),
                                               seq(min(Arth$Source_Year.c),
                                                   max(Arth$Source_Year.c),
                                                   1)), addx = T)

predict.Arth_timelagdf <- data.frame(Source_Year=seq(min(Arth$Source_Year),
                                     max(Arth$Source_Year),
                                     1),
                     fit=predict.Arth_timelag$pred,
                     upper=predict.Arth_timelag$ci.ub,
                     lower=predict.Arth_timelag$ci.lb,
                     stringsAsFactors=FALSE)

(Arth_timelagplot <- ggplot()+
    geom_hline(yintercept = 0, linetype="longdash", color = "grey", size = 1)+
    geom_point(data=Arth, aes(x=Source_Year, y=logRR, size = se), color="black", alpha = 0.1, shape=16)+
    geom_line(data = predict.Arth_timelagdf, aes(x = Source_Year, y = fit), color="#6D4092", size = 1)+
    geom_ribbon(data = predict.Arth_timelagdf, aes(x = Source_Year, y = fit, ymin = lower, ymax = upper), fill="#6D4092", alpha = 0.1)+
    ylab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    xlab("Publication year")+
    theme_classic(base_size = 30)+
    coord_cartesian(ylim = c(-9,9))+
    guides(size = "none")+
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          legend.title=element_blank()))

# Nematodes
# mean-centring year of publication to help with interpretation, particularly in S4.1.4.3.
Nema$Source_Year.c <- as.vector(scale(Nema$Source_Year, scale = F))

# Application of Equation 23 from the main manuscript
Nema_timelag <- rma.mv(logRR, 
                       V = NemaV, 
                       mods=~1 + se + Source_Year.c,
                       random=RandomEffects, 
                       method="REML",
                       test="t", # using t dist rather than z 
                       data=Nema)

# Check summary
summary(Nema_timelag)

# Check time-lag effects 
predict.Nema_timelag <- predict(Nema_timelag,
                                newmods= cbind(mean(Nema$se),
                                               seq(min(Nema$Source_Year.c),
                                                   max(Nema$Source_Year.c),
                                                   1)), addx = T)

predict.Nema_timelagdf <- data.frame(Source_Year=seq(min(Nema$Source_Year),
                                                     max(Nema$Source_Year),
                                                     1),
                                     fit=predict.Nema_timelag$pred,
                                     upper=predict.Nema_timelag$ci.ub,
                                     lower=predict.Nema_timelag$ci.lb,
                                     stringsAsFactors=FALSE)

(Nema_timelagplot <- ggplot()+
    geom_hline(yintercept = 0, linetype="longdash", color = "grey", size = 1)+
    geom_point(data=Nema, aes(x=Source_Year, y=logRR, size = se), color="black", alpha = 0.1, shape=16)+
    geom_line(data = predict.Nema_timelagdf, aes(x = Source_Year, y = fit), color="#E9AC56", size = 1)+
    geom_ribbon(data = predict.Nema_timelagdf, aes(x = Source_Year, y = fit, ymin = lower, ymax = upper), fill="#E9AC56", alpha = 0.1)+
    ylab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    xlab("Publication year")+
    theme_classic(base_size = 30)+
    guides(size = "none")+
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          legend.title=element_blank()))



# Richness ####
# Read data
datR <- read.csv("./Databases/RichnessDataset_ArthNema_imputed_Env.csv")

#Subset for the Arthropods with Order
ArthR <- datR %>% 
  dplyr::filter(Phylum %in% "Arthropoda") 

#Subset for the Arthropods with Order
NemaR <- datR %>% 
  dplyr::filter(Phylum %in% "Nematoda")


calc.vRRDelta <- function(x) {
  v = matrix((x$Control_SDCorrected[1]^2 / (x$Control_N[1] * x$Control_MeanCorrected[1]^2) + 
                0.5*(x$Control_SDCorrected[1]^4/(x$Control_N[1]^2*x$Control_MeanCorrected[1]^4))), nrow = nrow(x), ncol = nrow(x))
  diag(v) = x$VAR
  return(v)
} 

# Arrange based on Control ID to calculate properly the var-covar matrix
ArthR <- ArthR %>%
  arrange(Control_ID)

NemaR <- NemaR %>%
  arrange(Control_ID)

# For Arthropods
ArthVR <- bldiag(lapply(split(ArthR, ArthR$Control_ID), calc.vRRDelta))
is.positive.definite(ArthVR) # TRUE!

# For Nematodes
NemaVR <- bldiag(lapply(split(NemaR, NemaR$Control_ID), calc.vRRDelta))
is.positive.definite(NemaVR) # TRUE!

# Create a unique identifier for each effect size, which corresponds to the observational id (i.e. residual variance in the model)
ArthR$RowID <-  1:length(ArthR$Source)
NemaR$RowID <-  1:length(NemaR$Source)

# Set the Random Effects
RandomEffects <- list(~1|RowID, ~1|Source)

# Calculate seinv
ArthR$se <- sqrt(ArthR$VAR)
ArthR$invse <- 1/ArthR$se
ArthR$e_n <- with(ArthR, (4*(Control_N*Treatment_N)) / (Control_N + Treatment_N))
ArthR$inv_n_tilda <-  with(ArthR, (Control_N + Treatment_N)/(Control_N*Treatment_N))
ArthR$sqrt_inv_n_tilda <-  with(ArthR, sqrt(inv_n_tilda))

NemaR$se <- sqrt(NemaR$VAR)
NemaR$invse <- 1/NemaR$se
NemaR$e_n <- with(NemaR, (4*(Control_N*Treatment_N)) / (Control_N + Treatment_N))
NemaR$inv_n_tilda <-  with(NemaR, (Control_N + Treatment_N)/(Control_N*Treatment_N))
NemaR$sqrt_inv_n_tilda <-  with(NemaR, sqrt(inv_n_tilda))

# Egger test Richness ####
# Richness arthropods

ArthR_egger <- rma.mv(logRR, ArthVR,
                    mods= ~ 1 + sqrt_inv_n_tilda,
                    random=RandomEffects,
                    method="REML",
                    test = "t", 
                    data=ArthR)

# when p-val of mods is significant < 0.05 asimmetry (caused either by true heterogeneity or publication bias) is detected
summary(ArthR_egger)# pval > 0.05 for mods, good!

# Richness nematodes
NemaR_egger <- rma.mv(logRR, NemaVR,
                     mods= ~ 1 + sqrt_inv_n_tilda,
                     random=RandomEffects,
                     method="REML",
                     test = "t", 
                     data=NemaR)


# when p-val of mods is significant < 0.05 asimmetry (caused either by true heterogeneity or publication bias) is detected
summary(NemaR_egger)# pval > 0.05 for mods, good!

# plots
pdf("Results/ArthRichness_Eggerplot_EffectiveNsize.pdf", width = 7, height = 6)
funnel(ArthR$logRR, ArthR$VAR, ni = ArthR$e_n, yaxis="ni",
       #xlim = c(-3, 3),
       ylab = "Effective sample size",
       pch = 16,
       col = alpha('#0034A0', 0.35),
       shade = c("gray45", "gray65", "gray90") ,
       back =  "white",
       lty = c("blank","dashed"),
       xlab = "Effect size (lnRR)",
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex = 1.5)
dev.off()

pdf("Results/NemaRichness_Eggerplot_EffectiveNsize.pdf", width = 7, height = 6)
funnel(NemaR$logRR, NemaR$VAR, ni = NemaR$e_n, yaxis="ni",
       #xlim = c(-3, 3),
       ylab = "Effective sample size",
       pch = 16,
       col = alpha('#0034A0', 0.35),
       shade = c("gray45", "gray65", "gray90") ,
       back =  "white",
       lty = c("blank","dashed"),
       xlab = "Effect size (lnRR)",
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex = 1.5)
dev.off()


# Time-lag effects Richness ####
# Arthropods Richness
ArthR$Source_Year.c <- as.vector(scale(ArthR$Source_Year, scale = F))

# Application of Equation 23 from the main manuscript
ArthR_timelag <- rma.mv(logRR, 
                       V = ArthVR, 
                       mods=~1 + se + Source_Year.c,
                       random=RandomEffects, 
                       method="REML",
                       test="t", # using t dist rather than z 
                       data=ArthR)

# Check summary
summary(ArthR_timelag)

# Check time-lag effects 
predict.ArthR_timelag <- predict(ArthR_timelag,
                                newmods= cbind(mean(ArthR$se),
                                               seq(min(ArthR$Source_Year.c),
                                                   max(ArthR$Source_Year.c),
                                                   1)), addx = T)

predict.ArthR_timelagdf <- data.frame(Source_Year=seq(min(ArthR$Source_Year),
                                           max(ArthR$Source_Year),
                                           1),
                           fit=predict.ArthR_timelag$pred,
                           upper=predict.ArthR_timelag$ci.ub,
                           lower=predict.ArthR_timelag$ci.lb,
                           stringsAsFactors=FALSE)

(ArthR_timelagplot <- ggplot()+
    geom_hline(yintercept = 0, linetype="longdash", color = "grey", size = 1)+
    geom_point(data=ArthR, aes(x=Source_Year, y=logRR, size = se), color="black", alpha = 0.1, shape=16)+
    geom_line(data = predict.ArthR_timelagdf, aes(x = Source_Year, y = fit), color="#6D4092", size = 1)+
    geom_ribbon(data = predict.ArthR_timelagdf, aes(x = Source_Year, y = fit, ymin = lower, ymax = upper), fill="#6D4092", alpha = 0.1)+
    ylab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    xlab("Publication year")+
    coord_cartesian(ylim = c(-3,3))+
    theme_classic(base_size = 30)+
    guides(size = "none")+
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          legend.title=element_blank()))

# Nematodes
NemaR$Source_Year.c <- as.vector(scale(NemaR$Source_Year, scale = F))

# Application of Equation 23 from the main manuscript
NemaR_timelag <- rma.mv(logRR, 
                       V = NemaVR, 
                       mods=~1 + se + Source_Year.c,
                       random=RandomEffects, 
                       method="REML",
                       test="t", # using t dist rather than z 
                       data=NemaR)

# Check summary
summary(NemaR_timelag)

# Check time-lag effects 
predict.NemaR_timelag <- predict(NemaR_timelag,
                                newmods= cbind(mean(NemaR$se),
                                               seq(min(NemaR$Source_Year.c),
                                                   max(NemaR$Source_Year.c),
                                                   1)), addx = T)

newdat_NemaR <- data.frame(Source_Year=seq(min(NemaR$Source_Year),
                                          max(NemaR$Source_Year),
                                          1),
                          fit=predict.NemaR_timelag$pred,
                          upper=predict.NemaR_timelag$ci.ub,
                          lower=predict.NemaR_timelag$ci.lb,
                          stringsAsFactors=FALSE)

(NemaR_timelagplot <- ggplot()+
  geom_hline(yintercept = 0, linetype="longdash", color = "grey", size = 1)+
  geom_point(data=NemaR, aes(x=Source_Year, y=logRR, size = se), color="black", alpha = 0.1, shape=16)+
  geom_line(data = newdat_NemaR, aes(x = Source_Year, y = fit), color="#E9AC56", size = 1)+
  geom_ribbon(data = newdat_NemaR, aes(x = Source_Year, y = fit, ymin = lower, ymax = upper), fill="#E9AC56", alpha = 0.1)+
  ylab(expression(paste("Effect size (",lnRR^Delta, ")")))+
  xlab("Publication year")+
  coord_cartesian(xlim = c(2005,2021), ylim = c(-0.75,0.75))+
  scale_y_continuous(breaks = c(-0.5,0,0.5))+
  theme_classic(base_size = 30)+
  guides(size = "none")+
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        legend.title=element_blank()))

# Save multiplots ####
# Multiplot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

pdf("Results/Abundance_Timelagplots.pdf",  width = 15, height = 7)
multiplot(Arth_timelagplot, Nema_timelagplot,cols=2)
dev.off()

pdf("Results/Richness_Timelagplots.pdf",  width = 15, height = 7)
multiplot(ArthR_timelagplot, NemaR_timelagplot,cols=2)
dev.off()

