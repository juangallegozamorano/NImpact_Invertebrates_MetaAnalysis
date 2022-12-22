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

# load necessary functions
source("Scripts/00_Functions.R")

# Read data
dat <- read.csv("./Databases/AbundanceDataset_ArthNema_imputed_Env.csv")

#Subset for the Arthropods
Arth <- dat %>% 
  dplyr::filter(Phylum %in% "Arthropoda") 

#Subset for the Nematodes
Nema <- dat %>% 
  dplyr::filter(Phylum %in% "Nematoda")

# Function for var-covar matrix with Delta correction
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

# Var-covar Arthropods
ArthV <- bldiag(lapply(split(Arth, Arth$Control_ID), calc.vRRDelta))
is.positive.definite(ArthV) # TRUE!
# ArthV <- PDfunc(ArthV) # If negative PDfunc make the V positive

# Var-covar Nematodes
NemaV <- bldiag(lapply(split(Nema, Nema$Control_ID), calc.vRRDelta))
is.positive.definite(NemaV) # TRUE!

# Create a unique identifier for each effect size, which corresponds to the observational id (i.e. residual variance in the model)
Arth$RowID <-  1:length(Arth$Source)
Nema$RowID <-  1:length(Nema$Source)

# Set the Random Effects
RandomEffects <- list(~1|RowID, ~1|Source)

# RE model for Arthropods  ##############
Arth_REmodel <- rma.mv(logRR,
                       V = ArthV, 
                       data = Arth, 
                       random = RandomEffects)

summary(Arth_REmodel)

# RE model for Nematoda  ##############
Nema_REmodel <- rma.mv(logRR,
                       V = NemaV, 
                       data = Nema, 
                       random = RandomEffects)
summary(Nema_REmodel)

# Setting up the basic plot
# Arthropoda
forestArth <- Arth %>% 
  dplyr::select(logRR, VAR) %>% 
  arrange(logRR)

forestArth$plottingID <- 1:nrow(forestArth)

# Nematoda
forestNema <- Nema %>% 
  dplyr::select(logRR, VAR) %>% 
  arrange(logRR)

forestNema$plottingID <- 1:nrow(forestNema)


#Get standard errors from variances
forestArth$se <- sqrt(forestArth$VAR)

forestNema$se <- sqrt(forestNema$VAR)

#Calculate 95% CI values
forestArth$lowerci <- (-1.96*forestArth$se)+forestArth$logRR
forestArth$upperci <- (1.96*forestArth$se)+forestArth$logRR

forestNema$lowerci <- (-1.96*forestNema$se)+forestNema$logRR
forestNema$upperci <- (1.96*forestNema$se)+forestNema$logRR

# Get average values from the model
ArthMean <- as.data.frame(cbind(Arth_REmodel$b, 
                                Arth_REmodel$se^2,
                                -50,
                                Arth_REmodel$se,
                                Arth_REmodel$ci.lb, 
                                Arth_REmodel$ci.ub))
colnames(ArthMean) <- c("logRR", "VAR","plottingID","se","lowerci", "upperci")

ArthMeanpreds <-  predict(Arth_REmodel)
ArthMean$pi.lb <- ArthMeanpreds$pi.lb
ArthMean$pi.ub <- ArthMeanpreds$pi.ub

# percentage change
((exp(ArthMeanpreds$pred))-1)*100
((exp(ArthMeanpreds$ci.lb))-1)*100
((exp(ArthMeanpreds$ci.ub))-1)*100

((exp(ArthMeanpreds$pi.lb))-1)*100
((exp(ArthMeanpreds$pi.ub))-1)*100

# Get average values from the model
NemaMean <- as.data.frame(cbind(Nema_REmodel$b, 
                                Nema_REmodel$se^2,
                                -50,
                                Nema_REmodel$se,
                                Nema_REmodel$ci.lb, 
                                Nema_REmodel$ci.ub))
colnames(NemaMean) <- c("logRR", "VAR","plottingID","se","lowerci", "upperci")

NemaMeanpreds <-  predict(Nema_REmodel)
NemaMean$pi.lb <- NemaMeanpreds$pi.lb
NemaMean$pi.ub <- NemaMeanpreds$pi.ub

# percentage change
((exp(NemaMeanpreds$pred))-1)*100
((exp(NemaMeanpreds$ci.lb))-1)*100
((exp(NemaMeanpreds$ci.ub))-1)*100

((exp(NemaMeanpreds$pi.lb))-1)*100
((exp(NemaMeanpreds$pi.ub))-1)*100

# Forest plot for Arthropods   ##############
(arth_forestplot <- ggplot(data=forestArth, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci))+ 
    
    # add error bar below
    geom_linerange(size=1,color="gray")+
    
    #this adds the effect sizes to the plot
    geom_point(size = 1.5) +
    
    #adding a vertical line at the effect = 0 mark
    geom_vline(xintercept=0, color="black", linetype="longdash", size = 1)+
    
    # Add the mean estimate from the model and the line of the model
    geom_pointrange(data = ArthMean, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci), size=1, color="#6D4092",shape=18)+
    
    geom_vline(xintercept=ArthMean$logRR, color="#6D4092", linetype="solid", size = 1)+
    geom_vline(xintercept=ArthMean$lowerci, color="#6D4092", linetype="dashed")+
    geom_vline(xintercept=ArthMean$upperci, color="#6D4092", linetype="dashed")+
    
    geom_vline(xintercept=ArthMean$pi.lb, color="#6D4092", linetype="dotted", size = 2)+
    geom_vline(xintercept=ArthMean$pi.ub, color="#6D4092", linetype="dotted", size = 2)+
    
    # Add scale of the x axis
    coord_cartesian(xlim=c(-10,10)) +
    xlab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    #thematic stuff
    theme_classic(base_size =40)+ 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y = element_blank(),
          axis.line.y = element_blank()))
arth_forestplot

# Zoom plot 
(arth_forestplotzoom <- ggplot(data=forestArth, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci))+ 
    
    # add error bar below
    geom_linerange(size=1,color="gray")+
    
    #this adds the effect sizes to the plot
    geom_point(size = 1.5) +
    
    #adding a vertical line at the effect = 0 mark
    geom_vline(xintercept=0, color="black", linetype="longdash", size = 1)+
    
    # Add the mean estimate from the model and the line of the model
    geom_pointrange(data = ArthMean, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci), size=1, color="#6D4092",shape=18)+
    
    geom_vline(xintercept=ArthMean$logRR, color="#6D4092", linetype="solid", size = 1)+
    geom_vline(xintercept=ArthMean$lowerci, color="#6D4092", linetype="dashed", size = 1.5)+
    geom_vline(xintercept=ArthMean$upperci, color="#6D4092", linetype="dashed", size = 1.5)+
    
    geom_vline(xintercept=ArthMean$pi.lb, color="#6D4092", linetype="dotted", size = 2)+
    geom_vline(xintercept=ArthMean$pi.ub, color="#6D4092", linetype="dotted", size = 2)+
    
    # Add scale of the x axis
    coord_cartesian(xlim=c(-0.5,0.5), ylim = c(-60,-40)) +
    scale_x_continuous(breaks=c(-0.5,0,0.5))+
    xlab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    #thematic stuff
    theme_classic(base_size =40)+ 
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=2),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          axis.line.x = element_blank(),
          axis.line.y = element_blank()))

pdf("Results/Abundance_Arthzoom_Forestplot.pdf", width = 4, height = 3, bg = "transparent")
arth_forestplotzoom
dev.off()


# Forest for Nematodes  ##############
(nema_forestplot <- ggplot(data=forestNema, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci))+ 
  
    # add error bar below
    geom_linerange(size=1,color="gray")+
    
    #this adds the effect sizes to the plot
    geom_point(size = 1.5) +
    
    #adding a vertical line at the effect = 0 mark
    geom_vline(xintercept=0, color="black", linetype="longdash", size = 1)+
    xlab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    # Add the mean estimate from the model and the line of the model
    geom_pointrange(data = NemaMean, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci), size=1, color="#E9AC56",shape=18)+
    
    geom_point(data = NemaMean, aes(y=plottingID, x=logRR), color="#E9AC56", size=7,shape=18)+
    geom_vline(xintercept=NemaMean$logRR, color="#E9AC56", linetype="solid", size = 1)+
    geom_vline(xintercept=NemaMean$lowerci, color="#E9AC56", linetype="dashed")+
    geom_vline(xintercept=NemaMean$upperci, color="#E9AC56", linetype="dashed")+
    
    geom_vline(xintercept=NemaMean$pi.lb, color="#E9AC56", linetype="dotted", size = 2)+
    geom_vline(xintercept=NemaMean$pi.ub, color="#E9AC56", linetype="dotted", size = 2)+
  
  
    coord_cartesian(xlim=c(-10,10)) +
    
  #thematic stuff
  theme_classic(base_size =40)+ 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line.y = element_blank()))
nema_forestplot

(nema_forestplotzoom <- ggplot(data=forestNema, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci))+ 
    
    # add error bar below
    geom_linerange(size=1,color="gray")+
    
    #this adds the effect sizes to the plot
    geom_point(size = 1.5) +
    
    #adding a vertical line at the effect = 0 mark
    geom_vline(xintercept=0, color="black", linetype="longdash", size = 1)+
    xlab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    # Add the mean estimate from the model and the line of the model
    geom_pointrange(data = NemaMean, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci), size=1, color="#E9AC56",shape=18)+
    
    geom_point(data = NemaMean, aes(y=plottingID, x=logRR), color="#E9AC56", size=7,shape=18)+
    geom_vline(xintercept=NemaMean$logRR, color="#E9AC56", linetype="solid", size = 1)+
    geom_vline(xintercept=NemaMean$lowerci, color="#E9AC56", linetype="dashed", size = 1.5)+
    geom_vline(xintercept=NemaMean$upperci, color="#E9AC56", linetype="dashed", size = 1.5)+
    
    geom_vline(xintercept=NemaMean$pi.lb, color="#E9AC56", linetype="dotted", size = 2)+
    geom_vline(xintercept=NemaMean$pi.ub, color="#E9AC56", linetype="dotted", size = 2)+
    
    
    coord_cartesian(xlim=c(-0.5,0.5), ylim = c(-60,-40)) +
    scale_x_continuous(breaks=c(-0.5,0,0.5))+
    
    #thematic stuff
    theme_classic(base_size =40)+ 
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=2),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          axis.line.x = element_blank(),
          axis.line.y = element_blank()))

pdf("Results/Abundance_Nemazoom_Forestplot.pdf", width = 4, height = 3, bg = "transparent")
nema_forestplotzoom
dev.off()

# Multiplot  ##############
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

pdf("Results/Abundance_Forestplots.pdf", width = 20, height = 10)
multiplot(arth_forestplot, nema_forestplot,cols=2)
dev.off()

