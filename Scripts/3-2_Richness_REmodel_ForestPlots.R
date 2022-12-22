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
datR <- read.csv("./Databases/RichnessDataset_ArthNema_imputed_Env.csv")

#Subset for the Arthropods
ArthR <- datR %>% 
  dplyr::filter(Phylum %in% "Arthropoda") 

#Subset for the Nematodes
NemaR <- datR %>% 
  dplyr::filter(Phylum %in% "Nematoda")

# Function for var-covar matrix with Delta correction
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
ArthV <- bldiag(lapply(split(ArthR, ArthR$Control_ID), calc.vRRDelta))
is.positive.definite(ArthV) # TRUE!

# For Nematodes
NemaV <- bldiag(lapply(split(NemaR, NemaR$Control_ID), calc.vRRDelta))
is.positive.definite(NemaV) # TRUE!

# Create a unique identifier for each effect size, which corresponds to the observational id (i.e. residual variance in the model)
ArthR$RowID <-  1:length(ArthR$Source)
NemaR$RowID <-  1:length(NemaR$Source)

# Set the Random Effects
RandomEffects <- list(~1|RowID, ~1|Source)

# RE model for Arthropods   ##############
ArthR_REmodel <- rma.mv(logRR,
                       V = ArthV, 
                       data = ArthR, 
                       random = RandomEffects)


# RE model for Nematodes   ##############
NemaR_REmodel <- rma.mv(logRR,
                       V = NemaV, 
                       data = NemaR, 
                       random = RandomEffects)

# Setting up the basic plot
# Arthropoda
forestArthR <- ArthR %>% 
  dplyr::select(logRR, VAR) %>% 
  arrange(logRR)

forestArthR$plottingID <- 1:nrow(forestArthR)

# Nematoda
forestNemaR <- NemaR %>% 
  dplyr::select(logRR, VAR) %>% 
  arrange(logRR)

forestNemaR$plottingID <- 1:nrow(forestNemaR)


#Get standard errors from variances
forestArthR$se <- sqrt(forestArthR$VAR)

forestNemaR$se <- sqrt(forestNemaR$VAR)

#Calculate 95% CI values
forestArthR$lowerci <- (-1.96*forestArthR$se)+forestArthR$logRR
forestArthR$upperci <- (1.96*forestArthR$se)+forestArthR$logRR

forestNemaR$lowerci <- (-1.96*forestNemaR$se)+forestNemaR$logRR
forestNemaR$upperci <- (1.96*forestNemaR$se)+forestNemaR$logRR

# Get average values from the model
ArthMeanR <- as.data.frame(cbind(ArthR_REmodel$b, 
                                 ArthR_REmodel$se^2,
                                 -1,
                                 ArthR_REmodel$se,
                                 ArthR_REmodel$ci.lb, 
                                 ArthR_REmodel$ci.ub))
colnames(ArthMeanR) <- c("logRR", "VAR","plottingID","se","lowerci", "upperci")

ArthMeanRpreds <-  predict(ArthR_REmodel)
ArthMeanR$pi.lb <- ArthMeanRpreds$pi.lb
ArthMeanR$pi.ub <- ArthMeanRpreds$pi.ub

# percentage change
((exp(ArthMeanRpreds$pred))-1)*100
((exp(ArthMeanRpreds$ci.lb))-1)*100
((exp(ArthMeanRpreds$ci.ub))-1)*100

((exp(ArthMeanRpreds$pi.lb))-1)*100
((exp(ArthMeanRpreds$pi.ub))-1)*100

# Get average values from the model
NemaMeanR <- as.data.frame(cbind(NemaR_REmodel$b, 
                                 NemaR_REmodel$se^2,
                                 -1,
                                 NemaR_REmodel$se,
                                 NemaR_REmodel$ci.lb, 
                                 NemaR_REmodel$ci.ub))
colnames(NemaMeanR) <- c("logRR", "VAR","plottingID","se","lowerci", "upperci")

NemaMeanRpreds <-  predict(NemaR_REmodel)
NemaMeanR$pi.lb <- NemaMeanRpreds$pi.lb
NemaMeanR$pi.ub <- NemaMeanRpreds$pi.ub

# percentage change
((exp(NemaMeanRpreds$pred))-1)*100
((exp(NemaMeanRpreds$ci.lb))-1)*100
((exp(NemaMeanRpreds$ci.ub))-1)*100


((exp(NemaMeanRpreds$pi.lb))-1)*100
((exp(NemaMeanRpreds$pi.ub))-1)*100

# Forest for Arthropods ##############
(arthR_forestplot <- ggplot(data=forestArthR, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci))+ 
    
    # add error bar below
    geom_linerange(size=1,color="gray")+
    
    #this adds the effect sizes to the plot
    geom_point(size = 1.5) +
    
    #adding a vertical line at the effect = 0 mark
    geom_vline(xintercept=0, color="black", linetype="longdash", size = 1)+
    
    # Add the mean estimate from the model and the line of the model
    # geom_linerange(data = ArthMeanR, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci), size=1, color="#6D4092")+
    geom_pointrange(data = ArthMeanR, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci), size=1, color="#6D4092",shape=18)+
    
    geom_point(data = ArthMeanR, aes(y=plottingID, x=logRR), color="#6D4092", size=7,shape=18)+
    geom_vline(xintercept=ArthMeanR$logRR, color="#6D4092", linetype="solid", size = 1)+
    geom_vline(xintercept=ArthMeanR$lowerci, color="#6D4092", linetype="dashed")+
    geom_vline(xintercept=ArthMeanR$upperci, color="#6D4092", linetype="dashed")+
    
    geom_vline(xintercept=ArthMeanR$pi.lb, color="#6D4092", linetype="dotted", size = 2)+
    geom_vline(xintercept=ArthMeanR$pi.ub, color="#6D4092", linetype="dotted", size = 2)+
    
    # Adjust scale
    scale_x_continuous(limits=c(-3,3))+
    xlab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    #thematic stuff
    theme_classic(base_size =40)+ 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y = element_blank(),
          axis.line.y = element_blank()))

# Zoom plot
(arthR_forestplotzoom <- ggplot(data=forestArthR, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci))+ 
    
    # add error bar below
    geom_linerange(size=1,color="gray")+
    
    #this adds the effect sizes to the plot
    geom_point(size = 1.5) +
    
    #adding a vertical line at the effect = 0 mark
    geom_vline(xintercept=0, color="black", linetype="longdash", size = 1)+
    
    # Add the mean estimate from the model and the line of the model
    geom_pointrange(data = ArthMeanR, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci), size=1, color="#6D4092",shape=18)+
    
    geom_point(data = ArthMeanR, aes(y=plottingID, x=logRR), color="#6D4092", size=7,shape=18)+
    geom_vline(xintercept=ArthMeanR$logRR, color="#6D4092", linetype="solid", size = 1)+
    geom_vline(xintercept=ArthMeanR$lowerci, color="#6D4092", linetype="dashed", size = 1.5)+
    geom_vline(xintercept=ArthMeanR$upperci, color="#6D4092", linetype="dashed", size = 1.5)+
    
    geom_vline(xintercept=ArthMeanR$pi.lb, color="#6D4092", linetype="dotted", size = 2)+
    geom_vline(xintercept=ArthMeanR$pi.ub, color="#6D4092", linetype="dotted", size = 2)+
    
    coord_cartesian(xlim=c(-0.5,0.5), ylim = c(-2,0)) +
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

pdf("Results/Richness_Arthzoom_Forestplot.pdf", width = 4, height = 3, bg = "transparent")
arthR_forestplotzoom
dev.off()


# Forest for Nematodes   ##############
(nemaR_forestplot <- ggplot(data=forestNemaR, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci))+ 
    
    # add error bar below
    geom_linerange(size=1,color="gray")+
    
    #this adds the effect sizes to the plot
    geom_point(size = 1.5) +
    
    #adding a vertical line at the effect = 0 mark
    geom_vline(xintercept=0, color="black", linetype="longdash", size = 1)+
    
    # Add the mean estimate from the model and the line of the model
    geom_pointrange(data = NemaMeanR, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci), size=1, color="#E9AC56",shape=18)+
    
    geom_point(data = NemaMeanR, aes(y=plottingID, x=logRR), color="#E9AC56", size=7,shape=18)+
    geom_vline(xintercept=NemaMeanR$logRR, color="#E9AC56", linetype="solid", size = 1)+
    geom_vline(xintercept=NemaMeanR$lowerci, color="#E9AC56", linetype="dashed")+
    geom_vline(xintercept=NemaMeanR$upperci, color="#E9AC56", linetype="dashed")+
    
    geom_vline(xintercept=NemaMeanR$pi.lb, color="#E9AC56", linetype="dotted", size = 2)+
    geom_vline(xintercept=NemaMeanR$pi.ub, color="#E9AC56", linetype="dotted", size = 2)+
    
    
    # scale_x_continuous(limits=c(-1,1))+
    xlab(expression(paste("Effect size (",lnRR^Delta, ")")))+
    #thematic stuff
    theme_classic(base_size =40)+ 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y = element_blank(),
          axis.line.y = element_blank()))

# Zoom plot
(nemaR_forestplotzoom <- ggplot(data=forestNemaR, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci))+ 
    
    # add error bar below
    geom_linerange(size=1,color="gray")+
    
    #this adds the effect sizes to the plot
    geom_point(size = 1.5) +
    
    #adding a vertical line at the effect = 0 mark
    geom_vline(xintercept=0, color="black", linetype="longdash", size = 1)+
    
    # Add the mean estimate from the model and the line of the model
    geom_pointrange(data = NemaMeanR, aes(y=plottingID, x=logRR, xmin = lowerci, xmax = upperci), size=1, color="#E9AC56",shape=18)+
    
    geom_point(data = NemaMeanR, aes(y=plottingID, x=logRR), color="#E9AC56", size=7,shape=18)+
    geom_vline(xintercept=NemaMeanR$logRR, color="#E9AC56", linetype="solid", size = 1)+
    geom_vline(xintercept=NemaMeanR$lowerci, color="#E9AC56", linetype="dashed", size = 1.5)+
    geom_vline(xintercept=NemaMeanR$upperci, color="#E9AC56", linetype="dashed", size = 1.5)+
    
    geom_vline(xintercept=NemaMeanR$pi.lb, color="#E9AC56", linetype="dotted", size = 2)+
    geom_vline(xintercept=NemaMeanR$pi.ub, color="#E9AC56", linetype="dotted", size = 2)+
    
    
    coord_cartesian(xlim=c(-0.5,0.5), ylim = c(-2,0)) +
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

pdf("Results/Richness_Nemazoom_Forestplot.pdf", width = 4, height = 3, bg = "transparent")
nemaR_forestplotzoom
dev.off()

# Multiplot ##############
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

pdf("Results/Richness_Forestplots.pdf",  width = 20, height = 10)
multiplot(arthR_forestplot, nemaR_forestplot,cols=2)
dev.off()

