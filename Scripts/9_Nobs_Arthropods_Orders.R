# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com


############## Load libraries and dataset ##############

library(tidyverse)

# Read data
ArthO <- read.csv("./Databases/AbundanceDataset_Arth_CompleteOrders_imputed.csv")

# Number of orders per metamorphosis class ####
metaN <- ArthO %>% 
  group_by(Metamorphosis, Order) %>% 
  tally(name = "Nobs") %>% 
  arrange(Metamorphosis, -Nobs) %>% # Arrange so bigger columns are first
  ungroup() %>% 
  dplyr::mutate(id = row_number(),
                logNobs = log10(Nobs+1),
                Metamorphosis = as.factor(Metamorphosis)) # Add id

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
to_add <- data.frame(matrix(NA, empty_bar*nlevels(metaN$Metamorphosis), ncol(metaN)) )
colnames(to_add) <- colnames(metaN)
to_add$Metamorphosis <- rep(levels(metaN$Metamorphosis), each=empty_bar)
metaN <- rbind(metaN, to_add)
metaN <- metaN %>% arrange(Metamorphosis, -Nobs)
metaN$id <- seq(1, nrow(metaN))


# ----- This section prepare a dataframe for labels ---- #
# Get the name and the y position of each label
label_data <- metaN

# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)


# prepare a data frame for base lines
base_data <- metaN %>% 
  dplyr::group_by(Metamorphosis) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  dplyr::mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]



(metaNplot <- ggplot() +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0.477, xend = start, yend = 0.477), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1.477, xend = start, yend = 1.477), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 2, xend = start, yend = 2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 2.477, xend = start, yend = 2.477), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(metaN$id),5), y = c(0.5, 1, 1.5, 2,2.5), label = c("3", "10", "30", "100","300"),
           color="grey", size=4 , angle=0, fontface="bold", hjust=1.25) +
  
  
  # This add the bars
  geom_bar(data=metaN, aes(x=id, y=logNobs, fill=as.factor(Metamorphosis)), stat="identity", alpha = 0.8) +
  
  scale_fill_manual(values = c("#e7464c",
                               "#6F3F97",
                               "#FDBD63"), 
                    guide = guide_legend(override.aes = list(alpha = 1)))+
  
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-3,5) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=12),
    # plot.margin = unit(c(0.5,-2,0.5,-2), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=id, y=logNobs+0.5, label=Order, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=4, angle= label_data$angle, inherit.aes = FALSE )+
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -0.2, xend = end, yend = -0.2), 
               colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )+
    annotate("text", x = 0, y = -3, label = "Nº observations", size = 5))

pdf("Results/Arthropods/Plots/Nobs_Metamorphosis.pdf",  width = 8, height = 7)
metaNplot
dev.off()
