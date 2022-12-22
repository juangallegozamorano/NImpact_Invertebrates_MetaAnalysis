# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com


############## Load libraries and dataset ##############

library(rotl) # Interface to the 'Open Tree of Life' API, CRAN v3.0.11
library(tidyverse) # Easily Install and Load the 'Tidyverse', CRAN v1.3.1
library(ape) # Analyses of Phylogenetics and Evolution, CRAN v5.5
library(phytools) # Phylogenetic Tools for Comparative Biology (and Other Things), CRAN v0.7-47

load("./Results/Arthropods/Phylogeny/arth_tree_Orders_random_branchlength.Rdata")

# Plot tree
plot(phylo_branchOrders, label.offset =.01, no.margin = TRUE)

## Creating a vector of colours based on the Metamorphosis type per tip label
my_colours <- case_when(str_detect(phylo_branchOrders$tip.label, "(?i)Coleoptera|Diptera|Hymenoptera|Lepidoptera|Mecoptera|Neuroptera|Raphidioptera")~"#e7464c",
                        str_detect(phylo_branchOrders$tip.label, "(?i)Araneae|Blattodea|Dermaptera|Diplura|Geophilomorpha|Lithobiomorpha|Embioptera|Ephemeroptera|Hemiptera|Isopoda|Mesostigmata|Opiliones|Odonata|Orthoptera|Phasmatodea|Pseudoscorpiones|Psocoptera|Scolopendromorpha|Sarcoptiformes|Tetramerocerata|Thysanoptera|Trombidiformes|Polyxenida|Julida")~"#6F3F97",
                        str_detect(phylo_branchOrders$tip.label, "(?i)Entomobryomorpha|Neelipleona|Poduromorpha|Symphypleona|Zygentoma|Protura")~"#E9AC57")


# Tree with color tips
plot.phylo(phylo_branchOrders,
           # show.node.label = T,
           tip.color = my_colours, label.offset =.01)

pdf("Results/Arthropods/Plots/Phylogeny_MetamorphosisColor.pdf",  width = 8, height = 9)
plot.phylo(phylo_branchOrders,
     tip.color = my_colours, label.offset =.01)
axisPhylo(backward = F, col = "#686969", col.axis = "#686969")
dev.off()


# Test phylogenetic signal
## Creating a vector of Metamorphosis type per tip label (Trait)
trait <- case_when(str_detect(phylo_branchOrders$tip.label, "(?i)Coleoptera|Diptera|Hymenoptera|Lepidoptera|Mecoptera|Neuroptera|Raphidioptera")~"Complete",
                   str_detect(phylo_branchOrders$tip.label, "(?i)Araneae|Blattodea|Dermaptera|Diplura|Geophilomorpha|Lithobiomorpha|Embioptera|Ephemeroptera|Hemiptera|Isopoda|Mesostigmata|Opiliones|Odonata|Orthoptera|Phasmatodea|Pseudoscorpiones|Psocoptera|Scolopendromorpha|Sarcoptiformes|Tetramerocerata|Thysanoptera|Trombidiformes|Polyxenida|Julida")~"Incomplete-Gradual",
                   str_detect(phylo_branchOrders$tip.label, "(?i)Entomobryomorpha|Neelipleona|Poduromorpha|Symphypleona|Zygentoma|Protura")~"None")


names(trait)<-phylo_branchOrders$tip.label

traitNum <- case_when(str_detect(phylo_branchOrders$tip.label, "(?i)Coleoptera|Diptera|Hymenoptera|Lepidoptera|Mecoptera|Neuroptera|Raphidioptera")~1,
                      str_detect(phylo_branchOrders$tip.label, "(?i)Araneae|Blattodea|Dermaptera|Diplura|Geophilomorpha|Lithobiomorpha|Embioptera|Ephemeroptera|Hemiptera|Isopoda|Mesostigmata|Opiliones|Odonata|Orthoptera|Phasmatodea|Pseudoscorpiones|Psocoptera|Scolopendromorpha|Sarcoptiformes|Tetramerocerata|Thysanoptera|Trombidiformes|Polyxenida|Julida")~2,
                      str_detect(phylo_branchOrders$tip.label, "(?i)Entomobryomorpha|Neelipleona|Poduromorpha|Symphypleona|Zygentoma|Protura")~3)
names(traitNum)<-phylo_branchOrders$tip.label

# Phylosig Blomberg’s K
K.metamor <- phylosig(phylo_branchOrders, traitNum, method = "K", test = TRUE, nsim = 1000)
print(K.metamor)
plot(K.metamor)
