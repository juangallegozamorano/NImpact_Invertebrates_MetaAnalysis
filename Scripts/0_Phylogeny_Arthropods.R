# Author: Juan Gallego-Zamorano
# Date: 22-12-2022
# Department of Environmental Science, Radboud University Nijmegen
# e-mail: j.gallego.zamorano@gmail.com

# This code create a phylogeny of Arthropods at the family and order level
# However, only the one from Order was used in the analysis

############## Load libraries and dataset ##############

library(rotl)
library(tidyverse)
library(ape)

# Read data
dat <- read.csv("./Databases/AbundanceDataset_ArthNema_imputed_Env.csv")

#Subset for the Order
Arth <- dat %>% 
  filter(Phylum %in% "Arthropoda")

# obtaining dataframe listing the Open Tree identifiers potentially matching our list of species (be aware that this will take a few minutes, and you can load the data below)
ArthC <- Arth %>% 
  drop_na(Class)

ArthO <- Arth %>% 
  drop_na(Order)

# ArthF <- Arth %>% 
#   drop_na(Family)

classes <- sort(unique(ArthC$Class))
orders <- sort(unique(ArthO$Order))

Arthtaxa <- tnrs_match_names(names = classes)
ArthtaxaO <- tnrs_match_names(names = orders)# It seems doable!

# according to the `approximate_match` column, there might be 
# 0 typo in the species list
nrow(Arthtaxa[Arthtaxa$approximate_match==TRUE,])
Arthtaxa[Arthtaxa$approximate_match==TRUE,]  #### zero species to fix

nrow(ArthtaxaO[ArthtaxaO$approximate_match==TRUE,])
ArthtaxaO[ArthtaxaO$approximate_match==TRUE,]  #### zero species to fix

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
Arthtaxa[Arthtaxa$number_matches != 1,]  #0!

ArthtaxaO[ArthtaxaO$number_matches != 1,]  #3

inspect(ArthtaxaO, taxon_name = "Diplura")
inspect(ArthtaxaO, taxon_name = "Isopoda")
inspect(ArthtaxaO, taxon_name = "Sarcoptiformes")

# Check the rank of the taxax
info <- taxonomy_taxon_info(87288)
rotl::tax_rank(info)

# Get lineage of a given taxa
tax <- tnrs_match_names("Oribatida")
tax_lineage(taxonomy_taxon_info(ott_id(tax), include_lineage = TRUE))

# check synonyms and change name accordingly
Arthtaxa[Arthtaxa$is_synonym==TRUE,] # 0!
ArthtaxaO[ArthtaxaO$is_synonym==TRUE,] # 0!

# saving the taxonomic data created on the 26th of October to speed the process in the future
save(Arthtaxa, file = "./Results/Arthropods/Phylogeny/Arthropod_taxa_Classes_Open_Tree_of_Life.RData")
save(ArthtaxaO, file = "./Results/Arthropods/Phylogeny/Arthropod_taxa_Orders_Open_Tree_of_Life.RData")

# Getting the tree corresponding to our taxa
tree <- tol_induced_subtree(ott_ids = Arthtaxa$ott_id, label = "name")

plot(tree, no.margin = TRUE)

treeO <- tol_induced_subtree(ott_ids = ArthtaxaO$ott_id, label = "name")

plot(treeO, no.margin = TRUE)

# Check if all the orders are in the Synthetic Tree
in_tree <- is_in_tree(ott_id(ArthtaxaO))
in_tree # All TRUE! 

##############################################################
# Dealing with polytomies
##############################################################

# we can check for the existence of polytomies by running the 
# following code. If polytomies exist, the output will be 
# `FALSE`, and vice versa.

is.binary(tree) # there are some polytomies
is.binary(treeO) # there are some polytomies

# to take care of these polytomies, we are going to use a 
# randomization approach
set.seed(111) #making it replicable, at least for this version of R (i.e. v.4.0.2)
tree_random <- multi2di(tree,random=TRUE)
tree_randomOrders <- multi2di(treeO,random=TRUE)

is.binary(tree_random)
is.binary(tree_randomOrders)

##############################################################
# Final checks
##############################################################

# exploring whether our tree covers all the orders we wanted 
# it to include, and making sure that the orders names in our 
# database match those in the tree. We use the following code.

# Change tip labels
tree_random$tip.label
tree_random$tip.label[tree_random$tip.label=="mrcaott343ott948"]<-"Arachnida"


tree_randomOrders$tip.label
tree_randomOrders$tip.label[tree_randomOrders$tip.label=="Diplura_(order_in_Mandibulata)"]<-"Diplura"


classes[!classes %in% as.character(tree_random$tip.label)] #listed in our database but not in the tree, 0 classes!
tree_random$tip.label[!as.character(tree_random$tip.label) %in% classes] # listed in the tree but not in our database. 0 classes!

orders[!orders %in% as.character(tree_randomOrders$tip.label)] #listed in our database but not in the tree, 0 orders! (Oribatida should be Sarcoptiformes)
tree_randomOrders$tip.label[!as.character(tree_randomOrders$tip.label) %in% orders] # listed in the tree but not in our database. 0 orders!

# Final tree
plot(tree_random, label.offset =.1, no.margin = TRUE)

tiff("./Results/Arthropods/Phylogeny/arth_phylogenetic_Classes_tree.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, label.offset =.1, no.margin = TRUE)

dev.off()


# Final tree
plot(tree_randomOrders, label.offset =.1, no.margin = TRUE)

tiff("./Results/Arthropods/Phylogeny/arth_phylogenetic_Orders_tree.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_randomOrders, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
save(tree_random, file = "./Results/Arthropods/Phylogeny/arth_tree_random.Rdata")
save(tree_randomOrders, file = "./Results/Arthropods/Phylogeny/arth_tree_Orders_random.Rdata")

##############################################################
# Computing branch lengths
##############################################################

# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature
classes[!classes %in% as.character(tree_random$tip.label)] #listed in our database but not in the tree, 0 species!
tree_random$tip.label[!as.character(tree_random$tip.label) %in% classes] # listed in the tree but not in our database. 0 species!
# all good!

# compute branch lengths of tree
phylo_branch <- compute.brlen(tree_random, method = "Grafen", power = 1)


phylo_branchOrders <- compute.brlen(tree_randomOrders, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE
is.ultrametric(phylo_branchOrders) # TRUE

save(phylo_branch, file = "./Results/Arthropods/Phylogeny/arth_tree_random.Rdata")
save(phylo_branchOrders, file = "./Results/Arthropods/Phylogeny/arth_tree_Orders_random_branchlength.Rdata")


# 
##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
arth_phylo_cor <- vcv(phylo_branch, cor = T)
arth_phyloOrders_cor <- vcv(phylo_branchOrders, cor = T)

# finally, save matrix for future analyses
save(arth_phylo_cor, file = "./Results/Arthropods/Phylogeny/arth_phylo_cor.Rdata")
save(arth_phyloOrders_cor, file = "./Results/Arthropods/Phylogeny/arth_phylo_Orders_cor.Rdata")

# saving session information with all packages versions for reproducibility purposes
sink("./Results/Arthropods/Phylogeny/arth_phylogeny_R_session.txt")
sessionInfo()
sink()
