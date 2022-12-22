# NImpact_Invertebrates_MetaAnalysis
This repository contains the data and scripts to run the analyses for the paper:
Gallego-Zamorano, J., de Jonge, M. M. J., Runge, K., Huls, S.H., Wang, J., Huijbregts, M. A. J., Schipper, A. M. Context-dependent responses of terrestrial invertebrates to anthropogenic nitrogen enrichment: a meta-analysis. (In review).

## Data overview
The "Databases" folder contains: 
- Decision list of all papers full-text screened by the different reviewers: "SecondScreeningPapers_DecisionList.xlsx"
- Reference list for all the papers included in both datasets:  "References_AllPapersAdded.xlsx"
- The raw dataset for abundance and richness: "AbundanceDataset.xlsx" and "RichnessDataset.xlsx". And their csv versions detoned by "_csv"
- The processed dataset with the imputed SDs and the extracted environmetal variables: "AbundanceDataset_ArthNema_imputed_Env.csv" and "RichnessDataset_ArthNema_imputed_Env.cvs"
- The processed dataset for the Arthropods filtered for rows with Order information: "AbundanceDataset_Arth_CompleteOrders_imputed.csv"

## Scripts
The "Scripts" folder contains all scripts in order to reproduce the analysis. The scripts with a 0 at the beggining of the name are necessary to prepare the data for the final analysis.

## Results
The "Results" folder contains the final fitted models for Arthropods and Nematodes for both abundance and richness.
