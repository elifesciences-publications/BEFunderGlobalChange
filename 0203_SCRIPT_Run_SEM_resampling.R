####################################################################################################################
##
## Run SEM with data resampling for the four datasets (Stressors - Biodiversity, Stressors Abundance, Nutrients- Biodiversity and Nutrients - Abundance datasets)
##
###################################################################################################################


## Script to run SEMs and store the results 
# Iteratively for each of four datasets (takes a long time to run)

# calls the function "SEM_resampling.R" where models are iteratively fitted 1000 times
# to randomly resampled data accounting for multiple diversity or abundance data maching
# a single value of decomposition response in studies reporting e.g. different metrics of
# diversity or different taxonomic groups in the same litterbag.


# load data and functions
source("0201_LOAD_Data.R")

## CHEMICAL STRESSORS RESPONSES=========================================


## dataset biodiversity of decomposers & decomposition responses
df <- pol_BEF_es

## run the models
source("functions/SEM_resampling.R")

# Store Results
saveRDS(Res.metaregression, "data/03_Model_results/Resmetaregression_pol_bef.R")
write.csv(Res.Rmastat, "data/03_Model_results/ResRmastat_pol_bef.csv")
write.csv(Res.Predframe, "data/03_Model_results/ResPredframe_pol_bef.csv")
write.csv(Res.Noiter, "data/03_Model_results/ResNoiter_pol_bef.csv")
write.csv(Res.Resid, "data/03_Model_results/ResResid_pol_bef.csv")
write.csv(Res.Cstat, "data/03_Model_results/ResCstat_pol_bef.csv")
write.csv(Res.AIC, "data/03_Model_results/ResAIC_pol_bef.csv")
write.csv(Res.coefs, "data/03_Model_results/ResCoefs_pol_bef.csv")
write.csv(Res.AnovazLD, "data/03_Model_results/ResAnovazLD_pol_bef.csv")
write.csv(Res.AnovazB,  "data/03_Model_results/ResAnovazB_pol_bef.csv")
write.csv(Res.Rsq,  "data/03_Model_results/ResRsq_pol_bef.csv")


## Dataset of decomposer abundance and decomposition responses
df <- pol_DEF_es

source("functions/SEM_resampling.R")

# Store results
saveRDS(Res.metaregression, "data/03_Model_results/Resmetaregression_pol_def.R")
write.csv(Res.Rmastat, "data/03_Model_results/ResRmastat_pol_def.csv")
write.csv(Res.Predframe, "data/03_Model_results/ResPredframe_pol_def.csv")
write.csv(Res.Noiter, "data/03_Model_results/ResNoiter_pol_def.csv")
write.csv(Res.Resid, "data/03_Model_results/ResResid_pol_def.csv")
write.csv(Res.Cstat, "data/03_Model_results/ResCstat_pol_def.csv")
write.csv(Res.AIC, "data/03_Model_results/ResAIC_pol_def.csv")
write.csv(Res.coefs, "data/03_Model_results/ResCoefs_pol_def.csv")
write.csv(Res.AnovazLD, "data/03_Model_results/ResAnovazLD_pol_def.csv")
write.csv(Res.AnovazB,  "data/03_Model_results/ResAnovazB_pol_def.csv")
write.csv(Res.Rsq,  "data/03_Model_results/ResRsq_pol_def.csv")



## NUTRIENT ENRICHMENT RESPONSES=========================================


## dataset biodiversity of decomposers & decomposition responses
df <- nut_BEF_es

## run the models
source("functions/SEM_resampling.R")

# Store Results
saveRDS(Res.metaregression, "data/03_Model_results/Resmetaregression_nut_bef.R")
write.csv(Res.Rmastat, "data/03_Model_results/ResRmastat_nut_bef.csv")
write.csv(Res.Predframe, "data/03_Model_results/ResPredframe_nut_bef.csv")
write.csv(Res.Noiter, "data/03_Model_results/ResNoiter_nut_bef.csv")
write.csv(Res.Resid, "data/03_Model_results/ResResid_nut_bef.csv")
write.csv(Res.Cstat, "data/03_Model_results/ResCstat_nut_bef.csv")
write.csv(Res.AIC, "data/03_Model_results/ResAIC_nut_bef.csv")
write.csv(Res.coefs, "data/03_Model_results/ResCoefs_nut_bef.csv")
write.csv(Res.AnovazLD, "data/03_Model_results/ResAnovazLD_nut_bef.csv")
write.csv(Res.AnovazB,  "data/03_Model_results/ResAnovazB_nut_bef.csv")
write.csv(Res.Rsq,  "data/03_Model_results/ResRsq_nut_bef.csv")


## Dataset of decomposer abundance and decomposition responses
df <- nut_DEF_es

source("functions/SEM_resampling.R")

# Store results
saveRDS(Res.metaregression, "data/03_Model_results/Resmetaregression_nut_def.R")
write.csv(Res.Rmastat, "data/03_Model_results/ResRmastat_nut_def.csv")
write.csv(Res.Predframe, "data/03_Model_results/ResPredframe_nut_def.csv")
write.csv(Res.Noiter, "data/03_Model_results/ResNoiter_nut_def.csv")
write.csv(Res.Resid, "data/03_Model_results/ResResid_nut_def.csv")
write.csv(Res.Cstat, "data/03_Model_results/ResCstat_nut_def.csv")
write.csv(Res.AIC, "data/03_Model_results/ResAIC_nut_def.csv")
write.csv(Res.coefs, "data/03_Model_results/ResCoefs_nut_def.csv")
write.csv(Res.AnovazLD, "data/03_Model_results/ResAnovazLD_nut_def.csv")
write.csv(Res.AnovazB,  "data/03_Model_results/ResAnovazB_nut_def.csv")
write.csv(Res.Rsq,  "data/03_Model_results/ResRsq_nut_def.csv")