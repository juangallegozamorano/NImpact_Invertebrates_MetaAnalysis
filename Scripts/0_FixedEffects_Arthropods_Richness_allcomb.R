# Define combinations of variables
vars <- c("logNadd * logNdep", 
          "logNadd * Duration_Years", 
          "logNadd * Feeding_GuildUsed", 
          "logNadd * CRU_MAT_mean", 
          "logNadd * CRU_MAP_mean",
          "logNadd * CEC", 
          "logNadd * Habitat")

vars2 <- c("(logNadd + logNadd2) * logNdep", 
           "(logNadd + logNadd2) * Duration_Years", 
           "(logNadd + logNadd2) * Feeding_GuildUsed", 
           "(logNadd + logNadd2) * CRU_MAT_mean", 
           "(logNadd + logNadd2) * CRU_MAP_mean",
           "(logNadd + logNadd2) * CEC", 
           "(logNadd + logNadd2) * Habitat")


models <- list()
models2 <- list()

for (i in 1:(length(vars)-1)){
  varcomb <- combn(vars,i)
  for (j in 1:ncol(varcomb)){
    
    if("logNadd * CRU_MAT_mean" %in% varcomb[,j] & "logNadd * CRU_MAP_mean" %in% varcomb[,j]){
      next
    }
    
    
    if("logNadd * CRU_MAT_mean" %in% varcomb[,j] & "logNadd * logNdep" %in% varcomb[,j]){
      next
    }
    

    model <- as.formula(paste0("~", paste0(varcomb[,j], collapse = "+")))
    
    models <- c(models, model)
    
    
  }
}

for (i in 1:(length(vars2)-1)){
  varcomb2 <- combn(vars2,i)
  for (j in 1:ncol(varcomb2)){
    
    if("(logNadd + logNadd2) * CRU_MAT_mean" %in% varcomb2[,j] & "(logNadd + logNadd2) * CRU_MAP_mean" %in% varcomb2[,j]){
      next
    }
    
    if("(logNadd + logNadd2) * CRU_MAT_mean" %in% varcomb2[,j] & "(logNadd + logNadd2) * logNdep" %in% varcomb2[,j]){
      next
    }
    
    model2 <- as.formula(paste0("~", paste0(varcomb2[,j], collapse = "+")))
    
    models2 <- c(models2, model2)
    
    
  }
}

null <- as.formula("~ 1")
Nad <-  as.formula("~ logNadd")
Nad2 <-  as.formula("~ logNadd + logNadd2")

FixedEffects <- c(null, Nad, Nad2, models, models2)

saveRDS(FixedEffects,"./Scripts/FixedEffects_Arthropods_Richness.rds")
