####################################################################################################################
##
## Structural equation modelling (SEM) analysis with data resampling
##
####################################################################################################################


## This script contains the core SEM analysis to be applied to each dataset (Stressors-Abundance, Stressors-Biodiversity, Nutrients-Abundance and Nutrients-Biodiversity)

## This script has to be called from another script that defines what dataset 
# is used for the analysis and stores the results: 0203_SCRIPT_Run_SEM_resampling.R

## The code runs a loop that randomly resamples duplicated effect size on decomposition within case study at each out of 1000 iterations, 
## for each iteration, two SEM are fitted to the data, 
## the "null SEM" is the null model without path from biodiversity or from abundance responses to decomposition responses (biodiversity- or abundance-mediated path)
## the "full SEM" is the model with the biodiversity- or abundance-mediated path
## a bivariate meta-regression is also fitted to the data at each iteration (for plot purposes)
## The SEM combines two linear mixed effect models (= submodels = second level meta-analyses) one for decomposition, one for biodiversity or abundance
## model fits, statistics and main results across 1000 models are extracted and summarized






#==============================================================================================================================
set.seed(123)

# optimization method and lmeControl in nlme : to avoid convergence issues
ctrl <- lmeControl(sigma = 1, opt='optim', maxIter = 200, msMaxIter = 200)

## which columns contain the duplicates
## In this case, aFactor is a group of unique effect size on decomposition with different effect size on biodiversity (clusterID)
df$aFactor <- factor(df$clusterID_LD)
## and anInt is the repeatetd decomposition measures
df$anInt <- df$zcor.ECD.LD

## which columns are replicates
whichCols <- c("Case.study", "zcor.ECD.LD")


## Start a loop
n_times <- 1000 # How many times should we run the analysis with a random sample (1000)



#==============================================================================================================================
## Store model outputs in lists 

## Store meta-regression outputs in lists
metaregression <- list(NULL)

## Store residuals of each submodel of the SEM
resid.zLD <- list(NULL) ## submodel of decomposition responses
resid.zB <- list(NULL)  ## submodel of biodiversity or abundance responses

## Store fitted values of each submodel of the SEM
fitted.zLD <- list(NULL)
fitted.zB <- list(NULL)

# Fisher C statistics of SEM
fisherC0 <- list(NULL) ## null model
fisherC <- list(NULL)  ## full model (with b-ef path)

# AIC of SEM
AIC0 <- list(NULL)
AIC <- list(NULL)

# SEM path coefficients
coefs.sem <- list(NULL)

## R squared
rsquared <- list(NULL)

## ANOVA for individual models
Anova.zLD <- list(NULL)
Anova.zB <- list(NULL)



#==============================================================================================================================
## Start a loop

for(n in 1:n_times){ 
  
  tryCatch({# catch error messages (convergence problem) but don't stop the for loop
    
    df2 <- c() ## To put non-duplicated and sampled data into
    
    ## Go through each unique LD value individually 
    for (study in unique(df$aFactor)) {
      
      eachStudy <- df[df$aFactor == study, ]
      
      if (nrow(eachStudy) > 1) {
        # if there are duplicates of LD
        
        index <- sample(1:nrow(eachStudy), 1, replace = FALSE, prob = NULL)
        # randomly sample 1 of the duplicates
        
        df2 <- rbind(df2, eachStudy[index, ])
        
      } else{
        #if no duplicates, keep all the values
        
        df2 <- rbind(df2, eachStudy)
        
      }
    }
    
    rm(eachStudy, study, index)
      

    ### Perform the meta-regression analysis ===========================
    metaregression[[n]] <- rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                                  mods =~ zcor.ECD.B,
                                  random =~ 1|Case.study/ID,
                                  data = df2)
  
    ### Perform the SEM analysis============================================
  
  # for abundance data, 
    
  if(is.numeric(df2$zcor.ECD.D)){  # when data is for abundance,
    
    modelsem <- 
      psem(
        # sub model of LD response
        lme(
          zcor.ECD.LD ~  zcor.ECD.B + ECD_max + 
            study_type,
          random = ~ 1 | Case.study / ID,
          weights = varFixed( ~ var.zcor.ECD.LD),
          control = ctrl,
          data = df2),
        
        # sub model of B response (B stands for biodiversity or abundance to facilitate coding)
        lme(
          zcor.ECD.B ~ ECD_max +  
            study_type + taxonomic.group,
          random = ~ 1 | Case.study / ID,
          weights = varFixed( ~ var.zcor.ECD.B),
          control = ctrl,
          data = df2))
    
    
    model0 <- 
      psem(
        # sub model of LD response without b-ef path
        lme(
          zcor.ECD.LD ~  ECD_max + 
            study_type,
          random = ~ 1 | Case.study / ID,
          weights = varFixed( ~ var.zcor.ECD.LD),
          control = ctrl,
          data = df2),
        
        # sub model of B response
        lme(
          zcor.ECD.B ~ ECD_max +  
            study_type + taxonomic.group,
          random = ~ 1 | Case.study / ID,
          weights = varFixed( ~ var.zcor.ECD.B),
          control = ctrl,
          data = df2))
    
    
  } 
  else  { # when data is for decomposer diversity
    
    modelsem <- 
      psem(
        # sub model of LD response
        lme(
          zcor.ECD.LD ~  zcor.ECD.B + ECD_max + 
            study_type,
          random = ~ 1 | Case.study / ID,
          weights = varFixed( ~ var.zcor.ECD.LD),
          control = ctrl,
          data = df2),
        
        # sub model of B response
        lme(
          zcor.ECD.B ~ ECD_max +  
            study_type + taxonomic.group + B.metric,
          random = ~ 1 | Case.study / ID,
          weights = varFixed( ~ var.zcor.ECD.B),
          control = ctrl,
          data = df2))
  
  model0 <- 
    psem(
      # sub model of LD response without b-ef path
      lme(
        zcor.ECD.LD ~  ECD_max + 
          study_type,
        random = ~ 1 | Case.study / ID,
        weights = varFixed( ~ var.zcor.ECD.LD),
        control = ctrl,
        data = df2),
      
      # sub model of B response
      lme(
        zcor.ECD.B ~ ECD_max +  
          study_type + taxonomic.group + B.metric,
        random = ~ 1 | Case.study / ID,
        weights = varFixed( ~ var.zcor.ECD.B),
        control = ctrl,
        data = df2))
  }
  
  
  ## TEST OF MEDIATION : is B-EF path significant? 
  # run fisherC, extract C stat, pvalue
  fisherC0[[n]] <- fisherC(model0)
  fisherC[[n]] <- fisherC(modelsem)
  
  # compare AIC modelsem vs. model0
  AIC0[[n]] <- AIC(model0)
  AIC[[n]] <- AIC(modelsem)
  
  ## RESIDUALS and FITTED
  resid.zLD[[n]] <- residuals(modelsem, type = "normalized")[,1]
  fitted.zLD[[n]] <- fitted(modelsem[[1]])
  
  resid.zB[[n]] <- residuals(modelsem, type = "normalized")[,2]
  fitted.zB[[n]] <- fitted(modelsem[[2]])
  
  
  ## COEFFICIENTS
  coefs.sem[[n]] <- coefs(modelsem, intercepts = TRUE)

  
  ## ANOVA for individual models
  Anova.zLD[[n]] <- Anova(modelsem[[1]])
  Anova.zB[[n]] <- Anova(modelsem[[2]])
  
  
  ## R squared
  rsquared[[n]] <- rsquared(modelsem)

  
  }, error= function(e){cat("ERROR: ", conditionMessage(e), "\n")}) # catch error messages but don't stop for loop
}





#==============================================================================================================================
## OUTPUTS META-REGRESSION ===================================================================================================

# Save results from metaregressions
Res.metaregression <- metaregression[!sapply(metaregression,is.null)] # if model did not converge, don't include

# Extract intercept, slope, QM stat and sample sizes============================================================================
# calculate the mean estimates and stats across the 1000 models with data resampled
intercept <- mean(sapply(Res.metaregression, function(x){x$b[1]}))
slope <- mean(sapply(Res.metaregression, function(x){x$b[2]}))
QM_stat <- mean(sapply(Res.metaregression, function(x){x$QM}))
QM_Meanpvalue <- mean(sapply(Res.metaregression, function(x){x$QMp}))
QM_pvalue <- myfun_chisqpval(chi = QM_stat, df = 1) # recalculate the p-value associated with the mean QM
no.studies <- mean(sapply(Res.metaregression, function(x){x$s.nlevels[1]}))
no.obs <- mean(sapply(Res.metaregression, function(x){x$s.nlevels[2]}))

Res.Rmastat <- data.frame(intercept, slope, QM_stat, QM_Meanpvalue, QM_pvalue, no.studies, no.obs)



# Extract confidence intervals around the slope=================================================================================
# calculate confidence interval around the slope using function predict

res.Rmapredframe <- lapply(Res.metaregression, function(res){
  data.frame(
    zcor.ECD.LD = as.numeric(res$yi),
    var.zcor.ECD.LD = as.numeric(res$vi),
    preds = predict(res)$pred,
    lwr = predict(res)$ci.lb,
    upr = predict(res)$ci.ub)
}
)


# calculate the mean predictions and conf interval for each unique observation of z-LD

Predzcor.ECD.LD <- rowMeans(sapply(res.Rmapredframe, function(x) {x$zcor.ECD.LD} )) # mean is not useful here as all the zcor-LD are identical across models
Predvar.zcor.ECD.LD <- rowMeans(sapply(res.Rmapredframe, function(x) {x$var.zcor.ECD.LD} ))
Predpreds <-  rowMeans(sapply(res.Rmapredframe, function(x) {x$preds} ))
Predlwr <-  rowMeans(sapply(res.Rmapredframe, function(x) {x$lwr} ))
Predupr <-  rowMeans(sapply(res.Rmapredframe, function(x) {x$upr} ))

# gather mean predictions into a dataframe 
Res.Predframe <- data.frame(
  zcor.ECD.LD = Predzcor.ECD.LD,
  var.zcor.ECD.LD = Predvar.zcor.ECD.LD,
  preds = Predpreds,
  lwr = Predlwr,
  upr = Predupr)




## OUTPUTS SEM ===================================================================================================
# count no. models that did converge
Res.Noiter <- length(fisherC0[!sapply(fisherC0,is.null)])


## RESIDUALS==================================================================================================

# Remove models that did not converge
resid.zLD2 <- resid.zLD[!sapply(resid.zLD,is.null)]
resid.zB2<- resid.zB[!sapply(resid.zB,is.null)]
fitted.zLD2 <- fitted.zLD[!sapply(fitted.zLD,is.null)]
fitted.zB2 <- fitted.zB[!sapply(fitted.zB, is.null)]

 
# store residuals and fitted values for diagnostic plots
Res.Resid <- data.frame(
  
  # mean residuals across models
  residLD = rowMeans(sapply(resid.zLD2, c)),
  residB = rowMeans(sapply(resid.zB2, c)),
  
  # mean fitted values
  fittedLD = rowMeans(sapply(fitted.zLD2, c)),
  fittedB = rowMeans(sapply(fitted.zB2, c)))



## TEST OF MEDIATION===========================================================================================


# Remove models that did not converge
fisherC0 <- fisherC0[!sapply(fisherC0,is.null)]
fisherC <- fisherC[!sapply(fisherC,is.null)]


# Mean Fisher's C stat and AIC from null and full models

Cstat0 <- myfun_meanci(sapply(fisherC0, function(x){x$Fisher.C})) # means and ci
Cstat <- myfun_meanci(sapply(fisherC, function(x){x$Fisher.C}))

Cstat0$df <- sapply(fisherC0, function(x){x$df})[1] # degree of freedom
Cstat$df <- sapply(fisherC, function(x){x$df})[1]

## C stat follow a chi sq distrib with df = 2k if the causal model is true; k is the no. of independence claims in the sem
# Critical value of chi-square : how large would the statistic value need to be before we reject the causal model? 
# 1 - alpha = 0.95; df = degree of freedoms from sem 
Cstat0$chisq.theo <- qchisq(0.95, Cstat0$df)
Cstat$chisq.theo <- qchisq(0.95, Cstat$df)

# Reject the model if the C value is unlikely to have occurred by chance (i.e. is not in the 95% most likely values of Cstat with assoc. df)

# test the probability of the mean C-stat from resampled data vs. theo. distrib 
Cstat0$pval <- 1-pchisq(Cstat0$mean.est, Cstat0$df) # p < 0.05 = reject causal model
Cstat$pval <- 1-pchisq(Cstat$mean.est, Cstat$df)

# Create a summary table with Cstatistic of null and full models
Res.Cstat <- rbind(Cstat0, Cstat)


## AICs=======================================================================================================

# AIC from null and full SEMs
AIC.0 <- sapply(AIC0[!sapply(AIC0,is.null)], c)
AIC.full <- sapply(AIC[!sapply(AIC,is.null)], c)

# Table with AIC comparison
Res.AIC <- rbind(mod0 = myfun_meanci(AIC.0), modfull = myfun_meanci(AIC.full))# compare means and ci


## PATH COEFFICIENTS===========================================================================================

# Remove models that did not converge
coefs.sem <- coefs.sem[!sapply(coefs.sem,is.null)]

# Parameter estimates for each iteration
# unstandardized paths
Param_Estimates <- (sapply(coefs.sem, function(x) {as.numeric(x$Estimate)}))# Ignore NAs because they are for categorical variables
# standard error
Param_Std.Error <- (sapply(coefs.sem, function(x) {as.numeric(x$Std.Error)}))# Ignore NAs because they are for categorical variables
# standardized paths
Param_Std.Estimates <- (sapply(coefs.sem, function(x) {as.numeric(x$Std.Estimate)}))
# Critical value, df and p-value for each path
Param_Crit.Value <- (sapply(coefs.sem, function(x) {as.numeric(x$Crit.Value)}))
Param_DF <- (sapply(coefs.sem, function(x) {as.numeric(x$DF)}))
Param_P.Value <- (sapply(coefs.sem, function(x) {as.numeric(x$P.Value)}))


# Create a Summary Table of the results for the average model
Res.coefs <- data.frame(
  Response = coefs.sem[[1]]$Response, 
  Predictor = coefs.sem[[1]]$Predictor,
  
  Estimate = rowMeans(Param_Estimates),
  Std.Error = rowMeans(Param_Std.Error), # variance within model
  
  EstCi.lwr = apply(Param_Estimates, 1, function(x) quantile(x, 0.025, na.rm = TRUE)), # conf. int across models
  EstCi.upr = apply(Param_Estimates, 1, function(x) quantile(x, 0.975, na.rm = TRUE)),
  
  Crit.Value = rowMeans(Param_Crit.Value),
  DF = rowMeans(Param_DF),
  
  Std.Estimate = rowMeans(Param_Std.Estimates),
  Std.EstCi.lwr = apply(Param_Std.Estimates, 1, function(x) quantile(x, 0.025, na.rm = TRUE)), # conf. int across models
  Std.EstCi.upr = apply(Param_Std.Estimates, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
  )



#================================================================================
## OUTPUTS SUB-MODELS ================================================================================

# Remove models that did not converge
Anova.zLD <- Anova.zLD[!sapply(Anova.zLD,is.null)]
Anova.zB <- Anova.zB[!sapply(Anova.zB,is.null)]

# Recreate anova table for the based on means from data resampling: submodel decomposition
Res.AnovazLD <- data.frame(
  Predictor = rownames(Anova.zLD[[1]]),
  Chisq = rowMeans(sapply(Anova.zLD, function(x) {x$Chisq})),
  Df = rowMeans(sapply(Anova.zLD, function(x) {x$Df}))
  )

# calculate p-value associated with the mean Chi square
Res.AnovazLD$P.Value <- myfun_chisqpval(chi = Res.AnovazLD$Chisq, df = Res.AnovazLD$Df)

# Recreate anova table for the based on means from data resampling: submodel biodiversity or abundance
Res.AnovazB <- data.frame(
  Predictor = rownames(Anova.zB[[1]]),
  Chisq = rowMeans(sapply(Anova.zB, function(x) {x$Chisq})),
  Df = rowMeans(sapply(Anova.zB, function(x) {x$Df}))
)

# calculate p-value associated with the mean Chi square
Res.AnovazB$P.Value <- myfun_chisqpval(chi = Res.AnovazB$Chisq, df = Res.AnovazB$Df)


# R squares (although not sure about their relevance in a meta-analysis context)
Res.Rsq <- data.frame(
  Response = rsquared[[1]]$Response,
  MarginalR2 = rowMeans(sapply(rsquared[!sapply(rsquared,is.null)], function(x) x$Marginal)),
  ConditionalR2 = rowMeans(sapply(rsquared[!sapply(rsquared,is.null)], function(x) x$Conditional)))


