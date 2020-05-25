####################################################################################################################
##
## RESULTS from SEM analysis of stressors and nutrients effects on decomposers and decomposition
##
###################################################################################################################


## This scripts summarises the results of structural equation models based on resampled data
## Results are derived from 0203_SCRIPT_Run_SEM_resampling.R where models are iteratively fitted to
## randomly resampled data and results are stored in the data folder

source("0201_LOAD_Data.R")

# 1. Relationship between decomposition and decomposers responses-------------------------

# Results of meta-regressions between effect sizes on decomposition and effect sizes on 
# diversity or abundance of decomposers. Results are based on a data resampling approach. 

## Meta-regressions-------------------------------------------------------------

# bivariate metaregressions based on resampled data : slopes and stats
Res.rmastat_pol_bef <- read.csv("data/03_Model_results/ResRmastat_pol_bef.csv")
Res.rmastat_pol_def <- read.csv("data/03_Model_results/ResRmastat_pol_def.csv")
Res.rmastat_nut_bef <- read.csv("data/03_Model_results/ResRmastat_nut_bef.csv")
Res.rmastat_nut_def <- read.csv("data/03_Model_results/ResRmastat_nut_def.csv")


Res.rmastat_pol_bef
Res.rmastat_pol_def 
Res.rmastat_nut_bef 
Res.rmastat_nut_def 

## FIGURE BEF relationship---------------------------------------------------------------
## SCript to generate the figure showing relation between effect sizes of chemical stressors and 
## nutrient enrichment on decomposition and on decompsers Diversity


# confidence intervals around slope
Res.predframe_pol_bef <- read.csv("data/03_Model_results/ResPredframe_pol_bef.csv")
Res.predframe_pol_def <- read.csv("data/03_Model_results/ResPredframe_pol_def.csv")
Res.predframe_nut_bef <- read.csv("data/03_Model_results/ResPredframe_nut_bef.csv")
Res.predframe_nut_def <- read.csv("data/03_Model_results/ResPredframe_nut_def.csv")

## Load datasets ----------------------------
data.list <- mget(grep("EF_es", ls(), value = TRUE))



# set sizes of plotted elements (from Fig1)
sizetext <- 13
sizelegend <- 11
sizepoint <- 2.5
sizeline <- 0.8
sizeannotate <- 3.7

colo_lea <- c("#0072B2", "#D55E00")



## Function to create panel for each driver type 

myfun_bivariateBEFchem <- function(dat, slope, predframe, plottitle, xlimBorA, xlabBorA){
  
  ## create the group-mean data
  gd <- dat %>%
    group_by(clusterID_LD) %>%
    summarise(zcor.ECD.LD = mean(zcor.ECD.LD),
              var.zcor.ECD.LD  = mean(var.zcor.ECD.LD),
              zcor.ECD.B  = mean(zcor.ECD.B),
              ECD_max = mean(ECD_max),
              weight = 1/mean(var.zcor.ECD.LD),
              Case.study = mean(Case.study),
              ID = ID[1])
  
  
  # add zcorECDB to predframe
  predframe2 <- dplyr::left_join(gd, predframe, by = c("zcor.ECD.LD", "var.zcor.ECD.LD"))
  
  
  ## plot 
  ggplot(dat, aes(x=zcor.ECD.B, y=zcor.ECD.LD, colour = ECD_max, 
                  # size = weights(res.bivariate.pol))) +
                  size = weight)) +
    
    # alpha to render the details trasnparent
    geom_point(shape=18,  col = "gray", size = sizepoint/1.5)+
    
    # add the means in a bigger size
    geom_point(shape = 18, data = gd, size = sizepoint)+
    
    # change colour of points continuous scale
    scale_colour_gradient(high = "#132B43", low = "#56B1F7",
                          space = "Lab", na.value = "grey50", guide = "colourbar",
                          aesthetics = "colour", name = 'Stressor intensity') +
    
    # titles and axis labels
    ylab("Effect on decomposition")+
    # xlab("Effect on biodiversity")+
    xlab(xlabBorA)+
    ggtitle(plottitle) +
    
    # axis lenght
    ylim(c(-3.1,2.6))+
    # xlim(c(-3.2,4.1))+
    xlim(xlimBorA)+
    
    # mark the zeros lines
    geom_hline(yintercept = 0, color = "black")+
    geom_vline(xintercept = 0, color = "black")+
    
    ## add slopes from resampled data  
    geom_abline(slope = slope$slope, intercept = slope$intercept, col = colo_lea[1], size = sizeline,
                linetype = 1 + if(slope$QM_pvalue>0.05){1} else{0})+ # change lty according to p-value of the meta-regression
    geom_line(data = predframe2, aes(x = zcor.ECD.B, y = lwr), inherit.aes = FALSE,
              linetype = 1 + if(slope$QM_pvalue>0.05){1} else{0}, color = colo_lea[1], size = sizeline/1.5)+
    geom_line(data = predframe2, aes(x = zcor.ECD.B, y = upr), inherit.aes = FALSE,
              linetype = 1 + if(slope$QM_pvalue>0.05){1} else{0}, color = colo_lea[1], size = sizeline/1.5)+
    
    
    # theme stff
    theme_bw() +
    theme(
      axis.text.y = element_text(face = "bold", size = sizetext),
      axis.text.x = element_text(face = "bold", size = sizetext),
      axis.title.x = element_text(size = sizetext, face = "bold"),
      axis.title.y = element_text(size = sizetext, face = "bold"),
      
      #legend
      plot.title = element_text(size = sizetext),
      legend.title = element_text(size = sizelegend),
      legend.text = element_text(size = sizelegend))
}



myfun_bivariateBEFnut <- function(dat, slope, predframe, plottitle, xlimBorA, xlabBorA){
  
  ## create the group-mean data
  gd <- dat %>% 
    group_by(clusterID_LD) %>% 
    summarise(zcor.ECD.LD = mean(zcor.ECD.LD),
              var.zcor.ECD.LD  = mean(var.zcor.ECD.LD), 
              zcor.ECD.B  = mean(zcor.ECD.B), 
              ECD_max = mean(ECD_max), 
              weight = 1/mean(var.zcor.ECD.LD), 
              Case.study = mean(Case.study), 
              ID = ID[1])
  
  
  # add zcorECDB to predframe
  predframe2 <- dplyr::left_join(gd, predframe, by = c("zcor.ECD.LD", "var.zcor.ECD.LD"))
  
  # start the plot
  ggplot(dat, aes(x=zcor.ECD.B, y=zcor.ECD.LD, colour = ECD_max, 
                  # size = weights(res.bivariate.pol))) +
                  size = weight)) +
    
    # alpha to render the details trasnparent
    geom_point(shape=18,  col = "gray", size = sizepoint/1.5)+
    
    # add the means in a bigger size
    geom_point(shape = 18, data = gd, size = sizepoint)+
    
    scale_colour_gradient(high = "#990000", low = "#FFF666",
                          space = "Lab", na.value = "grey50", guide = "colourbar",
                          aesthetics = "colour", name = 'Nutrient intensity') +
    
    ylab("Effect on decomposition")+
    xlab(xlabBorA)+
    ggtitle(plottitle)+
    
    # axis lenght
    ylim(c(-3.1,2.6))+
    xlim(xlimBorA)+
    
    
    ## add slopes from resampled data  
    geom_abline(slope = slope$slope, intercept = slope$intercept, col = colo_lea[2], size = sizeline,
                linetype = 1 + if(slope$QM_pvalue>0.05){1} else{0})+ # change lty according to p-value of the meta-regression
    geom_line(data = predframe2, aes(x = zcor.ECD.B, y = lwr), inherit.aes = FALSE,
              linetype = 1 + if(slope$QM_pvalue>0.05){1} else{0}, color = colo_lea[2], size = sizeline/1.5)+
    geom_line(data = predframe2, aes(x = zcor.ECD.B, y = upr), inherit.aes = FALSE,
              linetype = 1 + if(slope$QM_pvalue>0.05){1} else{0}, color = colo_lea[2], size = sizeline/1.5)+
    
    
    # mark the zeros lines
    geom_hline(yintercept = 0, color = "black")+#"darkgray")+
    geom_vline(xintercept = 0, color = "black")+#"darkgray")+
    
    
    # theme stff
    theme_bw() +
    theme(
      axis.text.y = element_text(face = "bold", size = sizetext),
      axis.text.x = element_text(face = "bold", size = sizetext),
      axis.title.x = element_text(size = sizetext, face = "bold"),
      axis.title.y = element_text(size = sizetext, face = "bold"),
      
      #legend
      plot.title = element_text(size = sizetext),
      legend.title = element_text(size = sizelegend),
      legend.text = element_text(size = sizelegend))
  
}



## arrange 4 panels into a single figure
Fig2A <-  myfun_bivariateBEFchem(pol_BEF_es, 
                                 Res.rmastat_pol_bef, 
                                 Res.predframe_pol_bef,
                                 plottitle = 'Stressors - Diversity',
                                 xlimBorA = c(-3.2,3.2),
                                 xlabBorA = "Effect on biodiversity") +  theme(legend.position='none')+
  annotate('label', x= -1.6, y= - 3,  hjust = 0, label=paste("QM = 8.6; p = 0.003 (23; 100)"), 
           fill = "white", size = sizeannotate, label.size = NA)

Fig2B <- myfun_bivariateBEFchem(pol_DEF_es, 
                                Res.rmastat_pol_def, 
                                Res.predframe_pol_def,
                                plottitle = 'Stressors - Abundance',
                                xlimBorA = c(-3.2,4.1),
                                xlabBorA = "Effect on abundance") +
  annotate('label', x= -1.4, y= - 3,  hjust = 0, label=paste("QM = 9.2; p = 0.002 (28; 282)"), 
           fill = "white", size = sizeannotate, label.size = NA)

Fig2C <- myfun_bivariateBEFnut(nut_BEF_es, 
                               Res.rmastat_nut_bef, 
                               Res.predframe_nut_bef,
                               plottitle = 'Nutrients - Diversity',
                               xlimBorA = c(-3.2,3.2),
                               xlabBorA = "Effect on biodiversity")+  theme(legend.position='none')+
  annotate('label', x= -1.2, y= - 3,  hjust = 0, label=paste("QM = 0.4; p = 0.53 (26; 97)"), 
           fill = "white", size = sizeannotate, label.size = NA)
Fig2D <- myfun_bivariateBEFnut(nut_DEF_es, 
                               Res.rmastat_nut_def, 
                               Res.predframe_nut_def,
                               plottitle = 'Nutrients - Abundance',
                               xlimBorA = c(-3.2,4.1),
                               xlabBorA = "Effect on abundance") +
  annotate('label', x= -1.2, y= - 3,  hjust = 0, label=paste("QM = 1.3; p = 0.25 (39; 181)"), 
           fill = "white", size = sizeannotate, label.size = NA)


Figure2 <- Fig2A + Fig2B + Fig2C + Fig2D

## Save final Figure 4 of the Manuscript-------------------------------------------------------------
# save a png with high res
ppi <- 300 #final resolution 600 # resolution
w <- 21 # width in cm

png("figs/Fig4_bivariateBEF.png",
    width = w,
    height= w/1.3 ,
    units = "cm",
    res=ppi)

Figure2

dev.off()


## 2. SEM results-------------------------------------------------------------

## 2.1. Residual plots-------------------------------------------------------------

# Load datasets of residuals and fitted values of the submodels of piecewiseSEM
# Values are average across 1000 model iterations
Res.Resid_pol_bef <- read.csv("data/03_Model_results/ResResid_pol_bef.csv")
Res.Resid_pol_def <- read.csv("data/03_Model_results/ResResid_pol_def.csv")
Res.Resid_nut_bef <- read.csv("data/03_Model_results/ResResid_nut_bef.csv")
Res.Resid_nut_def <- read.csv("data/03_Model_results/ResResid_nut_def.csv")

# function for plotting residual diagnostic plots
myfun_residresampl <- function(resid, fitted, plottitle){
  par(mfrow = c(1,3))
  qqnorm(resid, main = plottitle)
  qqline(resid, col = "red")
  hist(resid, main = plottitle)
  plot(resid ~ fitted, main = paste("Residuals vs Fitted", plottitle)); abline(0,0)
}

# Run for each dataset
myfun_residresampl(Res.Resid_pol_bef$residLD, Res.Resid_pol_bef$fittedLD, "Stressors, Decomposition (div data)")
myfun_residresampl(Res.Resid_pol_bef$residB, Res.Resid_pol_bef$fittedB, "Stressors, Diversity")

myfun_residresampl(Res.Resid_pol_def$residLD, Res.Resid_pol_def$fittedLD, "Stressors, Decomposition (abdc data)")
myfun_residresampl(Res.Resid_pol_def$residB, Res.Resid_pol_def$fittedB, "Stressors, Abundance")

myfun_residresampl(Res.Resid_nut_bef$residLD, Res.Resid_nut_bef$fittedLD, "Nutrients, Decomposition (div data)")
myfun_residresampl(Res.Resid_nut_bef$residB, Res.Resid_nut_bef$fittedB, "Nutrients, Diversity")

myfun_residresampl(Res.Resid_nut_def$residLD, Res.Resid_nut_def$fittedLD, "Nutrients, Decomposition (abdc data)")
myfun_residresampl(Res.Resid_nut_def$residB, Res.Resid_nut_def$fittedB, "Nutrients, Abundance")

# residuals indicate poor fit especially for nutrient abundance data


## 2.2. MEDIATION TESTS -----------------------------------------------------------------

# Tests of mediation (C stat and AIC): comparing models with and without the biodiversity- or abundance-
# mediated paths


## C statistic of null model (no b-ef path)-----------------------------------
ResCstat_pol_bef <- read.csv("data/03_Model_results/ResCstat_pol_bef.csv")
ResCstat_pol_def <- read.csv("data/03_Model_results/ResCstat_pol_def.csv")
ResCstat_nut_bef <- read.csv("data/03_Model_results/ResCstat_nut_bef.csv")
ResCstat_nut_def <- read.csv("data/03_Model_results/ResCstat_nut_def.csv")

# C stat of models omitting the B- or A-path
ResCstatTable <- rbind(ResCstat_pol_bef[1,], 
                       ResCstat_pol_def[1,], 
                       ResCstat_nut_bef[1,], 
                       ResCstat_nut_def[1,])

# C stat of models with the B- or A-path (but there are no clear d-sep test in the full models)
ResCstatTablefull <- rbind(ResCstat_pol_bef[2,], 
                           ResCstat_pol_def[2,], 
                           ResCstat_nut_bef[2,], 
                           ResCstat_nut_def[2,])

ResCstatTable$X <- c("Stressors, Biodiv", "Stressors, Abdc", "Nutrient, Biodiv", "Nutrient, Abdc")

## AIC and delta AIC------------------------------------------------

ResAIC_pol_bef <- read.csv("data/03_Model_results/ResAIC_pol_bef.csv")
ResAIC_pol_def <- read.csv("data/03_Model_results/ResAIC_pol_def.csv")
ResAIC_nut_bef <- read.csv("data/03_Model_results/ResAIC_nut_bef.csv")
ResAIC_nut_def <- read.csv("data/03_Model_results/ResAIC_nut_def.csv")


# calculate delta AIC (AIC full model - AIC null model): negative values means the full model improved // null model
deltaAIC <- c(ResAIC_pol_bef$mean.est[2] - ResAIC_pol_bef$mean.est[1], 
              ResAIC_pol_def$mean.est[2] - ResAIC_pol_def$mean.est[1],
              ResAIC_nut_bef$mean.est[2] - ResAIC_nut_bef$mean.est[1],
              ResAIC_nut_def$mean.est[2] - ResAIC_nut_def$mean.est[1])


## add total n and no. studies------------------------------------
counts <- rbind(
  as.numeric(countstudies(pol_BEF_es)),
  as.numeric(countstudies(pol_DEF_es)),
  as.numeric(countstudies(nut_BEF_es)),
  as.numeric(countstudies(nut_DEF_es)))

colnames(counts) <- c("no. studies", "n")

## Final output table summarizing the results
ResCstatTable <- cbind(subset(ResCstatTable, select = c(X, mean.est, ci.lwr, ci.upr, df, pval)), deltaAIC, counts)
print(ResCstatTable)

# Save results in a table: values are reported on Figure 5
write.csv(ResCstatTable, "tables/ResCstatTable.csv")


## 3.3. Conditional and marginal R squares-----------------------------------
# those are derived from the piecewiseSEM package
ResRsq_pol_bef <- read.csv("data/03_Model_results/ResRsq_pol_bef.csv")
ResRsq_pol_def <- read.csv("data/03_Model_results/ResRsq_pol_def.csv")
ResRsq_nut_bef <- read.csv("data/03_Model_results/ResRsq_nut_bef.csv")
ResRsq_nut_def <- read.csv("data/03_Model_results/ResRsq_nut_def.csv")


ResRsqTable <- rbind(ResRsq_pol_bef, 
                     ResRsq_pol_def, 
                     ResRsq_nut_bef, 
                     ResRsq_nut_def)

ResRsqTable$X <- c(rep("Stressors, Diversity",2), rep("Stressors, Abundance",2), rep("Nutrients, Diversity",2), rep("Nutrients, Abundance", 2))

print(ResRsqTable)

# save table
write.csv(ResRsqTable, "tables/ResRsquaredSEMTable.csv")


## 3.4. Summary tables-------------------------------------------------------------------

# Summary tables showing the results from SEM with path coefficients and 
# mean per level of categorical moderators in the 4 models.

ResCoefs_pol_bef <- read.csv("data/03_Model_results/ResCoefs_pol_bef.csv")
ResCoefs_pol_def <- read.csv("data/03_Model_results/ResCoefs_pol_def.csv")
ResCoefs_nut_bef <- read.csv("data/03_Model_results/ResCoefs_nut_bef.csv")
ResCoefs_nut_def <- read.csv("data/03_Model_results/ResCoefs_nut_def.csv")

# Stressors - Biodiversity
ResCoefs_pol_bef

# Stressors - Abundance
ResCoefs_pol_def

# Nutrients - Biodiversity
ResCoefs_nut_bef

# Nutrients - Biodiversity
ResCoefs_nut_bef



## 3.5. TABLES OF UNSTANDARDIZED COEFFICIENTS----------------------------------------------

## TABLE with the standardized path coefs from SEMs

ResCoefs_pol_bef <- read.csv("data/03_Model_results/ResCoefs_pol_bef.csv")
ResCoefs_pol_def <- read.csv("data/03_Model_results/ResCoefs_pol_def.csv")
ResCoefs_nut_bef <- read.csv("data/03_Model_results/ResCoefs_nut_bef.csv")
ResCoefs_nut_def <- read.csv("data/03_Model_results/ResCoefs_nut_def.csv")

Reslist <-  mget(grep("ResCoefs_", ls(), value = TRUE))

# take only path coefficients (slopes)
PathCoefs <- 
  lapply(Reslist, function(x) {x[x$Predictor %in% c("zcor.ECD.B", "ECD_max"), ]})

# which column to keep in final table
WhichColumns <-  c("Estimate", "Std.Error", "Crit.Value", "DF")

# make table simplified 
ttable <-   do.call(rbind, lapply(PathCoefs, function(x) x[,WhichColumns]))

# calculate p-value associated with the critical value
myfun_tstat <- function(t, df) {2*pt(abs(t), df, lower.tail = FALSE)}# see e.g. https://stats.stackexchange.com/questions/5135/interpretation-of-rs-lm-output

ttable$P.Value <- myfun_tstat(t = ttable$Crit.Value, df = ttable$DF)

FinaltTable <- cbind(
  # Arrange final table with renamed columns
  data.frame(
    ECD.type=factor(ifelse(grepl("nut", rownames(ttable)), "Nutrients", "Stressors")),
    B.metric=factor(ifelse(grepl("bef",rownames(ttable)), "Diversity", "Abundance")),
    path=factor(rep(c("B -> LD", "ECD -> LD", "ECD -> B"), 4))),
  
  # Add coefficients, t value and P-value
  ttable
)

# reorder table with stressors first
FinaltTable$ECD.type <- factor(FinaltTable$ECD.type, c("Stressors", "Nutrients"))
FinaltTable <- with(FinaltTable, FinaltTable[order(ECD.type),])

print(FinaltTable)

# Save supplementary table: values are reported on Figure 5
write.csv(FinaltTable, "tables/ResPathCoefsResampledat.csv")


## 3.6. Direct and indirect effects------------------------------------------------------
# compare direct ECD -> LD path vs. indirect ECD->B->LD path

# standardized coefs
WhichColumns <-  c("Response","Predictor", "Std.Estimate")

# simplify the table
pathtable <-   lapply(PathCoefs, function(x) x[,WhichColumns])


DirectIndirect <- do.call(rbind, lapply(pathtable, function(x) c(x[1,3] * x[3,3],  #indirect effect = std coefs B-LD * ECD-B
                                                                 x[2,3])))         #direct effect = std coef ECD-LD
colnames(DirectIndirect) <-  c("Indirect_Bmediated", "Direct")

# table with ECD type and B metric
DirectIndirect <- cbind(
  # Arrange final table with renamed columns
  data.frame(
    ECD.type=factor(ifelse(grepl("nut", rownames(DirectIndirect)), "Nutrients", "Stressors")),
    B.metric=factor(ifelse(grepl("bef",rownames(DirectIndirect)), "Diversity", "Abundance")),
    
    # Add the coefficients
    DirectIndirect
  ))

# reorder table with stressors first
DirectIndirect$ECD.type <- factor(DirectIndirect$ECD.type, c("Stressors", "Nutrients"))
DirectIndirect <- with(DirectIndirect, DirectIndirect[order(ECD.type),])

# Store table, values are reported in Figure 5
print(DirectIndirect)
write.csv(DirectIndirect, "tables/DirectIndirectEffects.csv")




