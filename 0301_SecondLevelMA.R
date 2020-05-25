####################################################################################################################
##
## SECOND-LEVEL META-ANALYSIS
##
###################################################################################################################


# This script reproduces the results of the meta-analysis testing the effect of covariates on
# decomposers and decompositionr esponses to stressors and nutrients

## 1. LOAD Data---------------------------------------------------------------------------
source("0201_LOAD_Data.R")

## 2. MODELS------------------------------------------------------------------------------
## Run second-level meta-analyses of decomposers responses (diversity and abundance) 
# with metafor to derive confidence intervals, QM stats and p-value
# and using all the data (no need for data resampling in this univariate approach)
# these models correspond to the submodels for biodiversity and abundance responses in the piecewise SEMs

modelrma_polbef <-  rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                           mods = ~ ECD_max +  
                             study_type + taxonomic.group + B.metric,
                           random = ~ 1 | Case.study / ID, 
                           data = pol_BEF_es)
modelrma_poldef <-  rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                           mods = ~ ECD_max +  
                             study_type + taxonomic.group,
                           random = ~ 1 | Case.study / ID, 
                           data = pol_DEF_es)
modelrma_nutbef <-  rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                           mods = ~ ECD_max +  
                             study_type + taxonomic.group + B.metric,
                           random = ~ 1 | Case.study / ID, 
                           data = nut_BEF_es)
modelrma_nutdef <-  rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                           mods = ~ ECD_max +  
                             study_type + taxonomic.group,
                           random = ~ 1 | Case.study / ID, 
                           data = nut_DEF_es)


# as all factors have only two levels, the pvalues in the summary test the significance of the factor overall

# taxonomic.group effect
resanova_taxo <- data.frame(rbind(
  pol_bef = as.numeric(c(anova(modelrma_polbef, btt = 4)[1:2])),
  pol_def = as.numeric(c(anova(modelrma_poldef, btt = 4)[1:2])),
  nut_bef = as.numeric(c(anova(modelrma_nutbef, btt = 4)[1:2])),
  nut_def = as.numeric(c(anova(modelrma_nutdef, btt = 4)[1:2]))))
resanova_taxo$Predictor <- rep("Taxonomic group", 4)
resanova_taxo$Response <- rep(c("Diversity", "Abundance"), 2)
resanova_taxo$ECD.type <- c(rep("Stressors", 2), rep("Nutrients", 2))


# study type effect
resanova_study <- data.frame(rbind(
  pol_bef = as.numeric(c(anova(modelrma_polbef, btt = 3)[1:2])),
  pol_def = as.numeric(c(anova(modelrma_poldef, btt = 3)[1:2])),
  nut_bef = as.numeric(c(anova(modelrma_nutbef, btt = 3)[1:2])),
  nut_def = as.numeric(c(anova(modelrma_nutdef, btt = 3)[1:2]))))
resanova_study$Predictor <- rep("Study type", 4)
resanova_study$Response <- rep(c("Diversity", "Abundance"), 2)
resanova_study$ECD.type <- c(rep("Stressors", 2), rep("Nutrients", 2))



# B.metric effect
resanova_Bmetric <- data.frame(rbind(
  pol_bef = as.numeric(c(anova(modelrma_polbef, btt = 5)[1:2])),
  pol_def = rep(NA, 2),
  nut_bef = as.numeric(c(anova(modelrma_nutbef, btt = 5)[1:2])),
  nut_def = rep(NA, 2)))
resanova_Bmetric$Predictor <- rep("Diversity metric", 4)
resanova_Bmetric$Response <- rep(c("Diversity", "Abundance"), 2)
resanova_Bmetric$ECD.type <- c(rep("Stressors", 2), rep("Nutrients", 2))

names(resanova_taxo)[c(1,2)] <- names(resanova_study)[c(1,2)] <- names(resanova_Bmetric)[c(1,2)] <- c("QM", "P")

print(resanova_taxo)
print(resanova_study)
print(resanova_Bmetric)



# save results
write.csv(resanova_taxo, "tables/ResMATaxo.csv")
write.csv(resanova_study, "tables/ResMAStudy.csv")
write.csv(resanova_Bmetric, "tables/ResMABmetric.csv")


# second-level meta-analyses for decomposition: remove duplicates ES on decompo
modelrmald_polbef <-  rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                             mods = ~ zcor.ECD.B + ECD_max +  
                               study_type,
                             random = ~ 1 | Case.study / ID, 
                             data = pol_BEF_es[!duplicated(pol_BEF_es$clusterID_LD),])
modelrmald_poldef <-  rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                             mods = ~ zcor.ECD.B+ ECD_max +  
                               study_type,
                             random = ~ 1 | Case.study / ID, 
                             data = pol_DEF_es[!duplicated(pol_DEF_es$clusterID_LD),])
modelrmald_nutbef <-  rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                             mods = ~ zcor.ECD.B+ ECD_max +  
                               study_type,
                             random = ~ 1 | Case.study / ID, 
                             data = nut_BEF_es[!duplicated(nut_BEF_es$clusterID_LD),])
modelrmald_nutdef <-  rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                             mods = ~ zcor.ECD.B+ ECD_max +  
                               study_type,
                             random = ~ 1 | Case.study / ID, 
                             data = nut_DEF_es[!duplicated(nut_DEF_es$clusterID_LD),])


# as all factors have only two levels, the pvalues in the summary test the significance of the factor

# study type effect
resanovaLD_study <- data.frame(rbind(
  pol_bef = as.numeric(c(anova(modelrmald_polbef, btt = 4)[1:2])),
  pol_def = as.numeric(c(anova(modelrmald_poldef, btt = 4)[1:2])),
  nut_bef = as.numeric(c(anova(modelrmald_nutbef, btt = 4)[1:2])),
  nut_def = as.numeric(c(anova(modelrmald_nutdef, btt = 4)[1:2]))))
resanovaLD_study$Predictor <- rep("Study type", 4)
resanovaLD_study$Response <- rep(c("Decomposition (Div)", "Decomposition (Abd)"), 2)
resanovaLD_study$ECD.type <- c(rep("Stressors", 2), rep("Nutrients", 2))
names(resanovaLD_study)[c(1,2)] <- c("QM", "P")
print(resanovaLD_study)


write.csv(resanovaLD_study, "tables/ResMALDstudy.csv")


## 3. FIGURE Stressor and nutrient intensity effects-----------------------------------------

# This code creates Figure 6 of the manuscript

colo_lea <- c("#0072B2", "#D55E00")

# set sizes of plotted elements
sizetext <- 12
sizelegend <- 11
sizepoint <- 1
sizeline <- 0.8
sizeannotate <- 3.4


## Function to create panel for each outcome
# Panels LD
myfun_LD_ECD <- function(dat, plottitle, mycol){
  
  dat$weight = 1/dat$var.zcor.ECD.LD
  
  ## calculate slope B-EF with a meta-regression
  slope <- rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                  mods = ~ ECD_max + study_type,
                  random =  ~ 1 | Case.study/ID, 
                  data = dat)
  
  # calculate confidence interval around the slope
  predframe <- with(dat, 
                    data.frame(zcor.ECD.LD, ECD_max, 
                               preds = predict(slope)$pred, 
                               lwr = predict(slope)$ci.lb, 
                               upr = predict(slope)$ci.ub))
  
  # extract statistics
  QMstat <- round(anova(slope, btt=2)[1]$QM, 1)
  QMpval <- round(anova(slope, btt=2)[2]$QMp, 2)
  nstud <- slope$s.nlevels[1]
  nobs <- slope$k
  
  # bquote to annotate the plot with stats and sample sizes
  labelstatq <- bquote(QM[df == 1] == .(QMstat))
  labelstatp <- if(QMpval>0.001) {bquote(p ==.(QMpval))}else{ bquote(p < 0.001)}
  labelstatns <- bquote((.(nstud) ~ ";" ~ .(nobs)))
  labelstat <- bquote(list(.(labelstatq), .(labelstatp), .(labelstatns)))
  
  ## plot 
  ggplot(dat, aes(x=ECD_max, y=zcor.ECD.LD,
                  size = weight)) +
    
    geom_point(col = mycol, pch = 1)+ #, size = sizepoint) +
    
    # titles and axis labels
    ylab("Effect on decomposition")+
    xlab(if( dat$ECD.type[1]=="Stressors"){"Stressor intensity"} else{"Nutrient intensity"})+
    ggtitle(plottitle) +
    
    # axis lenght
    # ylim(c(-3.1,2.6))+
    xlim(c(-2.2,10.1))+
    
    # mark the zeros lines
    geom_hline(yintercept = 0, color = "black")+
    geom_vline(xintercept = 0, color = "black")+
    
    # annotate with stats and sample sizes
    annotate('label', x = 10.1, y = min(dat$zcor.ECD.LD)- 0.7, hjust = 1, 
             label=deparse(labelstat), parse = TRUE, size = sizeannotate, , fill = "white", label.size = NA)+
    
    # add slope if significant and conf intervals
    # add sloeps and conf intervals
    geom_abline(slope = slope$b[2], intercept = slope$b[1], col = mycol, size = sizeline,
                linetype = 1+ifelse(anova(slope, btt=2)[2]$QMp>0.05, 1, 0))+ # change lty according to p-value of the meta-regression
   
    
    # theme stff
    theme_bw() +
    theme(
      axis.text.y = element_text(face = "bold", size = sizetext),
      axis.text.x = element_text(face = "bold", size = sizetext),
      axis.title.x = element_text(size = sizetext, face = "bold"),
      axis.title.y = element_text(size = sizetext, face = "bold"),
      
      #legend
      plot.title = element_text(size = sizetext),
      
      legend.position = "none")
}



# Panels B
myfun_B_ECD <- function(dat, plottitle, mycol, xaxis){
  
  dat$weight = 1/dat$var.zcor.ECD.B
  
  ## calculate slope B-EF with a meta-regression
  if( plottitle=="Abundance"){
    slope <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                    mods = ~ ECD_max + study_type + taxonomic.group,
                    random =  ~ 1 | Case.study/ID, 
                    data = dat)}
  
  else{slope <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                       mods = ~ ECD_max + study_type + taxonomic.group + B.metric,
                       random =  ~ 1 | Case.study/ID, 
                       data = dat)}
  
  
  # calculate confidence interval around the slope
  predframe <- with(dat, 
                    data.frame(zcor.ECD.B, ECD_max, 
                               preds = predict(slope)$pred, 
                               lwr = predict(slope)$ci.lb, 
                               upr = predict(slope)$ci.ub))
  # extract statistics
  QMstat <- round(anova(slope, btt=2)[1]$QM, 1)
  QMpval <- round(anova(slope, btt=2)[2]$QMp, 2)
  nstud <- slope$s.nlevels[1]
  nobs <- slope$k
  
  # bquote to annotate the plot with stats and sample sizes
  labelstatq <- bquote(QM[df == 1] == .(QMstat))
  labelstatp <- if(QMpval>0.001) {bquote(p ==.(QMpval))}else{ bquote(p < 0.001)}
  labelstatns <- bquote((.(nstud) ~ ";" ~ .(nobs)))
  labelstat <- bquote(list(.(labelstatq), .(labelstatp), .(labelstatns)))
  
  ## plot 
  ggplot(dat, aes(x=ECD_max, y=zcor.ECD.B, size = weight)) +
    
    geom_point(col = mycol, pch = 1)+
    
    
    # titles and axis labels
    ylab("Effect on decomposers")+
    xlab(if( dat$ECD.type[1]=="Stressors"){"Stressor intensity"} else{"Nutrient intensity"})+
    ggtitle(plottitle) +
    
    # axis lenght
    # ylim(c(-3.1,2.6))+
    xlim(c(-2.2,10.1))+
    
    # mark the zeros lines
    geom_hline(yintercept = 0, color = "black")+
    geom_vline(xintercept = 0, color = "black")+
    
    # add sloeps and conf intervals
    geom_abline(slope = slope$b[2], intercept = slope$b[1], col = mycol, size = sizeline,
                linetype = 1+ifelse(anova(slope, btt=2)[2]$QMp>0.05, 1, 0))+ # change lty according to p-value of the meta-regression
    
    # annotate with stats and sample sizes
    annotate('label', x = 10.1, y = (min(dat$zcor.ECD.B)- 0.75), hjust = 1, 
             label=deparse(labelstat), parse = TRUE, size = sizeannotate, fill = "white", label.size = NA)+
    
    # theme stff
    theme_bw() +
    theme(
      axis.text.y = element_text(face = "bold", size = sizetext),
      axis.text.x = element_text(face = "bold", size = sizetext),
      axis.title.x = element_text(size = sizetext, face = "bold"),
      axis.title.y = element_text(size = sizetext, face = "bold"),
      
      #legend
      plot.title = element_text(size = sizetext),
      
      legend.position = "none")
}


Fig_Sup_ECDmax <-
  # upper panels for stressors
  myfun_B_ECD(pol_BEF_es, "Diversity", colo_lea[1])+
  myfun_B_ECD(pol_DEF_es, "Abundance", colo_lea[1])+
  myfun_LD_ECD(pol_DEF_es[!duplicated(pol_DEF_es$clusterID_LD),], "Decomposition", colo_lea[1]) +
  
  #lower panels for resources
  myfun_B_ECD(nut_BEF_es, "Diversity", colo_lea[2]) +
  myfun_B_ECD(nut_DEF_es, "Abundance", colo_lea[2]) +
  myfun_LD_ECD(nut_DEF_es[!duplicated(nut_DEF_es$clusterID_LD),], "Decomposition", colo_lea[2])

# Fig_Sup_ECDmax

# save a png with high res
ppi <- 300 #final: 600 # resolution
w <- 21 # width in cm

png("figs/Fig6_ECDintensity.png",
    width=w,
    height=w/1.5,
    units = "cm",
    res=ppi)

Fig_Sup_ECDmax

dev.off()

## 4. FIGURE Categorical moderators----------------------------------------------------------

## This script creates figure 7 - mean effect sizes on decomposer diversity and abundance 
# per level of categorical moderators

# For plotting purposes, separate meta-analyses are run to derive the mean effect sizes for the different datasets
# I used metafor package to derive proper confidence intervals for meta-analysis
# The complete dataset is used ( no need for data resampling in this univariate approach)


# a function to run meta-analysis on each categorical predictor and get the mean effect sizes, CI and stats
myfun_forestplot <- function(res){
  y<-res$b
  ci_l<-res$ci.lb
  ci_h<-res$ci.ub
  fgdf<-data.frame(cbind(y,ci_l,ci_h))
  colnames(fgdf)[1]<-"y"
  colnames(fgdf)[2]<-"ci_l"
  colnames(fgdf)[3]<-"ci_h"
  fgdf$catego_pred <- factor(rownames(res$b))
  return(fgdf)
}

## 4.1. STUDY TYPE
# extract mean effect sizes and CI for levels of study type (observational versus experimental studies)
study_pol_bef <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, # Stressors - Biodiv dataset
                        mods = ~ study_type - 1,
                        random = ~ 1 | Case.study / ID, 
                        data = pol_BEF_es)
res.study_pol_bef <- myfun_forestplot(study_pol_bef)
res.study_pol_bef$no.stu <- as.numeric(countstudies(pol_BEF_es, study_type)$no.stu)
res.study_pol_bef$no.obs <- as.numeric(countstudies(pol_BEF_es, study_type)$no.obs)

study_pol_def <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, # Stressors - Abdc dataset
                        mods = ~ study_type - 1,
                        random = ~ 1 | Case.study / ID, 
                        data = pol_DEF_es)
res.study_pol_def <- myfun_forestplot(study_pol_def)
res.study_pol_def$no.stu <- as.numeric(countstudies(pol_DEF_es, study_type)$no.stu)
res.study_pol_def$no.obs <- as.numeric(countstudies(pol_DEF_es, study_type)$no.obs)

study_nut_bef <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, # Nutrients - Biodiv dataset
                        mods = ~ study_type - 1,
                        random = ~ 1 | Case.study / ID, 
                        data = nut_BEF_es)
res.study_nut_bef <- myfun_forestplot(study_nut_bef)
res.study_nut_bef$no.stu <- as.numeric(countstudies(nut_BEF_es, study_type)$no.stu)
res.study_nut_bef$no.obs <- as.numeric(countstudies(nut_BEF_es, study_type)$no.obs)

study_nut_def <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, # Nutrients - Abdc dataset
                        mods = ~ study_type - 1,
                        random = ~ 1 | Case.study / ID, 
                        data = nut_DEF_es)
res.study_nut_def <- myfun_forestplot(study_nut_def)
res.study_nut_def$no.stu <- as.numeric(countstudies(nut_DEF_es, study_type)$no.stu)
res.study_nut_def$no.obs <- as.numeric(countstudies(nut_DEF_es, study_type)$no.obs)



## 4.2. TAXO GROUP
# extract mean effect sizes and CI for levels of taxonomic group (animal versus microbial decomposers)
taxo_pol_bef <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                       mods = ~ taxonomic.group - 1,
                       random = ~ 1 | Case.study / ID, 
                       data = pol_BEF_es)
res.taxo_pol_bef <- myfun_forestplot(taxo_pol_bef)
res.taxo_pol_bef$no.stu <- as.numeric(countstudies(pol_BEF_es, taxonomic.group)$no.stu)
res.taxo_pol_bef$no.obs <- as.numeric(countstudies(pol_BEF_es, taxonomic.group)$no.obs)

taxo_pol_def <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                       mods = ~ taxonomic.group - 1,
                       random = ~ 1 | Case.study / ID, 
                       data = pol_DEF_es)
res.taxo_pol_def <- myfun_forestplot(taxo_pol_def)
res.taxo_pol_def$no.stu <- as.numeric(countstudies(pol_DEF_es, taxonomic.group)$no.stu)
res.taxo_pol_def$no.obs <- as.numeric(countstudies(pol_DEF_es, taxonomic.group)$no.obs)

taxo_nut_bef <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                       mods = ~ taxonomic.group - 1,
                       random = ~ 1 | Case.study / ID, 
                       data = nut_BEF_es)
res.taxo_nut_bef <- myfun_forestplot(taxo_nut_bef)
res.taxo_nut_bef$no.stu <- as.numeric(countstudies(nut_BEF_es, taxonomic.group)$no.stu)
res.taxo_nut_bef$no.obs <- as.numeric(countstudies(nut_BEF_es, taxonomic.group)$no.obs)

taxo_nut_def <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                       mods = ~ taxonomic.group - 1,
                       random = ~ 1 | Case.study / ID, 
                       data = nut_DEF_es)
res.taxo_nut_def <- myfun_forestplot(taxo_nut_def)
res.taxo_nut_def$no.stu <- as.numeric(countstudies(nut_DEF_es, taxonomic.group)$no.stu)
res.taxo_nut_def$no.obs <- as.numeric(countstudies(nut_DEF_es, taxonomic.group)$no.obs)



## Store the results for each dataset
Res.forest_pol_bef <- rbind(res.study_pol_bef, res.taxo_pol_bef)
Res.forest_pol_bef$categomods <- factor(ifelse(grepl("expe", rownames(Res.forest_pol_bef)), "Experimental", 
                                               ifelse(grepl("obser", rownames(Res.forest_pol_bef)), "Observational",
                                                      ifelse(grepl("groupA", rownames(Res.forest_pol_bef)), "Animals", "Microbes"))))

Res.forest_pol_def <- rbind(res.study_pol_def, res.taxo_pol_def)
Res.forest_pol_def$categomods <- factor(ifelse(grepl("expe", rownames(Res.forest_pol_def)), "Experimental", 
                                               ifelse(grepl("obser", rownames(Res.forest_pol_def)), "Observational",
                                                      ifelse(grepl("groupA", rownames(Res.forest_pol_def)), "Animals", "Microbes"))))

Res.forest_nut_bef <- rbind(res.study_nut_bef, res.taxo_nut_bef)
Res.forest_nut_bef$categomods <- factor(ifelse(grepl("expe", rownames(Res.forest_nut_bef)), "Experimental", 
                                               ifelse(grepl("obser", rownames(Res.forest_nut_bef)), "Observational",
                                                      ifelse(grepl("groupA", rownames(Res.forest_nut_bef)), "Animals", "Microbes"))))

Res.forest_nut_def <- rbind(res.study_nut_def, res.taxo_nut_def)
Res.forest_nut_def$categomods <- factor(ifelse(grepl("expe", rownames(Res.forest_nut_def)), "Experimental", 
                                               ifelse(grepl("obser", rownames(Res.forest_nut_def)), "Observational",
                                                      ifelse(grepl("groupA", rownames(Res.forest_nut_def)), "Animals", "Microbes"))))



## FIGURE - make a Forest plots
# pick colors
colo_lea <- c("#0072B2", "#D55E00")

# set sizes of plotted elements (from Fig1)
sizetext <- 13
sizepoint <- 3
widtherrorbar <- 0.1
sizeerrorbar <- 0.4


## function to make a forest ggplot
myfun_Forestggplot_B2 <- function(df,  plottitle, mycol){
  
  # reorder factor levels
  df$categomods2 <- factor(df$categomods, c("Observational", "Experimental",  "Microbes", "Animals"))
  
  
  # make plot
  ggplot(df, aes(x=categomods2, y=y, shape = categomods2))+ 
    
    # error bars are conf intervals 95%
    geom_errorbar(width=widtherrorbar, 
                  size = sizeerrorbar, 
                  aes(ymin = df$ci_l,
                      ymax = df$ci_h),           # confidence intervals from Std err of models
                  col = mycol) + 
    
    # points shape and colors
    geom_point(size= sizepoint, col = mycol, fill = mycol)+
    
    # change their shape (pch)
    scale_shape_manual(values=c(17,2,19,1))+  # Use a hollow circle and triangle
    
    # axis
    ylim(-1.2, 1.2)+
    ylab("Effect size")+
    xlab(" ")+
    
    # flip the coordinates to make forest plot
    coord_flip()+
    
    # add lines 
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 2.5, lty= 2)+
    
    
    # add no. studies and observation to axis labels
    scale_x_discrete(breaks=c("Observational", "Experimental", "Microbes", "Animals"),
                     labels=c(paste("Obs. (",df$no.stu[2], "; ", df$no.obs[2], ")", sep = ""),
                              paste("Expe. (", df$no.stu[1], "; ", df$no.obs[1],")", sep = ""),
                              paste("Microbes (", df$no.stu[4], "; ", df$no.obs[4],")", sep = ""),
                              paste("Animals (", df$no.stu[3], "; ", df$no.obs[3],")", sep = ""))) +
    
    # theme and design
    theme_bw() +
    ggtitle(plottitle) +
    
    
    
    # theme stff
    theme(axis.text.y=element_text(face = "bold", size = sizetext), 
          axis.text.x=element_text(face = "bold", size = sizetext),
          axis.title.x = element_text(size=sizetext, face = "bold"),
          axis.title.y = element_text(size=sizetext, face = "bold"),
          plot.title = element_text(size = sizetext), 
          legend.title = element_text(size = sizetext),
          legend.position = "none",
          legend.text = element_text(size = sizetext))
  
  
  
  
}

# Plots for the 4 datasets
Fig_catmods <-
  myfun_Forestggplot_B2(df = Res.forest_pol_bef,
                        plottitle = "Stressors - Biodiversity", 
                        mycol = colo_lea[1]  )+
  
  myfun_Forestggplot_B2(df = Res.forest_pol_def, 
                        plottitle = "Stressors - Abundance", 
                        mycol = colo_lea[1]  )+
  myfun_Forestggplot_B2(df = Res.forest_nut_bef,
                        plottitle = "Nutrients - Biodiversity", 
                        mycol = colo_lea[2]  )+
  
  myfun_Forestggplot_B2(df = Res.forest_nut_def, 
                        plottitle = "Nutrients - Abundance", 
                        mycol = colo_lea[2]  )

# Fig_catmods

# save a png with high res
ppi <- 300 # 600 final resolution
w <- 21 # width in cm

png("figs/Fig7_CategoMods.png",
    width=w,
    height=w/1.3,
    units = "cm",
    res=ppi)

Fig_catmods

dev.off()
