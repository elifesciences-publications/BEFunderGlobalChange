####################################################################################################################
##
## LOAD the data for Structural equation models (SEM) of biodiversity mediated effects of stressors and nutrients on decomposition
##
###################################################################################################################


## This scripts loads the effect size data necessary to recreate the SEM results 

rm(list=ls())

## 1. LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------

library(dplyr)
library(ggplot2)
library(metafor)
library(broom)
library(knitr)
library(tidyverse)
library(devtools)
library(nlme)
library(car)
library(patchwork) #devtools::install_github("thomasp85/patchwork")
library(boot) # required to do bootstrap
library("ggsci") # nice color palettes : https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
library(lattice)
library(reshape2)
library("gridExtra") # for grid.arrange
library(piecewiseSEM)


## notin operator
'%notin%' <- function(x,y)!('%in%'(x,y))

# count no. studies and observations per level of categorical moderator
countstudies <- function(dat, ...){
  return(dat %>% 
           group_by_(.dots = lazyeval::lazy_dots(...)) %>% 
           summarize(no.stu = countlevels(Case.study), no.obs = n()))
}

## count levels of a non factor variable
countlevels <- function(x){
  return(length(levels(as.factor(x))))
}

## dplyr::bind rows transfo characters back to factors
## https://stackoverflow.com/questions/42278020/bind-rows-of-data-frames-with-some-factor-columns#42278444
bind_rows_keep_factors <- function(...) {
  ## Identify all factors
  factors <- unique(unlist(
    map(list(...), ~ select_if(..., is.factor) %>% names())
  ))
  ## Bind dataframes, convert characters back to factors
  suppressWarnings(bind_rows(...)) %>% 
    mutate_at(vars(one_of(factors)), factor)  
}

## Function to check the residuals of rma.mv visually
function_resid <- function(res, plottitle){
  par(mfrow = c(1,3))
  qqnorm(residuals(res,type="pearson"),main=plottitle)
  qqline(residuals(res,type="pearson"),col="red")
  hist(residuals(res, type = "pearson"), main =plottitle)
  plot(residuals(res, type = "pearson") ~ fitted(res), main = paste("Residuals vs Fitted", plottitle)); abline(0,0)
}

# check collinearity with pairs function
## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


## Function for mean and CI
myfun_meanci <- function(x){
  data.frame(mean.est = mean(x, na.rm = TRUE),
             ci.lwr = quantile(x, 0.025, na.rm = TRUE),
             ci.upr = quantile(x, 0.975, na.rm = TRUE))
  
}

# calculate pvalue 0.05 associated with an observed chisq value
myfun_chisqpval <- function(chi,df){ 1-pchisq(chi, df)} 



## 2. LOAD DATA---------------------------------------------------------------------------


## This script loads the effect size datasets, and prepare the data for subsequent meta-analysis


## Load the datasets abundance-decomposition (DandEF1) and biodiv-decomposition (BandEF1)
effect_size_BandEF1 = read.csv("data/effect_size_BandEF1.csv", h = TRUE, sep = ",")
effect_size_DandEF1 = read.csv("data/effect_size_DandEF1.csv", h = TRUE, sep = ",")

# rename effect sizes_D into B to facilitate scripts: 
effect_size_DandEF1$zcor.ECD.B <- effect_size_DandEF1$zcor.ECD.D
effect_size_DandEF1$zcor.ECD.LD <- effect_size_DandEF1$zcor.ECD.LDd

# rename variances
effect_size_DandEF1$var.zcor.ECD.B <- effect_size_DandEF1$var.zcor.ECD.D
effect_size_DandEF1$var.zcor.ECD.LD <- effect_size_DandEF1$var.zcor.ECD.LDd

# short name ECD max 
effect_size_BandEF1$ECD_max <- effect_size_BandEF1$ECD_intensity_norm_max
effect_size_DandEF1$ECD_max <- effect_size_DandEF1$ECD_intensity_norm_max

# remove NA in DandEF dataset
effect_size_DandEF1 <- effect_size_DandEF1 %>% filter(!is.na(ECD_max))# remove NA for ECD max (study 395, 1 obs)
effect_size_DandEF1 <- effect_size_DandEF1 %>% filter(!is.na(var.zcor.ECD.B))# remove NA for ECD max (study 395, 1 obs)


## Prep data -----------------------------------------------------------------------------

# Type of environmental change driver (ECD)
# Create a categorical variable for type of ECD: Resources stands for nutrient enrichment studies and Stressors for chemical stressors
summary(effect_size_BandEF1$ECD_subcategory)
effect_size_BandEF1 = dplyr::mutate(effect_size_BandEF1, ECD.type = factor(ifelse(
  grepl("eutroph|nitro|nutri|wastewater", ECD_subcategory), "Resources", "Stressors")))

summary(effect_size_DandEF1$ECD_subcategory)
effect_size_DandEF1 = dplyr::mutate(effect_size_DandEF1, ECD.type = factor(ifelse(
  grepl("eutroph|nitro|nutri|wastewater", ECD_subcategory), "Resources", "Stressors")))

## Observation ID
effect_size_BandEF1$studyID <- factor(effect_size_BandEF1$studyID)
effect_size_DandEF1$studyID <- factor(effect_size_DandEF1$studyID)

# assign unique ID to each observation
effect_size_BandEF1$ID <-  factor(1 : nrow(effect_size_BandEF1))
effect_size_DandEF1$ID <-  factor(1 : nrow(effect_size_DandEF1))

## Duplicates ID
# create identifier for duplicated effect sizes on decomposition in each Case.study: cluster ID.
# observations with the same clusterID are duplicates
effect_size_DandEF1 <- transform(
  
  effect_size_DandEF1, 
  
  clusterID_LD = factor(
    interaction(Case.study, zcor.ECD.LD, drop=TRUE)), 
  
  clusterID_B = factor(
    interaction(Case.study, zcor.ECD.B, drop = TRUE)))

# create identifier for duplicated effect sizes on decomposition in each Case.study: cluster ID.
effect_size_BandEF1 <- transform(
  
  effect_size_BandEF1, 
  
  clusterID_LD = factor(
    interaction(Case.study, zcor.ECD.LD, drop=TRUE)), 
  
  clusterID_B = factor(
    interaction(Case.study, zcor.ECD.B, drop = TRUE)))

# check that clusterID identifies the duplicated effect sizes on LD
all.equal(
  subset(effect_size_BandEF1, !duplicated(zcor.ECD.LD)), ## identical zLD values
  subset(effect_size_BandEF1, !duplicated(clusterID_LD)) ## clusterID_LD identifies duplicates
)

all.equal(
  subset(effect_size_DandEF1, !duplicated(zcor.ECD.LD)), # and for abundance data
  subset(effect_size_DandEF1, !duplicated(clusterID_LD))
)


# Biodiversity metric
# there are not enough studies to test shannon and evenness separately in the meta-analysis
# create a variable separating biodiversity observations in terms of taxa richness (richness) and of diversity indices (shannon and evenness)
effect_size_BandEF1 <- effect_size_BandEF1 %>% mutate(
  B.metric = factor(
    ifelse(B.metric == "taxa richness", "richness", "div_indices")))


# dataset for stressors-------------------------------------------------------------------

## Biodiversity data

pol_BEF_es <- effect_size_BandEF1 %>% filter(ECD.type == "Stressors")
pol_BEF_es <- droplevels(pol_BEF_es)
nrow(pol_BEF_es) # 100 effect sizes
countlevels(pol_BEF_es$studyID) # 21 papers
countlevels(pol_BEF_es$Case.study) # 23 studies

## Abundance data
pol_DEF_es <- effect_size_DandEF1 %>% filter(ECD.type == "Stressors")
pol_DEF_es <- droplevels(pol_DEF_es)
nrow(pol_DEF_es) # 282 effect sizes
countlevels(pol_DEF_es$studyID) # 25 studies
countlevels(pol_DEF_es$Case.study) # 28 studies

# rename D.metric variable to facilitate coding
pol_DEF_es$B.metric <- pol_DEF_es$D.metric


# dataset for nutrients----------------------------------------------------------------
# biodiv data
nut_BEF_es <- effect_size_BandEF1 %>% filter(ECD.type == "Resources")
nut_BEF_es <- droplevels(nut_BEF_es)

nrow(nut_BEF_es) # 97 effect sizes
countlevels(nut_BEF_es$studyID) # 21 studies
countlevels(nut_BEF_es$Case.study) # 26 cases

# density data
nut_DEF_es <- effect_size_DandEF1 %>% filter(ECD.type == "Resources")
nut_DEF_es <- droplevels(nut_DEF_es)

nrow(nut_DEF_es) # 181 effect sizes
countlevels(nut_DEF_es$studyID) # 32 studies
countlevels(nut_DEF_es$Case.study) # 39 cases


