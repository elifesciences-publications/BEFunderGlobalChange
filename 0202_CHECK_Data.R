####################################################################################################################
##
## CHECK the data : publication bias, independence of covariates etc.
##
###################################################################################################################


## This scripts loads and check the effect size data for the SEM analysis. 
# It contains test for publication bias, extreme values, and independence between covariates
# It creates supplementary tables and figures summarising the publication bias results.

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

## 3. CHECK DATA ---------------------------------------------------------------------

# Check for publication bias, outliers and severe non-normality, and independence of covariates.

# Main covariates: 
# * study type : experimental or observational studies 
# * system type : terrestrial or aquatic ecosystems
# * ECD intensity and magnitude : difference between the normalized chemical stressors intensity of the highest and the lowest impacted sites (i.e. size of the gradient or contrast in a common scale)
# * litter quality : labile litter (forbs and woody deciduous species (except oaks)) versus recalcitrant (graminoids and grass-like monocots; woody evergreen species; any woody material)
# * taxonomic group : animals or microbes
# * biodiversity metric : taxa richness or diversity indices (Shannon and evenness indices together because only one study reported evenness in this dataset)

# 3.1 STRESSORS------------------------------------------------------------------------

# 3.1.1. Publication bias--------------------------------------------------------------

## FIGURE Funnel plots-------------------------------------------------------------

# including ECD magnitude as covariate 
## models for diversity data
res.B1 <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                 random = ~ 1|Case.study/ID, 
                 mods = ~ ECD_max,
                 data = pol_BEF_es)
res.LD1 <- rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                  random = ~ 1|Case.study/ID, 
                  mods = ~ ECD_max,
                  data = subset(pol_BEF_es, !duplicated(clusterID_LD)))  ## remove duplicates 

## models for abundance data
res.B2 <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                 random = ~ 1|Case.study/ID, 
                 mods = ~ ECD_max,
                 data = pol_DEF_es)
res.LD2 <- rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                  random = ~ 1|Case.study/ID, 
                  mods = ~ ECD_max,
                  data =  subset(pol_DEF_es, !duplicated(clusterID_LD)))  ## remove duplicates 

### FIGURE - funnel plots
par(mfrow=c(2,2))
funnel(res.B1, main = "Stressors - Diversity")
funnel(res.LD1, main = "Stressors - Decomposition (Div)")

funnel(res.B2, main = "Stressors - Abundance ")
funnel(res.LD2, main = "Stressors - Decomposition (Abd)")

## Save for SUpplementary figs

# save a png with high res
ppi <- 600 # resolution
w <- 21 # width in cm

png("figs/FigS1_FunnelPlots_pol.png",
    width=w,
    height=w/1.3,
    units = "cm",
    res=ppi)

par(mfrow=c(2,2))
funnel(res.B1, main = "Stressors - Diversity")
funnel(res.LD1, main = "Stressors - Decomposition (Div)")

funnel(res.B2, main = "Stressors - Abundance ")
funnel(res.LD2, main = "Stressors - Decomposition (Abd)")

dev.off()


## TABLES Egger's regressions------------------------------------------------------------
# test for funnel plot assymetry

#### Biodiv dataset
dat <- pol_BEF_es

# standardized effect sizes (E*)
dat$std.es.LD <- dat$zcor.ECD.LD/sqrt(dat$var.zcor.ECD.LD)
dat$std.es.B <- dat$zcor.ECD.B/sqrt(dat$var.zcor.ECD.B)

# standard errors (SE)
dat$se.es.LD <- sqrt(dat$var.zcor.ECD.LD)
dat$se.es.B <-  sqrt(dat$var.zcor.ECD.B)

# precision (1/SE)
dat$precision.es.LD <- 1/dat$se.es.LD
dat$precision.es.B <- 1/dat$se.es.B

# Egger's regressions : E* vs. 1/SE and incoporating the main covariate ECD magnitude.
res.eggerLD <- rma.mv(std.es.LD, var.zcor.ECD.LD, 
                      random =  ~ 1 | Case.study/ID,
                      data = subset(dat, !duplicated(clusterID_LD)),
                      mods =~ ECD_max + precision.es.LD)

res.eggerB <- rma.mv(std.es.B, var.zcor.ECD.B, 
                     random =  ~ 1 | Case.study/ID,
                     data = dat,
                     mods =~ ECD_max + precision.es.B)
rm(dat)

#### Abundance data
dat <- pol_DEF_es

# standardized effect sizes (E*)
dat$std.es.LD <- dat$zcor.ECD.LD/sqrt(dat$var.zcor.ECD.LD)
dat$std.es.B <- dat$zcor.ECD.B/sqrt(dat$var.zcor.ECD.B)

# standard errors (SE)
dat$se.es.LD <- sqrt(dat$var.zcor.ECD.LD)
dat$se.es.B <-  sqrt(dat$var.zcor.ECD.B)

# precision (1/SE)
dat$precision.es.LD <- 1/dat$se.es.LD
dat$precision.es.B <- 1/dat$se.es.B

# Egger's regressions : E* vs. 1/SE and incoporating the main covariate ECD magnitude.
res.eggerLD2 <- rma.mv(std.es.LD, var.zcor.ECD.LD, 
                       random =  ~ 1 | Case.study/ID,
                       data = subset(dat, !duplicated(clusterID_LD)),
                       mods =~ ECD_max + precision.es.LD)

res.eggerB2 <- rma.mv(std.es.B, var.zcor.ECD.B, 
                      random =  ~ 1 | Case.study/ID,
                      data = dat,
                      mods =~ ECD_max + precision.es.B)
rm(dat)


# Get Results :
# extract pvalues
pub.bias.pol <- c(res.eggerB$pval[1],res.eggerLD$pval[1],
                  res.eggerB2$pval[1],res.eggerLD2$pval[1])

# results: If there is publication bias, intercept of ES ~ sd(ES) should be different than 0 


df.pubbias <- data.frame(dataset = c(rep("BEF", 2), rep("DEF",2)),
                         variable = rep(c("B", "LD"), 2),
                         PublicationBias_Pvalue = pub.bias.pol) %>% 
  mutate(PublicationBias = ifelse(pub.bias.pol > 0.05, "no", "pub. bias"))



# Save table supplementary information
intercept.pub.bias <- c(res.eggerB$b[1],res.eggerLD$b[1],
                        res.eggerB2$b[1],res.eggerLD2$b[1])
se.pub.bias <- c(res.eggerB$se[1],res.eggerLD$se[1],
                 res.eggerB2$se[1],res.eggerLD2$se[1])

TableS0 <- data.frame(df.pubbias, Intercept = intercept.pub.bias, SE = se.pub.bias)

TableS0

write.csv(TableS0, "tables/TableS0_EggersRegressions_pol.csv")

## 3.1.2. Distribution : outliers and normality------------------------------------------

### Boxplots and histograms

### Identify extreme values of effect sizes and of stressor intensity (ECD_max) in the biodiversity dataset
par(mfrow = c(3,2))
boxplot(pol_BEF_es$zcor.ECD.B, main = "ES Div")
hist(pol_BEF_es$zcor.ECD.B, main = "ES Div")

boxplot(subset(pol_BEF_es, !duplicated(clusterID_LD))$zcor.ECD.LD, main = "ES decomposition")
hist(subset(pol_BEF_es, !duplicated(clusterID_LD))$zcor.ECD.LD, main = "ES decomposition")

boxplot(subset(pol_BEF_es, !duplicated(clusterID_LD))$ECD_max, main = "ECD intensity")
hist(subset(pol_BEF_es, !duplicated(clusterID_LD))$ECD_max, main = "ECD intensity")

### Identify extreme values of effect sizes and of stressor intensity (ECD_max) in the abundance dataset
par(mfrow = c(3,2))
boxplot(pol_DEF_es$zcor.ECD.B, main = "ES Abundance")
hist(pol_DEF_es$zcor.ECD.B, main = "ES Abundance")

boxplot(subset(pol_DEF_es, !duplicated(clusterID_LD))$zcor.ECD.LD, main = "ES decompo")
hist(subset(pol_DEF_es, !duplicated(clusterID_LD))$zcor.ECD.LD, main = "ES decompo")

boxplot(subset(pol_DEF_es, !duplicated(clusterID_LD))$ECD_max, main = "ECD intensity")
hist(subset(pol_DEF_es, !duplicated(clusterID_LD))$ECD_max, main = "ECD intensity")


## function to identify extreme values (i.e. values beyond whiskers of boxplots)
identify_outliers <- function(dat, es){
  return(subset(dat, es %in% boxplot.stats(es)$out))
}

# count no. outlier per dataset
df.out <- data.frame(df.pubbias,
                     n.outliers = c(
                       nrow(identify_outliers(pol_BEF_es, pol_BEF_es$zcor.ECD.B)),
                       nrow(identify_outliers(pol_BEF_es, pol_BEF_es$zcor.ECD.LD)), 
                       nrow(identify_outliers(pol_DEF_es, pol_DEF_es$zcor.ECD.B)), 
                       nrow(identify_outliers(pol_DEF_es, pol_DEF_es$zcor.ECD.LD)))
)

df.out

list(
  # inpsect observations that have extreme values of effect sizes in the biodiversiy dataset
  BEF_B = select(
    identify_outliers(pol_BEF_es, pol_BEF_es$zcor.ECD.B), # on biodiversity
    Case.study, zcor.ECD.B, var.zcor.ECD.B, study_type, system, ECD_subcategory),
  
  
  BEF_LD =  select(
    identify_outliers(pol_BEF_es, pol_BEF_es$zcor.ECD.LD), # on decomposition
    Case.study, zcor.ECD.LD, var.zcor.ECD.LD, study_type, system, ECD_subcategory),
  
  
  # and abundance dataset
  DEF_D = select(identify_outliers(pol_DEF_es, pol_DEF_es$zcor.ECD.B), # effects on abundance
                 Case.study, zcor.ECD.B, var.zcor.ECD.B, study_type, system, ECD_subcategory),
  
  DEF_LD = select(
    identify_outliers(pol_DEF_es, pol_DEF_es$zcor.ECD.LD) , # extreme effects on decomposition
    Case.study, zcor.ECD.LD, var.zcor.ECD.LD, study_type, system, ECD_subcategory)
)


## 3.1.3. Independence of covariates-----------------------------------------------------

### Number of studies and observations per level of categorical moderators

## Biodiv dataset
dat <- subset(pol_BEF_es, !duplicated(clusterID_LD))# remove duplicates
countstudies(dat, study_type)
countstudies(dat, system)
countstudies(dat, taxonomic.group)
countstudies(dat, B.metric)

## Abundance  dataset
dat <- subset(pol_DEF_es, !duplicated(clusterID_LD))
countstudies(dat, study_type)
countstudies(dat, system)
countstudies(dat, taxonomic.group)
countstudies(dat, D.metric)

### Number of studies and observations per level of crossed categorical moderators

## Biodiv dataset
dat <- subset(pol_BEF_es, !duplicated(clusterID_LD))# remove duplicates
countstudies(dat, study_type, system)
countstudies(dat, study_type, taxonomic.group)
countstudies(dat, study_type, B.metric)

countstudies(dat, system, taxonomic.group)
countstudies(dat, system, B.metric)

countstudies(dat, taxonomic.group, B.metric)

## Abundance  dataset
dat <- subset(pol_DEF_es, !duplicated(clusterID_LD))
countstudies(dat, study_type, system)
countstudies(dat, study_type, taxonomic.group)
countstudies(dat, study_type, D.metric)

countstudies(dat, system, taxonomic.group)
countstudies(dat, system, D.metric)

countstudies(dat, taxonomic.group, D.metric)

### 3.1.4. Relationships-----------------------------------------------------------------
dat.pairplot.bef <- pol_BEF_es %>%  select(
  zcor.ECD.LD, zcor.ECD.B,
  ECD_magnitude, ECD_max,
  study_type, system, taxonomic.group, B.metric,
  ECD_mix)
pairs(
  dat.pairplot.bef,
  upper.panel = panel.smooth,
  lower.panel = panel.cor,
  diag.panel = panel.hist,
  main = "Biodiversity dataset"
)


## 3.2. NUTRIENTS-------------------------------------------------------

# 3.2.1. Publication bias------------------------------------------------------

## FIGURE Funnel plots------------------------------------------------------

# models including ECD magnitude as covariate 
res.B1 <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                 random = ~ 1|Case.study/ID, 
                 mods = ~ ECD_max,
                 data = nut_BEF_es)
res.LD1 <- rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                  random = ~ 1|Case.study/ID, 
                  mods = ~ ECD_max,
                  data = subset(nut_BEF_es, !duplicated(clusterID_LD)))  ## remove duplicates 

res.B2 <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                 random = ~ 1|Case.study/ID, 
                 mods = ~ ECD_max,
                 data = nut_DEF_es)
res.LD2 <- rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                  random = ~ 1|Case.study/ID, 
                  mods = ~ ECD_max,
                  data =  subset(nut_DEF_es, !duplicated(clusterID_LD)))  ## remove duplicates 

### FIGURE - funnel plots

par(mfrow=c(2,2))
funnel(res.B1, main = "Nutrients - Diversity")
funnel(res.LD1, main = "Nutrients - Decomposition (Div)")

funnel(res.B2, main = "Nutrients - Abundance ")
funnel(res.LD2, main = "Nutrients - Decomposition (Abd)")


## Save for SUpplementary figs

# save a png with high res
ppi <- 600 # resolution
w <- 21 # width in cm

png("figs/FigS1_FunnelPlots_nut.png",
    width=w,
    height=w/1.3,
    units = "cm",
    res=ppi)

par(mfrow=c(2,2))
funnel(res.B1, main = "Nutrients - Diversity")
funnel(res.LD1, main = "Nutrients - Decomposition (Div)")

funnel(res.B2, main = "Nutrients - Abundance ")
funnel(res.LD2, main = "Nutrients - Decomposition (Abd)")

dev.off()


### TABLES Egger's regression------------------------------------------------------


#### Biodiv dataset
# test for funnel plot assymetry

dat <- nut_BEF_es


# standardized effect sizes (E*)
dat$std.es.LD <- dat$zcor.ECD.LD/sqrt(dat$var.zcor.ECD.LD)
dat$std.es.B <- dat$zcor.ECD.B/sqrt(dat$var.zcor.ECD.B)

# standard errors (SE)
dat$se.es.LD <- sqrt(dat$var.zcor.ECD.LD)
dat$se.es.B <-  sqrt(dat$var.zcor.ECD.B)

# precision (1/SE)
dat$precision.es.LD <- 1/dat$se.es.LD
dat$precision.es.B <- 1/dat$se.es.B

# Egger's regressions : E* vs. 1/SE and incoporating the main covariate ECD magnitude.
res.eggerLD <- rma.mv(std.es.LD, var.zcor.ECD.LD, 
                      random =  ~ 1 | Case.study/ID,
                      data = subset(dat, !duplicated(clusterID_LD)),
                      mods =~ ECD_max + precision.es.LD)

res.eggerB <- rma.mv(std.es.B, var.zcor.ECD.B, 
                     random =  ~ 1 | Case.study/ID,
                     data = dat,
                     mods =~ ECD_max + precision.es.B)
rm(dat)

#### Abundance dataset

dat <- nut_DEF_es

# standardized effect sizes (E*)
dat$std.es.LD <- dat$zcor.ECD.LD/sqrt(dat$var.zcor.ECD.LD)
dat$std.es.B <- dat$zcor.ECD.B/sqrt(dat$var.zcor.ECD.B)

# standard errors (SE)
dat$se.es.LD <- sqrt(dat$var.zcor.ECD.LD)
dat$se.es.B <-  sqrt(dat$var.zcor.ECD.B)

# precision (1/SE)
dat$precision.es.LD <- 1/dat$se.es.LD
dat$precision.es.B <- 1/dat$se.es.B

# Egger's regressions : E* vs. 1/SE and incoporating the main covariate ECD magnitude.
res.eggerLD2 <- rma.mv(std.es.LD, var.zcor.ECD.LD, 
                       random =  ~ 1 | Case.study/ID,
                       data = subset(dat, !duplicated(clusterID_LD)),
                       mods =~ ECD_max + precision.es.LD)

res.eggerB2 <- rma.mv(std.es.B, var.zcor.ECD.B, 
                      random =  ~ 1 | Case.study/ID,
                      data = dat,
                      mods =~ ECD_max + precision.es.B)
rm(dat)


# Get Results :
# extract pvalues
pub.bias.nut <- c(res.eggerB$pval[1],res.eggerLD$pval[1],
                  res.eggerB2$pval[1],res.eggerLD2$pval[1])

# results: If there is publication bias, intercept of ES ~ sd(ES) should be different than 0 


df.pubbias <- data.frame(dataset = c(rep("BEF", 2), rep("DEF",2)),
                         variable = rep(c("B", "LD"), 2),
                         PublicationBias_Pvalue = pub.bias.nut) %>% 
  mutate(PublicationBias = ifelse(pub.bias.nut > 0.05, "no", "pub. bias"))



# Save table supplementary information
intercept.pub.bias <- c(res.eggerB$b[1],res.eggerLD$b[1],
                        res.eggerB2$b[1],res.eggerLD2$b[1])
se.pub.bias <- c(res.eggerB$se[1],res.eggerLD$se[1],
                 res.eggerB2$se[1],res.eggerLD2$se[1])

TableS0 <- data.frame(df.pubbias, Intercept = intercept.pub.bias, SE = se.pub.bias)

TableS0

write.csv(TableS0, "tables/TableS0_EggersRegressions_nut.csv")


## 3.2.2. Distribution : outliers and normality

### Boxplots and histograms

### Identify extreme values of effect sizes and of stressor intensity (ECD_max) in the biodiversity dataset
par(mfrow = c(3,2))
boxplot(nut_BEF_es$zcor.ECD.B, main = "ES Div")
hist(nut_BEF_es$zcor.ECD.B, main = "ES Div")

boxplot(subset(nut_BEF_es, !duplicated(clusterID_LD))$zcor.ECD.LD, main = "ES decomposition")
hist(subset(nut_BEF_es, !duplicated(clusterID_LD))$zcor.ECD.LD, main = "ES decomposition")

boxplot(subset(nut_BEF_es, !duplicated(clusterID_LD))$ECD_max, main = "ECD intensity")
hist(subset(nut_BEF_es, !duplicated(clusterID_LD))$ECD_max, main = "ECD intensity")



### Identify extreme values of effect sizes and of stressor intensity (ECD_max) in the abundance dataset
par(mfrow = c(3,2))
boxplot(nut_DEF_es$zcor.ECD.B, main = "ES Abundance")
hist(nut_DEF_es$zcor.ECD.B, main = "ES Abundance")

boxplot(subset(nut_DEF_es, !duplicated(clusterID_LD))$zcor.ECD.LD, main = "ES decompo")
hist(subset(nut_DEF_es, !duplicated(clusterID_LD))$zcor.ECD.LD, main = "ES decompo")

boxplot(subset(nut_DEF_es, !duplicated(clusterID_LD))$ECD_max, main = "ECD intensity")
hist(subset(nut_DEF_es, !duplicated(clusterID_LD))$ECD_max, main = "ECD intensity")

### Number of extreme values
# count no. outlier per dataset
df.out <- data.frame(df.pubbias,
                     n.outliers = c(
                       nrow(identify_outliers(nut_BEF_es, nut_BEF_es$zcor.ECD.B)),
                       nrow(identify_outliers(nut_BEF_es, nut_BEF_es$zcor.ECD.LD)), 
                       nrow(identify_outliers(nut_DEF_es, nut_DEF_es$zcor.ECD.B)), 
                       nrow(identify_outliers(nut_DEF_es, nut_DEF_es$zcor.ECD.LD)))
)

df.out


### Inspection of extreme values

list(
  # inpsect observations that have extreme values of effect sizes in the biodiversiy dataset
  BEF_B = select(
    identify_outliers(nut_BEF_es, nut_BEF_es$zcor.ECD.B), # on biodiversity
    Case.study, zcor.ECD.B, var.zcor.ECD.B, study_type, system, ECD_subcategory),
  
  
  BEF_LD =  select(
    identify_outliers(nut_BEF_es, nut_BEF_es$zcor.ECD.LD), # on decomposition
    Case.study, zcor.ECD.LD, var.zcor.ECD.LD, study_type, system, ECD_subcategory),
  
  
  # and abundance dataset
  DEF_D = select(identify_outliers(nut_DEF_es, nut_DEF_es$zcor.ECD.B), # effects on abundance
                 Case.study, zcor.ECD.B, var.zcor.ECD.B, study_type, system, ECD_subcategory),
  
  DEF_LD = select(
    identify_outliers(nut_DEF_es, nut_DEF_es$zcor.ECD.LD) , # extreme effects on decomposition
    Case.study, zcor.ECD.LD, var.zcor.ECD.LD, study_type, system, ECD_subcategory)
)



## 3.2.3. Independence of covariates------------------------------------------------------

### Number of studies and observations per level of categorical moderators

## Biodiv dataset
dat <- subset(nut_BEF_es, !duplicated(clusterID_LD))# remove duplicates
countstudies(dat, study_type)
countstudies(dat, system)
countstudies(dat, taxonomic.group)
countstudies(dat, B.metric)


## Abundance  dataset
dat <- subset(nut_DEF_es, !duplicated(clusterID_LD))
countstudies(dat, study_type)
countstudies(dat, system)
countstudies(dat, taxonomic.group)
countstudies(dat, D.metric)

### Number of studies and observations per level of crossed categorical moderators

## Biodiv dataset
dat <- subset(nut_BEF_es, !duplicated(clusterID_LD))# remove duplicates
countstudies(dat, study_type, system)
countstudies(dat, study_type, taxonomic.group)
countstudies(dat, study_type, B.metric)

countstudies(dat, system, taxonomic.group)
countstudies(dat, system, B.metric)

countstudies(dat, taxonomic.group, B.metric)



## Abundance  dataset
dat <- subset(nut_DEF_es, !duplicated(clusterID_LD))
countstudies(dat, study_type, system)
countstudies(dat, study_type, taxonomic.group)
countstudies(dat, study_type, D.metric)

countstudies(dat, system, taxonomic.group)
countstudies(dat, system, D.metric)

countstudies(dat, taxonomic.group, D.metric)



### 3.2.4. Relationships------------------------------------------------------------------

dat.pairplot.bef <- nut_BEF_es %>%  select(
  zcor.ECD.LD, zcor.ECD.B,
  ECD_magnitude, ECD_max,
  study_type, system, taxonomic.group, B.metric)
pairs(
  dat.pairplot.bef,
  upper.panel = panel.smooth,
  lower.panel = panel.cor,
  diag.panel = panel.hist,
  main = "Biodiversity dataset"
)

dat.pairplot.def <- nut_DEF_es %>%  select(
  zcor.ECD.LD, zcor.ECD.D,
  ECD_magnitude, ECD_max,
  study_type, system, taxonomic.group, D.metric)
pairs(
  dat.pairplot.def,
  upper.panel = panel.smooth,
  lower.panel = panel.cor,
  diag.panel = panel.hist,
  main = "Abundance dataset"
)
