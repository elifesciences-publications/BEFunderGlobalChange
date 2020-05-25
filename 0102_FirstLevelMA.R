####################################################################################################################
##
## First level meta-analysis : Overall effects of stressors and nutrients on decomposers and decomposition
##
###################################################################################################################

## This script runs the first-level meta-analysis that compares the overall mean effect of two types of drivers of global change 
# on the three response variables: decomposer diversity, decomposer abundance or biomass and litter decomposition rates

# First load datasets of effect sizes,
# Then performs the meta-analysis
# And returns and save a figure showing the grand mean effect sizes of stressors and nutrients on the three response variables
# The figure summarizes the results from the first level meta-analysis (Figure 3 of the manuscript)


## 1. LOAD FUNCTIONS AND PACKAGES--------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(ggplot2)
library(metafor)
library(broom)
library(knitr)
library(tidyverse)
library(devtools)
library(piecewiseSEM)
library(nlme)
library(car)
library(patchwork) #devtools::install_github("thomasp85/patchwork")

## notin operator
'%notin%' <- function(x,y)!('%in%'(x,y))

## count levels of a non factor variable
countlevels <- function(x){
  return(length(levels(as.factor(x))))
}


# count no. studies and observations per level of categorical moderator
countstudies <- function(dat, ...){
  return(dat %>% 
           group_by_(.dots = lazyeval::lazy_dots(...)) %>% 
           summarize(no.stu = countlevels(Case.study), no.obs = n()))
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

## Function to check the residuals visually
function_resid <- function(res){
  par(mfrow = c(1,3))
  qqnorm(residuals(res,type="pearson"),main="QQ plot: residuals")
  qqline(residuals(res,type="pearson"),col="red")
  hist(residuals(res, type = "pearson"))
  plot(residuals(res, type = "rstandard") ~ predict(res)$pred, main = "Residuals vs Fitted"); abline(0,0)
}


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

## 2. LOAD DATA---------------------------------------------------------------------------

## Load the datasets abundance-decomposition and biodiv-decomposition means and sd for gradients and control-treatment studies

effect_size_B1 = read.csv("data/effect_size_B1.csv", h = TRUE, sep = ",")
effect_size_D1 = read.csv("data/effect_size_D1.csv", h = TRUE, sep = ",")
effect_size_LD1 = read.csv("data/effect_size_LD1.csv", h = TRUE, sep = ",")

## Prep data -----------------------------------------------------------------------------

# ECD type
summary(effect_size_B1$ECD_subcategory)
effect_size_B1 = dplyr::mutate(effect_size_B1, ECD.type = factor(ifelse(grepl("eutroph|nitro|nutri|wastewater", ECD_subcategory), "Resources", "Stressors")))

summary(effect_size_D1$ECD_subcategory)
effect_size_D1 = dplyr::mutate(effect_size_D1, ECD.type = factor(ifelse(grepl("eutroph|nitro|nutri|wastewater", ECD_subcategory), "Resources", "Stressors")))

summary(effect_size_LD1$ECD_subcategory)
effect_size_LD1 = dplyr::mutate(effect_size_LD1, ECD.type = factor(ifelse(grepl("eutroph|nitro|nutri|wastewater", ECD_subcategory), "Resources", "Stressors")))


## 3. ANALYSIS------------------------------------------------------------------

## 3.1. BIODIVERSITY RESPONSE TO NUTRIENTS VS. STRESSORS---------------------------------

# size of the dataset
# no. observation in total
nrow(effect_size_B1)

# no. observation per driver type
with(effect_size_B1, table(ECD.type))

# no. studies per ECD type
countstudies(effect_size_B1, ECD.type)

# effect size distribution ----------------------------------------------------------------------------------------------
hist(effect_size_B1$zcor.ECD.B)
shapiro.test(effect_size_B1$zcor.ECD.B)

summary(effect_size_B1$zcor.ECD.B)

# outliers ----------------------------------------------------------------------------------------------
dotchart(sort(effect_size_B1$zcor.ECD.B))

# homogeneity of variances
boxplot(effect_size_B1$zcor.ECD.B ~ effect_size_B1$ECD.type)

# publication bias (without covariates)-------------------------------------------------------------------
funnel(x = effect_size_B1$zcor.ECD.B, vi = effect_size_B1$var.zcor.ECD.B)

# random effect at the observation level
effect_size_B1$ID <- 1 : nrow(effect_size_B1)

## meta-analysis diversity overall------------------------------------------------------------------
res.B <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                random =  ~ 1 | Case.study/ID, 
                data = effect_size_B1)
res.B # mean effect size negative : **

## meta-analysis diversity as function of ECD type ------------------------------------------------
res.B.mod_ecd <- rma.mv(zcor.ECD.B, var.zcor.ECD.B, 
                        random =  ~ 1 | Case.study/ID, 
                        data = effect_size_B1,
                        mods =~ ECD.type-1) #removing the intercept to have the mean effects per ECD type
res.B.mod_ecd 

# Q statistic : test difference between mean effect sizes for Res vs Stress.
anova(res.B.mod_ecd, btt = c(1:2))

# test if mean effect sizes differ from 0 : Stress. only 
res.B.mod_ecd$pval

# Heterogeneity and I2 value: see 
W <- diag(1/effect_size_B1$var.zcor.ECD.B)
X <- model.matrix(res.B.mod_ecd)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * res.B.mod_ecd$sigma2 / (sum(res.B.mod_ecd$sigma2) + (res.B.mod_ecd$k-res.B.mod_ecd$p)/sum(diag(P)))
# About 75 % of variability is estimated to be due to between cluster heterogeneity (i.e. Case study)
# 16 % of variability remaining explained by within cluster heterogeneity (i.e. Observations within case study)
# remaining 10 % are the sampling variance


## FIG Residuals ---------------------------------------------------------------------------------------
function_resid(res.B.mod_ecd)

## FIG Funnel plots ECD effect on B -----------------------------------------------------------------
par(mfrow = c(1,1))
funnel(res.B.mod_ecd, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)# symetric round mean effect size, but no pattern of increasing variance at low sample size


## 3.2. ABUNDANCE RESPONSE TO NUTRIENTS VS. STRESSORS---------------------------------

# size of the dataset
# no. observation in total
nrow(effect_size_D1)

# no. observation per ECD type
with(effect_size_D1, table(ECD.type))

# no. studies per ECD type
countstudies(effect_size_D1, ECD.type)

# effect size distribution ----------------------------------------------------------------------------------------------
hist(effect_size_D1$zcor.ECD.D)
shapiro.test(effect_size_D1$zcor.ECD.D)

summary(effect_size_D1$zcor.ECD.D)


## meta-analysis abundance overall------------------------------------------------------------------
effect_size_D1$ID <- 1:nrow(effect_size_D1)

res.D <- rma.mv(zcor.ECD.D, var.zcor.ECD.D, 
                random =  ~ 1 | Case.study/ID, 
                data = effect_size_D1)
res.D # mean effect size ci overlap zero

## meta-analysis abundance as function of ECD type ------------------------------------------------
res.D.mod_ecd <- rma.mv(zcor.ECD.D, var.zcor.ECD.D, 
                        random =  ~ 1 | Case.study/ID,
                        data = effect_size_D1,
                        mods =~ ECD.type-1) #removing the intercept to have the mean effects per ECD type

res.D.mod_ecd 

# test difference between type of ECD : *
anova(res.D.mod_ecd, btt = c(1:2))

# do mean effect sizes differ from 0 : Stress. only **
res.D.mod_ecd$pval  


## FIG Residuals ---------------------------------------------------------------------------------------
# model with ECD type
function_resid(res.D.mod_ecd)

## FIG Funnel plots ECD effect on B -----------------------------------------------------------------
# no sign of asymetry
par(mfrow = c(1,2))
funnel(res.D.mod_ecd, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)# symetric round mean effect size, but no pattern of increasing variance at low sample size

## 3.3. DECOMPOSITION RESPONSE TO NUTRIENTS VS. STRESSORS---------------------------------

# size of the dataset
# no. observation in total
nrow(effect_size_LD1)

# no. observations having variance estimate
nrow(dplyr::filter(effect_size_LD1, !is.na(var.zcor.ECD.LD)))

# no. observation per ECD type
with(effect_size_LD1, table(ECD.type))

# no. studies per ECD type
countstudies(effect_size_LD1, ECD.type) # studies
effect_size_LD1 %>% group_by(ECD.type) %>% summarise(no.studies = countlevels(studyID)) # papers

# effect size distribution ----------------------------------------------------------------------------------------------
hist(effect_size_LD1$zcor.ECD.LD)
shapiro.test(effect_size_LD1$zcor.ECD.LD)

summary(effect_size_LD1$zcor.ECD.LD)


## meta-analysis diversity overall------------------------------------------------------------------
effect_size_LD1$ID <- 1:nrow(effect_size_LD1)

res.LD <- rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                 random =  ~ 1 | Case.study/ID, 
                 
                 data = effect_size_LD1)
res.LD # n.s

## meta-analysis diversity as function of ECD type ------------------------------------------------
res.LD.mod_ecd <- rma.mv(zcor.ECD.LD, var.zcor.ECD.LD, 
                         random =  ~ 1 | Case.study/ID, 
                         data = effect_size_LD1,
                         mods =~ ECD.type-1) #removing the intercept to have the mean effects per ECD type

res.LD.mod_ecd 

# test difference between type of ECD : ***
anova(res.LD.mod_ecd, btt = c(1:2))

# do mean effect sizes differ from 0 : Negative effect of stress. *** (p = 0.0004) vs. positive of resources (p = 0.06), but marginally signif
res.LD.mod_ecd$pval


## FIG Residuals ---------------------------------------------------------------------------------------
# model with ECD type : residuals don't look good at the tails
function_resid(res.LD.mod_ecd)

## FIG Funnel plots ECD effect on B -----------------------------------------------------------------
# no sign of asymetry
par(mfrow = c(1,1))
funnel(res.LD.mod_ecd, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)# symetric round mean effect size, but no pattern of increasing variance at low sample size



## 4. FIGURES----------------------------------------------------------------------------

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Mean effect sizes from MA models  -------------------------------------------------------
# fgb dataset for ggplot
y<-summary(res.B.mod_ecd)$b
ci_l<-summary(res.B.mod_ecd)$ci.lb
ci_h<-summary(res.B.mod_ecd)$ci.ub
fgb<-data.frame(cbind(y,ci_l,ci_h))
colnames(fgb)[1]<-"y"
colnames(fgb)[2]<-"ci_l"
colnames(fgb)[3]<-"ci_h"
fgb$ECD<-c("Nutrients","Stressors")
fgb$ECD<-as.factor(fgb$ECD)
fgb$no.obs <- c(with(effect_size_B1, table(ECD.type))[1], with(effect_size_B1, table(ECD.type))[2])
no.stu.b <- effect_size_B1 %>% group_by(ECD.type) %>% summarise(no.studies = countlevels(Case.study))
fgb$no.stu <- c(no.stu.b$no.studies[1], no.stu.b$no.studies[2])


# fgd dataset for ggplot
y<-summary(res.D.mod_ecd)$b
ci_l<-summary(res.D.mod_ecd)$ci.lb
ci_h<-summary(res.D.mod_ecd)$ci.ub
fgd<-data.frame(cbind(y,ci_l,ci_h))
colnames(fgd)[1]<-"y"
colnames(fgd)[2]<-"ci_l"
colnames(fgd)[3]<-"ci_h"
fgd$ECD<-c("Nutrients","Stressors")
fgd$ECD<-as.factor(fgd$ECD)
fgd$no.obs <- c(with(effect_size_D1, table(ECD.type))[1], with(effect_size_D1, table(ECD.type))[2])
no.stu.d <- effect_size_D1 %>% group_by(ECD.type) %>% summarise(no.studies = countlevels(Case.study))
fgd$no.stu <- c(no.stu.d$no.studies[1], no.stu.d$no.studies[2])

# fgld dataset for ggplot
y<-summary(res.LD.mod_ecd)$b
ci_l<-summary(res.LD.mod_ecd)$ci.lb
ci_h<-summary(res.LD.mod_ecd)$ci.ub
fgld<-data.frame(cbind(y,ci_l,ci_h))
colnames(fgld)[1]<-"y"
colnames(fgld)[2]<-"ci_l"
colnames(fgld)[3]<-"ci_h"
fgld$ECD<-c("Nutrients","Stressors")
fgld$ECD<-as.factor(fgld$ECD)
fgld$no.obs <- c(with(effect_size_LD1, table(ECD.type))[1], with(effect_size_LD1, table(ECD.type))[2])
no.stu.ld <- effect_size_LD1 %>% group_by(ECD.type) %>% summarise(no.studies = countlevels(Case.study))
fgld$no.stu <- c(no.stu.ld$no.studies[1], no.stu.ld$no.studies[2])


sizetext <- 13
sizeyaxis <- 11
sizepoint <- 3
widtherrorbar <- 0.1
sizeerrorbar <- 0.4


# Plot of effect sizes of B as function of type of ECD (nut vs. chem. stress.)
plot_B <-
  ggplot(fgb, aes(x=ECD, y=y)) +
  geom_errorbar(width=widtherrorbar, aes(ymin=ci_l, ymax=ci_h), col = c("#D55E00", "#0072B2"), size = sizeerrorbar) +
  geom_point(shape=21, size = sizepoint, col=c("#D55E00", "#0072B2"), fill = c("#D55E00", "#0072B2")) +
  ylim(-1.5,1.5)+
  ylab("Effect size")+
  xlab(" ")+
  coord_flip()+
  geom_hline(yintercept = 0)+
  
  scale_x_discrete(breaks=c("Nutrients", "Stressors"),
                   labels=c(paste("Nutrients (",fgb$no.stu[1], "; ", fgb$no.obs[1], ")", sep = ""), 
                            paste("Stressors (", fgb$no.stu[2], "; ", fgb$no.obs[2],")", sep = ""))) +
  annotate(geom="text", x=2.1, y=fgb$y[2], vjust = 0.5, label=paste("***"), colour="#0072B2", size=0.3*sizetext)+
  theme_bw() +
  ggtitle("Biodiversity") +
  theme(axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size= sizeyaxis, face = "bold"),
        plot.title = element_text(size = sizetext))



# Plot of effect sizes of B as function of type of ECD (nut vs. chem. stress.)
plot_D <- 
  ggplot(fgd, aes(x=ECD, y=y)) +
  geom_errorbar(width=widtherrorbar, aes(ymin=ci_l, ymax=ci_h), col = c("#D55E00", "#0072B2"), size = sizeerrorbar) +
  geom_point(shape=21, size = sizepoint, col=c("#D55E00", "#0072B2"), fill = c("#D55E00", "#0072B2")) +
  ylim(-1.6,1.6)+
  ylab("Effect size")+
  xlab(" ")+
  coord_flip()+
  geom_hline(yintercept = 0)+
  scale_x_discrete(breaks=c("Nutrients", "Stressors"),
                   labels=c(paste("(",fgd$no.stu[1], "; ", fgd$no.obs[1], ")", sep = ""), 
                            paste("(", fgd$no.stu[2], "; ", fgd$no.obs[2],")", sep = ""))) +
  annotate(geom="text", x=2.1, y=fgd$y[2], vjust = 0.5, label=paste("*"), colour="#0072B2", size=0.3*sizetext)+
  theme_bw() +
  ggtitle("Abundance") +
  theme(axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size=sizeyaxis, face = "bold"),
        plot.title = element_text(size = sizetext))



# Plot of effect sizes of B as function of type of ECD (nut vs. chem. stress.)
plot_LD <- 
  ggplot(fgld, aes(x=ECD, y=y, group=1)) +
  geom_errorbar(width=widtherrorbar, aes(ymin=ci_l, ymax=ci_h), col = c("#D55E00", "#0072B2"), size = sizeerrorbar) +
  geom_point(shape=21, size = sizepoint, col=c("#D55E00", "#0072B2"), fill = c("#D55E00", "#0072B2")) +
  ylim(-1.6,1.6)+
  ylab("Effect size")+
  xlab(" ")+
  coord_flip()+
  geom_hline(yintercept = 0)+
  scale_x_discrete(breaks=c("Nutrients", "Stressors"),
                   labels=c(paste("(",fgld$no.stu[1], "; ", fgld$no.obs[1], ")", sep = ""), 
                            paste("(", fgld$no.stu[2], "; ", fgld$no.obs[2],")", sep = ""))) +
  annotate(geom="text", x=2.1, y=fgld$y[2], vjust = 0.5, label=paste("***"), colour="#0072B2", size=0.3*sizetext)+
  # annotate(geom="text", x=1.12, y=fgld$ci_h[1] + 0.2, label=paste("."), colour="#D55E00", size=0.7*sizetext)+
  #geom_point(aes(2,fgld$ci_l[2]),shape=8) +
  theme_bw() +
  ggtitle("Litter decomposition") +
  theme(axis.text.y=element_text(face = "bold", size = sizeyaxis), 
        axis.text.x=element_text(size = sizeyaxis),
        axis.title.x = element_text(size=sizeyaxis, face = "bold"),
        plot.title = element_text(size = sizetext))


# save a png with high res
ppi <- 300 # resolution
w <- 21 # width in cm

png("figs/Fig3_GrandMeans.png", 
    width=w, 
    height=0.3*w, 
    units = "cm",
    res=ppi)

plot_B + plot_D + plot_LD


dev.off()

## 5. SUPPLEMENTARY TABLE----------------------------------------------------------------

## Creating Table S1 - Wald type test for stressors vs. nutrients overall effects

## ## Table with stats : omnibus test of moderators=================================================
table_Qmstat <- data.frame(Response = c("Diversity", "Abundance", "Litter decomposition"),
                           rbind(as.numeric(anova(res.B.mod_ecd)[c(1,5,4, 2)]), 
                                 as.numeric(anova(res.D.mod_ecd)[c(1,5,4,2)]), 
                                 as.numeric(anova(res.LD.mod_ecd)[c(1,5,4,2)])))
colnames(table_Qmstat) <- c("Response", "QM", "df","n", "pVal")

table_Qmstat

write.csv(table_Qmstat, "tables/TableS1_QMstat_polvsnut.csv")

