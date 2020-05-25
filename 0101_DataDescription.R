####################################################################################################################
##
## ## Data description
##
###################################################################################################################

## This script provide a description of the data: count number of studies and observations
# for different categorical factor levels and makes Figure 1 - map and barplots of no. studies
rm(list=ls())

## 0. Load packages-----------------------------
library(dplyr)
library(ggplot2)
library(patchwork) #devtools::install_github("thomasp85/patchwork")


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


## 1. Load datasets ----------------------------

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


## 2. Prep data -----------------------------------------------------------------------------

# Type of environmental change driver (ECD)
# Create a categorical variable for type of ECD: Resources stands for nutrient enrichment studies and Stressors for chemical stressors
summary(effect_size_BandEF1$ECD_subcategory)
effect_size_BandEF1 = dplyr::mutate(effect_size_BandEF1, ECD.type = factor(ifelse(
  grepl("eutroph|nitro|nutri|wastewater", ECD_subcategory), "Resources", "Stressors")))

summary(effect_size_DandEF1$ECD_subcategory)
effect_size_DandEF1 = dplyr::mutate(effect_size_DandEF1, ECD.type = factor(ifelse(
  grepl("eutroph|nitro|nutri|wastewater", ECD_subcategory), "Resources", "Stressors")))


# Biodiversity metric
# there are not enough studies to test shannon and evenness separately in the meta-analysis
# create a variable separating biodiversity observations in terms of taxa richness (richness) and of diversity indices (shannon and evenness)
effect_size_BandEF1 <- effect_size_BandEF1 %>% mutate(
  B.metric = factor(
    ifelse(B.metric == "taxa richness", "richness", "div_indices")))


## Subset variables of interest
whichvar <- c("studyID", "Case.study", "ECD.type" ,
              "zcor.ECD.LD", "var.zcor.ECD.LD",
              "zcor.ECD.B", "var.zcor.ECD.B",
              "ECD_magnitude", "ECD_max",
              "country",
              "study_type", "system", "taxonomic.group", 
              "community.nature", "consumer.taxa", "ES_origin", "meshsize..mm.", 
              "litter.type", "ECD_subcategory", "ECD_identity")

## Diversity dataset
bef <- subset(effect_size_BandEF1, select = whichvar)
bef$B.metric <- factor(rep("Biodiversity", nrow(bef)))
bef$Bindices <- effect_size_BandEF1$B.metric

## Abundance dataset
def <- subset(effect_size_DandEF1, select = whichvar)
def$B.metric <- factor(rep("Abundance", nrow(def)))
def$Bindices <- effect_size_DandEF1$D.metric

## Overall
overalldat <- bind_rows_keep_factors(bef, def)


## 3. Count studies (main text of results)-------------------------------------------------
## Count total no. studies and observations
countstudies(overalldat) # no. case studies and obs
length(levels(overalldat$studyID)) # no. publications


## geographical coverage (continent-wide)
library(countrycode)

overalldat$continent <- countrycode(sourcevar = overalldat$country,
                                    origin = "country.name",
                                    destination = "continent")
countstudies(overalldat, continent)


## system and study types
countstudies(overalldat, system)
countstudies(overalldat, study_type)


## taxonomic coverage
countstudies(overalldat, taxonomic.group)


## Diversity vs. Abundance data
countstudies(overalldat, B.metric)
countstudies(overalldat, Bindices)


## Env change drivers
countstudies(overalldat, ECD.type, ECD_subcategory)


## Identity of drivers
summary(overalldat$ECD_identity)


## litterbag meshsizes
overalldat$meshsize <- factor(ifelse(overalldat$meshsize..mm. %in% c(".2", "0.25", "0.3", "0.5", "fine"), "fine", 
                                     ifelse(is.na(overalldat$meshsize..mm.), "NA", "coarse")))

countstudies(overalldat, meshsize)

## decomposition, litter types
summary(overalldat$litter.type)


## no. studies reporting evenness
effect_size_BandEF1 = read.csv("data/effect_size_BandEF1.csv", h = TRUE, sep = ",")
countstudies(effect_size_BandEF1, B.metric)




## 4. FIGURE DATA DESCRIPTION-------------------------------------------------
# set sizes of text in plots
sizetext <- 12
sizelegend <- 11

## 4.1. Figure map================
library(ggmap)
library(rnaturalearth)
library(viridis)

stud <- countstudies(overalldat, country)

# country codes 
stud <- stud %>% 
  mutate(iso =  countrycode(stud$country, origin = "country.name", destination = "iso3c")) %>% 
  dplyr::group_by(iso) %>% 
  summarize(no.obs = sum(no.obs), no.stu = sum(no.stu))

wm <- ne_countries(scale = 110, returnclass = "sf") %>%
  left_join(stud, by = c("adm0_a3" = "iso"))

wm <- sf::st_transform(wm,"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

colorpalettemap <- c("#874b03","#bd7719", "#edc74a", "#60dbc8",  "#018c80") # brighter than colorbrewer

map <- ggplot(wm)+
  geom_sf(aes(fill = log(no.obs)))+
  scale_fill_gradientn(colours = colorpalettemap, na.value = "black",
                       name = "log(no.obs)",
                       breaks = c(0,log(10),log(25),log(100)),
                       labels = c(1, 10,25,100))+
  theme_bw()+
  theme(legend.position = c(0,0),
        legend.justification = c(-0.2, -0.1),
        legend.title = element_text(size=sizelegend),
        legend.text = element_text(size=sizelegend))


## 4.2. Figure barplot system types ~ driver*taxonomic group============

datplot <- countstudies(overalldat, ECD.type, system, taxonomic.group)

# names and order of factor levels
levels(datplot$ECD.type) <- c("Nutrient Enrichment", "Chemical Stressors")
datplot$ECD.type <- factor(datplot$ECD.type, levels = c("Chemical Stressors", "Nutrient Enrichment"))
levels(datplot$taxonomic.group) <-  c("Animal", "Microbial")
levels(datplot$system) <-  c("Aquatic", "Terrestrial")

# colors
pale_aquaterr <- c("white", "black")

Fig_systemtaxo <- ggplot(data=datplot, aes(x=taxonomic.group, y=no.stu, fill=system)) +
  geom_bar(stat = "identity", col = "black")+
  scale_fill_manual(values=pale_aquaterr)+
  
  ylab("no. studies")+
  xlab(" ")+
  facet_grid(~ECD.type)+
  
  # theme stuff
  theme_bw() +
  theme(axis.text.y=element_text(face = "bold", size = sizetext), 
        axis.text.x=element_text(size = sizetext),
        axis.title.y = element_text(size=sizetext, face = "bold"),
        
        #legend
        legend.title = element_blank() ,
        legend.text = element_text(size=sizelegend), 
        legend.position = "bottom",
        legend.spacing.x = unit(0.25, 'cm'),
        
        #for the facets
        strip.text.x = element_text(size=sizelegend, face = "bold"),
        strip.background = element_rect(colour="black", fill="white"))



## 4.3 Figure barplot : B and D metrics================

datplot2 <- countstudies(overalldat,  Bindices, B.metric)
levels(datplot2$Bindices) <-  c("Abundance", "Biomass", "Diversity indices", "Taxa richness")


# colors
pale_divind <- c("#cccccc",  "#f7f7f7", "#525252", "#969696")#, "#f6e8c3") # GRAYS

Fig_bindices <- ggplot(data=datplot2, aes(x=B.metric, y=no.stu, fill = Bindices)) +
  geom_bar(stat = "identity", position = "dodge", col = "black")+
  
  scale_fill_manual(values=pale_divind)+
  
  ylab("no. studies")+
  xlab(" ")+
  
  # theme stuff
  theme_bw() +
  theme(axis.text.y=element_text(face = "bold", size = sizetext), 
        axis.text.x=element_text(size = sizetext),
        axis.title.y = element_text(size=sizetext, face = "bold"),
        
        #legend
        legend.title = element_blank() ,
        legend.text = element_text(size=10), 
        legend.position = "bottom",
        legend.spacing.x = unit(0.15, 'cm'),
        
        
        #for the facets
        strip.text.x = element_text(size=sizelegend, face = "bold"),
        strip.background = element_rect(colour="black", fill="white"))


## 5. Final figure------------------------------------------------------------------------
ppi <- 600# final: 600 # resolution
w <- 20 # width in cm

png("figs/Fig2_DataDescription.png",
    width=w,
    height=w,
    units = "cm",
    res=ppi)

map + labs(tag = "A") + {
  Fig_bindices +labs(tag = "B") + Fig_systemtaxo +labs(tag = "C")} + plot_layout(ncol = 1, heights = c(2,1))

dev.off()
