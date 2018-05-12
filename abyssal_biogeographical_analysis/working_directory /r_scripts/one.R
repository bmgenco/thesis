### Import CeDAMar and GBIf biological data sets and fix spatial and typographical errors ###

rm(list=ls())

#### Set directories #### 

#setwd("/home/brandon/minerva/abyssal_biogeographical_analysis/") #! For running from absolute location
setwd("..") #! For running script from its own directory.
wd<-(path="./working_directory/") 
data<-(path="./data/")
out<-(path="./output/")
rs<-(path="./r_scripts/")
#setwd(wd) 

#### Install/load packages ####

f.ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("sp")
f.ipak(packages)

#### load CeDAMar & GBIF ####
c<-read.csv("data/DatabaseSyndeepCedamar.csv")
g<-read.delim("data/GBIF_DARWIN.txt")

# spatial exclusion of bad spatial points in g
g<-subset(g, g$hasGeospatialIssues == "false" & g$hasCoordinate =="true")
# remove empty columns
g<- Filter(function(x)!all(is.na(x)), g)

# save as r objects. # @ future @

saveRDS(g, file="./output/gbif.RDS")
saveRDS(c, file="./output/cedamar.RDS")

print("Finished script!")

