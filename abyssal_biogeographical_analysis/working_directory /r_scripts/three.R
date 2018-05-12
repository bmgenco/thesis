#### chi-sqaured gear ####

rm(list=ls())
script_start <- Sys.time()
print("Calculating Ch-sqaured test on species records for all taxonomic classes by gear type: CeDAMar dataset")

#### Set directories: #### 


#setwd("/home/brandon/minerva/abyssal_biogeographical_analysis/") #! For running from absolute location
setwd("..") #! For running script from its own directory. Up one directory

wd<-(path="./working_directory/") 
data<-(path="./data/")
out<-(path="./output/")
rs<-(path="./r_scripts/")
#setwd(wd)

#### Install/load packages: ####

f.ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("psych", "ca", "vcd", "plyr", "ryouready")
f.ipak(packages) 
rm(f.ipak, packages)

#### Load data: ####

c<-readRDS("./output/cedamar.RDS") # % previous %

### data subsetting ####

# animals only
c.k<-c[which(c$WORMS_KINGDOM== c$WORMS_KINGDOM[1]) , ]
# categories of interest
cg<-c.k[c("WORMS_CLASS", "GEAR_CATEGORY")]
# contigency table - drop unused factors
test<-table(cg$WORMS_CLASS, cg$GEAR_CATEGORY, exclude=c("Other/Unknown", "Err:512"))
test<-test[-1,]
# remove classes not represented
new<-test[which(rowSums(test) > 0),]
#convert to datqframe to make further operations easier


## remove white spaces individually  ##
# *improve code*
r<-rownames(new)
r[8]<-"Ascidiacea"
r[13]<-"Caudofoveata"
r[36]<-"Maxillopoda"
r[40]<-"Nematoda_incertae_sedis"
rownames(new)<-r


## sum duplicates/spelling mistakes than remove them individually  ###
# *improve code*


new[7,]<-new[7,]+new[8,]
new[12,]<-new[12,]+new[13,]
new[35,]<-new[35,]+new[36,]+new[37]
new[44,]<-new[43,]+new[44,]


final=new[c(-8,-13,-36, -37, -43),]



rm(c, c.k, cg, new, test)
CeDAMar_class_and_gear<-as.table(final)

# remove gear types with low rperesnation. Note this does not give the warning when computing chi sqaured, yet returns same p-value
#v2<-CeDAMar_class_and_gear[c(-5,-8,-1)]
#b<-chisq.test(v2)

# only major gear types:
#v3<-CeDAMar_class_and_gear[c(2,3,4,7,9)]
#c<-chisq.test(v3)

### chi-squared test ####
a<-chisq.test(CeDAMar_class_and_gear)
print(a)

print("Rejection of independence. Species/lowest taxanomic identifer collected per class is dependent on gear type")

print("Script Finished")
script_end<-Sys.time()
script_end - script_start
