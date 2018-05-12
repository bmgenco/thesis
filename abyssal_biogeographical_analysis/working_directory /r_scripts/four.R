### Main MDE Anlaysis CeDAMar ###

rm(list=ls())
print("run script four and four_b, or run four twice following instructions in script header before second run")

# Instruction for running script second time:
# Search for " ^ alter ^ ". Uncomment "HUNDRED ONLY" lines. Comment out "GEOMETRIC MEAN".
# Select full function use CTRL + SHIFT + C.
# Five locations: @ Rarefaction step. @ Rename variables 1 & 2 steps. @ MDE step. @ Save step.
# Need to run both time for full analysis and graphing.


#### Set directories: #### 

script_start <- Sys.time()
print("calculating MDE CeDAMar")
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
packages <- c("tidyr", "permute", "lattice", "vegan", "rich", "stringr","sp") 
f.ipak(packages) 
rm(packages)

#### Load Data ####

c<-readRDS("./output/cedamar.RDS") # % Previous %
provinces<-readRDS("./output/provinces.RDS") # % Previous %
final_centroids<-readRDS("./output/final_centroids.RDS") # % Previous % 

#### Clean data ####

c[c==""]  <- NA
c<-c[which(c$WORMS_KINGDOM== c$WORMS_KINGDOM[1]) , ]
c<-subset(c, !is.na(c$SPECIES))
c<-subset(c, !is.na(c$COUNT))
c<-subset(c, !is.na(c$WORMS_CLASS))
c<-c[ , c("COUNT", "SPECIES", "WORMS_GENUS", "WORMS_CLASS", "LON_DEG", "LAT_DEG")]

tmp<-c$WORMS_CLASS
t1<-str_trim(tmp)
c$WORMS_CLASS<-t1
c$WORMS_CLASS<-as.character(c$WORMS_CLASS)
c$WORMS_CLASS<- gsub("Maxilopoda","Maxillopoda", c$WORMS_CLASS)
c$WORMS_CLASS<- gsub("Phascolosomati","Phascolosomatidea", c$WORMS_CLASS)
c$WORMS_CLASS<- gsub("Phascolosomatideadea","Phascolosomatidea", c$WORMS_CLASS)
c$WORMS_CLASS<-as.factor(c$WORMS_CLASS)
c<-droplevels(c)

#### Subset CeDAMar by province ####

coords<-cbind(c$LON_DEG, c$LAT_DEG)
cd = SpatialPointsDataFrame(coords,c)
proj4string(cd)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
fix1<-("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m +lon_wrap=80")
fix2<-("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m +lon_wrap=180")

cd.f1<-spTransform(cd, CRS(fix1))
cd.f2<-spTransform(cd, CRS(fix2))

cd1<-cd.f1[provinces[[1]],]

cd2<-cd[provinces[[2]],]
cd3<-cd[provinces[[3]],] # no samples
cd4<-cd[provinces[[4]],]
cd5<-cd[provinces[[5]],]
cd6<-cd[provinces[[6]],]

cd7<-cd.f2[provinces[[7]],]

cd8<-cd[provinces[[8]],]
cd9<-cd[provinces[[9]],] # no samples

cd10<-cd.f2[provinces[[10]],]
cd11<-cd.f2[provinces[[11]],]
cd12<-cd.f2[provinces[[12]],]
cd13<-cd.f2[provinces[[13]],] # no samples

cd14<-cd[provinces[[14]],] # no samples

rm(cd13, cd14, cd9, cd3) # *improve code*
# change to condtional if empty dataframe
rm(c, cd, coords, cd.f1, cd.f2, provinces, fix1, fix2)

#### Create a list object of datasets ----

all<-list(cd1, cd2, cd4, cd5, cd6, cd7, cd8, cd10, cd11, cd12)
rm(cd1, cd2, cd4, cd5, cd6, cd7, cd8, cd10, cd11, cd12)
temp<-lapply(all, function(x) {x<- transform(x, site= as.numeric(interaction(LON_DEG, LAT_DEG, drop=TRUE)))}) # assign unique sites
names(temp)<-c("p1", "p2", "p4", "p5", "p6", "p7", "p8", "p10", "p11", "p12")

f.1<-function(x){x<-x[,c(1:6,10)]}
a<-lapply(temp, f.1)
rm(all, temp, f.1)

#### Data rearrangement ####

# ninety nine problems
f.1<-function(x){x[x$COUNT < 9999,]} 
b<-lapply(a, f.1)
rm(a, f.1)

# class issues
f.2<-function(x){split(x, x$WORMS_CLASS, drop=T)} 
a<-lapply(b, f.2)
rm(b, f.2)

# list made up of a list (per province) of class dataframes
f.1<-function(x){lapply(x, function(x)split(x, x$site, drop =T))} # site splitter 
b<-lapply(a, f.1)
rm(a, f.1)

# fix duplicate species columns
f.2<-function(x){lapply(x, function(x){lapply(x, function(x)cbind((aggregate(COUNT~SPECIES, data=x, FUN =sum, drop=T)),(x[1,4:7])))})} 
options(warn=-1) # function returns warnings b/c of removing duplicates
a<-lapply(b, f.2)
options(warn = 0)
rm(b, f.2)

# rearrange rarefication preperation #
f.1<-function(x){lapply(x, function(x){lapply(x, function(x)cbind((spread(x[,1:2], SPECIES, COUNT, fill = F, convert = F, drop = T)),(x[1,3:6])))})}
b<-lapply(a, f.1)
rm(a, f.1)

detach("package:tidyr", unload = T)
#detach("package:dplyr", unload = T)
f.ipak("plyr")

# rearranges column order $WORM_CLASS, SITE, lat long, not placed at end....
f.2<-function(x){lapply(x, rbind.fill)} 
a<-lapply(b, f.2)
rm(b, f.2)

detach("package:plyr", unload=T)
f.ipak("dplyr")

f.1<-function(x){x %>% select(-WORMS_CLASS, -LON_DEG, -LAT_DEG, -site, everything())} 
f.2<-function(x){lapply(x, f.1)}
b<-lapply(a, f.2)
rm(a, f.1, f.2)

f.1<-function(x){
  t<-x[,1:((ncol(x))-4)]
  s<-x[, ((ncol(x))-3):(ncol(x))]
  t[is.na(t)]<-0
  cbind(t,s)
}
f.2<-function(x){lapply(x, f.1)}
a<-lapply(b, f.2)
rm(b, f.1, f.2)

# remove all instances only represented by one species:
f.1<-function(x)(lapply(x, function(x) if(length(x) > 5) {x}))
b<-lapply(a, f.1)
rm(a, f.1)

f.2<-function(x) x[!sapply(x, is.null)] 
rare_prep<-lapply(b, f.2)
rm(b, f.2)

#### Rarefaction ####

detach("package:dplyr", unload = T)
library(plyr)
options(warn=-1)


f.H_rare<-function(x){
  # function f.RARE()
  # requires package "plyr"
  # set waring options to: options(warn=-1) ....
  # b/c rareficaton sub sampling sizes = medain plus one std, and sum of all individials => larger than minium per site.
  # remove bootstraping b/c of returned errors
  
  splta<-x[,1:((ncol(x))-4)]
  splta<-round(splta) #counts must be integers 
  spltb<-x[, ((ncol(x))-3):(ncol(x))]
  f1<-cbind(splta, spltb)
  
  rm(splta, spltb)
  
  fa<-f1[,1:((ncol(f1))-4)]
  f1<-f1[which(rowSums(fa) >= 1), ]
  
  a<-f1[,1:((ncol(f1))-4)]
  b<-f1[, ((ncol(f1))-3):(ncol(f1))]
  
  #choose rarfication 
  inds<-vector()
  for(i in 1:nrow(a)) inds[i]<-sum(a[i,])
  #samp_size<-round_any((median(inds)+1),1)
  inds<-sort(inds)
  #samp_size<-inds[1]
  hun<-100
  #max<-sum(inds)
  
  # calculations
  #c<-rarefy(a, samp_size)
  #d<-rich(a, nrandom=1000)
  #e<-rarefy(a, max)
  f<-rarefy(a, hun)
  height<-nrow(a)
  #mean<-rep(d$mr, height)
  #boot<-rep(d$bootCR$cr.boot, height)
  #number_rare<-rep(samp_size, height) 
  #z<-cbind(mean,boot,e,c,number_rare,b)
  #names(z)<-c("mean_richness", "cumulative_boot", "max_estimated", "estimated_richness", "sub_sample_size", "WORMS_CLASS", "LON_DEG", "LAT_DEG", "site")                       
  sp<-rep(length(a), height)
  count<-rep(sum(inds), height)
  #geometric mean
  #m<-cbind(c,e,f)
  #gm<-apply(m,1, function(x) exp(mean(log(x))))
  z<-cbind(f, count,sp,b)
  names(z)<-c("100_richness",  "total_individuals",   "total_species", "WORMS_CLASS", "LON_DEG", "LAT_DEG", "site") 
  x<-z
  
} #HUNDRED ONLY
f.rare<-function(x){
  # function f.RARE()
  # requires package "plyr"
  # set waring options to: options(warn=-1) ....
  # b/c rareficaton sub sampling sizes = medain plus one std, and sum of all individials => larger than minium per site.
  # remove bootstraping b/c of returned errors
  
  splta<-x[,1:((ncol(x))-4)]
  splta<-round(splta) #counts must be integers 
  spltb<-x[, ((ncol(x))-3):(ncol(x))]
  f1<-cbind(splta, spltb)
  
  rm(splta, spltb)
  
  fa<-f1[,1:((ncol(f1))-4)]
  f1<-f1[which(rowSums(fa) >= 1), ]
  
  a<-f1[,1:((ncol(f1))-4)]
  b<-f1[, ((ncol(f1))-3):(ncol(f1))]
  
  #choose rarfication 
  inds<-vector()
  for(i in 1:nrow(a)) inds[i]<-sum(a[i,])
  samp_size<-round_any((median(inds)+1),1)
  #inds<-sort(inds)
  #samp_size<-inds[1]
  hun<-100
  max<-sum(inds)
  
  # calculations
  c<-rarefy(a, samp_size)
  #d<-rich(a, nrandom=1000)
  e<-rarefy(a, max)
  f<-rarefy(a, hun)
  height<-nrow(a)
  #mean<-rep(d$mr, height)
  #boot<-rep(d$bootCR$cr.boot, height)
  number_rare<-rep(samp_size, height) 
  #z<-cbind(mean,boot,e,c,number_rare,b)
  #names(z)<-c("richness", "cumulative_boot", "max_estimated", "estimated_richness", "sub_sample_size", "WORMS_CLASS", "LON_DEG", "LAT_DEG", "site")                       
  sp<-rep(length(a), height)
  count<-rep(sum(inds), height)
  #geometric mean
  m<-cbind(c,e,f)
  gm<-apply(m,1, function(x) exp(mean(log(x))))
  z<-cbind(gm,e,f,c,count,number_rare,sp,b)
  names(z)<-c("richness", "max_inds_richness", "hundred_richness","sub_sample_richness", "total_individuals",  "sub_sample_size_median", "total_species", "WORMS_CLASS", "LON_DEG", "LAT_DEG", "site") 
  x<-z
   } #GEOMETRIC MEAN

# ^ alter ^ 
# Comment out/in appropriate f.1! 

# f.1<-function(x){lapply(x, f.H_rare)} #HUNDRED ONLY
# OR
f.1<-function(x){lapply(x, f.rare)} #GEOMETRIC MEAN

a<-lapply(rare_prep, f.1)
rarefied<-lapply(a, function(x)(x<-do.call("rbind", x)))
names(rarefied)<-(names(a))
rm(a, rare_prep, f.1, f.rare, f.H_rare)

#### Add centroids ####

#distance prep
f.1<-function(x){
  coords<-cbind(x$LON_DEG, x$LAT_DEG)
  z<-SpatialPointsDataFrame(coords,x)
  proj4string(z)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  x<-z
}
b<-lapply(rarefied, f.1)
names(b)<-names(rarefied)

rm(rarefied, f.1)
b$p1@data$distance<-spDistsN1(b$p1, final_centroids[1,], longlat=T)
b$p2@data$distance<-spDistsN1(b$p2, final_centroids[2,], longlat=T)
b$p4@data$distance<-spDistsN1(b$p4, final_centroids[4,], longlat=T)
b$p5@data$distance<-spDistsN1(b$p5, final_centroids[5,], longlat=T)
b$p6@data$distance<-spDistsN1(b$p6, final_centroids[6,], longlat=T)
b$p7@data$distance<-spDistsN1(b$p7, final_centroids[7,], longlat=T)
b$p8@data$distance<-spDistsN1(b$p8, final_centroids[8,], longlat=T)
b$p10@data$distance<-spDistsN1(b$p10, final_centroids[10,], longlat=T)
b$p11@data$distance<-spDistsN1(b$p11, final_centroids[11,], longlat=T)
b$p12@data$distance<-spDistsN1(b$p12, final_centroids[12,], longlat=T)
rm(final_centroids)
#### Rename variables 1----

# ^ alter ^ 
# Comment out/in appropriate f.2! 

# f.2=function(x){
#   rownames(x)<-NULL
#   p<-as.data.frame(x)
#   rownames(p)<-NULL
#   x<-p[,c(1:8)]
# } #HUNDRED ONLY
# OR
f.2=function(x){
  rownames(x)<-NULL
  p<-as.data.frame(x)
  rownames(p)<-NULL
  x<-p[,c(1:12)]} #GEOMETRIC MEAN

#### Centroids cont. ####

b[[1]]<-f.2(b[[1]])
b[[2]]<-f.2(b[[2]])
b[[3]]<-f.2(b[[3]])
b[[4]]<-f.2(b[[4]])
b[[5]]<-f.2(b[[5]])
b[[6]]<-f.2(b[[6]])
b[[7]]<-f.2(b[[7]])
b[[8]]<-f.2(b[[8]])
b[[9]]<-f.2(b[[9]])
b[[10]]<-f.2(b[[10]])

b$p1$province<-"1"
b$p2$province<-"2"
b$p4$province<-"4"
b$p5$province<-"5"
b$p6$province<-"6"
b$p7$province<-"7"
b$p8$province<-"8"
b$p10$province<-"10"
b$p11$province<-"11"
b$p12$province<-"12"

#### Rename variables 2----

# ^ alter ^ 
# Comment out/in appropriate f.3!

# f.3<-function(x){
#   z<-x
#   names(z)<-c("richness",  "total_inds",  "total_species", "WORMS_CLASS", "LON_DEG", "LAT_DEG", "site", "distance", "province")
#   x<-z
# } #HUNDRED ONLY
# # OR
f.3<-function(x){
  z<-x
  names(z)<-c("richness",  "max_inds_richness_",  "100_richness", "median_sample_richness", "total_inds", "sub_sample_size", "total_species", "WORMS_CLASS", "LON_DEG", "LAT_DEG", "site", "distance", "province")
  x<-z
} #GEOMETRIC MEAN

a<-lapply(b, f.3)
rm(f.1, f.2, f.3, b)

#### remove single sites ----

# class issues
f.2<-function(x){split(x, x$WORMS_CLASS, drop=T)} 
b<-lapply(a, f.2)
rm(a, f.2)

# remove classes with only one site
f.1<-function(x)(lapply(x, function(x) if(nrow(x) > 1) {x}))
a<-lapply(b, f.1)
rm(b, f.1)

f.2<-function(x) x[!sapply(x, is.null)] 
mde_prep<-lapply(a, f.2)
rm(f.1, f.2, a)

#### MDE: distance ~ rarefied diversity ####

f.H_mde<-function(x){
  distance<-x$distance
  richness<-x$richness
  reg<-lm(richness ~ distance)
  s<-summary(reg)
  m<-coefficients(reg)[["distance"]]
  b<-coefficients(reg)[["(Intercept)"]]
  p_m<-s$coefficients[2,4]
  p_b<-s$coefficients[1,4]
  cl<-(x[1,4])
  pr<-x[1,9]
  site_n<-nrow(x)
  #distance stats
  r<-(max(x$distance)-min(x$distance))
  std<-sd(x$distance)
  med<-median(x$distance)
  z<-cbind(m,b,s$r.squared,p_m,p_b,as.character(cl),site_n,r,med,std,pr)
  z<-as.data.frame(z)
  names(z)<-c("slope", "intercept","r_squared", "p_slope", "p_intercept", "class", "sites", "site_range", "median_distance", "standard deviation",  "province")
  x<-z
} #HUNDRED ONLY
f.mde<-function(x){
  distance<-x$distance
  richness<-x$richness
  reg<-lm(richness ~ distance)
  s<-summary(reg)
  m<-coefficients(reg)[["distance"]]
  b<-coefficients(reg)[["(Intercept)"]]
  p_m<-s$coefficients[2,4]
  p_b<-s$coefficients[1,4]
  cl<-(x[1,8])
  pr<-x[1,13]
  site_n<-nrow(x)
  #distance stats
  r<-(max(x$distance)-min(x$distance))
  std<-sd(x$distance)
  med<-median(x$distance)
  z<-cbind(m,b,s$r.squared,p_m,p_b,as.character(cl),site_n,r,med,std,pr)
  z<-as.data.frame(z)
  names(z)<-c("slope", "intercept","r_squared", "p_slope", "p_intercept", "class", "sites", "site_range", "median_distance", "standard deviation",  "province")
  x<-z
} #GEOMETRIC MEAN


f.1<-function(x)(lapply (x, f.H_mde))
f.2<-function(x)(lapply (x, f.mde))

# ^ alter ^ 
# Comment out/in appropriate stats! 

#stats<-lapply(mde_prep, f.1) #HUNDRED ONLY
# OR
stats<-lapply(mde_prep, f.2) #GEOMETRIC MEAN

#rm(mde_prep)
rm(f.1, f.2, f.mde, f.H_mde, mde_prep)

#### Stats rearrangement ####

detach("package:plyr", unload = T)
library(tidyr) 

f.1=function(x)(x<-do.call("rbind", x))
a<-lapply(stats, f.1)
names(a)<-names(stats)
b<-f.1(a)
rownames(b)<-NULL
rm(a, f.1)

b[] <- lapply(b, as.character)
b[,c(1:5, 7:10)]<-lapply(b[,c(1:5, 7:10)], as.numeric)
b$class<-as.factor(b$class)
b$province<-as.factor(b$province)
stats<-b
rm(b)

# basic subsetting of data 
all_sig<-subset(stats, p_slope < 0.1)
a<-nrow(all_sig)
non_sig<-subset(stats, p_slope >=0.1)
b<-nrow(non_sig)
mde<-subset(all_sig, slope <0)
c<-nrow(mde)
antimde<-subset(all_sig, slope >=0)
d<-nrow(antimde)
sample<-rbind.data.frame(all_sig, non_sig)

print("CeDAMar MDE")
print(paste0("Total significant (alpha = 0.1) class replicates = ", a))
print(paste0("Total non significant (alpha = 0.1) class replicates = ", b))
print(paste0("significant class replicates supporting MDE (negative regression slope values) = ", c))
print(paste0("significant class replicates opposing MDE (positive or zero regression slope values) = ", d ))

#### Save results ####

print("saving results to output folder")

# ^ alter ^ 
# Comment out/in appropriate save

#saveRDS(stats, file="./output/hundred_mde_cedamar_stats.RDS") # @ future @
# OR
saveRDS(stats, file="./output/mde_cedamar_stats.RDS") # @ future @


print("Script Finished")
script_end<-Sys.time()
script_end - script_start



