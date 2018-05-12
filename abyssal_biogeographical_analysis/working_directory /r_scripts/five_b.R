### Main MDE Anlaysis GBIF ###
# body of script copied from script four. Differences in data rearrangement and rarefication analysis

rm(list=ls())

# Instruction for running script second time:
# Search for " ^ alter ^ ". Uncomment "HUNDRED ONLY" lines. Comment out "GEOMETRIC MEAN".
# Select full function use CTRL + SHIFT + C.
# Four locations: @ Rarefactipon step. @ Rename variables 1 step. @ MDE step. @ Save step.
# Need to run both time for full analysis and graphing.

#### Set directories: #### 
script_start <- Sys.time()

print("Calculating MDE GBIF hundred Rarefaction only")
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

g<-readRDS("./output/gbif.RDS") # % Previous %
provinces<-readRDS("./output/provinces.RDS") # % Previous %
final_centroids<-readRDS("./output/final_centroids.RDS") # % Previous % 

#### Clean data ####

g<-g[which(g$kingdom== g$kingdom[1]) , ]
g<-droplevels(g)
a<-cbind.data.frame(g[,1], g$decimalLongitude, g$decimalLatitude, g$class, g$classKey, g$order, g$orderKey, g$family, g$familyKey, g$genus, g$genusKey, g$species, g$speciesKey)
names(a)<-c("id", "longitude", "latitude", "class", "classKey", "order", "orderKey", "family", "familyKey", "genus", "genusKey", "species", "speciesKey")
b<-a[ , c("longitude", "latitude", "class", "genus", "species")]
g<-b
rm(a,b)

#### Subset GBIF by province ####

coords<-cbind(g$longitude, g$latitude)
cd<-SpatialPointsDataFrame(coords,g)
proj4string(cd)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
fix1<-("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m +lon_wrap=80")
fix2<-("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m +lon_wrap=180")

cd.f1<-spTransform(cd, CRS(fix1))
cd.f2<-spTransform(cd, CRS(fix2))

cd1<-cd.f1[provinces[[1]],]

cd2<-cd[provinces[[2]],]
cd3<-cd[provinces[[3]],] 
cd4<-cd[provinces[[4]],]
cd5<-cd[provinces[[5]],]
cd6<-cd[provinces[[6]],]

cd7<-cd.f2[provinces[[7]],]

cd8<-cd[provinces[[8]],]
cd9<-cd[provinces[[9]],]

cd10<-cd.f2[provinces[[10]],]
cd11<-cd.f2[provinces[[11]],]
cd12<-cd.f2[provinces[[12]],]
cd13<-cd.f2[provinces[[13]],] 
cd14<-cd[provinces[[14]],] 

rm(g, cd, coords, cd.f1, cd.f2, provinces, fix1, fix2)
#### Create a list object of datasets ----

a<-list(cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10, cd11, cd12, cd13, cd14)

# add sites
all<-lapply(a, function(x) {x<- transform(x, site= as.numeric(interaction(longitude, latitude, drop=TRUE)))})
names(all)<-c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10", "p11", "p12", "p13", "p14")
rm(a,cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10, cd11, cd12, cd13, cd14)

#### Data rearrangement ####

# class issues
f.1<-function(x){split(x, x$class, drop=T)}
a<-lapply(all, f.1)

# list made up of a list (per province) of class dataframes
f.2<-function(x){lapply(x, function(x)split(x, x$site, drop =T))}
b<-lapply(a, f.2)
rm(a, all, f.1, f.2)

#split by genus
f.1<-function(x){lapply(x, function(x)(lapply(x, function(x)split(x, x$genus, drop = T)))) }
a<-lapply(b, f.1)
rm(b,f.1)

f.2<-function(x) {
  sp<-nrow(x)
  w<-x[1,c(1:4,9)]
  z<-cbind(w, sp)
  names(z)<-c("longitude", "latitude", "class", "genus", "site", "species_count")
  x<-z
}
f.3<-function(x){lapply(x, f.2)}
f.4<-function(x){lapply(x, function(x)(lapply(x, f.3)))}

b<-lapply(a, f.4)
rm(a, f.2, f.3, f.4)

# recombine
f.1<-function(x)(x<-do.call("rbind", x))
f.2<-function(x){lapply(x, f.1)}
f.3<-function(x){lapply(x, f.2)}
a<-lapply(b, f.3)
rm(b, f.1, f.2, f.3)

# rearrange
f.1<-function(x){lapply(x, function(x){lapply(x, function(x)cbind((spread(x[,c(4,6)], genus, species_count, fill = F, convert = F, drop = T)),(x[1,c(1:3,5)])))})}
b<-lapply(a, f.1)
rm(a, f.1)

detach("package:tidyr", unload = T)
f.ipak("plyr")

f.2<-function(x){lapply(x, rbind.fill)} # rearranges column order
a<-lapply(b, f.2)
rm(b, f.2)

detach("package:plyr", unload=T)
f.ipak("dplyr")

f.1<-function(x){x %>% select(-class, -longitude, -latitude, -site, everything())}
f.2<-function(x){lapply(x, f.1)}
b<-lapply(a, f.2)
rm(a, f.1, f.2)

f.1<-function(x){
  a<-x[,1:((ncol(x))-4)]
  b<-x[, ((ncol(x))-3):(ncol(x))]
  a[is.na(a)]<-0
  cbind(a,b)
}
f.2<-function(x){lapply(x, f.1)}
a<-lapply(b, f.2)
rm(b, f.1, f.2)

# remove all only represent by one genus and nans:
f.1<-function(x)(lapply(x, function(x) if(length(x) > 5) {x}))
b<-lapply(a, f.1)

f.2<-function(x) {x[!sapply(x, is.null)]} 
rare_prep<-lapply(b,f.2)
rm(a, b, f.1, f.2)

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
  names(z)<-c("richness",  "total_species",   "total_genera","longitude", "latitude", "class", "site" ) 
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
  names(z)<-c("richness", "max_sp_richness_", "hundred_richness","sub_sample_richness", "total_species",  "sub_sample_size", "total_genera",  "longitude", "latitude", "class", "site") 
  x<-z
  
}

# ^ alter ^ 
# Comment out/in appropriate f.1! 

f.1<-function(x){lapply(x, f.H_rare)} #HUNDRED ONLY
# OR
# f.1<-function(x){lapply(x, f.rare)} #GEOMETRIC MEAN


a<-lapply(rare_prep, f.1)
rarefied<-lapply(a, function(x)(x<-do.call("rbind", x)))
names(rarefied)<-(names(a))
rm(a, rare_prep, f.1, f.rare)

#### Add centroids ####

f.1<-function(x){
  coords<-cbind(x$longitude, x$latitude)
  z<-SpatialPointsDataFrame(coords,x)
  proj4string(z)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  x<-z
}
b<-lapply(rarefied, f.1)
names(b)<-names(rarefied)
rm(rarefied, f.1)

b$p1@data$distance<-spDistsN1(b$p1, final_centroids[1,], longlat=T)
b$p2@data$distance<-spDistsN1(b$p2, final_centroids[2,], longlat=T)
b$p3@data$distance<-spDistsN1(b$p3, final_centroids[3,], longlat=T)
b$p4@data$distance<-spDistsN1(b$p4, final_centroids[4,], longlat=T)
b$p5@data$distance<-spDistsN1(b$p5, final_centroids[5,], longlat=T)
b$p6@data$distance<-spDistsN1(b$p6, final_centroids[6,], longlat=T)
b$p7@data$distance<-spDistsN1(b$p7, final_centroids[7,], longlat=T)
b$p8@data$distance<-spDistsN1(b$p8, final_centroids[8,], longlat=T)
b$p9@data$distance<-spDistsN1(b$p9, final_centroids[9,], longlat=T)
b$p10@data$distance<-spDistsN1(b$p10, final_centroids[10,], longlat=T)
b$p11@data$distance<-spDistsN1(b$p11, final_centroids[11,], longlat=T)
b$p12@data$distance<-spDistsN1(b$p12, final_centroids[12,], longlat=T)
b$p13@data$distance<-spDistsN1(b$p13, final_centroids[13,], longlat=T)
b$p14@data$distance<-spDistsN1(b$p14, final_centroids[14,], longlat=T)
rm(final_centroids)

#### Rename variables 1----

# ^ alter ^ 
# Comment out/in appropriate f.2! 

f.2=function(x){
  rownames(x)<-NULL
  p<-as.data.frame(x)
  rownames(p)<-NULL
  x<-p[,c(1:8)]
} #HUNDRED ONLY
# OR
# f.2=function(x){
#   rownames(x)<-NULL
#   p<-as.data.frame(x)
#   rownames(p)<-NULL
#   x<-p[,c(1:12)]
# } #GEOMETRIC MEAN

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
b[[11]]<-f.2(b[[11]])
b[[12]]<-f.2(b[[12]])
b[[13]]<-f.2(b[[13]])
b[[14]]<-f.2(b[[14]])

b$p1$province<-"1"
b$p2$province<-"2"
b$p3$province<-"3"
b$p4$province<-"4"
b$p5$province<-"5"
b$p6$province<-"6"
b$p7$province<-"7"
b$p8$province<-"8"
b$p9$province<-"9"
b$p10$province<-"10"
b$p11$province<-"11"
b$p12$province<-"12"
b$p13$province<-"13"
b$p14$province<-"14"
rm(f.2)

#### remove single sites ----

# class issues
f.1<-function(x){split(x, x$class, drop=T)} 
a<-lapply(b, f.1)
rm(b, f.2)

# remove classes with only one site
f.2<-function(x)(lapply(x, function(x) if(nrow(x) > 1) {x}))
b<-lapply(a, f.2)
rm(a, f.2)

f.1<-function(x) x[!sapply(x, is.null)] 
mde_prep<-lapply(b, f.1)
rm(f.1, f.2, b)

#### MDE distance ~ rarfeied diversity ####

f.H_mde<-function(x){
  distance<-x$distance
  richness<-x$richness
  reg<-lm(richness ~ distance)
  s<-summary(reg)
  m<-coefficients(reg)[["distance"]]
  b<-coefficients(reg)[["(Intercept)"]]
  p_m<-s$coefficients[2,4]
  p_b<-s$coefficients[1,4]
  cl<-(x[1,6])
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
  cl<-(x[1,10])
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

stats<-lapply(mde_prep, f.1) #HUNDRED ONLY
# OR
#stats<-lapply(mde_prep, f.2) #GEOMETRIC MEAN

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

print("GBIF hundred MDE")
print(paste0("Total significant (alpha = 0.1) class replicates = ", a))
print(paste0("Total non significant (alpha = 0.1) class replicates = ", b))
print(paste0("significant class replicates supporting MDE (negative regression slope values) = ", c))
print(paste0("significant class replicates opposing MDE (positive or zero regression slope values) = ", d ))

#### Save results ####

print("saving results to output folder")

# ^ alter ^ 
# Comment out/in appropriate save

saveRDS(stats, file="./output/hundred_mde_gbif_stats.RDS") # @ future @
# OR
#saveRDS(stats, file="./output/mde_gbif_stats.RDS") # @ future @



print("Script Finished")
script_end<-Sys.time()
script_end - script_start








