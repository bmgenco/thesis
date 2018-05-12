##  SUBPROVINCE SCALE & GASTROPOD DIVERSITY cASE STUDY: NORTH ATLANTIC PROVINCE  ##

rm(list=ls())

#### Set directories: #### 
script_start <- Sys.time()
print("Calculating Gastpod analysis subprovince scale")

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
packages <- c("lattice", "marmap", "sp", "geoR", "rgdal", "rgeos", "raster", "tidyr", "dplyr", "ggplot2" ) # "ncdf4"
f.ipak(packages) 
rm(f.ipak, packages)

#### Conversion of netcdf to r object   ####
# r object is provided in data folder, negating this step
# Coordinates New England/ NW Atlantic75/-60/31.5/42.5 #
# nc<-raster("/ne.grd", varname="z")
# nenc<-as.bathy(nc)
# # saveRDS(nenc, file="/nenc.rda")

#### Load Data ####

c<-read.csv("data/DatabaseSyndeepCedamar.csv")
NE<-readRDS(file.path(data, "nenc.rda"))

#### Select NW Atl sites ####

osites<-unique(c[c("LON_DEG", "LAT_DEG")])
oloc<-cbind(osites$LON_DEG, osites$LAT_DEG)
opts<-SpatialPoints(oloc)
proj4string(opts)<-CRS("+proj=longlat +datum=WGS84")

h<-c[c$LON_DEG> -75 & c$LON_DEG< -60 & c$LAT_DEG> 31.5 & c$LAT_DEG< 42.5,]

NEbound<-as.SpatialGridDataFrame(NE)
neopts<-opts[NEbound,] # @ future @ used in figures
rm(c, opts, oloc, osites, NEbound)

#Date conversion!!
sites<-unique(h$EnvSiteID)
h$START_DATE<-as.Date(h$START_DATE, "%d/%m/%Y")

#### Lat long-> UTM; utm zone 19n for use in geoR  #####

### shelf points

d<-h[!duplicated(h$EnvSiteID), c("EnvSiteID", "LON_DEG", "LAT_DEG", "START_DATE", "MAX_DEPTH")]
names(d)<-c("EnvSiteID",  "Longitude", "Latitude", "Date" , "Depth")
d$Depth<-as.numeric(paste(d$Depth))
shelf<-d[d$Depth> 3500 & d$Depth< 4000,]
shloc<-cbind(shelf$Longitude, shelf$Latitude)
shelfpts<-SpatialPoints(shloc)
proj4string(shelfpts)<-CRS("+proj=longlat +datum=WGS84")

hloc<-cbind(h$LON_DEG, h$LAT_DEG)
hpts<-SpatialPoints(hloc)
proj4string(hpts)<-CRS("+proj=longlat +datum=WGS84")
res <- spTransform(hpts, CRS("+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
h.utm<-as.data.frame(res)
h$LON_DEG<-h.utm$coords.x1
h$LAT_DEG<-h.utm$coords.x2
h$MAX_DEPTH<-as.numeric(paste(h$MAX_DEPTH))
rm(res, h.utm, hloc, d, shelf, shloc)

#### Alpha diversity: ####

id<-unique(h$EnvSiteID) 
a<-(1:n_distinct(h$EnvSiteID))
for (i in  1:n_distinct(h$EnvSiteID)){
a[i]<-sum(h$EnvSiteID==id[i])}
b<-cbind(id,a)
rm(a, id, i)
b<-as.data.frame(b)
names(b)<-c("EnvSiteID", "Diversity") #b is alpha diversity per site

d<-h[!duplicated(h$EnvSiteID), c("EnvSiteID", "LON_DEG", "LAT_DEG", "START_DATE", "MAX_DEPTH")]
names(d)<-c("EnvSiteID",  "Longitude", "Latitude", "Date" , "Depth")
# site 1652 had date error, listed as 04/31/1966, there are only 30 days in april.. produces NAs 
# diversity site + depth, coords, date:
h.div<-merge(b, d, by="EnvSiteID")
rm(b,d)
h.div$EnvSiteID<-(1:24)
names(h.div)<-c("Site", "Diversity", "Longitude", "Latitude", "Date" , "Depth")

#### Abundance: ####

#*improve code*
#species selection (replace this with random & set seed):
#sample(1:81, 3, replace=F)
#gives: 12, 39, 27

#h$SPECIES[12]
#h$SPECIES[39]
#h$SPECIES[27]

#*improve code*
# creates 0 counts for 3 unique data sets based on each species number
COUNT<-rep(0,24)
hsites<-unique(h[c("LON_DEG", "LAT_DEG")])

a<-cbind(hsites$LON_DEG, hsites$LAT_DEG, sites, COUNT)
a<-as.data.frame(a)
names(a)<-c("LON_DEG", "LAT_DEG","EnvSiteID", "COUNT")
a$EnvSiteID<-as.integer(a$EnvSiteID)

#Aceteon melampoides
h.12<-h[which(h$SPECIES== h$SPECIES[12]) , ]
h.12<-h.12[, c(24,23,1,20)]
b<-setdiff(a[,(1:3)], h.12[,(1:3)])
b$COUNT<-(rep(0,length(b$LON_DEG)))
h.12<-union(b, h.12)

#Benthonella tenella
h.39<-h[which(h$SPECIES== h$SPECIES[39]) , ]
h.39<-h.39[, c(24,23,1,20)]
b<-setdiff(a[,(1:3)], h.39[,(1:3)])
b$COUNT<-(rep(0,length(b$LON_DEG)))
h.39<-union(b, h.39)

#Theta lyronuclea
h.27<-h[which(h$SPECIES== h$SPECIES[27]) , ]
h.27<-h.27[, c(24,23,1,20)]
b<-setdiff(a[,(1:3)], h.27[,(1:3)])
b$COUNT<-(rep(0,length(b$LON_DEG)))
h.27<-union(b, h.27)

rm(a, hsites, COUNT, b, hpts)

#### Non-Spatial Analysis: ####
div_by_dept<-lm(h.div$Diversity ~ h.div$Depth)

#### Geostats NE: ####
#### Diversity: ----
d<-as.geodata(h.div, coords.col =3:4, data.col = 2 )
d.var<-(d$data)
d.v<-variog(d)
d.v4<-variog4(d)

#### Abundances: ----

#species a: Aceteon melampoides
a<-as.geodata(h.12, coords.col=1:2,data.col=4)
a.var<-(a$data)
a.v <- variog(a)
a.v4 <- variog4(a)
wlsm.a <- variofit(a.v, ini= c(35,2*10^5), cov.model="matern", nugget=20)
reml.a <- likfit(a, ini.cov.pars = c(30,2*10^5), nug = 20, cov.model="exponential", lik.method = "RML")
#summary(wlsm.a)

#species b: Benthonella tenella
b<-as.geodata(h.39, coords.col=1:2,data.col=4)
b.var<-(b$data)
b.v<-variog(b)
b.v4 <- variog4(b)
wlsm.b <- variofit(b.v, ini= c(600,4*10^5), cov.model="matern", nugget=500)
reml.b <- likfit(b, ini.cov.pars = c(30,2*10^5), nug = 20, cov.model="exponential",lik.method = "RML")
#summary(wlsm.b)

#species Theta lyronuclea
t<-as.geodata(h.27, coords.col=1:2,data.col=4)
summary(t)
t.var<-(t$data)
t.v<-variog(t)
t.v4 <- variog4(t)
wlsm.t <- variofit(t.v, ini= c(0.3,5*10^5), cov.model="matern", nugget=1.7)
#summary(wlsm.t)

####   Krigging: ----
#set area to predict, fit to
# not included in thesis
# xr<-max(h.27$LON_DEG)-min(h.27$LON_DEG)
# yr<-max(h.27$LAT_DEG)-min(h.27$LAT_DEG)
# 
# xo<-(min(h.27$LON_DEG)-2000)
# xn<-(max(h.27$LON_DEG)+2000)
# yo<-(min(h.27$LAT_DEG)-2000)
# yn<-(max(h.27$LAT_DEG)+2000)
# 
# #resolution?
# x<-seq(xo,xn,250)
# y<-seq(yo,yn,250)
# 
# d1 <- expand.grid(x=x,y=y)
# t.conv <- krige.conv(t, loc=d1, krige=krige.control(obj.m = wlsm.t))
# image(t.conv)


#### Save for plotting ####
print("saving results to output folder")
save(h, h.div, div_by_dept, d, d.v, d.v4, a, a.v, a.v4, wlsm.a, b, b.v, b.v4, wlsm.b, t, t.v, t.v4, wlsm.t, neopts, shelfpts, file ="./output/data_from_six.RData")
print("Script Finished")
script_end<-Sys.time()
script_end - script_start
