## Centroid Calculation (Figure 3.1)  ##

rm(list=ls())

#### Set directories: #### 
script_start <- Sys.time()

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
packages <- c("sp", "rgdal", "rgeos") # add parallel package
f.ipak(packages) 
rm(f.ipak, packages)

#### Load data: ####

shp<-readOGR((file.path(data, "Provinces/Abyssal_shape/")), "GOODSprovinces_abyssal")
#load("./output/all_centroids.RData")

#### Split provinces  ####
# method for dealing with province spatial objects

for (i in unique(shp@data$Province)){
nam<- paste("p", i, sep =".")
a<-assign(nam, shp[shp@data$Province==i,])
} # loop which creates new SpatialPolygonDataFrame object per province

provinces<-c(p.1, p.2, p.3, p.4, p.5, p.6, p.7, p.8, p.9, p.10, p.11, p.12, p.13, p.14) # @ future @ 
rm(nam, a, i, shp) 

#### Save as KMLs ####
# This section is for visually determining province extent to reproject provinces that are split in
# normal -180 - 180 projections. This causes errors in centroid calculations.
# Reprojections are coded into centroids calculations in following section. Left this section in code to indicate method.

#plot(shp)

#b<-lapply(provinces, function(x){
  #x<-x[1]
  #a<-bbox(x)
  #return(a)}) # review bbox extent to determine which provinces need to be reprojected

#setwd(out)

#writeOGR(p.1, dsn="province_1.kml", layer="Individual Polygons", driver="KML")
#writeOGR(p.7, dsn="province_7.kml", layer="Individual Polygons", driver="KML")
#writeOGR(p.11, dsn="province_11.kml", layer="Individual Polygons", driver="KML")
#writeOGR(p.12, dsn="province_12.kml", layer="Individual Polygons", driver="KML")
#writeOGR(p.13, dsn="province_13.kml", layer="Individual Polygons", driver="KML")

#setwd("..")

#### Reproject certain provinces ####

fix1<-("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m +lon_wrap=80")
fix2<-("+proj=longlat +ellps=WGS84 +datum=WGS84 +units=m +lon_wrap=180")

newp.1<-spTransform(p.1, CRS(fix1))
newp.7<-spTransform(p.7, CRS(fix2))
newp.10<-spTransform(p.10, CRS(fix2))
newp.11<-spTransform(p.11, CRS(fix2))
newp.12<-spTransform(p.12, CRS(fix2))
newp.13<-spTransform(p.13, CRS(fix2))

provinces[[1]]<-newp.1
provinces[[7]]<-newp.7
provinces[[10]]<-newp.10
provinces[[11]]<-newp.11
provinces[[12]]<-newp.12
provinces[[13]]<-newp.13

rm(p.1, p.2, p.3, p.4, p.5, p.6, p.7, p.8, p.9, p.10, p.11, p.12, p.13, p.14, newp.1, newp.7, newp.10, newp.11, newp.12, newp.13, fix2, fix1)

#### Centroids: ####
# See Table 3.3: Province centroids, and pg 54 of thesis for text decription of steps
print("Centroids being calculated, CPU intensive and time consuming: an hour +")

f.centroids<-function(x){
  # Function calculates centroids per province on orignal mercator projection. Uses that point for reprojection to Azimuthal Equidistant. 
  # The process is repreated for projection on lamebert eqaul area for final centroid.   
  
  x<-x[1]
  
  cent.merc<-gCentroid(x, byid=F)
  
  nlong<-cent.merc@coords[1]
  nlat<-cent.merc@coords[2]
  aeqd<-paste("+proj=aeqd +lat_0=", nlat, " +lon_0=", nlong, " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", sep = "")
  p_aeqd<-spTransform(x, CRS(aeqd))
  
  cent.aeqd<-gCentroid(p_aeqd, byid=F)
  
  nnlong<-cent.aeqd@coords[1]
  nnlat<-cent.aeqd@coords[2]
  laea<-paste("+proj=laea +lat_0=", nlat, " +lon_0=", nlong, " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", sep = "")
  p_laea<-spTransform(x, CRS(laea))
  
  cent.laea<-gCentroid(p_laea, byid=F)
  
  cent.final<-spTransform(cent.laea, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))
  multi<-c(cent.merc, cent.aeqd, cent.laea, cent.final)
  return(multi)}

cent_start<-Sys.time()
all_centroids<-lapply(provinces, f.centroids) 
final_centroids<-sapply(all_centroids, function(x){
  a<-x[4]
  return(a)
}) # @ future @ 
end_time <- Sys.time()
end_time - cent_start

#### Convert final centroids to dataframe as opposed to a list ####
# *improve code*

m<-matrix(nrow=14, ncol=2)
a<-final_centroids

m[1,1]<-a[[1]]@coords[1]
m[1,2]<-a[[1]]@coords[2]
m[2,1]<-a[[2]]@coords[1]
m[2,2]<-a[[2]]@coords[2]
m[3,1]<-a[[3]]@coords[1]
m[3,2]<-a[[3]]@coords[2]
m[4,1]<-a[[4]]@coords[1]
m[4,2]<-a[[4]]@coords[2]
m[5,1]<-a[[5]]@coords[1]
m[5,2]<-a[[5]]@coords[2]
m[6,1]<-a[[6]]@coords[1]
m[6,2]<-a[[6]]@coords[2]
m[7,1]<-a[[7]]@coords[1]
m[7,2]<-a[[7]]@coords[2]
m[8,1]<-a[[8]]@coords[1]
m[8,2]<-a[[8]]@coords[2]
m[9,1]<-a[[9]]@coords[1]
m[9,2]<-a[[9]]@coords[2]
m[10,1]<-a[[10]]@coords[1]
m[10,2]<-a[[10]]@coords[2]
m[11,1]<-a[[11]]@coords[1]
m[11,2]<-a[[11]]@coords[2]
m[12,1]<-a[[12]]@coords[1]
m[12,2]<-a[[12]]@coords[2]
m[13,1]<-a[[13]]@coords[1]
m[13,2]<-a[[13]]@coords[2]
m[14,1]<-a[[14]]@coords[1]
m[14,2]<-a[[14]]@coords[2]

b<-as.data.frame(m)
tmp<-SpatialPointsDataFrame(coords = b, proj4string = a[[1]]@proj4string, data = b)
names(tmp)<-c("x","y")

final_centroids<-tmp #! comment out if desired.

#### Save results ####

print("saving results to output folder")

#save(all_centroids, file = "./output/all_centroids_v2.RData") 
saveRDS(final_centroids, file = "./output/final_centroids.RDS") # @ future @
saveRDS(provinces, file="./output/provinces.RDS") # @ future @

print("Script Finished")
script_end<-Sys.time()
script_end - script_start







