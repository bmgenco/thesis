### loads data and makes figures sections corresponding to numbers and titles within thesis  ###
# Should be similar, but in some cases may vary (plot labels, different centroid locations, post figure editing, etc.)
# useful as a general template
# figure 3.2 & 3.3 require shape files which are too large to upload to git hub. 
print("If script fails to run may be due to font issues")

rm(list=ls())
script_start <- Sys.time()
print("making figures..")
#### Set directories #### 

#setwd("/home/brandon/minerva/abyssal_biogeographical_analysis/") #! For running from absolute location
setwd("..") #! For running script from its own directory. Up one directory

#wd<-(path="./working_directory/") 
data<-(path="./data/")
out<-(path="./output/")
fig<-(path="./output/figures/")
rs<-(path="./r_scripts/")
#setwd(wd)

#### Install/load packages ####

f.ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("sp", "rgdal", "extrafont", "maptools", "maps", "mapproj", "rgeos", "ggplot2", "raster", "marmap", "geoR") #, "colorspace")
f.ipak(packages)

#### Color & Fonts ####

#for figure 3.1 & others:
#package: "colorspace", choose_palette()
#used to produce "prov_colors", sequential color ramps, which tonally should be 
# accessible to more common types of colorblindness.
# see http://colorbrewer2.org/
#pie(rep(1,14), col=prov_colors, labels=prov_colors)
# rearranged color scheme for provinces. Object "prov_colors"

#font_import(pattern = "Arial")
#may need to install fonts first        

#### CHAPTER 1 ####
#### Figure 1.1: Global abyssal depths ####
# created using program "gmt tools"
#### CHAPTER  2 ####
# all figures created in "Mathematica 10"
#### CHAPTER 3  ####
#requires scripts 1, 2, 4, 5 to be run
#### Load objects for chapter 3 tables and graphs ####

c<-read.csv("data/DatabaseSyndeepCedamar.csv")
g<-read.delim("data/GBIF_DARWIN.txt")
ocean<-readOGR((file.path(data,"/ocean_NE_10/")), layer="ne_10m_ocean")
land50<-readOGR((file.path(data,"land_ne_50/shp/")),"ne_50m_land")
land<-readOGR((file.path(data,"land_NE_110/shp/")),"ne_110m_land")
shp<-readOGR((file.path(data, "Provinces/Abyssal_shape/")), "GOODSprovinces_abyssal")
final_centroids<-readRDS((file.path(out, "final_centroids.RDS"))) 
#newshp<-readOGR((file.path(data,"equidistant/")),"abyssal_equidistant_rproj") #created in arcgis from original
#newland<-readOGR((file.path(data,"/equidistant/")),"land_equidistant_reproj") #created in arcgis from original

c.mde<-readRDS((file.path(out, "mde_cedamar_stats.RDS")))
g.mde<-readRDS((file.path(out, "mde_gbif_stats.RDS")))
c.hun<-readRDS((file.path(out, "hundred_mde_cedamar_stats.RDS")))
g.hun<-readRDS((file.path(out, "hundred_mde_gbif_stats.RDS")))

#### 3.1: Abyssal provinces ABS 3.0 ####
# photoshop was used in final figure for thesis to add province names/legend

prov_colors<-c("#007495", "#2D3184", "#00A89A", "#00989A", "#A8DB9E", "#61C499", "#87D09B", "#006190", "#E2EDAF", "#32B69A", "#008798", "#004B8A", "#C7E5A5", "#F9F2BD")

postscript("output/figures/3_1 Abyssal provinces ABS 3.0.eps", family = "ArialMT")
plot(land, col="black", bg="grey40", border=FALSE)
plot(shp, col=(prov_colors[shp$Province]), border=FALSE,  add=TRUE)
plot(final_centroids, add=TRUE, col="grey90", pch=21, bg=prov_colors, cex=2.1 )
text(final_centroids$x, final_centroids$y, labels=c(1:14), col="grey5", cex=0.69)
dev.off()
rm(prov_colors)

#### 3.2: Province area ####
# new.shp<-spTransform(shp, CRS("+proj=laea +towgs84=0,0,0"))
# new.land<-spTransform(land, CRS("+proj=laea +towgs84=0,0,0"))
# 
# pdf("output/figures/3_2 Province area.pdf",  height=6, width=6)
# plot(new.shp, border="gray70", bg="grey40", col="white", lty=0)
# plot(new.land, border="black", bg="transparent", lty=2, col="transparent", lwd=0.5, add=T )
# dev.off()
# rm(new.land, new.shp, land)

#### 3.3: Equal area & Equidistant abyssal projections ####

# areashp<-spTransform(newshp, CRS("+proj=laea +lat_0=-27.75 +lon_0=-163.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# arealand<-spTransform(newland, CRS("+proj=laea +lat_0=-27.75 +lon_0=-163.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# # *improve code* grey fill issues
# # requires photoshop editing to remove extranous polygon lines (e.g., Antartica)
# 
# pdf("output/figures/3_3 Equal area & Equidistant abyssal projections", height=8, width=16)
# par(mfrow=c(1,2))
# plot(areashp, border="gray70", bg="grey40", col="white", lty=0, main="A: Equal Area")
# plot(arealand, border="black", bg="transparent", lty=2, col="transparent", lwd=0.5, add=T )
# 
# plot(newshp, border="gray70", bg="grey40", col="white", lty=0, main="B: Equal Distance")
# plot(newland, border="black", bg="transparent", lty=2, col="transparent", lwd=0.5, add=T )
# dev.off()
# rm(newshp, newland, arealand, areashp)

#### 3.4: Biological databases ocean sites ####

g<-subset(g, g$hasGeospatialIssues == "false" & g$hasCoordinate =="true")
g<- Filter(function(x)!all(is.na(x)), g)

gsites<-unique(g[c("decimalLongitude", "decimalLatitude")])
gloc<-cbind(gsites$decimalLongitude, gsites$decimalLatitude)
gpts<-SpatialPoints(gloc)

csites<-unique(c[c("LON_DEG", "LAT_DEG")])
cloc<-cbind(csites$LON_DEG, csites$LAT_DEG)
cpts<-SpatialPoints(cloc)

proj4string(cpts)<-CRS("+proj=longlat +datum=WGS84")
proj4string(gpts)<-CRS("+proj=longlat +datum=WGS84")
cpts<-spTransform (cpts, proj4string(shp))
gpts<-spTransform (gpts, proj4string(shp))

o.gpts<-gpts[ocean,]
o.cpts<-cpts[ocean,]

m.land<-SpatialPolygons2map(land50, namefield=NULL)

setEPS()
postscript("output/figures/3_4 Biological databases ocean sites.eps")
map(m.land, col = "grey55", fill = T, border=FALSE)
map.axes(cex.axis=0.8)
plot(o.gpts, pch=20, cex=0.15, add=TRUE, col="green")
plot(o.cpts, pch=20, cex=0.15, add=TRUE, col="magenta")
dev.off()
rm(c, cloc, csites, g, gloc, gsites, land50, ocean, o.gpts, o.cpts)

#### 3.5: Biological databases abyssal sites ####
# run required lines from previous figure
# *improve code* Fill issues in province, requires post editing to fill color issues (e.g., North Atl provinnce Mid ocean ridge)

a.gpts<-gpts[shp,]
a.cpts<-cpts[shp,]
m.shp<-SpatialPolygons2map(shp)

setEPS()
postscript("output/figures/3_5 Biological databases abyssal sites.eps")
map(m.shp, col="white", fill=T, lwd=.4, bg="grey68", border="grey68")
map(m.land, col = "grey50", fill = T, border=FALSE, add=T)
map.axes(cex.axis=0.7)
plot(a.gpts, pch=20, cex=0.15, add=TRUE, col="green")
plot(a.cpts, pch=20, cex=0.15, add=TRUE, col="magenta")
dev.off()
rm(a.cpts, a.gpts, cpts, gpts, m.land, m.shp, shp)

#### 3.6: Class replicates adjusted species slopes per province - CeDAMar ####

all_sig<-subset(c.mde, p_slope < 0.1)
nonsig<-subset(c.mde, p_slope >=0.1)
nulls<-nonsig
nulls$slope<-0

sample<-rbind.data.frame(all_sig, nonsig)
sample$province<-as.factor(sample$province)
sample_changed<-rbind.data.frame(all_sig, nulls)
sample_changed$province<-as.factor(sample$province)

a<-names(sample_changed)
a[1]<-("regression slope diversity to centroid distance km")
names(sample_changed)<-a
b<-names(sample_changed)
b[1]<-("Slope")
names(sample_changed)<-b
rm(all_sig, nonsig, nulls, sample)

postscript("output/figures/3_6 Class replicates adjusted species slopes per province - CeDAMar.eps", family = "ArialMT")
par(mar=c(8, 5, 4.1, 5))
t<-ggplot(sample_changed, aes(x = class, y =Slope)) + geom_point() + facet_grid(~province)
t + facet_wrap(~province, nrow=2, scale="free", dir="h") + 
  labs(y="Plots per province: Class regression values = alpha diversity by km from centroid") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
rm(a,b, sample_changed, t, c.mde)

#### 3.7: Box-plots per class replicates - CeDAMar ####
# may need to change orientation of final figure

all_sig<-subset(c.hun, p_slope < 0.1)
all_sig<-droplevels(all_sig)

postscript("output/figures/3_7 Box-plots per class replicates - CeDAMar.eps ", family = "ArialMT")
par(mar=c(8, 5.5, 4.1, 2.1))
par(mgp=c(4.5,1,0))
plot(all_sig$class,  all_sig$slope, las=2, ylab="Regression slope: Alpha diversity by km from centroid")
abline(h=0, col="red", lty=2)
dev.off()
rm(all_sig, c.hun)

#### 3.8: Class replicates genus slopes per province - GBIF ####

all_sig<-subset(g.mde, p_slope < 0.1)
all_sig$province<-as.factor(all_sig$province)

postscript("output/figures/3_8 Class replicates genus slopes per province - GBIF.eps", family = "ArialMT")
par(mar=c(8, 5, 4.1, 5.2))
t<-ggplot(all_sig, aes(x = class, y =slope)) + geom_point() + facet_grid(~province)
t + facet_wrap(~province, nrow=2, scale="free", dir="h") + 
  labs(y="Plots per province: Signifigant class regression values = Genera diversity by km from centroid") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

#### 3.9: Box-plots per class replicates - GBIF ####

all_sig<-subset(g.hun, p_slope < 0.1)
all_sig<-droplevels(all_sig)

postscript("output/figures/3_9 Box-plots per class replicates - GBIF.eps ", family = "ArialMT")
par(mar=c(8, 5.5, 4.1, 2.1))
par(mgp=c(4.5,1,0))
plot(all_sig$class,  all_sig$slope, las=2, ylab="Regression slope: Genera diversity by km from centroid")
abline(h=0, col="red", lty=2)
dev.off()
rm(all_sig, g.hun)

#### CHAPTER 4  ####
#requires scripts 2 & 6 to be run
#### Load object for chapter 4 tables and graphs ####

provinces<-readRDS((file.path(out, "provinces.RDS")))
load("./output/data_from_six.RData")
NE<-readRDS(file.path(data, "nenc.rda"))
land<-readOGR((file.path(data,"land_NE_110/shp/")),"ne_110m_land")

#### 4.1: Case study sites and North Atlantic province ####
# *improve code* Fill issues in province, requires post editing to fill color issues (province 2, Mid ocean ridge)

p2<-provinces[[2]]
coast<-crop(land, extent(-75, 0, 20, 50))
print("cpu intensive step")
p2cropped<-crop(p2, extent(-75, 0, 20, 50)) # time consuming
print("cpu intensive step over")

m.coast<-SpatialPolygons2map(coast, namefield=NULL)
m.crop<-SpatialPolygons2map(p2cropped, namefield=NULL)

pdf("output/figures/4_1 Case study sites and North Atlantic province.pdf")
map(m.coast, fill =T, col="grey50", bg="grey68")
map(m.crop, fill=T, col="white", add=T, border="grey30")
map.axes(cex.axis=0.7)
plot(neopts, add=T, pch=21, bg ="orangered", col="midnightblue", cex=1.2)
dev.off()
rm(m.coast, m.crop, p2cropped, coast, p2)

#### 4.2: Study Site: NW Atlantic - 24 sites. ####

ramp<-c("#2D3184", "#254289", "#1C518F", "#115F96", "#046C9C", "#0078A2", "#0084A7",
"#068FAB", "#189AAF", "#28A4B3", "#45B6B8", "#54BEBA", "#62C5BC", "#70CCBD", "#7DD2BF", 
"#8BD8C0", "#97DDC1", "#A3E2C3", "#AFE5C4", "#BAE9C6", "#C4EBC8", "#CDEECA", "#D6F0CD",
"#DEF1D0", "#E4F2D3", "#EAF2D6", "#EEF2DA", "#F2F2DE", "#F3F1E4")

bmp("output/figures/4_2 Study Site: NW Atlantic - 24 sites.bmp", res=600, height=4800, width=4800)
plot(NE, bpal=ramp,  image=T, deep=c(-6500,0), shallow=c(-50,0), step=c(200,0), 
     lwd=c(0.4,0.8), lty=c(1,1))
scaleBathy(NE, deg=2, y=42.5,x=-75)
points(neopts, pch=21, col="orange",bg=col2alpha("red",.4),cex=.8)
dev.off()
rm(ramp)

#### 4.3: Diversity and abundance by depth ####

h.lm<-lm(h$COUNT~h$MAX_DEPTH)
tiff("output/figures/4._3 Diversity and abundance by depth.tiff",  res=600, height=4800, width=9600)
par(mfrow=c(1,2))
plot(h.div$Diversity ~ h.div$Depth, xlab=" Depth - m", ylab="Species per Site", main="A: Diversity",
     pch=25, col='blue', bg='blue',cex=4, cex.lab=1.5, cex.axis=1.5)
abline(div_by_dept, col="red")

plot(h$COUNT~ sort(h$MAX_DEPTH), xlab= "Depth - m", ylab="Individuals per Species",  main="B: Abundance (all species)", pch=21, bg="orange", col="black", cex.lab=1.5, cex.axis=1.5)
abline(h.lm, col="red")
dev.off()

#### 4.4: Diversity semivariograms ####
# *improve code* does not correctly provide a title fo Panel A

setEPS()
postscript("output/figures/4_4 Diversity semivariograms.eps", width=8, height=4.5)
par(mfrow=c(1,2))
plot(d.v4, xlab="Distance - meters", main="A: Omndirectional semivariogram",  lwd=2)
abline(h=var(d$data), col="orange")
plot(d.v, type="b", xlab="Distance - meters", main="B: Directional semivariograms",  lwd=2, cex.main=0.75)
abline(h=var(d$data), col="orange")
dev.off()


#### 4.5: Three species semivariograms ####

setEPS()
postscript("output/figures/4_5 Three species semivariograms.eps", width=9, height=8)
par(mfrow=c(2,3))
plot(t.v, pch=19, col="red", cex=2, cex.axis=1,  cex.sub=1.5, main="T. lyronuclea", xlab="Distance - meters")
abline(h=var(t$data), col="red")
lines(wlsm.t, lwd=2)

plot(a.v, pch=24, col="black", bg="orange", cex=2,  cex.axis=1,  cex.sub=1.5, main="A. melampoides", xlab="Distance - meters")
abline(h=var(a$data), col="orange")
lines(wlsm.a, lwd=2)

plot(b.v, pch=15, col="blue", cex=2,  cex.axis=1,  cex.sub=1.5, main="B. tenella", xlab="Distance - meters")
abline(h=var(b$data), col="blue")
lines(wlsm.b, lwd=2)

plot(t.v4, xlab="Distance - meters", cex=2, cex.lab=1.5, cex.axis=1,  cex.sub=1.5, lwd=2)
abline(h=var(t$data), col="red")  

plot(a.v4, xlab="Distance - meters",cex=2, cex.lab=1.5, cex.axis=1,  cex.sub=1.5, lwd=2)
abline(h=var(a$data), col="orange")  

plot(b.v4, xlab="Distance - meters", cex=2, cex.lab=1.5, cex.axis=1,  cex.sub=1.5, lwd=2)
abline(h=var(b$data), col="blue") 
dev.off()

#### 4.6: Continental rise sites ####

ramp<-c("#2D3184", "#254289", "#1C518F", "#115F96", "#046C9C", "#0078A2", "#0084A7",
        "#068FAB", "#189AAF", "#28A4B3", "#45B6B8", "#54BEBA", "#62C5BC", "#70CCBD", "#7DD2BF", 
        "#8BD8C0", "#97DDC1", "#A3E2C3", "#AFE5C4", "#BAE9C6", "#C4EBC8", "#CDEECA", "#D6F0CD",
        "#DEF1D0", "#E4F2D3", "#EAF2D6", "#EEF2DA", "#F2F2DE", "#F3F1E4")



bmp("output/figures/4_6 Continental rise sites.bmp", res=600, height=4800, width=4800)
plot(NE, bpal=ramp,  image=T, deep=c(-6500,0), shallow=c(-50,0), step=c(200,0), 
     lwd=c(0.4,0.8), lty=c(1,1))
scaleBathy(NE, deg=2, y=42.5,x=-75)
points(shelfpts, pch=21, col="orange",bg=col2alpha("red",.4),cex=.8)
dev.off()



rm(h, h.div, div_by_dept, d, d.v, d.v4, a, a.v, a.v4, wlsm.a, b, b.v, b.v4, wlsm.b, t, t.v, t.v4, wlsm.t, neopts, shelfpts, NE)
print("plots saved in figures directory")
print("Script Finished")
script_end<-Sys.time()
script_end - script_start
