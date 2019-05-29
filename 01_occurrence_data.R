# Load --------------------------------------------------------------------

# ATTENTION: The next line should set the working directory to the source file location. This only works in RStudio. If you use RGUI or experience problems please replace with full path!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(raster)
library(maptools)
library(sp)
library(dismo)
library(rgdal)
library(rgeos)

# data can be downloaded from https://figshare.com/collections/A_high-resolution_pan-European_tree_occurrence_dataset/3288407
EUForestspecies <- read.csv("EUForestspecies.csv")
# available from NUTS or Natural Earth (needs pre-processing via union() of countries)
eu_shp <- readOGR(dsn=paste0(getwd(),"/europe_coastline.shp"))
# shapefile of adminstrative borders of Germany is available from GADM
ger_shp <- readOGR(dsn=paste0(getwd(),"/DEU_adm0.shp"))


# Retrieve species ------------------------------------------------------------------

aalba_eu <- EUForestspecies[EUForestspecies$`SPECIES.NAME`=="Abies alba",]
fsylva_eu <- EUForestspecies[EUForestspecies$`SPECIES.NAME`=="Fagus sylvatica",]

# mixed forest
posmix <- match(interaction(aalba_eu$X,aalba_eu$Y),interaction(fsylva_eu$X,fsylva_eu$Y))
posmix <- na.omit(posmix)
mix_eu <- fsylva_eu[posmix,]



# Convert coordinates -----------------------------------------------------

fsylva_sp <- fsylva_eu[,c(1,2)]
coordinates(fsylva_sp) <- ~X+Y
crs(fsylva_sp) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
fsylva_transf <- spTransform(fsylva_sp,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
euf_fsylva_eu <- as.data.frame(fsylva_transf)


aalba_sp <- aalba_eu[,c(1,2)]
coordinates(aalba_sp) <- ~X+Y
crs(aalba_sp) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
aalba_transf <- spTransform(aalba_sp,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
euf_aalba_eu <- as.data.frame(aalba_transf)

mix_sp <- mix_eu[,c(1,2)]
coordinates(mix_sp) <- ~X+Y
crs(mix_sp) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
mix_transf <- spTransform(mix_sp,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
euf_mix_eu <- as.data.frame(mix_transf)


# Display occurrences -----------------------------------------------------------------

data("wrld_simpl")
plot(wrld_simpl, main="Species distribution ranges",xlim=c(-10,28), ylim=c(37,60), axes=T, col="light yellow")
legend(30,45,legend=c("beech","fir","beech-fir"),col=c("lightgreen","darkgreen","orange"),cex=1, pch=20)
points(euf_fsylva_eu$X,euf_fsylva_eu$Y,pch=20,cex=0.5,col="lightgreen")
points(euf_aalba_eu$X,euf_aalba_eu$Y,pch=20,cex=0.5,col="darkgreen")
points(euf_mix_eu$X,euf_mix_eu$Y,pch=20,cex=0.5,col="orange")
box()



# [SKIP FOR NOW] Select pseudo-absences for Europe--------------------------------------------------

# create circles-model (buffer) around presences
c_a_eu <- circles(aalba_sp, d=50000, lonlat=TRUE)
c_f_eu <- circles(fsylva_sp, d=50000, lonlat=TRUE)
c_m_eu <- circles(mix_sp, d=50000, lonlat=TRUE)

pol_a_eu <- polygons(c_a_eu)
pol_f_eu <- polygons(c_f_eu)
pol_m_eu <- polygons(c_m_eu)

# clip with EU coastline to exclude offshore pixels
pol_a_eu_clip <- gIntersection(pol_a_eu,eu_shp,byid=T,drop_lower_td = TRUE)
pol_f_eu_clip <- gIntersection(pol_f_eu,eu_shp,byid=T,drop_lower_td = TRUE)
pol_m_eu_clip <- gIntersection(pol_m_eu,eu_shp,byid=T,drop_lower_td = TRUE)




# sample randomly from all circles, prevalence arbitrarily set to 0.5
set.seed(1992)
samp_a_eu <- spsample(pol_a_eu_clip, 20000, type='random', iter=10)
samp_f_eu <- spsample(pol_f_eu_clip, 40000, type='random', iter=10)
samp_m_eu <- spsample(pol_m_eu_clip, 10000, type='random', iter=10)
mask <- raster(paste0(getwd(),"/eu_AnnMeanT.tif"))

# get unique cells
cells_a_eu <- cellFromXY(mask, samp_a_eu)
cells_a_eu <- unique(cells_a_eu)
xy_a_eu <- xyFromCell(mask, cells_a_eu)

cells_f_eu <- cellFromXY(mask, samp_f_eu)
cells_f_eu <- unique(cells_f_eu)
xy_f_eu <- xyFromCell(mask, cells_f_eu)

cells_m_eu <- cellFromXY(mask, samp_m_eu)
cells_m_eu <- unique(cells_m_eu)
xy_m_eu <- xyFromCell(mask, cells_m_eu)

# check results on map
data("wrld_simpl")
par(mfrow=c(2,2))
plot(wrld_simpl, main="silver fir",xlim=c(-10,28), ylim=c(37,60), axes=T, col="light yellow")
points(xy_a_eu, cex=0.2, pch=20, col='blue')
plot(aalba_sp,add=T,pch=20,col="red",cex=0.2)

plot(wrld_simpl, main="beech",xlim=c(-10,28), ylim=c(37,60), axes=T, col="light yellow")
points(xy_f_eu, cex=0.2, pch=20, col='blue')
plot(fsylva_sp,add=T,pch=20,col="red",cex=0.2)

plot(wrld_simpl, main="mixed",xlim=c(-10,28), ylim=c(37,60), axes=T, col="light yellow")
points(xy_m_eu, cex=0.2, pch=20, col='blue')
plot(mix_sp,add=T,pch=20,col="red",cex=0.2)


euf_aalba_eu_abs <- cbind(PRESENCE=0,xy_a)
colnames(euf_aalba_eu_abs) <- c("PRESENCE","X","Y")
euf_fsylva_eu_abs <- cbind(PRESENCE=0,xy_f)
colnames(euf_fsylva_eu_abs) <- c("PRESENCE","X","Y")
euf_mix_eu_abs <- cbind(PRESENCE=0,xy_m)
colnames(euf_mix_eu_abs) <- c("PRESENCE","X","Y")

euf_aalba_eu <- rbind(as.data.frame(aalba_eu),euf_aalba_eu_abs)
euf_fsylva_eu <- rbind(as.data.frame(fsylva_eu),euf_fsylva_eu_abs)
euf_mix_eu <- rbind(as.data.frame(mix_eu),euf_mix_eu_abs)

save(euf_aalba_eu, file="euf_aalba_eu")
save(euf_fsylva_eu, file="euf_fsylva_eu")
save(euf_mix_eu, file="euf_mix_eu")


# Delimiting the geographical background ----------------------------------

ger_proj <- spTransform(ger_shp,CRS("+init=epsg:25832"))
# set 50 km buffer around Germany, so to include Vogese Mountains and Bohemian Forest
ger_buf50 <- gBuffer(ger_proj,byid=T,width=50000)
ger_buf50 <- spTransform(ger_buf50,CRS(as.character(crs(ger_shp))))
ger_buf50_clip <- gIntersection(ger_buf50,eu_shp,byid=T,drop_lower_td = TRUE)


crs(aalba_transf) <- crs(ger_buf50_clip)
crs(fsylva_transf) <- crs(ger_buf50_clip)
crs(mix_transf) <- crs(ger_buf50_clip)
aalba_ger <- aalba_transf[ger_buf50_clip,]
fsylva_ger <- fsylva_transf[ger_buf50_clip,]
mix_ger <- mix_transf[ger_buf50_clip,]



# Select pseudo-absences for Germany --------------------------------------

# create circles-model (buffer) around presences
x_a_ger <- circles(aalba_ger, d=10000, lonlat=TRUE)
x_f_ger <- circles(fsylva_ger, d=10000, lonlat=TRUE)
x_m_ger <- circles(mix_ger, d=10000, lonlat=TRUE)

pol_a_ger <- polygons(x_a_ger)
pol_f_ger <- polygons(x_f_ger)
pol_m_ger <- polygons(x_m_ger)

pol_a_ger_clip <- gIntersection(pol_a_ger,ger_buf50_clip,byid=T,drop_lower_td = TRUE)
pol_f_ger_clip <- gIntersection(pol_f_ger,ger_buf50_clip,byid=T,drop_lower_td = TRUE)
pol_m_ger_clip <- gIntersection(pol_m_ger,ger_buf50_clip,byid=T,drop_lower_td = TRUE)

# sample randomly from all circles
set.seed(1992)
samp_a_ger <- spsample(pol_a_ger_clip, 7000, type='random', iter=10)
samp_f_ger <- spsample(pol_f_ger_clip, 25000, type='random', iter=10)
samp_m_ger <- spsample(pol_m_ger_clip, 5000, type='random', iter=10)
mask <- raster(paste0(getwd),"/eu_AnnMeanT.tif")

# get unique cells
cells_a_ger <- cellFromXY(mask, samp_a_ger)
cells_a_ger <- unique(cells_a_ger)
xy_a_ger <- xyFromCell(mask, cells_a_ger)

cells_f_ger <- cellFromXY(mask, samp_f_ger)
cells_f_ger <- unique(cells_f_ger)
xy_f_ger <- xyFromCell(mask, cells_f_ger)

cells_m_ger <- cellFromXY(mask, samp_m_ger)
cells_m_ger <- unique(cells_m_ger)
xy_m_ger <- xyFromCell(mask, cells_m_ger)

# check results on map
data("wrld_simpl")
par(mfrow=c(2,2))
plot(wrld_simpl, main="silver fir",ylim=c(47,56), xlim=c(9,15), axes=T, col="light yellow")
points(xy_a_ger, cex=0.2, pch=20, col='blue')
plot(aalba_ger,add=T,pch=20,col="red",cex=0.2)
plot(wrld_simpl, main="beech",ylim=c(47,56), xlim=c(9,15), axes=T, col="light yellow")
points(xy_f_ger, cex=0.2, pch=20, col='blue')
plot(fsylva_ger,add=T,pch=20,col="red",cex=0.2)
plot(wrld_simpl, main="mixed",ylim=c(47,56), xlim=c(9,15), axes=T, col="light yellow")
points(xy_m_ger, cex=0.2, pch=20, col='blue')
plot(mix_ger,add=T,pch=20,col="red",cex=0.2)
plot.new()
par(xpd=TRUE)
legend("center",legend=c("presence","(pseudo-) absence"),col=c("red","blue"),pch=20,cex=1)
par(xpd=FALSE)

# (pseudo-) absences
euf_aalba_ger_abs <- cbind(PRESENCE=0,xy_a_ger)
colnames(euf_aalba_ger_abs) <- c("PRESENCE","X","Y")
euf_fsylva_ger_abs <- cbind(PRESENCE=0,xy_f_ger)
colnames(euf_fsylva_ger_abs) <- c("PRESENCE","X","Y")
euf_mix_ger_abs <- cbind(PRESENCE=0,xy_m_ger)
colnames(euf_mix_ger_abs) <- c("PRESENCE","X","Y")
# presences
euf_aalba_ger_pres <- cbind(PRESENCE=1,as.data.frame(aalba_ger))
colnames(euf_aalba_ger_pres) <- c("PRESENCE","X","Y")
euf_fsylva_ger_pres <- cbind(PRESENCE=1,as.data.frame(fsylva_ger))
colnames(euf_fsylva_ger_pres) <- c("PRESENCE","X","Y")
euf_mix_ger_pres <- cbind(PRESENCE=1,as.data.frame(mix_ger))
colnames(euf_mix_ger_pres) <- c("PRESENCE","X","Y")

# bind presences with pseudo absences
euf_aalba_ger <- rbind(euf_aalba_ger_pres,euf_aalba_ger_abs)
euf_fsylva_ger <- rbind(euf_fsylva_ger_pres,euf_fsylva_ger_abs)
euf_mix_ger <- rbind(euf_mix_ger_pres,euf_mix_ger_abs)


# Save results ------------------------------------------------------------


save(euf_aalba_ger, file="euf_aalba_ger")
save(euf_fsylva_ger, file="euf_fsylva_ger")
save(euf_mix_ger, file="euf_mix_ger")

