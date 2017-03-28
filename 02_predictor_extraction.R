# Load data and libraries -------------------------------------------------

# ATTENTION: The next line should set the working directory to the source file location. This only works in RStudio. If you use RGUI or experience problems please replace with full path!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
library(raster)
library(maptools)
library(rgdal)
library(gdalUtils)
library(fmsb)
library(usdm)

# EU-Forest occurrences
load(paste0(getwd(),"/euf_aalba_ger"))
load(paste0(getwd(),"/euf_fsylva_ger"))
load(paste0(getwd(),"/euf_mix_ger"))



# Stack predictor files ---------------------------------------------------

### NOTE: preprocessing of environmental data excluded from code. Load ready results

files_eu <- list.files(getwd(),pattern = ".tif",full.names=T)
pred_eu_sel <- stack(files_eu)


# Predictor data statistics -----------------------------------------------

### VIF
vif_eu_sel <- vif(pred_eu_sel)

### Moran's I
# CAUTION: heavy computation
mori_eu <- matrix(data=NA, nrow=10, ncol=2)
mori_eu[,1] <- names(pred_eu_sel)
for(l in 1:10){
  mori_eu[l,2] <- Moran(pred_eu_sel[[l]])
}

### Correlation matrix
corr_eu_sel <- layerStats(pred_eu_sel,'pearson', na.rm=T)
corrmat_eu_sel <- corr_eu_sel$'pearson correlation coefficient'
library(corrplot)
corrplot(corrmat_eu_sel,method="number")


# Extract environmental values to occurrences ----------------------------------

euf_aalba_ger_xy <- cbind(euf_aalba_ger$X,euf_aalba_ger$Y)
euf_aalba_ger_env <- extract(pred_eu_sel,euf_aalba_ger_xy,df=TRUE,sp=TRUE)
euf_aalba_ger_env <- euf_aalba_ger_env[,-1] # remove automatically added ID
euf_aalba_ger_env <- cbind(euf_aalba_ger,euf_aalba_ger_env) # add coordinates again

euf_fsylva_ger_xy <- cbind(euf_fsylva_ger$X,euf_fsylva_ger$Y)
euf_fsylva_ger_env <- extract(pred_eu_sel,euf_fsylva_ger_xy,df=TRUE,sp=TRUE)
euf_fsylva_ger_env <- euf_fsylva_ger_env[,-1]
euf_fsylva_ger_env <- cbind(euf_fsylva_ger,euf_fsylva_ger_env) 

euf_mix_ger_xy <- cbind(euf_mix_ger$X,euf_mix_ger$Y)
euf_mix_ger_env <- extract(pred_eu_sel,euf_mix_ger_xy,df=TRUE,sp=TRUE)
euf_mix_ger_env <- euf_mix_ger_env[,-1]
euf_mix_ger_env <- cbind(euf_mix_ger,euf_mix_ger_env) 



# Extracted values statistics ---------------------------------------------


### Boxplots
bionames <- c("Aridity Index", "Aspect","Bulk density","Elevation","Isothermality", "Max Temperature of Warmest Month","Min Temperature of Coldest Month","pH","Precipitation of Warmest Quarter","Slope (degree)")

# Background of Abies alba
par(mfrow=c(5,2),mar=c(2,1,2,1))
for(n in 4:13){
  boxplot(euf_aalba_ger_env[,n],horizontal=TRUE, main=bionames[n-3],las=1,ylim=c(minValue(pred_eu_sel[[n-3]]),maxValue(pred_eu_sel[[n-3]])))
  #segments(minValue(pred_eu_sel[[n-2]]),1,maxValue(pred_eu_sel[[n-2]]))
  print(n)
}

# Background of Fagus sylvatica 
for(n in 4:13){
  boxplot(euf_fsylva_ger_env[,n],horizontal=TRUE, main=bionames[n-3],las=1,ylim=c(minValue(pred_eu_sel[[n-3]]),maxValue(pred_eu_sel[[n-3]])))
  #segments(minValue(pred_eu_sel[[n-2]]),1,maxValue(pred_eu_sel[[n-2]]))
  print(n)
}

# Background of Mixed forests
for(n in 4:13){
  boxplot(euf_mix_ger_env[,n],horizontal=TRUE, main=bionames[n-3],las=1,ylim=c(minValue(pred_eu_sel[[n-3]]),maxValue(pred_eu_sel[[n-3]])))
  #segments(minValue(pred_eu_sel[[n-2]]),1,maxValue(pred_eu_sel[[n-2]]))
  print(n)
}


### Scatterplot matrices
euf_aalba_ger_env_pres <- euf_aalba_ger_env[euf_aalba_ger_env$PRESENCE==1,]
euf_aalba_ger_env_pres$species <- "abies"
euf_fsylva_ger_env_pres <- euf_fsylva_ger_env[euf_fsylva_ger_env$PRESENCE==1,]
euf_fsylva_ger_env_pres$species <- "fagus"
species_ger_env <- rbind(euf_fsylva_ger_env_pres,euf_aalba_ger_env_pres)
species_ger_env$species <- as.factor(species_ger_env$species)

par(mfrow=c(1,1))
pairs(species_ger_env[,4:13],main="Species scatterplot (red= Fagus, black= Abies)",col=species_ger_env$species)



# Save results ------------------------------------------------------------

save(euf_aalba_ger_env,file="euf_aalba_ger_env")
save(euf_fsylva_ger_env,file="euf_fsylva_ger_env")
save(euf_mix_ger_env,file="euf_mix_ger_env")
