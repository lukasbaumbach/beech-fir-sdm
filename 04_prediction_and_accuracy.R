
# load libraries, data and results --------------------------------------

# ATTENTION: The next line should set the working directory to the source file location. This only works in RStudio. If you use RGUI or experience problems please replace with full path!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
library(PresenceAbsence)
library(sdm)
library(raster)

# occurrences
load("D:/Masterarbeit/Data/r/occurrences/euf_aalba_ger")
load("D:/Masterarbeit/Data/r/occurrences/euf_fsylva_ger")
load("D:/Masterarbeit/Data/r/occurrences/euf_mix_ger")

# models
load(paste0(getwd(),"/aa_ger_all"))
load(paste0(getwd(),"/fs_ger_all"))
load(paste0(getwd(),"/mix_ger_all"))


# environmental data
files_ger <- list.files(getwd(),pattern = ".tif",full.names=T)
pred_ger_sel <- stack(files_ger)

names(pred_ger_sel) <- c("eu_ai","eu_aspect_res","eu_bulk_0cm","eu_dem_res","eu_Isotherm","eu_MaxTWarmM","eu_MinTColdM","eu_ph_0cm","eu_PrecWarmQ","eu_slope_res")


# predict -----------------------------------------------------------------

# aggregate environmental data for first coarse overview
pred_ger_agg <- aggregate(pred_ger_sel,fact=10,fun=mean)

# predict
aa_ger_all_predtest <- predict(aa_ger_all,pred_ger_agg,mean=T)
fs_ger_all_predtest <- predict(fs_ger_all,pred_ger_agg,mean=T)
mix_ger_all_predtest <- predict(mix_ger_all,pred_ger_agg,mean=T)

names(aa_ger_all_predtest) <- c("GLM","GAM","BRT","SVM","RF","MaxEnt")

# plot predictions
labbrks <- seq(0,1,0.2)
plot(aa_ger_all_predtest, lab.breaks=labbrks, zlim=c(0,1)) 

plot(fs_ger_all_predtest, lab.breaks=labbrks, zlim=c(0,1)) 

plot(mix_ger_all_predtest, lab.breaks=labbrks, zlim=c(0,1)) 


# prepare data for calibration plots --------------------------------------

aalba_ger_sp <- euf_aalba_ger
coordinates(aalba_ger_sp) <- ~X+Y
fsylva_ger_sp <- euf_fsylva_ger
coordinates(fsylva_ger_sp) <- ~X+Y
mix_ger_sp <- euf_mix_ger
coordinates(mix_ger_sp) <- ~X+Y


aa_all_ext <- extract(aa_ger_all_predtest,aalba_ger_sp)
id <- matrix(data=seq(1,length(aalba_ger_sp[,1]),1), nrow=length(aalba_ger_sp[,1]), ncol=1)
aa_all_ext <- cbind(id,aalba_ger_sp$PRESENCE,aa_all_ext)
aa_all_ext <- data.frame(aa_all_ext)
aa_all_ext$V1 <- as.character(aa_all_ext$V1)
colnames(aa_all_ext) <- c("id","observed","glm","gam","brt","svm","rf","maxent")

fs_all_ext <- extract(fs_ger_all_predtest,fsylva_ger_sp)
id <- matrix(data=seq(1,length(fsylva_ger_sp[,1]),1), nrow=length(fsylva_ger_sp[,1]), ncol=1)
fs_all_ext <- cbind(id,fsylva_ger_sp$PRESENCE,fs_all_ext)
fs_all_ext <- data.frame(fs_all_ext)
fs_all_ext$V1 <- as.character(fs_all_ext$V1)
colnames(fs_all_ext) <- c("id","observed","glm","gam","brt","svm","rf","maxent")

mix_all_ext <- extract(mix_ger_all_predtest,mix_ger_sp)
id <- matrix(data=seq(1,length(mix_ger_sp[,1]),1), nrow=length(mix_ger_sp[,1]), ncol=1)
mix_all_ext <- cbind(id,mix_ger_sp$PRESENCE,mix_all_ext)
mix_all_ext <- data.frame(mix_all_ext)
mix_all_ext$V1 <- as.character(mix_all_ext$V1)
colnames(mix_all_ext) <- c("id","observed","glm","gam","brt","svm","rf","maxent")

# goodness of fit ---------------------------------------------------------

### calibration plots
par(mfrow=c(2,3))
for(i in 1:5){
  calibration.plot(aa_all_ext,na.rm=T,which=i)
}
par(mfrow=c(2,3))
for(i in 1:5){
  calibration.plot(fs_all_ext,na.rm=T,which=i)
}
par(mfrow=c(2,3))
for(i in 1:5){
  calibration.plot(mix_all_ext,na.rm=T,which=i)
}

### ROC plots (also calculates AUC)
par(mfrow=c(1,1))
auc.roc.plot(aa_all_ext,na.rm=T,legend.cex = 0.5,color=TRUE)
auc.roc.plot(fs_all_ext,na.rm=T,legend.cex = 0.5,color=TRUE)
auc.roc.plot(mix_all_ext,na.rm=T,legend.cex = 0.5,color=TRUE)



