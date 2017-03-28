
# Load data and libraries -------------------------------------------------

# ATTENTION: The next line should set the working directory to the source file location. This only works in RStudio. If you use RGUI or experience problems please replace with full path!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
library(sdm)
library(dismo)
library(rJava)
library(maptools)
library(raster)

# Occurrence data with extracted environmental values
load(paste0(getwd(),"/euf_aalba_ger_env"))
load(paste0(getwd(),"/euf_fsylva_ger_env"))
load(paste0(getwd(),"/euf_mix_ger_env"))

# Predictor data
files_eu <- list.files(getwd(),pattern = ".tif",full.names=T)
pred_eu_sel <- stack(files_eu)

# Preparation -------------------------------------------------------------

# convert to SpatialPointsDataFrame
train_ger_aa <- euf_aalba_ger_env
coordinates(train_ger_aa) <- ~X+Y
crs(train_ger_aa) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
train_ger_fs <- euf_fsylva_ger_env
coordinates(train_ger_fs) <- ~X+Y
crs(train_ger_fs) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
train_ger_mix <- euf_mix_ger_env
coordinates(train_ger_mix) <- ~X+Y
crs(train_ger_mix) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

# build formula
fmla <- as.formula(paste("PRESENCE ~ ",paste(names(pred_eu_sel),collapse="+")))


# Build model(s) -------------------------------------------------------------

# create sdmData objects
d_ger_aa <- sdmData(formula=fmla, train=train_ger_aa)
d_ger_fs <- sdmData(formula=fmla, train=train_ger_fs)
d_ger_mix <- sdmData(formula=fmla, train=train_ger_mix)

# store locations of presence/absence within sdmData
p_ger_aa <- d_ger_aa@species$PRESENCE@presence
a_ger_aa <- d_ger_aa@species$PRESENCE@absence
p_ger_fs <- d_ger_fs@species$PRESENCE@presence
a_ger_fs <- d_ger_fs@species$PRESENCE@absence
p_ger_mix <- d_ger_mix@species$PRESENCE@presence
a_ger_mix <- d_ger_mix@species$PRESENCE@absence



### GENERAL SETTINGS for all models: 
### prevalence = 0.5 (500 pres/1000 abs)
### 20% test data
### 100 iterations

# NOTE: Running all models together may take several hours! Change the 'methods' argument to pick only the desired model(s). MaxEnt is known to cause problems from time to time, as it runs externally with Java. Try running a MaxEnt model separately before doing an all-model-run!
# Abies alba
aa_maxent_test <- sdm(formula = fmla,data=d_ger_aa[c(sample(p_ger_aa,50,replace = F),sample(a_ger_aa,100,replace = F)),drop=F],methods='maxent',test.percent=20, n=10)

aa_ger_sdm_all <- sdm(formula = fmla,data=d_ger_aa[c(sample(p_ger_aa,500,replace = F),sample(a_ger_aa,1000,replace = F)),drop=F],methods=c('glm','gam','brt','svm','rf','maxent'),test.percent=20, n=100)

# Fagus sylvatica
fs_maxent_test <- sdm(formula = fmla,data=d_ger_fs[c(sample(p_ger_fs,50,replace = F),sample(a_ger_fs,100,replace = F)),drop=F],methods='maxent',test.percent=20, n=10)

fs_ger_sdm_all <- sdm(formula = fmla,data=d_ger_fs[c(sample(p_ger_fs,500,replace = F),sample(a_ger_fs,1000,replace = F)),drop=F],methods=c('glm','gam','brt','svm','rf','maxent'),test.percent=20, n=100)

# Mixed forest
mix_maxent_test <- sdm(formula = fmla,data=d_ger_mix[c(sample(p_ger_mix,50,replace = F),sample(a_ger_mix,100,replace = F)),drop=F],methods='maxent',test.percent=20, n=10)

mix_ger_sdm_all <- sdm(formula = fmla,data=d_ger_mix[c(sample(p_ger_mix,500,replace = F),sample(a_ger_mix,1000,replace = F)),drop=F],methods=c('glm','gam','brt','svm','rf','maxent'),test.percent=20, n=100)

# Model overview (run as is)
aa_ger_sdm_all
fs_ger_sdm_all
mix_ger_sdm_all

# Save results ------------------------------------------------------------

save(aa_ger_sdm_all,file="mix_ger_sdm_all")
save(fs_ger_sdm_all,file="mix_ger_sdm_all")
save(mix_ger_sdm_all,file="mix_ger_sdm_all")


