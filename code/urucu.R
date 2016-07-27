# title: "Soil Spatial Modelling"
# author: "Alessandro Samuel Rosa"
# 
# Spatial exploratory data analysis were performed using QGIS and SAGA GIS. The results are described in 
# \code\exploratory.Rmd. The GRASS GIS database was created using \code\dbGRASS.R. All covariates were loaded
# to the database, as well as all existing vector (polygon) data.

spgrass6::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = "data/GRASS", location = "urucu", mapset = "database",
  override = TRUE, pid = Sys.getpid())

# Prepare calibration points ##################################################################################

# Load necessary data
target_soil_map <- spgrass6::readVECT("target_soil_map")
access_limit <- spgrass6::readVECT("access_limit")
cal_profiles <- raster::shapefile("data/vector/Perfis.shp")
cal_boreholes <- raster::shapefile("data/vector/Tradagens.shp")
val_boreholes <- raster::shapefile("data/vector/Tradagens_Validacao.shp")

# Merge field observations, removing those that fall outlide the accessible area.
# Here we create the field calibration dataset.
out <- which(is.na(sp::over(cal_boreholes, access_limit)$OBJECTID))
cal_field <- 
  as.data.frame(
    rbind(sp::coordinates(cal_profiles)[, ], sp::coordinates(cal_boreholes)[-out, ], 
          sp::coordinates(val_boreholes)[, ]))
head(cal_field)
sp::coordinates(cal_field) <- ~ coords.x1 + coords.x2
sp::proj4string(cal_field) <- sp::proj4string(cal_profiles)
colnames(cal_field@data) <- "UM"
cal_field$UM <- sp::over(cal_field, target_soil_map)$UM
str(cal_field)
rm(cal_boreholes, cal_profiles, val_boreholes, out)
gc()

# Select a stratified random sample (n ~ 2000).
# Here we create the probabilistic calibration dataset.
set.seed(2001)
cal_random <- sp::spsample(target_soil_map, n = 2000, type = "stratified")
tmp <- sp::over(cal_random, target_soil_map)$UM
cal_random <- sp::SpatialPointsDataFrame(cal_random, data = data.frame(UM = tmp))
str(cal_random)

# Load expert points.
# Here we create the expect calibration dataset.
cal_expert <- raster::shapefile("data/vector/Trein_Classes.shp")
cal_expert$UM <- as.factor(cal_expert$MDS)
levels(cal_expert$UM) <- levels(cal_random$UM)[c(1, 2, 2, 3, 4)]
cal_expert@data <- data.frame(UM = cal_expert$UM)
str(cal_expert)
rm(tmp)

# Prepare calibration data ####################################################################################

# Load calibration points to GRASS DB
spgrass6::writeVECT(SDF = cal_field, vname = "cal_field")
spgrass6::writeVECT(SDF = cal_random, vname = "cal_random")
spgrass6::writeVECT(SDF = cal_expert, vname = "cal_expert")

# Setup database of calibration points
maps <- c("elevation", "flow_up", "ln_flow_up", "curvature", "prof_curv", "plan_curv", "slope", "flow_down",
          "flow_accum", "ln_flow_accum", "twi")
cols <- paste(maps, "DOUBLE PRECISION")
cols <- paste(cols, collapse = ",")
cols_samp <- substring(maps, first = 1, last = 10)

# Field calibration
pts <- "cal_field"
system(paste("v.info -c ", pts))
cmd <- paste("v.db.addcol map=", pts, " columns='", cols, "'", sep = "")
system(cmd)
system(paste("v.info -c ", pts))
cmd <- paste("v.what.rast vect=", pts, " raster=", maps, " column=", cols_samp, sep = "")
lapply(cmd, system)
system(paste("v.info -c ", pts))
rm(cmd, pts)

# Random calibration
pts <- "cal_random"
system(paste("v.info -c ", pts))
cmd <- paste("v.db.addcol map=", pts, " columns='", cols, "'", sep = "")
system(cmd)
system(paste("v.info -c ", pts))
cmd <- paste("v.what.rast vect=", pts, " raster=", maps, " column=", cols_samp, sep = "")
lapply(cmd, system)
system(paste("v.info -c ", pts))
rm(cmd, pts)

# Expert calibration
pts <- "cal_expert"
system(paste("v.info -c ", pts))
cmd <- paste("v.db.addcol map=", pts, " columns='", cols, "'", sep = "")
system(cmd)
system(paste("v.info -c ", pts))
cmd <- paste("v.what.rast vect=", pts, " raster=", maps, " column=", cols_samp, sep = "")
lapply(cmd, system)
system(paste("v.info -c ", pts))
rm(cmd, pts)

rm(cols, maps)

# Calibrate soil prediction models ############################################################################

form <- formula(paste("UM ~ ", paste(cols_samp[c(1, 2, 4, 5, 6, 7, 8, 9, 11)], collapse = " + ")))
form <- formula(paste("UM ~ ", paste(cols_samp[c(1, 4, 6, 7, 8, 10)], collapse = " + "))) # log and cor

# Field calibration points
cal_field <- spgrass6::readVECT("cal_field")
cal_field$UM <- as.factor(cal_field$UM)

str(cal_field)
head(cal_field@data)
nrow(na.omit(cal_field@data)) - nrow(cal_field@data)
which(sapply(lapply(cal_field@data, is.na), any))

fit_field <- MASS::lda(form, data = cal_field@data, na.action = na.omit)
fit_field_fitted <- predict(fit_field, newdata = na.omit(cal_field@data), prior = fit_field$prior)
sum(diag(table(data.frame(obs = na.omit(cal_field@data)$UM, fit = fit_field_fitted$class)))) / 
  nrow(na.omit(cal_field@data))

# Random calibration points
cal_random <- spgrass6::readVECT("cal_random")
cal_random$UM <- as.factor(cal_random$UM)

str(cal_random)
head(cal_random@data)
nrow(na.omit(cal_random@data)) - nrow(cal_random@data)
which(sapply(lapply(cal_random@data, is.na), any))

fit_random <- MASS::lda(form, data = cal_random@data, na.action = na.omit)
fit_random_fitted <- predict(fit_random, newdata = na.omit(cal_random@data), prior = fit_random$prior)
sum(diag(table(data.frame(obs = na.omit(cal_random@data)$UM, fit = fit_random_fitted$class)))) / 
  length(cal_random$UM)


tmp <- subselect::ldaHmat(form, na.omit(cal_field@data))
subselect::improve(tmp$mat, 1, 5)
round(cor(na.omit(cal_field@data[, -c(1:2)])), 2)





# Expert calibration points
cal_expert <- spgrass6::readVECT("cal_expert")
cal_expert$UM <- as.factor(cal_expert$UM)

str(cal_expert)
head(cal_expert@data)
nrow(na.omit(cal_expert@data)) - nrow(cal_expert@data)
which(sapply(lapply(cal_expert@data, is.na), any))

fit_expert <- MASS::lda(form, data = cal_expert@data, na.action = na.omit)
fit_expert_fitted <- predict(fit_expert, newdata = na.omit(cal_expert@data), prior = fit_expert$prior)
sum(diag(table(data.frame(obs = na.omit(cal_expert@data)$UM, fit = fit_expert_fitted$class)))) / 
  nrow(na.omit(cal_expert@data))
