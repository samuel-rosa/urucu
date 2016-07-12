# title: "Soil Spatial Modelling"
# author: "Alessandro Samuel Rosa"
# 
# Spatial exploratory data analysis were performed using QGIS and SAGA GIS. The results are described in 
# \code\exploratory.Rmd. The Grass GIS database was created using \code\dbGRASS.R. All covariates were loaded
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

# Covariates
# There are p = 11 covariates:
# - elevation,  slope,  curvature, plan curvatura, profile curvature, downslope frow length, upslope flow
#   length, flow accumulation, flow direction, aspect and topographic wetness index.

# Build Grass GIS database
