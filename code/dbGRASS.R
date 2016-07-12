# title: "Create GRASS GIS database"
# author: "Alessandro Samuel Rosa"

rm(list = ls())

# CREATE GRASS GIS DATABASE ###################################################################################

# Load elevation raster surface.
elevation <- rgdal::readGDAL("data/grid/mdehc5x5.tif")
str(elevation)

# Set GRASS GIS DATABASE
dbGRASS <- "data/GRASS"
spgrass6::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS, SG = elevation, location = "urucu", mapset = "database",
  override = TRUE, pid = Sys.getpid())

# Set PERMANENT coordinate reference system
spgrass6::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS, location = "urucu", mapset = "PERMANENT",
  override = TRUE, pid = Sys.getpid())
system("g.proj -c epsg=32720")

# Update projection information in 'database'
spgrass6::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS, location = "urucu", mapset = "database",
  override = TRUE, pid = Sys.getpid())
system("g.region -d")

# Load elevation raster surface and set it as default
spgrass6::writeRAST(elevation, "elevation")
system("g.region rast=elevation")
spgrass6::gmeta()
rm(elevation)
gc()

# LOAD RASTER SURFACES INTO DATABASE ##########################################################################
dir <- "/home/lgcs-mds/projects/urucu/data/grid/"

# Aspect
system(paste("r.in.gdal input=", dir, "aspect.tif output=aspect", sep = ""))
system("r.info map=aspect")

# Upslope flow path length
system(paste("r.in.gdal input=", dir, "comdrenmon.tif output=upslope", sep = ""))
system("r.info map=upslope")
system(paste("r.in.gdal input=", dir, "ln_comdrenmon.tif output=ln_upslope", sep = ""))
system("r.info map=ln_upslope")

# Curvature
system(paste("r.in.gdal input=", dir, "curvatura.tif output=curvature", sep = ""))
system("r.info map=curvature")

# Profile curvature
system(paste("r.in.gdal input=", dir, "curvperf.tif output=profile_curvature", sep = ""))
system("r.info map=profile_curvature")

# Plan curvature
system(paste("r.in.gdal input=", dir, "curvplan.tif output=plan_curvature", sep = ""))
system("r.info map=plan_curvature")

# Slope
system(paste("r.in.gdal input=", dir, "declgrauhc.tif output=slope", sep = ""))
system("r.info map=slope")

# Flow direction
system(paste("r.in.gdal input=", dir, "dirfluxcort2.tif output=flow_direction", sep = ""))
system("r.info map=flow_direction")

# Downslope flow path length
system(paste("r.in.gdal input=", dir, "flowlend_n.tif output=downslope", sep = ""))
system("r.info map=downslope")

# Flow accumulation
system(paste("r.in.gdal input=", dir, "fluxacc5x5.tif output=accumulation", sep = ""))
system("r.info map=accumulation")
system(paste("r.in.gdal input=", dir, "ln_flow_accumulation.tif output=ln_accumulation", sep = ""))
system("r.info map=ln_accumulation")

# Topographic wetness index
system(paste("r.in.gdal input=", dir, "itu5x5.tif output=twi", sep = ""))
system("r.info map=twi")

# Northerness
system(paste("r.in.gdal input=", dir, "northerness.tif output=northerness", sep = ""))
system("r.info map=northerness")


















# Load boundary vector to GRASS and use it to create a MASK.
boundary <- raster::shapefile("data/QGIS/boundary.shp")
spgrass6::writeVECT(boundary, "boundary")
system("v.to.rast -h")
system("v.to.rast input=boundary output=boundary use=val")
system("r.mask input=boundary")
system("d.mon start=x0")
system("d.rast past_landuse")
system("d.mon stop=x0")

# CURRENT LAND USE ############################################################################################
# Load current land use data to GRASS
current_landuse <- raster::shapefile("data/QGIS/landuse.shp")
current_landuse@data$landuse <- as.factor(current_landuse@data$landuse)
current_landuse@data$id <- as.integer(c(1, 1, 0))
spgrass6::writeVECT(current_landuse, "current_landuse", v.in.ogr_flags = "overwrite")
system("v.to.rast -h")
system("v.to.rast --o input=current_landuse output=current_landuse column=id")

system("d.mon x0")
system("d.rast current_landuse")
system("d.mon stop=x0")

# FIGURES #####################################################################################################

# Current land use
load(file = "data/R/pointData.rda")
pointData <- pointData[seq(1, nrow(pointData), 5), c("x", "y")]
map <- 
  sp::spplot(
    current_landuse, "landuse", col.regions = c("lemonchiffon", "darkolivegreen1"), col = c("lightgray"),
    scales = list(draw = TRUE), 
    panel = function(x, y, ...) {
      sp::panel.polygonsplot(x, y, ...)
      lattice::panel.points(
        x = pointData$x, y = pointData$y, pch = 21, fill = "darkolivegreen3", col.symbol = "dimgrey")})
dev.off()
png("res/fig/current-landuse.png", width = 700, height = 700, res = 150)
map
dev.off()
rm(map)

# Past land use
past_landuse <- spgrass6::readRAST("past_landuse")
past_landuse$past_landuse <- (past_landuse$past_landuse - min(past_landuse$past_landuse, na.rm = TRUE)) / 
  (max(past_landuse$past_landuse, na.rm = TRUE) - min(past_landuse$past_landuse, na.rm = TRUE))
past_landuse$past_landuse <- round(past_landuse$past_landuse, 4)
past_landuse <- as(past_landuse, "SpatialPixelsDataFrame")
map <- 
  sp::spplot(
    past_landuse, "past_landuse", col.regions = rev(terrain.colors(100)), scales = list(draw = TRUE),
    at = seq(0, 1, 0.01),
    panel = function(x, y, ...) {
      lattice::panel.levelplot(x, y, ...)
      lattice::panel.points(
        x = pointData$x, y = pointData$y, pch = 21, fill = "darkolivegreen3", col.symbol = "dimgrey")})

dev.off()
png("res/fig/past-landuse.png", width = 700, height = 700, res = 150)
map
dev.off()
rm(map)
