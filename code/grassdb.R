# title: "Create GRASS GIS database on Linux machine"
# author: "Alessandro Samuel Rosa"

# Clean up and source user defined helper functions ###########################################################
rm(list = ls())
source("code/helper.R")

# CREATE GRASS GIS DATABASE ###################################################################################
# Set GRASS GIS DATABASE
dbGRASS <- "data/GRASS"
spgrass7::initGRASS(
  gisBase = "/usr/lib/grass70/", gisDbase = dbGRASS, location = "urucu", mapset = "database",
  override = TRUE, pid = Sys.getpid())

# Set PERMANENT coordinate reference system
spgrass7::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS, location = "urucu", mapset = "PERMANENT",
  override = TRUE, pid = Sys.getpid())
grassGis("g.proj -c epsg=29190")

# Update projection information in 'database'
spgrass7::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS, location = "urucu", mapset = "database",
  override = TRUE, pid = Sys.getpid())
grassGis("g.region -d")

# PROCESS ELEVATION DATA ######################################################################################
# Processing elevation data is necessary because the DEM used by Villela (2013) has noise resulting from the
# use of contour lines, as well as many sinks on hilltops and flat areas. All processing steps are performed
# using QGIS, GRASS and SAGA.
# Afterwards, all covariates used by Villela (2013) are recomputed using ArcGIS.

# Load raw elevation data
dir <- path.expand("~/projects/urucu/data/grid/old_dem/")
cmd <- paste("r.in.gdal input=", dir, "mdehc5x5.tif output=mdehc5x5", sep = "")
grassGis(cmd)

# Remove stair-like artifacts using thee window sizes
size <- c(3, 7, 15)
cmd <- paste("r.neighbors --o in=mdehc5x5 out=mdehc5x5.smooth.", size, " method=average size=", size, sep = "")
lapply(cmd, grassGis)

# Derive a profile curvature raster surface from the original DEM and each of the processed DEMs
# Open data in QGIS to check which was the most efficient filter
cmd <- paste("r.slope.aspect elevation=mdehc5x5.smooth.", size, " pcurv=pcurv.", size, sep = "")
parallel::mclapply(cmd, grassGis, mc.cores = 2)
cmd <- paste("r.slope.aspect elevation=mdehc5x5 pcurv=pcurv", sep = "")
grassGis(cmd)

# Use the DEM filtered using a window size of 3 cells

# Fill sinks repeatedly to guarantee the quality of results
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.smooth.", size[1], " elevation=mdehc5x5.fill.", size[1], 
  " direction=mdehc5x5.dir.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, grassGis, mc.cores = 2)
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.fill.", size[1], " elevation=mdehc5x5.filled.", size[1], 
  " direction=mdehc5x5.direc.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, grassGis, mc.cores = 2)
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.filled.", size[1], " elevation=mdehc5x5.filled2.", size[1], 
  " direction=mdehc5x5.direc2.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, grassGis, mc.cores = 2)
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.filled2.", size[1], " elevation=mdehc5x5.filled3.", size[1], 
  " direction=mdehc5x5.direc3.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, grassGis, mc.cores = 2)
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.filled3.", size[1], " elevation=mdehc5x5.filled4.", size[1], 
  " direction=mdehc5x5.direc4.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, grassGis, mc.cores = 2)
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.filled4.", size[1], " elevation=mdehc5x5.filled5.", size[1], 
  " direction=mdehc5x5.direc5.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, grassGis, mc.cores = 2)

# Export saga grid
cmd <- paste(
  "r.out.gdal input=mdehc5x5.filled5.", size[1], " output=/home/lgcs-mds/tmp/tmp_mde_filled.sgrd", sep = "")
grassGis(cmd)

# Look for flat areas
# Open in QGIS to check results
cmd <- paste("saga_cmd ta_preprocessor 0 -DEM /home/lgcs-mds/tmp/tmp_mde_filled.sgrd",
             "-FLATS /home/lgcs-mds/tmp/flats.sgrd")
grassGis(cmd)

# Derive a profile curvature raster surface from the selected dem
# Open data in QGIS to check the efficiency of the algorithm. Compare with location of 
# flat surfaces.
cmd <- paste("r.slope.aspect elevation=mdehc5x5.filled5.", size[1], " pcurv=pcurv.filled5.", size[1], sep = "")
parallel::mclapply(cmd, grassGis, mc.cores = 2)

# Remove flat surfaces
RSAGA::rsaga.get.usage(lib = "ta_morphometry", module = 4)
cmd <- paste(
  "saga_cmd ta_preprocessor 4 -ELEV /home/lgcs-mds/tmp/tmp_mde_filled.sgrd ", 
  "-FILLED /home/lgcs-mds/tmp/tmp_mde_sinkless.sgrd ",
  "-WSHED /home/lgcs-mds/tmp/tmp_mde_watershed.sgrd ", sep = "")
grassGis(cmd)
cmd <- paste("saga_cmd ta_preprocessor 0 -DEM /home/lgcs-mds/tmp/tmp_mde_sinkless.sgrd ",
             "-FLATS /home/lgcs-mds/tmp/flats.sgrd", sep = "")
grassGis(cmd)
# It now seems that all flats were removed.

# Compute the difference between original and processed DEM
cmd <- paste("r.mapcalc 'mdehc5x5.diff.", size[1], "=mdehc5x5.filled5.", size[1], "-mdehc5x5'", sep = "")
grassGis(cmd)

# Round values to 2 digits
cmd <- paste(
  "r.mapcalc 'elevation=(double(round(mdehc5x5.filled5.", size[1], "*100))/100)'", sep = "")
parallel::mclapply(cmd, grassGis, mc.cores = 2)
cmd <- paste("r.out.gdal input=elevation output=/home/lgcs-mds/projects/urucu/data/grid/elevation.tif")
grassGis(cmd)

# LOAD RASTER SURFACES INTO DATABASE ##########################################################################
# All primary covariates were derived using ArcGIS after sinks and noise had been removed.
# Set parameters of the GRASS GIS region based on the elevation raster surface.
dir <- "/home/lgcs-mds/projects/urucu/data/grid/"
grassGis("d.mon x0")

# Elevation
grassGis(paste("r.in.gdal input=", dir, "elevation.tif output=elevation", sep = ""))
grassGis("r.info elevation")
grassGis("g.region rast=elevation")
spgrass7::gmeta()
grassGis("d.rast.leg elevation")
grassGis("d.histogram elevation")

# Upslope flow path length
grassGis(paste("r.in.gdal input=", dir, "flow_up.tif output=flow_up", sep = ""))
grassGis("r.info flow_up")
grassGis("d.rast.leg flow_up")
grassGis("d.histogram flow_up")
# This covariate has a very assymetric distribution. It would be wise to reduce the assymetry to yielding
# stronger correlations with the target soil variable. One solution is to use a logarithmic transform (but we
# could also create categories). Round values to the second decimal place.
cmd <- c("r.mapcalc 'ln_flow_up=log(flow_up + 1)'")
grassGis(cmd)
cmd <- paste("r.mapcalc 'ln_flow_up = double(round(ln_flow_up * 100)) * 0.01'")
grassGis(cmd)
grassGis("r.info ln_flow_up")
grassGis("d.rast.leg ln_flow_up")
grassGis("d.histogram ln_flow_up")

# Curvature
# This covariate has a very peaky distribution, with most values around zero. A solution might be to create
# categories.
# Round values to the second decimal place.
grassGis(paste("r.in.gdal input=", dir, "curvature.tif output=curvature", sep = ""))
cmd <- paste("r.mapcalc 'curvature = double(round(curvature * 100)) * 0.01'")
grassGis(cmd)
grassGis("r.info curvature")
grassGis("d.rast.leg curvature")
grassGis("d.histogram curvature")

# Profile curvature
# This covariate has a very peaky distribution, with most values around zero. A solution might be to create
# categories.
# Round values to the second decimal place.
grassGis(paste("r.in.gdal input=", dir, "prof_curv.tif output=prof_curv", sep = ""))
cmd <- paste("r.mapcalc 'prof_curv = double(round(prof_curv * 100)) * 0.01'")
grassGis(cmd)
grassGis("r.info prof_curv")
grassGis("d.rast.leg prof_curv")
grassGis("d.histogram prof_curv")

# Plan curvature
# This covariate has a very peaky distribution, with most values around zero. A solution might be to create
# categories.
# Round values to the second decimal place.
grassGis(paste("r.in.gdal input=", dir, "plan_curv.tif output=plan_curv", sep = ""))
cmd <- paste("r.mapcalc 'plan_curv = double(round(plan_curv * 100)) * 0.01'")
grassGis(cmd)
grassGis("r.info plan_curv")
grassGis("d.rast.leg plan_curv")
grassGis("d.histogram plan_curv")

# Slope
# The empirical distribution is very assymetric, with most values close to zero. Larger values are found close
# to the drainage network.
# Correct flat surfaces adding 0.1º to the data. Then round values to the second decimal place.
grassGis(paste("r.in.gdal input=", dir, "slope.tif output=slope", sep = ""))
cmd <- paste("r.mapcalc 'slope =slope + 0.1'")
grassGis(cmd)
cmd <- paste("r.mapcalc 'slope = double(round(slope * 100)) * 0.01'")
grassGis(cmd)
grassGis("r.info slope")
grassGis("d.rast.leg slope")
grassGis("d.histogram slope")

# Downslope flow path length
grassGis(paste("r.in.gdal input=", dir, "flow_down.tif output=flow_down", sep = ""))
grassGis("r.info flow_down")
grassGis("d.rast.leg flow_down")
grassGis("d.histogram flow_down")

# Flow accumulation
# Like flow direction, this covariate has erroneous values outside the study area. This will likely not be a
# problem because we use a masking variable.
# Compute catchment area using the cell size (25 m²)
grassGis(paste("r.in.gdal input=", dir, "flow_accum.tif output=flow_accum", sep = ""))
cmd <- paste("r.mapcalc 'flow_accum = (flow_accum + 1) * 25'")
grassGis(cmd)
grassGis("r.info flow_accum")
grassGis("d.rast.leg flow_accum")
grassGis("d.histogram flow_accum")
# This covariate has an assymetric distribution. It could be helpful to transform it using a logarithmic
# function. Then round values to the second decimal place.
cmd <- paste("r.mapcalc 'ln_flow_accum = log(flow_accum)'")
grassGis(cmd)
cmd <- paste("r.mapcalc 'ln_flow_accum = double(round(ln_flow_accum * 100)) * 0.01'")
grassGis(cmd)
grassGis("r.info ln_flow_accum")
grassGis("d.rast.leg ln_flow_accum")
grassGis("d.histogram ln_flow_accum")

# Topographic wetness index
# Reference: Wilson and Gallant (2000)
# Round values to the second decimal place.
cmd <- paste("r.mapcalc 'twi = log((flow_accum / 5) / tan(slope))'", sep = "")
grassGis(cmd)
cmd <- paste("r.mapcalc 'twi = double(round(twi * 100)) * 0.01'")
grassGis(cmd)
grassGis("r.info twi")
grassGis("d.rast.leg twi")
grassGis("d.histogram twi")

grassGis("d.mon stop=x0")
rm(dir)

# LOAD VECTOR DATA INTO DATABASE ##############################################################################

dir <- paste(getwd(), "/data/vector/", sep = "")

# Target soil map
# Clean data.frame and convert target soil map to raster
tmp <- raster::shapefile(paste(dir, "target_soil_map.shp", sep = ""))
tmp@data <- data.frame(UM = tmp$UM)
levels(tmp$UM)
tmp$UM <- as.integer(tmp$UM)
spgrass7::writeVECT(tmp, "target_soil_map", v.in.ogr_flags = "overwrite")
grassGis("v.info -c target_soil_map")
cmd <- paste("v.to.rast --o input=target_soil_map output=target_soil_map column=UM")
grassGis(cmd)
grassGis("d.mon x0")
grassGis("d.rast.leg target_soil_map")

# Outer limit
grassGis(paste("v.in.ogr dsn=", dir, "outer_area.shp output=outer_limit", sep = ""))
grassGis("v.info map=outer_limit")

# Accesible area
grassGis(paste("v.in.ogr dsn=", dir, "access_area.shp output=access_limit", sep = ""))
grassGis("v.info map=access_limit")
grassGis("v.to.rast input=access_limit output=access_limit use=val")

# Non-accessible area
grassGis(paste("v.in.ogr dsn=", dir, "non_access_area.shp output=non_access_limit", sep = ""))
grassGis("v.info map=non_access_limit")

# Transform vector to raster to create a MASK.
grassGis("v.to.rast input=non_access_limit output=non_access_limit use=val")
grassGis("r.mask --o input=non_access_limit")

# Map inset
# Polygon defining the small area where spatial predictions are made using all spatial predictions models to
# compare their performance visually. The area presents a high topographic variation, with all four map units
# of the target variable present.
# Produce a perfect rectangle using pedometrics::bbox2sp.
# We use the vector data to create raster data which will be used as a raster mask.
tmp <- raster::shapefile(paste(dir, "map_inset.shp", sep = ""))
tmp <- pedometrics::bbox2sp(tmp, sp = "SpatialPolygonsDataFrame")
spgrass7::writeVECT(tmp, "map_inset", v.in.ogr_flags = "overwrite")
# cmd <- paste("v.in.ogr input=", dir, "map_inset.shp output=map_inset", sep = "")
# grassGis(cmd)
grassGis("v.to.rast --overwrite input=map_inset output=map_inset use=val")

# Geology
grassGis(paste("v.in.ogr --o input=", dir, "geology.shp output=geology snap=1", sep = ""))
grassGis("v.info map=geology")

rm(dir)
