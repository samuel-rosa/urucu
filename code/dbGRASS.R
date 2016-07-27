# title: "Create GRASS GIS database"
# author: "Alessandro Samuel Rosa"

rm(list = ls())

# CREATE GRASS GIS DATABASE ###################################################################################

# Set GRASS GIS DATABASE
dbGRASS <- "data/GRASS"
spgrass6::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS, location = "urucu", mapset = "database",
  override = TRUE, pid = Sys.getpid())

# Set PERMANENT coordinate reference system
spgrass6::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS, location = "urucu", mapset = "PERMANENT",
  override = TRUE, pid = Sys.getpid())
system("g.proj -c epsg=29190")

# Update projection information in 'database'
spgrass6::initGRASS(
  gisBase = "/usr/lib/grass64/", gisDbase = dbGRASS, location = "urucu", mapset = "database",
  override = TRUE, pid = Sys.getpid())
system("g.region -d")

# PROCESS ELEVATION DATA ######################################################################################
# Processing elevation data is necessary because the DEM used by Villela (2013) has noise resulting from the
# use of contour lines, as well as many sinks on hilltops and flat areas. All processing steps are performed
# using QGIS, GRASS and SAGA.
# Afterwards, all covariates used by Villela (2013) are recomputed using ArcGIS.

# Load raw elevation data
dir <- path.expand("~/projects/urucu/data/grid/old_dem/")
cmd <- paste("r.in.gdal input=", dir, "mdehc5x5.tif output=mdehc5x5", sep = "")
system(cmd)

# Remove stair-like artifacts using thee window sizes
size <- c(3, 7, 15)
cmd <- paste("r.neighbors --o in=mdehc5x5 out=mdehc5x5.smooth.", size, " method=average size=", size, sep = "")
lapply(cmd, system)

# Derive a profile curvature raster surface from the original DEM and each of the processed DEMs
# Open data in QGIS to check which was the most efficient filter
cmd <- paste("r.slope.aspect elevation=mdehc5x5.smooth.", size, " pcurv=pcurv.", size, sep = "")
parallel::mclapply(cmd, system, mc.cores = 2)
cmd <- paste("r.slope.aspect elevation=mdehc5x5 pcurv=pcurv", sep = "")
system(cmd)

# Use the DEM filtered using a window size of 3 cells

# Fill sinks repeatedly to guarantee the quality of results
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.smooth.", size[1], " elevation=mdehc5x5.fill.", size[1], 
  " direction=mdehc5x5.dir.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, system, mc.cores = 2)
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.fill.", size[1], " elevation=mdehc5x5.filled.", size[1], 
  " direction=mdehc5x5.direc.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, system, mc.cores = 2)
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.filled.", size[1], " elevation=mdehc5x5.filled2.", size[1], 
  " direction=mdehc5x5.direc2.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, system, mc.cores = 2)
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.filled2.", size[1], " elevation=mdehc5x5.filled3.", size[1], 
  " direction=mdehc5x5.direc3.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, system, mc.cores = 2)
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.filled3.", size[1], " elevation=mdehc5x5.filled4.", size[1], 
  " direction=mdehc5x5.direc4.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, system, mc.cores = 2)
cmd <- paste(
  "r.fill.dir --o input=mdehc5x5.filled4.", size[1], " elevation=mdehc5x5.filled5.", size[1], 
  " direction=mdehc5x5.direc5.", size[1], " type=answers", sep = "")
parallel::mclapply(cmd, system, mc.cores = 2)

# Export saga grid
cmd <- paste(
  "r.out.gdal input=mdehc5x5.filled5.", size[1], " output=/home/lgcs-mds/tmp/tmp_mde_filled.sgrd", sep = "")
system(cmd)

# Look for flat areas
# Open in QGIS to check results
cmd <- paste("saga_cmd ta_preprocessor 0 -DEM /home/lgcs-mds/tmp/tmp_mde_filled.sgrd",
             "-FLATS /home/lgcs-mds/tmp/flats.sgrd")
system(cmd)

# Derive a profile curvature raster surface from the selected dem
# Open data in QGIS to check the efficiency of the algorithm. Compare with location of 
# flat surfaces.
cmd <- paste("r.slope.aspect elevation=mdehc5x5.filled5.", size[1], " pcurv=pcurv.filled5.", size[1], sep = "")
parallel::mclapply(cmd, system, mc.cores = 2)

# Remove flat surfaces
RSAGA::rsaga.get.usage(lib = "ta_morphometry", module = 4)
cmd <- paste(
  "saga_cmd ta_preprocessor 4 -ELEV /home/lgcs-mds/tmp/tmp_mde_filled.sgrd ", 
  "-FILLED /home/lgcs-mds/tmp/tmp_mde_sinkless.sgrd ",
  "-WSHED /home/lgcs-mds/tmp/tmp_mde_watershed.sgrd ", sep = "")
system(cmd)
cmd <- paste("saga_cmd ta_preprocessor 0 -DEM /home/lgcs-mds/tmp/tmp_mde_sinkless.sgrd ",
             "-FLATS /home/lgcs-mds/tmp/flats.sgrd", sep = "")
system(cmd)
# It now seems that all flats were removed.

# Compute the difference between original and processed DEM
cmd <- paste("r.mapcalc 'mdehc5x5.diff.", size[1], "=mdehc5x5.filled5.", size[1], "-mdehc5x5'", sep = "")
system(cmd)

# Round values to 2 digits
cmd <- paste(
  "r.mapcalc 'elevation=(double(round(mdehc5x5.filled5.", size[1], "*100))/100)'", sep = "")
parallel::mclapply(cmd, system, mc.cores = 2)
cmd <- paste("r.out.gdal input=elevation output=/home/lgcs-mds/projects/urucu/data/grid/elevation.tif")
system(cmd)

# LOAD RASTER SURFACES INTO DATABASE ##########################################################################
# All primary covariates were derived using ArcGIS after sinks and noise had been removed.
# Set parameters of the GRASS GIS region based on the elevation raster surface.
dir <- "/home/lgcs-mds/projects/urucu/data/grid/"
system("d.mon x0")

# Elevation
system(paste("r.in.gdal input=", dir, "elevation.tif output=elevation", sep = ""))
system("r.info elevation")
system("g.region rast=elevation")
spgrass6::gmeta()
system("d.rast.leg elevation")
system("d.histogram elevation")

# Upslope flow path length
system(paste("r.in.gdal input=", dir, "flow_up.tif output=flow_up", sep = ""))
system("r.info flow_up")
system("d.rast.leg flow_up")
system("d.histogram flow_up")
# This covariate has a very assymetric distribution. It would be wise to reduce the assymetry to yielding
# stronger correlations with the target soil variable. One solution is to use a logarithmic transform (but we
# could also create categories). Round values to the second decimal place.
cmd <- c("r.mapcalc 'ln_flow_up=log(flow_up + 1)'")
system(cmd)
cmd <- paste("r.mapcalc 'ln_flow_up = double(round(ln_flow_up * 100)) * 0.01'")
system(cmd)
system("r.info ln_flow_up")
system("d.rast.leg ln_flow_up")
system("d.histogram ln_flow_up")

# Curvature
# This covariate has a very peaky distribution, with most values around zero. A solution might be to create
# categories.
# Round values to the second decimal place.
system(paste("r.in.gdal input=", dir, "curvature.tif output=curvature", sep = ""))
cmd <- paste("r.mapcalc 'curvature = double(round(curvature * 100)) * 0.01'")
system(cmd)
system("r.info curvature")
system("d.rast.leg curvature")
system("d.histogram curvature")

# Profile curvature
# This covariate has a very peaky distribution, with most values around zero. A solution might be to create
# categories.
# Round values to the second decimal place.
system(paste("r.in.gdal input=", dir, "prof_curv.tif output=prof_curv", sep = ""))
cmd <- paste("r.mapcalc 'prof_curv = double(round(prof_curv * 100)) * 0.01'")
system(cmd)
system("r.info prof_curv")
system("d.rast.leg prof_curv")
system("d.histogram prof_curv")

# Plan curvature
# This covariate has a very peaky distribution, with most values around zero. A solution might be to create
# categories.
# Round values to the second decimal place.
system(paste("r.in.gdal input=", dir, "plan_curv.tif output=plan_curv", sep = ""))
cmd <- paste("r.mapcalc 'plan_curv = double(round(plan_curv * 100)) * 0.01'")
system(cmd)
system("r.info plan_curv")
system("d.rast.leg plan_curv")
system("d.histogram plan_curv")

# Slope
# The empirical distribution is very assymetric, with most values close to zero. Larger values are found close
# to the drainage network.
# Correct flat surfaces adding 0.1º to the data. Then round values to the second decimal place.
system(paste("r.in.gdal input=", dir, "slope.tif output=slope", sep = ""))
cmd <- paste("r.mapcalc 'slope =slope + 0.1'")
system(cmd)
cmd <- paste("r.mapcalc 'slope = double(round(slope * 100)) * 0.01'")
system(cmd)
system("r.info slope")
system("d.rast.leg slope")
system("d.histogram slope")

# Downslope flow path length
system(paste("r.in.gdal input=", dir, "flow_down.tif output=flow_down", sep = ""))
system("r.info flow_down")
system("d.rast.leg flow_down")
system("d.histogram flow_down")

# Flow accumulation
# Like flow direction, this covariate has erroneous values outside the study area. This will likely not be a
# problem because we use a masking variable.
# Compute catchment area using the cell size (25 m²)
system(paste("r.in.gdal input=", dir, "flow_accum.tif output=flow_accum", sep = ""))
cmd <- paste("r.mapcalc 'flow_accum = (flow_accum + 1) * 25'")
system(cmd)
system("r.info flow_accum")
system("d.rast.leg flow_accum")
system("d.histogram flow_accum")
# This covariate has an assymetric distribution. It could be helpful to transform it using a logarithmic
# function. Then round values to the second decimal place.
cmd <- paste("r.mapcalc 'ln_flow_accum = log(flow_accum)'")
system(cmd)
cmd <- paste("r.mapcalc 'ln_flow_accum = double(round(ln_flow_accum * 100)) * 0.01'")
system(cmd)
system("r.info ln_flow_accum")
system("d.rast.leg ln_flow_accum")
system("d.histogram ln_flow_accum")

# Topographic wetness index
# Reference: Wilson and Gallant (2000)
# Round values to the second decimal place.
cmd <- paste("r.mapcalc 'twi = log((flow_accum / 5) / tan(slope))'", sep = "")
system(cmd)
cmd <- paste("r.mapcalc 'twi = double(round(twi * 100)) * 0.01'")
system(cmd)
system("r.info twi")
system("d.rast.leg twi")
system("d.histogram twi")

system("d.mon stop=x0")
rm(dir)

# LOAD VECTOR DATA INTO DATABASE ##############################################################################

dir <- "/home/lgcs-mds/projects/urucu/data/vector/"

# Target soil map

system(paste("v.in.ogr dsn=", dir, "target_soil_map.shp output=target_soil_map", sep = ""))
system("v.info map=target_soil_map")

# Outer limit
system(paste("v.in.ogr dsn=", dir, "outer_area.shp output=outer_limit", sep = ""))
system("v.info map=outer_limit")

# Accesible area
system(paste("v.in.ogr dsn=", dir, "access_area.shp output=access_limit", sep = ""))
system("v.info map=access_limit")

# Non-accessible area
system(paste("v.in.ogr dsn=", dir, "non_access_area.shp output=non_access_limit", sep = ""))
system("v.info map=non_access_limit")

# Transform vector to raster to create a MASK.
system("v.to.rast input=non_access_limit output=non_access_limit use=val")
system("r.mask --o input=non_access_limit")

rm(dir)
