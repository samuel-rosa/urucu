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
system("r.info map=elevation")
system("g.region rast=elevation")
spgrass6::gmeta()
system("d.rast.leg elevation")
system("d.histogram elevation")
system("r.mask elevation")

# Aspect
# ArcGIS help: Cells in the input raster that are flat—with zero slope—are assigned an aspect of -1.
# Given the conditions of the study area, it seems reasonable to set -1 as 0º. But let us first decide if 
# aspect will be used as a covariate or will be dropped from the dataset.
system(paste("r.in.gdal input=", dir, "aspect.tif output=aspect", sep = ""))
system("r.info map=aspect")
system("d.rast.leg aspect")
system("d.histogram aspect")

# Upslope flow path length
system(paste("r.in.gdal input=", dir, "flow_up.tif output=flow_up", sep = ""))
system("r.info flow_up")
system("d.rast.leg flow_up")
system("d.histogram flow_up")
# This covariate has a very assymetric distribution. It would be wise to reduce the assymetry to yielding
# stronger correlations with the target soil variable. One solution is to use a logarithmic transform (but we
# could also create categories).
cmd <- c("r.mapcalc 'ln_flow_up=log(flow_up + 1)'")
system(cmd)
system("r.info ln_flow_up")
system("d.rast.leg ln_flow_up")
system("d.histogram ln_flow_up")

# Curvature
# This covariate has a very peaky distribution, with most values around zero. Besides, it may be necessary to
# round values to the second decimal place. A solution might be to create categories.
system(paste("r.in.gdal input=", dir, "curvature.tif output=curvature", sep = ""))
system("r.info curvature")
system("d.rast.leg curvature")
system("d.histogram curvature")

# Profile curvature
# This covariate has a very peaky distribution, with most values around zero. Besides, it may be necessary to
# round values to the second decimal place. A solution might be to create categories.
system(paste("r.in.gdal input=", dir, "prof_curv.tif output=prof_curv", sep = ""))
system("r.info prof_curv")
system("d.rast.leg prof_curv")
system("d.histogram prof_curv")

# Plan curvature
# This covariate has a very peaky distribution, with most values around zero. Besides, it may be necessary to
# round values to the second decimal place. A solution might be to create categories.
system(paste("r.in.gdal input=", dir, "plan_curv.tif output=plan_curv", sep = ""))
system("r.info plan_curv")
system("d.rast.leg plan_curv")
system("d.histogram plan_curv")

# Slope
# It may be necessary to round values to the second decimal place.
# The empirical distribution is very assymetric, with most values close to zero. Larger values are found close
# to the drainage network.
system(paste("r.in.gdal input=", dir, "slope.tif output=slope", sep = ""))
system("r.info slope")
system("d.rast.leg slope")
system("d.histogram slope")

# Flow direction
# This covariate has erroneous values outside the study area. This will likely not be a problem because we use 
# a masking variable.
system(paste("r.in.gdal input=", dir, "flow_dir.tif output=flow_dir", sep = ""))
system("r.info flow_dir")
system("d.rast.leg flow_dir")
system("d.histogram flow_dir")

# Downslope flow path length
system(paste("r.in.gdal input=", dir, "flow_down.tif output=flow_down", sep = ""))
system("r.info flow_down")
system("d.rast.leg flow_down")
system("d.histogram flow_down")

# Flow accumulation
# Like flow direction, this covariate has erroneous values outside the study area. This will likely not be a
# problem because we use a masking variable.
system(paste("r.in.gdal input=", dir, "flow_accum.tif output=flow_accum", sep = ""))
system("r.info flow_accum")
system("d.rast.leg flow_accum")
system("d.histogram flow_accum")
# This covariate has an assymetric distribution. It would could be helpful to transform it using a logarithmic
# function.
cmd <- paste("r.mapcalc 'ln_flow_accum = log(flow_accum + 1)'")
system(cmd)
system("r.info ln_flow_accum")
system("d.rast.leg ln_flow_accum")
system("d.histogram ln_flow_accum")







# Topographic wetness index
system(paste("r.in.gdal input=", dir, "itu5x5.tif output=twi", sep = ""))
system("r.info map=twi")

# Northerness
system(paste("r.in.gdal input=", dir, "northerness.tif output=northerness", sep = ""))
system("r.info map=northerness")

# Check raster surfaces
system("d.mon start=x0")
system("d.rast twi")


system("d.mon stop=x0")
rm(dir)

# LOAD VECTOR DATA INTO DATABASE ##############################################################################

dir <- "/home/lgcs-mds/projects/urucu/data/vector/"

# Target soil map
system(paste("v.in.ogr dsn=", dir, "MCS_Simplif_Poligon_correct.shp output=target_soil_map", sep = ""))
system("v.info map=target_soil_map")

# Outer limit
system(paste("v.in.ogr dsn=", dir, "Lim50000Semidet.shp output=outer_limit", sep = ""))
system("v.info map=outer_limit")

# Accesible area
system(paste("v.in.ogr dsn=", dir, "Lim10000Detal.shp output=access_limit", sep = ""))
system("v.info map=access_limit")

# Non-accessible area
# system(paste("v.in.ogr dsn=", dir, "Lim25000Semidet.shp output=non_access_limit", sep = ""))
# system("v.info map=non_access_limit")

# # Load boundary vector to GRASS and use it to create a MASK.
# boundary <- raster::shapefile("data/QGIS/boundary.shp")
# spgrass6::writeVECT(boundary, "boundary")
# system("v.to.rast -h")
# system("v.to.rast input=boundary output=boundary use=val")
# system("r.mask input=boundary")

rm(dir)

