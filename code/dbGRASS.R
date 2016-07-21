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
# Processing elevation data is necessary because the existing DEM has noise resulting from the use of 
# contour lines, as well as many sinks on hilltops and flat areas. All processing steps are performed using
# QGIS, GRASS and SAGA.

# Load raw elevation data
dir <- path.expand("~/projects/urucu/data/grid/")
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
dir <- "/home/lgcs-mds/projects/urucu/data/grid/"

# Elevation
system(paste("r.in.gdal input=", dir, "mdehc5x5.tif output=elevation", sep = ""))
system("r.info map=elevation")
system("g.region rast=elevation")
spgrass6::gmeta()

# Aspect
system(paste("r.in.gdal input=", dir, "aspect.tif output=aspect", sep = ""))
system("r.info map=aspect")

# Upslope flow path length
system(paste("r.in.gdal input=", dir, "comdrenmon.tif output=upslope", sep = ""))
system("r.info map=upslope")
cmd <- paste("r.in.gdal --overwrite input=", dir, "ln_comdrenmon.tif output=ln_upslope", sep = "")
system(cmd)
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

