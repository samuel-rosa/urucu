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
