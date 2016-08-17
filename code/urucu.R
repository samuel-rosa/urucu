# title: "Soil Spatial Modelling"
# author: "Alessandro Samuel Rosa"
# 
# Spatial exploratory data analysis were performed using QGIS and SAGA GIS. The results are described in 
# \code\exploratory.Rmd. The GRASS GIS database was created using \code\dbGRASS.R. All covariates were loaded
# to the database, as well as all existing vector (polygon) data.

rm(list = ls())

# Start OS dependent GRASS GIS ################################################################################

if (.Platform$OS.type == "unix") {
  gisBase <- "/usr/lib/grass64/"
} else {
  gisBase <- "C:/Program Files (x86)/GRASS GIS 6.4.4"
}

grassGis <- function (cmd) {
  if (.Platform$OS.type == "unix") {
    system(cmd)
  } else {
    shell(cmd)
  }
}

spgrass6::initGRASS(
  gisBase = gisBase,
  # gisBase = "/usr/lib/grass64/", 
  gisDbase = "data/GRASS", location = "urucu", mapset = "database",
  override = TRUE, pid = Sys.getpid())
system("r.mask -o target_soil_map")

# User defined functions ######################################################################################

# Purity
purity <- 
  function (obs, fit) {
    tab <- table(data.frame(fit, obs))
    diagonal <- diag(tab)
    overall_purity <- sum(diag(tab)) / length(obs)
    map_unit_purity <- diagonal / rowSums(tab)
    class_representation <- diagonal / colSums(tab)
    res <- 
      list(overall_purity = overall_purity, 
           by_class = data.frame(mu_purity = map_unit_purity, class_rep = class_representation))
    return(res)
  }

# Prepare calibration data ####################################################################################

# Load necessary data
target_soil_map <- spgrass6::readVECT("target_soil_map")
access_limit <- spgrass6::readVECT("access_limit")
non_access_limit <- spgrass6::readVECT("non_access_limit")
cal_profiles <- raster::shapefile("data/vector/Perfis.shp")
cal_boreholes <- raster::shapefile("data/vector/Tradagens.shp")
val_boreholes <- raster::shapefile("data/vector/Tradagens_Validacao.shp")

# Merge field observations, removing those that fall outside the accessible area.
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
cal_field$UM <- as.factor(sp::over(cal_field, target_soil_map)$UM)
levels(cal_field$UM) <- 
  levels(raster::shapefile("data/vector/target_soil_map.shp", stringsAsFactors = TRUE)$UM)
str(cal_field)
rm(cal_boreholes, cal_profiles, val_boreholes, out)
gc()

# Load expert points, removing those that fall outside the non-accessible area.
# Here we create the expect calibration dataset.
cal_expert <- raster::shapefile("data/vector/Trein_Classes.shp")
out <- sp::over(cal_expert, non_access_limit)
out <- which(is.na(sp::over(cal_expert, non_access_limit)$cat))
cal_expert <- cal_expert[-out, ]
cal_expert$UM <- as.factor(cal_expert$MDS)
levels(cal_expert$UM) <- levels(cal_field$UM)[c(1, 2, 2, 3, 4)]
cal_expert@data <- data.frame(UM = cal_expert$UM)
str(cal_expert)

# Select three probability samples of sizes equal to n = 'cal_field', n = 'cal_expert' and n = 2000.
# Here we create the probabilistic calibration datasets. First we need to load the covariate data to create
# our balanced sample. Due to the large size of the study region, it is wiser to export a csv file omiting
# cells with NA values to a temporary file.
id <- c("elevation", "curvature", "plan_curv", "slope", "flow_down", "twi")
cmd <- paste("r.out.xyz input=",  paste(id, collapse = ","), " output=data/tmp/access_covars.csv", sep = "")
system(cmd)
covars <- read.table("data/tmp/access_covars.csv", sep = "|")
colnames(covars) <- c("x", "y", id)
head(covars)
# Now we produce the random field balanced samples
set.seed(2001)
cal_random_field <- BalancedSampling::cube(
  prob = rep(length(cal_field) / nrow(covars), nrow(covars)), Xbal = as.matrix(covars[, -c(1:2)]))
cal_random_field <- covars[cal_random_field, ]
sp::coordinates(cal_random_field) <- ~ x + y
sp::proj4string(cal_random_field) <- sp::proj4string(target_soil_map)
cal_random_field$UM <- as.factor(sp::over(cal_random_field[, 1:2], target_soil_map)$UM)
levels(cal_random_field$UM) <- 
  levels(raster::shapefile("data/vector/target_soil_map.shp", stringsAsFactors = TRUE)$UM)
str(cal_random_field)
# Now we produce the random expert balanced samples
set.seed(2001)
cal_random_expert <- BalancedSampling::cube(
  prob = rep(length(cal_expert) / nrow(covars), nrow(covars)), Xbal = as.matrix(covars[, -c(1:2)]))
cal_random_expert <- covars[cal_random_expert, ]
sp::coordinates(cal_random_expert) <- ~ x + y
sp::proj4string(cal_random_expert) <- sp::proj4string(target_soil_map)
cal_random_expert$UM <- as.factor(sp::over(cal_random_expert[, 1:2], target_soil_map)$UM)
levels(cal_random_expert$UM) <- 
  levels(raster::shapefile("data/vector/target_soil_map.shp", stringsAsFactors = TRUE)$UM)
str(cal_random_expert)
# Now we produce the random large balanced samples
set.seed(2001)
cal_random_large <- BalancedSampling::cube(
  prob = rep(2000 / nrow(covars), nrow(covars)), Xbal = as.matrix(covars[, -c(1:2)]))
cal_random_large <- covars[cal_random_large, ]
sp::coordinates(cal_random_large) <- ~ x + y
sp::proj4string(cal_random_large) <- sp::proj4string(target_soil_map)
cal_random_large$UM <- as.factor(sp::over(cal_random_large[, 1:2], target_soil_map)$UM)
levels(cal_random_large$UM) <- 
  levels(raster::shapefile("data/vector/target_soil_map.shp", stringsAsFactors = TRUE)$UM)
str(cal_random_large)

# Load field and expert calibration points to GRASS GIS
spgrass6::writeVECT(SDF = cal_field, vname = "cal_field", v.in.ogr_flags = c("overwrite"))
spgrass6::writeVECT(SDF = cal_expert, vname = "cal_expert", v.in.ogr_flags = c("overwrite"))

# Setup database of calibration points
system("r.mask -o non_access_limit")
cols <- paste(id, "DOUBLE PRECISION")
cols <- paste(cols, collapse = ",")
cols_samp <- substring(id, first = 1, last = 10)

# Field calibration
pts <- "cal_field"
system(paste("v.info -c ", pts))
cmd <- paste("v.db.addcol map=", pts, " columns='", cols, "'", sep = "")
system(cmd)
system(paste("v.info -c ", pts))
cmd <- paste("v.what.rast vect=", pts, " raster=", id, " column=", cols_samp, sep = "")
lapply(cmd, system)
system(paste("v.info -c ", pts))
rm(cmd, pts)

# Expert calibration
pts <- "cal_expert"
system(paste("v.info -c ", pts))
cmd <- paste("v.db.addcol map=", pts, " columns='", cols, "'", sep = "")
system(cmd)
system(paste("v.info -c ", pts))
cmd <- paste("v.what.rast vect=", pts, " raster=", id, " column=", cols_samp, sep = "")
lapply(cmd, system)
system(paste("v.info -c ", pts))
rm(cmd, pts)

rm(cols)

# Calibrate soil prediction models ############################################################################

# Set formula
form <- formula(paste("UM ~ ", paste(cols_samp, collapse = " + ")))

# Field calibration points
cal_field <- spgrass6::readVECT("cal_field")
cal_field$UM <- as.factor(cal_field$UM)
# Look for NAs
str(cal_field)
head(cal_field@data)
nrow(na.omit(cal_field@data)) - nrow(cal_field@data)
which(sapply(lapply(cal_field@data, is.na), any))
# Calibrate LDA model
fit_field_lda <- MASS::lda(form, data = cal_field@data, na.action = na.omit)
fit_field_lda$fitted <- predict(fit_field_lda, newdata = na.omit(cal_field@data), prior = fit_field_lda$prior)
# Calibrate RF model
set.seed(2001)
fit_field_rf <- randomForest::randomForest(form, data = cal_field@data, na.action = na.omit)

# Expert calibration points
cal_expert <- spgrass6::readVECT("cal_expert")
cal_expert$UM <- as.factor(cal_expert$UM)
# Look for NAs
str(cal_expert)
head(cal_expert@data)
nrow(na.omit(cal_expert@data)) - nrow(cal_expert@data)
which(sapply(lapply(cal_expert@data, is.na), any))
# Calibrate LDA model
fit_expert_lda <- MASS::lda(form, data = cal_expert@data, na.action = na.omit)
fit_expert_lda$fitted <- 
  predict(fit_expert_lda, newdata = na.omit(cal_expert@data), prior = fit_expert_lda$prior)
# Calibrate RF model
set.seed(2001)
fit_expert_rf <- randomForest::randomForest(form, data = cal_expert@data, na.action = na.omit)

# Random field calibration points
# Look for NAs
str(cal_random_field)
head(cal_random_field@data)
nrow(na.omit(cal_random_field@data)) - nrow(cal_random_field@data)
which(sapply(lapply(cal_random_field@data, is.na), any))
# Calibrate LDA model
fit_random_field_lda <- MASS::lda(form, data = cal_random_field@data)
fit_random_field_lda$fitted <- 
  predict(fit_random_field_lda, newdata = cal_random_field@data, prior = fit_random_field_lda$prior)
# Calibrate RF model
set.seed(2001)
fit_random_field_rf <- randomForest::randomForest(form, data = cal_random_field@data)

# Random expert calibration points
# Look for NAs
str(cal_random_expert)
head(cal_random_expert@data)
nrow(na.omit(cal_random_expert@data)) - nrow(cal_random_expert@data)
which(sapply(lapply(cal_random_expert@data, is.na), any))
# Calibrate LDA model
fit_random_expert_lda <- MASS::lda(form, data = cal_random_expert@data)
fit_random_expert_lda$fitted <- 
  predict(fit_random_expert_lda, newdata = cal_random_expert@data, prior = fit_random_expert_lda$prior)
# Calibrate RF model
set.seed(2001)
fit_random_expert_rf <- randomForest::randomForest(form, data = cal_random_expert@data)

# Random large calibration points
# Look for NAs
str(cal_random_large)
head(cal_random_large@data)
nrow(na.omit(cal_random_large@data)) - nrow(cal_random_large@data)
which(sapply(lapply(cal_random_large@data, is.na), any))
# Calibrate LDA model
fit_random_large_lda <- MASS::lda(form, data = cal_random_large@data)
fit_random_large_lda$fitted <-
  predict(fit_random_large_lda, newdata = cal_random_large@data, prior = fit_random_large_lda$prior)
# Calibrate RF model
set.seed(2001)
fit_random_large_rf <- randomForest::randomForest(form, data = cal_random_large@data)

# Prepare validation data #####################################################################################

# Load required data
system("r.mask -o access_limit")
cmd <- paste("r.out.xyz input=target_soil_map output=data/tmp/target_soil_map.csv", sep = "")
system(cmd)
soil_map <- read.table("data/tmp/target_soil_map.csv", sep = "|")
soil_map$id <- 1:nrow(soil_map)
soil_map <- soil_map[order(soil_map$V3), ]

# Get stratified simple random sample (proportional to area)
size <- round(2000 * (summary(as.factor(soil_map$V3)) / length(soil_map$V3)))
area <- list(total = (length(soil_map$V3) * 25), strata = (size * 25))
set.seed(1984)
val_sample <- sampling::strata(
  data = soil_map[order(soil_map$V3), ], stratanames = "V3", size = size, method = "srswor")
str(val_sample)
summary(as.factor(val_sample$V3))

# Prepare validation data
soil_map <- soil_map[val_sample$ID_unit, ]
val_sample <- list(sample = val_sample, data = cbind(soil_map[, -c(1:2)], covars[soil_map$id, ]))
val_sample$data$V3 <- as.factor(val_sample$data$V3)
levels(val_sample$data$V3) <- 
  levels(raster::shapefile("data/vector/target_soil_map.shp", stringsAsFactors = TRUE)$UM)
colnames(val_sample$data) <- c("UM", colnames(val_sample$data)[-1])
str(val_sample$data)

# Field data
val_field_lda <- purity(obs = val_sample$data$UM, fit = predict(fit_field_lda, val_sample$data)$class)
val_field_rf <- purity(obs = val_sample$data$UM, fit = predict(fit_field_rf, val_sample$data))

# Expert data
val_expert_lda <- purity(obs = val_sample$data$UM, fit = predict(fit_expert_lda, val_sample$data)$class)
val_expert_rf <- purity(obs = val_sample$data$UM, fit = predict(fit_expert_rf, val_sample$data))

# Random data (field)
val_random_field_lda <- 
  purity(obs = val_sample$data$UM, fit = predict(fit_random_field_lda, val_sample$data)$class)
val_random_field_rf <- purity(obs = val_sample$data$UM, fit = predict(fit_random_field_rf, val_sample$data))

# Random data (expert)
val_random_expert_lda <- 
  purity(obs = val_sample$data$UM, fit = predict(fit_random_expert_lda, val_sample$data)$class)
val_random_expert_rf <- purity(obs = val_sample$data$UM, fit = predict(fit_random_expert_rf, val_sample$data))

# Random data (large)
val_random_large_lda <- 
  purity(obs = val_sample$data$UM, fit = predict(fit_random_large_lda, val_sample$data)$class)
val_random_large_rf <- purity(obs = val_sample$data$UM, fit = predict(fit_random_large_rf, val_sample$data))
