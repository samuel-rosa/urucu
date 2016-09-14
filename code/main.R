# title: "Soil Spatial Modelling"
# author: "Alessandro Samuel Rosa"
# 
# Spatial exploratory data analysis were performed using QGIS and SAGA GIS. The results are described in 
# \code\exploratory.Rmd. The GRASS GIS database was created using \code\dbGRASS.R. All covariates were loaded
# to the database, as well as all existing vector (polygon) data.

# Clean up and source user defined helper functions ###########################################################
rm(list = ls())
source("code/helper.R")

# Start GRASS GIS #############################################################################################
spgrass7::initGRASS(
  gisBase = gisBase, gisDbase = "data/GRASS", location = "urucu", mapset = "database",
  override = TRUE, pid = Sys.getpid())
grassGis("r.mask --o target_soil_map")

# Convert all GRASS 6 vector maps to GRASS 7
# Source: https://grasswiki.osgeo.org/wiki/Convert_all_GRASS_6_vector_maps_to_GRASS_7
# grassGis("v.build.all")
# grassGis("db.connect -d")
# grassGis("db.connect -p")
# grassGis("v.db.reconnect.all -cd")

# Load core packages ##########################################################################################
require(MASS)
require(randomForest)

# Prepare calibration data ####################################################################################

# Load necessary data
target_soil_map <- spgrass7::readVECT("target_soil_map")
access_limit <- spgrass7::readVECT("access_limit")
non_access_limit <- spgrass7::readVECT("non_access_limit")
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
plot(cal_field@coords)
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
plot(cal_expert@coords)

# Select three probability samples of sizes equal to n = 'cal_field', n = 'cal_expert' and n = 2000.
# Here we create the probabilistic calibration datasets. First we need to load the covariate data to create
# our balanced sample. Due to the large size of the study region, it is wiser to export a csv file omiting
# cells with NA values to a temporary file.
id_covars <- c("elevation", "curvature", "plan_curv", "slope", "flow_down", "twi")
cmd <- paste(
  "r.out.xyz input=", paste(id_covars, collapse = ","), " output=data/tmp/access_covars.csv", sep = "")
grassGis(cmd)
covars <- read.table("data/tmp/access_covars.csv", sep = "|")
colnames(covars) <- c("x", "y", id_covars)
head(covars)

# Start with the field-based random balanced samples, i.e. the probability sample with as many observations
# as in the field calibration sample.
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
plot(cal_random_field@coords)

# Now we produce the expert-based random balanced sample, i.e. the probability sample with as many 
# observations as in the expert calibration sample.
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
plot(cal_random_expert@coords)

# Finally, we produce the large random balanced sample, i.e. the probability sample with as many observations
# as we can get (n = 2000).
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
plot(cal_random_large@coords)

# Load field and expert calibration points to GRASS GIS
spgrass7::writeVECT(SDF = cal_field, vname = "cal_field", v.in.ogr_flags = c("overwrite"))
spgrass7::writeVECT(SDF = cal_expert, vname = "cal_expert", v.in.ogr_flags = c("overwrite"))

# Setup database of calibration points (sample data from covariates)
grassGis("r.mask -o non_access_limit")
cols <- paste(id, "DOUBLE PRECISION")
cols <- paste(cols, collapse = ",")
cols_samp <- substring(id, first = 1, last = 10)

# Field calibration
pts <- "cal_field"
grassGis(paste("v.info -c ", pts))
cmd <- paste("v.db.addcol map=", pts, " columns='", cols, "'", sep = "")
grassGis(cmd)
grassGis(paste("v.info -c ", pts))
cmd <- paste("v.what.rast vect=", pts, " raster=", id, " column=", cols_samp, sep = "")
lapply(cmd, grassGis)
grassGis(paste("v.info -c ", pts))
rm(cmd, pts)
cal_field <- spgrass7::readVECT("cal_field")
cal_field$UM <- as.factor(cal_field$UM)

# Expert calibration
pts <- "cal_expert"
grassGis(paste("v.info -c ", pts))
cmd <- paste("v.db.addcol map=", pts, " columns='", cols, "'", sep = "")
grassGis(cmd)
grassGis(paste("v.info -c ", pts))
cmd <- paste("v.what.rast vect=", pts, " raster=", id, " column=", cols_samp, sep = "")
lapply(cmd, grassGis)
grassGis(paste("v.info -c ", pts))
rm(cmd, pts)
cal_expert <- spgrass7::readVECT("cal_expert")
cal_expert$UM <- as.factor(cal_expert$UM)

rm(cols)

# Save calibration points
save(cal_field, cal_expert, cal_random_field, cal_random_expert, cal_random_large, 
     file = "data/R/calibration_points.rda")

# Prepare figure with calibration observations
# Transform coordinates to kilometres to improve plotting
# This is to see the spatial distribution of the sample point in the accessible area.
n <- sapply(list(cal_field, cal_expert, cal_random_field, cal_random_expert, cal_random_large), nrow)
id <- c("Field", "Expert", rep("Random", 3))
id <- paste(id, " (n = ", n, ")", sep = "")
xy <- list(cal_field, cal_expert, cal_random_field, cal_random_expert, cal_random_large)
xy <- lapply(seq_along(xy), function (i) {
  res <- as.data.frame(sp::coordinates(xy[[i]]))
  colnames(res) <- c("x", "y")
  res$cal <- id[i]
  return (res)
})
xy <- do.call(rbind, xy)
xy[, 1:2] <- xy[, 1:2] / 1000
map <- lattice::xyplot(
  y ~ x | cal, xy, xlab = "Easting (km)", ylab = "Northing (km)", aspect = "iso", layout = c(3, 2),
  par.settings = list(fontsize = list(text = 8, points = 6)),
  panel = function (x, y, ...) {
    lattice::panel.grid(v = -1, h = -1)
    lattice::panel.polygon(
      x = access_limit@polygons[[1]]@Polygons[[1]]@coords[, 1] / 1000, 
      y = access_limit@polygons[[1]]@Polygons[[1]]@coords[, 2] / 1000, 
      col = "lightcoral", border = "lightcoral")
    lattice::panel.xyplot(x, y, col = "black", pch = 20, cex = 0.1)
  })
map$index.cond[[1]] <- c(4, 5, 3, 2, 1)
dev.off()
png("res/fig/calibration_points.png", width = 480*3, height = 480*2, res = 300)
map
dev.off()

rm(map, xy)
gc()

# Prepare figure with the distribution of points per category
xy <- list(cal_field, cal_expert, cal_random_field, cal_random_expert, cal_random_large)
xy <- lapply(xy, function (x) x$UM)
names(xy) <- id
xy <- lapply(xy, function (x) summary(x) / sum(summary(x))) 
xy <- cbind(stack(xy), um = names(xy[[1]]))
col <- c("lightsalmon", "cornflowerblue", "azure2", "lightgoldenrod1")
p <- lattice::barchart(
  values ~ um | ind, data = xy, col = col, 
  scales = list(x = list(draw = FALSE)), ylab = "Proportion",
  key = list(text = list(as.character(unique(xy$um))), rectangles = list(col = col)),
  par.settings = list(fontsize = list(text = 8, points = 6)),
  panel = function (...) {
    lattice::panel.grid(v = -1, h = -1)
    lattice::panel.barchart(...)
  })
p$index.cond[[1]] <- c(4, 5, 3, 2, 1)
names(p$legend) <- "bottom"
dev.off()
png("res/fig/cal_points_prop.png", width = 480*3, height = 480*2, res = 300)
p
dev.off()

rm(p, col, xy)

# Check how well the random samples representat the covariate data
tmp <- rbind(
  cbind(stack(covars[, 3:ncol(covars)]), id = "pop"), 
  cbind(stack(cal_random_field@data[, -ncol(cal_random_field@data)]), id = "n = 383"), 
  cbind(stack(cal_random_expert@data[, -ncol(cal_random_expert@data)]), id = "n = 838"),
  cbind(stack(cal_random_large@data[, -ncol(cal_random_large@data)]), id = "n = 2003"))
p <- lattice::bwplot(
  values ~ id | ind, data = tmp, scales = list(y = list(relation = "free")),
  panel = function (...) {
    lattice::panel.grid(v = -1, h = -1)
    lattice::panel.bwplot(...)
  })
p$par.settings <- list(plot.symbol = list(cex = 0.2), box.dot = list(cex = 0.5))
dev.off()
png("res/fig/balanced_sampling.png", width = 480*4, height = 480*2, res = 250)
p
dev.off()

# Calibrate soil prediction models ############################################################################

# Set formula
form <- formula(paste("UM ~ ", paste(cols_samp, collapse = " + ")))

# Field calibration points
str(cal_field)
head(cal_field@data)
nrow(na.omit(cal_field@data)) - nrow(cal_field@data)
which(sapply(lapply(cal_field@data, is.na), any))
# Calibrate LDA model
fit_field_lda <- MASS::lda(form, data = cal_field@data, na.action = na.omit)
fit_field_lda$predicted <- 
  predict(fit_field_lda, newdata = na.omit(cal_field@data), prior = fit_field_lda$prior)$class
# Calibrate RF model
set.seed(2001)
fit_field_rf <- randomForest::randomForest(form, data = cal_field@data, na.action = na.omit)

# Expert calibration points
# Look for NAs
str(cal_expert)
head(cal_expert@data)
nrow(na.omit(cal_expert@data)) - nrow(cal_expert@data)
which(sapply(lapply(cal_expert@data, is.na), any))
# Calibrate LDA model
fit_expert_lda <- MASS::lda(form, data = cal_expert@data, na.action = na.omit)
fit_expert_lda$predicted <- 
  predict(fit_expert_lda, newdata = na.omit(cal_expert@data), prior = fit_expert_lda$prior)$class
# Calibrate RF model
set.seed(2001)
fit_expert_rf <- randomForest::randomForest(form, data = cal_expert@data, na.action = na.omit)

# Field-based random calibration points
# Look for NAs
str(cal_random_field)
head(cal_random_field@data)
nrow(na.omit(cal_random_field@data)) - nrow(cal_random_field@data)
which(sapply(lapply(cal_random_field@data, is.na), any))
# Calibrate LDA model
fit_random_field_lda <- MASS::lda(form, data = cal_random_field@data)
fit_random_field_lda$predicted <- 
  predict(fit_random_field_lda, newdata = cal_random_field@data, prior = fit_random_field_lda$prior)$class
# Calibrate RF model
set.seed(2001)
fit_random_field_rf <- randomForest::randomForest(form, data = cal_random_field@data)

# Expert-based random calibration points
# Look for NAs
str(cal_random_expert)
head(cal_random_expert@data)
nrow(na.omit(cal_random_expert@data)) - nrow(cal_random_expert@data)
which(sapply(lapply(cal_random_expert@data, is.na), any))
# Calibrate LDA model
fit_random_expert_lda <- MASS::lda(form, data = cal_random_expert@data)
fit_random_expert_lda$predicted <- 
  predict(fit_random_expert_lda, newdata = cal_random_expert@data, prior = fit_random_expert_lda$prior)$class
# Calibrate RF model
set.seed(2001)
fit_random_expert_rf <- randomForest::randomForest(form, data = cal_random_expert@data)

# Large set of random calibration points
# Look for NAs
str(cal_random_large)
head(cal_random_large@data)
nrow(na.omit(cal_random_large@data)) - nrow(cal_random_large@data)
which(sapply(lapply(cal_random_large@data, is.na), any))
# Calibrate LDA model
fit_random_large_lda <- MASS::lda(form, data = cal_random_large@data)
fit_random_large_lda$predicted <-
  predict(fit_random_large_lda, newdata = cal_random_large@data, prior = fit_random_large_lda$prior)$class
# Calibrate RF model
set.seed(2001)
fit_random_large_rf <- randomForest::randomForest(form, data = cal_random_large@data)

# Deterministic landform classification algorithm
fit_lca <- fit_land_classification(cal_field@data)
overallPurity(obs = cal_field$UM, fit = fit_lca)

# Save calibrated models
save(fit_land_classification, fit_field_rf, fit_field_lda, fit_expert_rf, fit_expert_lda, fit_random_field_rf,
     fit_random_field_lda, fit_random_expert_rf, fit_random_expert_lda, fit_random_large_rf, 
     fit_random_large_lda, file = "data/R/calibrated_models.rda")

# Compute calibration purity and save as a csv file
cal_purity <- 
  list(
    field = sapply(list(fit_field_lda, fit_field_rf), 
                   function (x) overallPurity(obs = cal_field$UM, fit = x$predicted)),
    expert = sapply(list(fit_expert_lda, fit_expert_rf), 
                    function (x) overallPurity(obs = cal_expert$UM, fit = x$predicted)),
    random_field = sapply(list(fit_random_field_lda, fit_random_field_rf), 
                          function (x) overallPurity(obs = cal_random_field$UM, fit = x$predicted)),
    random_expert = sapply(list(fit_random_expert_lda, fit_random_expert_rf), 
                           function (x) overallPurity(obs = cal_random_expert$UM, fit = x$predicted)),
    random_large = sapply(list(fit_random_large_lda, fit_random_large_rf), 
                          function (x) overallPurity(obs = cal_random_large$UM, fit = x$predicted)))
cal_purity$field <- c(cal_purity$field, overallPurity(obs = cal_field$UM, fit = fit_lca))
cal_purity <- stack(cal_purity)
cal_purity$model <- c("FLD", "BRF", "LCA", rep(c("FLD", "BRF"), 4))

# Save csv file with calibration purity
tmp <- data.frame(x = round(cal_purity$values * 100, 2))
rownames(tmp) <- paste(cal_purity$ind, "_", cal_purity$model, sep = "")
write.csv(tmp, file = "res/tab/calibration_purity.csv")
rm(tmp)

# n <- sapply(list(cal_field, cal_expert, cal_random_field, cal_random_expert, cal_random_large), nrow)
# id <- c("Field", "Expert", rep("Random", 3))
# id <- paste(id, "\n(n = ", n, ")", sep = "")
# good_fit$ind <- factor(rep(id, each = 2), levels = id)
# col <- c("darkseagreen1", "lightgoldenrod1")
# p <- lattice::barchart(
#   values ~ ind, group = model, good_fit, col = col, ylab = "Calibration agreement", 
#   xlab = "Calibration dataset", 
#   key = list(text = list(unique(good_fit$model)), rectangles = list(col = col)),
#   par.settings = list(fontsize = list(text = 8, points = 6))
#   )
# names(p$legend) <- "inside"
# p$legend$inside$x <- 0.6
# p$legend$inside$y <- 0.875
# dev.off()
# png("res/fig/fit_accuracy.png", width = 8, height = 8, units = "cm", res = 600)
# p
# dev.off()
# 
# rm(p, good_fit, n, id)
# gc()

# Prepare validation data #####################################################################################

# Load required data
grassGis("r.mask -o access_limit")
cmd <- paste("r.out.xyz input=target_soil_map output=data/tmp/target_soil_map.csv", sep = "")
grassGis(cmd)
soil_map <- read.table("data/tmp/target_soil_map.csv", sep = "|")
soil_map$id <- 1:nrow(soil_map)
soil_map <- soil_map[order(soil_map$V3), ]

# Get stratified simple random sample (proportional to area) using the reference soil map
size <- round(2000 * (summary(as.factor(soil_map$V3)) / length(soil_map$V3)))
area <- list(total = length(soil_map$V3), strata = summary(as.factor(soil_map$V3)))
set.seed(1984)
val_sample <- sampling::strata(data = soil_map, stratanames = "V3", size = size, method = "srswor")
str(val_sample)
summary(as.factor(val_sample$V3))

# Prepare validation data
# Get IDs of selected units
soil_map <- soil_map[val_sample$ID_unit, ]
val_sample <- list(sample = val_sample, data = cbind(soil_map[, -c(1:2)], covars[soil_map$id, ]))
val_sample$data$V3 <- as.factor(val_sample$data$V3)
levels(val_sample$data$V3) <- 
  levels(raster::shapefile("data/vector/target_soil_map.shp", stringsAsFactors = TRUE)$UM)
colnames(val_sample$data) <- c("UM", colnames(val_sample$data)[-1])
str(val_sample$data)

# save data
save(val_sample, file = "data/R/validation_points.rda")

# Validation ##################################################################################################
# We compute the overall actual purity for each sampling strata, which are used to compute a weigthed average
# that represents the overall actual purity of the soil map. The weights are the relative size of each stratum.
# Because the strata coincide with the categories of the target variable, the strata-based overall actual 
# purity is equivalent to the class representation (producer's accuracy), i.e. the proportion of the area
# occupied with a given class that is mapped as that class.

validation <- data.frame(
  observed = val_sample$data$UM,
  strata = as.factor(val_sample$sample$Stratum),
  field_lda = predict(fit_field_lda, val_sample$data)$class,
  field_rf = predict(fit_field_rf, val_sample$data),
  field_lca = fit_land_classification(val_sample$data),
  expert_lda = predict(fit_expert_lda, val_sample$data)$class,
  expert_rf = predict(fit_expert_rf, val_sample$data),
  random_field_lda = predict(fit_random_field_lda, val_sample$data)$class,
  random_field_rf = predict(fit_random_field_rf, val_sample$data),
  random_expert_lda = predict(fit_random_expert_lda, val_sample$data)$class,
  random_expert_rf = predict(fit_random_expert_rf, val_sample$data),
  random_large_lda = predict(fit_random_large_lda, val_sample$data)$class,
  random_large_rf = predict(fit_random_large_rf, val_sample$data)
)

# # Landform classification algorithm
# cal_profiles <- raster::shapefile("data/vector/Perfis.shp")
# cal_boreholes <- raster::shapefile("data/vector/Tradagens.shp")
# val_boreholes <- raster::shapefile("data/vector/Tradagens_Validacao.shp")
# tmp1 <- cal_profiles@data[, c(30, 31, 32, 40, 41, 39, 36, 37)]
# tmp2 <- cal_boreholes@data[, c(18, 19, 20, 27, 28, 25, 35, 24)]
# tmp3 <- val_boreholes@data[, c(26, 27, 28, 36, 37, 34, 32, 33)]
# colnames(tmp1) <- c("x", "y", id_covars)
# colnames(tmp2) <- c("x", "y", id_covars)
# colnames(tmp3) <- c("x", "y", id_covars)
# tmp <- rbind(tmp1, tmp2, tmp3)
# sp::coordinates(tmp) <- ~ x + y
# sp::proj4string(tmp) <- sp::proj4string(target_soil_map)
# tmp@data$UM <- sp::over(tmp, target_soil_map)$UM
# tmp <- tmp[-which(is.na(tmp@data$UM)), ]
# tmp2 <- tmp[sample(1:nrow(tmp), 96), ]
# 
# # Compare DEM data
# fit1 <- fit_land_classification(cal_field@data) # processed dem data
# fit2 <- fit_land_classification(tmp2@data) # original dem data
# overallPurity(obs = cal_field$UM, fit = fit1)
# overallPurity(obs = tmp2$UM, fit = as.factor(fit2))
# 
# # Simulate validation datasets
# a <- vector()
# i <- 1
# for (i in 1:1000) {
#   tmp2 <- tmp[sample(1:nrow(tmp), 96), ]
#   fit2 <- fit_land_classification(tmp2@data)
#   a[i] <- overallPurity(obs = tmp2$UM, fit = as.factor(fit2))
# }

# Compute actual class representation and save data in a csv file
class_representation <- list(
  mean = sapply(3:ncol(validation), function (i) {
    c(by(validation, validation$strata, function (x) overallPurity(obs = x[, 1], fit = x[, i])))
  }))
colnames(class_representation$mean) <- names(validation)[-c(1:2)]
rownames(class_representation$mean) <- levels(val_sample$data$UM)
class_representation$se <- sqrt((class_representation$mean * (1 - class_representation$mean)) / (size - 1))

tmp <- lapply(class_representation, function (x) round(x * 100, 2))
tmp <- matrix(paste(tmp$mean, " (", tmp$se, ")", sep = ""), nrow = nrow(tmp$mean))
colnames(tmp) <- colnames(class_representation$mean)
rownames(tmp) <- rownames(class_representation$mean)
tmp <- t(tmp)
write.csv(tmp, file = "res/tab/class_representation.csv")
rm(tmp)

# Compute overall actual purity and save csv file 
val_purity <- 
  list(mean = apply(class_representation$mean, 2, function (x) sum(x * (area$strata / area$total))))
w2 <- (area$strata / area$total) ^ 2
val_purity$se <- 
  sqrt(apply(class_representation$mean, 2, function (x) sum(w2 * ((x * (1 - x)) / (size - 1)))))
val_purity <- data.frame(val_purity)
tmp <- apply(val_purity, 1, function (x) 
  paste(round(x[1] * 100, 2), " (", round(x[2] * 100, 2), ")", sep = ""))
write.csv(tmp, file = "res/tab/validation_purity.csv")
rm(tmp)

# Compute the probability of occurence of each category
occurence_prob <- list(
  field_lda = predict(fit_field_lda, val_sample$data)$posterior,
  field_rf = predict(fit_field_rf, val_sample$data, type = "prob"),
  expert_lda = predict(fit_expert_lda, val_sample$data)$posterior,
  expert_rf = predict(fit_expert_rf, val_sample$data, type = "prob"),
  random_field_lda = predict(fit_random_field_lda, val_sample$data)$posterior,
  random_field_rf = predict(fit_random_field_rf, val_sample$data, type = "prob"),
  random_expert_lda = predict(fit_random_expert_lda, val_sample$data)$posterior,
  random_expert_rf = predict(fit_random_expert_rf, val_sample$data, type = "prob"),
  random_large_lda = predict(fit_random_large_lda, val_sample$data)$posterior,
  random_large_rf = predict(fit_random_large_rf, val_sample$data, type = "prob")
)

# Compute theoretical purity and save data in a csv file
# We compute the overall mean, median and standard error. The median is compared to the mean to see 
# if the theoretical purity is a skewed variable. The overall theoretical purity is a left skewed 
# variable (the median is always larger than the mean). Besides, the difference between the median and the
# mean increases with the mean.
tmp <- sapply(occurence_prob, function (x) apply(x, 1, max))
theo_purity <- list(
  mean = sapply(1:ncol(tmp), function (i) {
    c(by(tmp, as.factor(val_sample$sample$Stratum), function (x) mean(x[, i])))
  }),
  median = sapply(1:ncol(tmp), function (i) {
    c(by(tmp, as.factor(val_sample$sample$Stratum), function (x) median(x[, i])))
  }),
  var = sapply(1:ncol(tmp), function (i) {
    c(by(tmp, as.factor(val_sample$sample$Stratum), function (x) var(x[, i])))
  })
)
tmp <- 
  data.frame(mean = apply(theo_purity$mean, 2, function (x) sum(x * (area$strata / area$total))),
             median = apply(theo_purity$median, 2, function (x) sum(x * (area$strata / area$total))),
             se = sqrt(apply(theo_purity$var, 2, function (x) sum(w2 * x))))
rownames(tmp) <- rownames(val_purity)[-3]
tmp
tmp <- apply(tmp, 1, function (x) 
  paste(round(x[1] * 100, 2), " (", round(x[3] * 100, 2), ")", sep = ""))
write.csv(tmp, file = "res/tab/theoretical_purity.csv")
rm(tmp)

# Compute validation entropy and save data in a csv file
# The validation entropy seems to be less skewd than the theoretical purity. However, it still is left skewed,
# which means that the median says we are more uncertain than the mean says. This relation changes as the 
# entropy decreases, i.e. the median says we are more certain than the mean says. In other words, the entropy
# becomes right skewed with decreasing values. I guess these are common characteristics of variables bounded 
# to an interval.
val_entropy <- sapply(occurence_prob, function (x) apply(x, 1, entropy))
val_entropy <- list(
  mean = sapply(1:ncol(val_entropy), function (i) {
    c(by(val_entropy, as.factor(val_sample$sample$Stratum), function (x) mean(x[, i])))
  }),
  median = sapply(1:ncol(val_entropy), function (i) {
    c(by(val_entropy, as.factor(val_sample$sample$Stratum), function (x) median(x[, i])))
  }),
  var = sapply(1:ncol(val_entropy), function (i) {
    c(by(val_entropy, as.factor(val_sample$sample$Stratum), function (x) var(x[, i])))
  })
  )
val_entropy$overal <-
  data.frame(mean = apply(val_entropy$mean, 2, function (x) sum(x * (area$strata / area$total))),
             median = apply(val_entropy$median, 2, function (x) sum(x * (area$strata / area$total))),
             se = sqrt(apply(val_entropy$var, 2, function (x) sum(w2 * x))))
rownames(val_entropy$overal) <- rownames(val_purity)[-3]
val_entropy$overal

tmp <- apply(val_entropy$overal[, -2], 1, function (x) 
  paste(round(x[1] * 100, 2), " (", round(x[2] * 100, 2), ")", sep = ""))
write.csv(tmp, file = "res/tab/validation_entropy.csv")
rm(tmp)

# Prepare table with model performance and save as a csv file
id <- c("calibration_purity", "theoretical_purity", "validation_entropy", "validation_purity")
tmp <- lapply(id, function (x) {
  res <- read.csv(paste("res/tab/", x, ".csv", sep = ""))
  colnames(res) <- c("model", x)
  res
})
tmp[[2]] <- rbind(tmp[[2]][1:2, ], NA, tmp[[2]][3:nrow(tmp[[2]]), ])
tmp[[3]] <- rbind(tmp[[3]][1:2, ], NA, tmp[[3]][3:nrow(tmp[[3]]), ])
tmp <- do.call(cbind, tmp)[, c("model", id)]
write.csv(tmp, "res/tab/overall_model_performance.csv")
rm(tmp)

# Spatial predictions #########################################################################################

# Load spatial prediction models
load("data/R/calibrated_models.rda")

# Prepare covariates
# We first work in a small area to check the perfomance of all prediction models.
# We export a text file containing only the data from inside the prediction region (NAs are not exported).
grassGis("r.mask --o map_inset")
id_covars <- c("elevation", "curvature", "plan_curv", "slope", "flow_down", "twi")
cmd <- paste(
  "r.out.xyz input=", paste(id_covars, collapse = ","), " output=data/tmp/inset_covars.csv", sep = "")
grassGis(cmd)
inset_covars <- read.table("data/tmp/inset_covars.csv", sep = "|")
colnames(inset_covars) <- c("x", "y", id_covars)
str(inset_covars)

# Field calibration data (n = 383)
pred_field_lda <- spPredict(fit_field_lda, inset_covars)
pred_field_rf <- spPredict(fit_field_rf, inset_covars)
fit_land_classification(tool = "grass", vname = "inset_land_class")
cmd <- paste("r.out.xyz --overwrite input=inset_land_class output=data/tmp/inset_land_class.csv", sep = "")
grassGis(cmd)
pred_field_lca <- read.table("data/tmp/inset_land_class.csv", sep = "|")
colnames(pred_field_lca) <- c("x", "y", "UM")
pred_field_lca$UM <- as.factor(pred_field_lca$UM)
levels(pred_field_lca$UM) <- levels(pred_expert_lda$UM)
sp::gridded(pred_field_lca) <- ~ x + y

# Expert calibration data (n = 837)
pred_expert_lda <- spPredict(fit_expert_lda, inset_covars)
pred_expert_rf <- spPredict(fit_expert_rf, inset_covars)

# Random calibration data (n = 383)
pred_random_field_lda <- spPredict(fit_random_field_lda, inset_covars)
pred_random_field_rf <- spPredict(fit_random_field_rf, inset_covars)

# Random calibration data (n = 837)
pred_random_expert_lda <- spPredict(fit_random_expert_lda, inset_covars)
pred_random_expert_rf <- spPredict(fit_random_expert_rf, inset_covars)

# Random calibration data (n = 2000)
pred_random_large_lda <- spPredict(fit_random_large_lda, inset_covars)
pred_random_large_rf <- spPredict(fit_random_large_rf, inset_covars)



sp::spplot(pred_field_lca, "UM", col.regions = soil.colors)








# We export a text file containing only the data from inside the prediction region (NAs are not exported).
# Attention: processing the full grid at once requires a lot of RAM (16 Gb).
grassGis("r.mask --o non_access_limit")


grassGis(cmd)
covars <- read.table("data/tmp/non_access_covars.csv", sep = "|")
colnames(covars) <- c("x", "y", id_covars)
head(covars)
str(covars)

# Field calibration data
pred_field_lda <- predict(fit_field_lda, newdata = covars)$posterior
save(pred_field_lda, file = "data/R/pred_field_lda.rda")
rm(pred_field_lda)
gc()

pred_field_rf <- predict(fit_field_rf, newdata = covars, type = "prob")
save(pred_field_rf, file = "data/R/pred_field_rf.rda")
rm(pred_field_rf)
gc()

coords <- covars[, 1:2]
rm(covars)
gc()

load("data/R/pred_field_rf.rda")
coords <- cbind(coords, pred_field_rf)
write.csv(coords, file = "data/tmp/pred_rf.csv")
rm(coords)
gc()
