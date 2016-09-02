# Start OS dependent GRASS GIS ----
if (.Platform$OS.type == "unix") {
  # gisBase <- "/usr/lib/grass64/"
  gisBase <- "/usr/lib/grass70/"
} else {
  gisBase <- "C:/Program Files (x86)/GRASS GIS 6.4.4"
}

# OS dependent function to run GRASS features ----
grassGis <- function (cmd) {
  if (.Platform$OS.type == "unix") {
    system(cmd)
  } else {
    shell(cmd)
  }
}

# Overall purity ----
overallPurity <- 
  function (obs, fit, weights) {
    tab <- table(data.frame(fit, obs))
    diagonal <- diag(tab)
    if (missing(weights)) {
      res <- sum(diag(tab)) / length(obs)
    } else {
      res <- sum(diagonal * weights)
    }
    # map_unit_purity <- diagonal / rowSums(tab)
    # class_representation <- diagonal / colSums(tab)
    # res <- 
    # list(overall_purity = overall_purity, 
    # by_class = data.frame(mu_purity = map_unit_purity, class_rep = class_representation))
    return (res)
  }

# Return the columns with the maximum value ----
maxCol <- 
  function (x, colnames = TRUE, as.factor = TRUE, ...) {
    if (colnames) {
      res <- colnames(x)[max.col(x, ...)]
    } else {
      res <- max.col(x, ...)
    }
    
    if (as.factor) {
      res <- as.factor(res)
    }
    return (res)
  }

# Shannon entropy ----
entropy <-
  function (x) {
    - sum(x * log(x, base = length(x)), na.rm = TRUE)
  }

# Deterministic landform classification algorithm ----
fit_land_classification <-
  function (x) {
    res <- vector()
    for (i in 1:nrow(x)) {
      if (x$slope[i] > 3.50) {
        # if (x$slope[i] > 3.60) {
        res[i] <- "CXal+PVA+PVal"
      } else {
        if (x$curvature[i] < -0.05 || x$slope[i] < 1.70 && x$elevation[i] < 60.00) {
          # if (x$curvature[i] < -0.04 || x$slope[i] < 3.6 && x$elevation[i] < 62) {
          res[i] <- "GXvd+RYq+CXbd+ESkg"
        } else {
          if (x$flow_down[i] > 15000.00 && x$twi[i] > 8.00) {
            # if (x$flow_down[i] > 15000.00 && x$twi[i] > 5.2) {
            res[i] <- "PACd"
          } else {
            res[i] <- "PAd+PAal+PAa"
          }
        }
      }
    }
    return (as.factor(res))
  }

# Colour ramps ----
uncertainty.colors <- 
  colorRampPalette(c("olivedrab", "khaki", "maroon1"))
soil.colors <- 
  c("lightsalmon", "cadetblue1", "azure2", "lightgoldenrod1")
