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

# Deterministic landform classification algorithm for R ----
fit_land_classification <-
  function (x, tool = "r", vname) {
    
    param <- c(3.50, -0.05, 1.70, 60.00, 15000.00, 8.00)
    # param <- c(3.60, -0.04, 3.60, 62.00, 15000.00, 5.20)
    
    if (tool == "r") {
      pb <- txtProgressBar(min = 0, max = nrow(x), style = 3)
      res <- vector()
      for (i in 1:nrow(x)) {
        if (x$slope[i] > param[1]) {
          res[i] <- "CXal+PVA+PVal"
        } else {
          if (x$curvature[i] < param[2] || x$slope[i] < param[3] && x$elevation[i] < param[4]) {
            res[i] <- "GXvd+RYq+CXbd+ESkg"
          } else {
            if (x$flow_down[i] > param[5] && x$twi[i] > param[6]) {
              res[i] <- "PACd"
            } else {
              res[i] <- "PAd+PAal+PAa"
            }
          }
        }
        setTxtProgressBar(pb, value = i)
      }
      return (as.factor(res))
    } else {
      if (tool == "grass") {
        cmd <- 
          paste("r.mapcalc --overwrite '", vname, " = ",
                "if(slope > ", param[1], ", 1, ", 
                "if(curvature < ", param[2], " || slope < ", param[3], " && elevation < ", param[4], ", 2, ",
                "if(flow_down > ", param[5], " && twi > ", param[6], ", 3, 4)))'", sep = "")
        grassGis(cmd)
      }
    }
  }

# Colour ramps ----
uncertainty.colors <- 
  colorRampPalette(c("olivedrab", "khaki", "maroon1"))
soil.colors <- 
  # c("lightsalmon", "cadetblue1", "azure2", "lightgoldenrod1")
  c("lightsalmon", "cornflowerblue", "azure2", "lightgoldenrod1")

# Spatial prediction -----
spPredict <-
  function (model, covariates, sp = TRUE, probs = FALSE) {
    if (inherits(model, "lda")) {
      res <- data.frame(predict(model, newdata = covariates)$posterior)
    } else {
      res <- data.frame(predict(model, covariates, type = "prob"))
    }
    
    if (probs) {
      res$UM <- maxCol(res)
      res$entropy <- apply(res[, 1:4], 1, entropy)  
    } else {
      res <- data.frame(
        UM = maxCol(res), entropy = apply(res[, 1:4], 1, entropy))
    }
    
    if (sp) {
      res <- cbind(covariates[, 1:2], res)
      sp::gridded(res) <- ~ x + y
    }
    return (res)
  }
