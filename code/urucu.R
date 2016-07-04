

# Spatial exploratory data analysis
tmp <- raster::shapefile("data/vector/Trein_Classes.shp")
str(tmp@data)
plot(as.factor(tmp@data$MDS))
