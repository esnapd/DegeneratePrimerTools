#' Function to draw results on geographic map.
#'
#' @param geoData Required. Data to be displayed on the map (dataframe)
#' 
#' @importFrom maps map
#' @export
draw_geomap_results <- function(geoData){
  
  map(database="world", xlim = c(-175,-50), ylim = c(25,75), col = "grey80", fill = TRUE)
  lon <- geoData$Longitude
  lat <- geoData$Latitude
  #coords <- mapproject(lon, lat, proj="Molecules")
  points(lon, lat, pch=20, col="red")
  #map("state", col = "grey80", fill = TRUE)
  #points(lon, lat, pch=20, col="red")
  
}
