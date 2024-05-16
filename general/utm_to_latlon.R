library(sf)

#Converts a matrix of eastings and northings (eastings first column, northings second column) to UTM
#Inputs:
#	EastNorths: [N x 2 matrix] of eastings (col 1) and northings (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#Outputs:
#	LonLats or LatLons: [N x 2 matrix] of longitudes and latitudes - lons are first column 
utm_to_latlon <- function(EastNorths, utm.zone = 34,southern_hemisphere=TRUE){
  utms <- data.frame(X=EastNorths[,1],Y=EastNorths[,2])
  if(southern_hemisphere){
    crs_string <- paste0("+proj=utm +datum=WGS84 +units=m +no_defs +south +zone=", utm.zone)
  } else{
    crs_string <- paste0("+proj=utm +datum=WGS84 +units=m +no_defs +north +zone=", utm.zone)
  }
  utms_sf <- utms %>% sf::st_as_sf(coords = c('X','Y'), crs = crs_string)
  lonlat_sf <- utms_sf %>% sf::st_transform(crs = 4326)
  LonsLats <- sf::st_coordinates(lonlat_sf)
  return(LonsLats)
  
}
