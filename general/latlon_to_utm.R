library(sf)

#LAT/LON TO UTM CONVERSIONS (AND VICE VERSA)
#Converts a matrix of lons and lats (lons first column, lats second column) to UTM
#Inputs:
#	LonsLats: [N x 2 matrix] of lons (col 1) and lats (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#Outputs:
#	EastNorths: [N x 2 matrix] of Eastings and Northings - eastings are first column
latlon_to_utm <- function(LonsLats,utm.zone='34',southern_hemisphere=TRUE){
  lonlat <- data.frame(lon = LonsLats[,1], lat = LonsLats[,2])
  lonlat_sf <- lonlat %>% sf::st_as_sf(coords = c('lon','lat'), crs = 4326)
  if(southern_hemisphere){
    crs_string <- paste0("+proj=utm +datum=WGS84 +units=m +no_defs +south +zone=", utm.zone)
  } else{
    crs_string <- paste0("+proj=utm +datum=WGS84 +units=m +no_defs +north +zone=", utm.zone)
  }
  xy_sf <- lonlat_sf %>% sf::st_transform(crs = crs_string)
  EastsNorths <- sf::st_coordinates(xy_sf)
  return(EastsNorths)
  
}