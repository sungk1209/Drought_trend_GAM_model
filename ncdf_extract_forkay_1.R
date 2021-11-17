

###########################################################################
###  Load functions
###########################################################################
require(tidyverse)

require(sf)
require(mapproj)
#require(scales)
require(rgdal)

require(ncdf4)

###########################################################################
###  Set up projection
###########################################################################

### Modis meter sinusoidal
modis <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
### WGS84 This is the default lat lon
wgs84 <- "+init=epsg:4326"  ### This is what cliwoc uses

###########################################################################
###  Set up point for grid
###########################################################################
### Single point
p1 <- st_sf(a = "SanJuan,CO", geom=st_sfc(st_point(c(-106.75,37.75)), crs = 4326))

###########################################################################
###  Read in, check position
###########################################################################

### Choose a single NCDF file
ncdf_info <- "https://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-CM4/ssp585/r1i1p1f1/day/pr/gr1/v20180701/pr_day_GFDL-CM4_ssp585_r1i1p1f1_gr1_20150101-20341231.nc"

### Open the NCDF file
nc_file <- nc_open(ncdf_info)

### Extract lat, lon, and time info
lat_list <- ncvar_get(nc_file, "lat")
lon_list <- ncvar_get(nc_file, "lon")
date_list <- ncvar_get(nc_file, "time") + as.Date("1850-01-01")

### Recreate the Lat Lon grid
lat_lon_df <- expand.grid(lat = lat_list, lon=lon_list)
lat_lon_sf <- st_as_sf(lat_lon_df, coords = c("lon", "lat"),  crs = 4326)

### Project the lat lon to meters
lat_lon_subset <- lat_lon_sf %>%
	st_transform(crs = modis) 

### Find nearest point
nearest_test <- p1 %>%
	st_transform(crs = modis) %>%
	st_nearest_feature(lat_lon_subset)

### Extract lat lon
nearest_pt <-lat_lon_df[nearest_test,]
nearest_pt$a <- p1$a

### Find the grid cell that corresponds
lat_i <- match(nearest_pt$lat, lat_list)
lon_i <- match(nearest_pt$lon, lon_list)
nearest_pt$lat_index <-lat_i
nearest_pt$lon_index <- lon_i

### Extract variable data
var_data_j <-  ncvar_get(nc_file, "pr", start=c(nearest_pt$lon_index, nearest_pt$lat_index, 1), count=c(1,1,-1)) 

### Extract variable
var_df <- tibble(site = nearest_pt$a, date = date_list, variable = "precip", value = var_data_j)

### Close the nc file
nc_close(nc_file)

p1
nearest_pt
head(var_df)

