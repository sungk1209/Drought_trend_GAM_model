# *------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: .R
# | DATE: 
# | CREATED BY:  Jim Stagge -Kay edited        
# *----------------------------------------------------------------
# | PURPOSE:  download GCM with CMIP6 scenario and construct 
# |          3 months precipitation 
# | 
# |
# *------------------------------------------------------------------

###########################################################################
## Set the Paths
###########################################################################

### Path for Data and Output	
data_path <- "../data/"
output_path <- "../output/"

### Set up output folders
write_output_path <- output_path
#write_output_path <- file.path(output_path, "flow_model")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
#write_figures_path <- file.path(write_output_path, "figures")
#dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)


###########################################################################
###  Load functions
###########################################################################
require(tidyverse)
require(here)
require(lubridate)
require(zoo)
#require(dataRetrieval)
require(ggplot2)
#require(ggrepel)

require(sf)
require(mapproj)
#require(scales)
require(rgdal)

require(ncdf4)

select <- dplyr::select


###########################################################################
###  Read in Data
###########################################################################
load(file.path(write_output_path, "gridmet_results.RData"))

ls()

###################################################################
###  Define Projections
###################################################################
#nad83 <- "+init=epsg:4269"f
wgs84 <- "+init=epsg:4326"  ### This is what cliwoc uses

### Modis meter sinusoidal
modis <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

###################################################################
###  Create points for gauges - use these to download climate data
###################################################################

lon_pt <- 254
lat_pt <- 38

### Single point
p1 <-  st_sf(a = "SanJuan", geom=st_sfc(st_point(c(lon_pt,lat_pt)), crs = 4326))

### Multiple gauges

plot(p1, axes=TRUE, graticule=TRUE)

### Project point
p1_proj <-  p1 %>%
	st_transform(crs = modis) 
plot(p1, axes=TRUE, graticule=TRUE)

###################################################################
###  Create a buffer
###################################################################

### Create a buffer around DC
buffer_zone <- p1 %>%
	st_transform(crs = modis) %>% ### Needs to be projected to m for proper distances
	st_union() %>%
	st_buffer(dist = 300000) %>% # 300 km 
	st_transform(crs = wgs84)

plot(buffer_zone, axes=TRUE, graticule=TRUE)
plot(st_geometry(us_states), add=TRUE)

plot(st_geometry(us_states),  axes=TRUE, graticule=TRUE)
plot(buffer_zone, add=TRUE)

###########################################################################
###  Prepare a list of CMIP6 data to download
###########################################################################
### Create an object to hold the data links
ncdf_historical <- tibble(short_name = "pr", 
	standard_name= "precipitation_flux", 
	begin_date = paste0(seq(1850, 2010, 20), "0101"), 
	end_date = c(paste0(seq(1869, 2009, 20), "1231"), "20141231"), 
	scenario = "historical", 
	model="GFDL-CM4")
### Add URL
ncdf_historical <- ncdf_historical %>% 
	mutate(url = paste0('https://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/day/pr/gr1/v20180701/pr_day_GFDL-CM4_historical_r1i1p1f1_gr1_', begin_date, '-', end_date, '.nc#fillmismatch'))

ncdf_historical


### Create an object to hold the data links
ncdf_ssp585 <- tibble(short_name = "pr", 
	standard_name= "precipitation_flux", 
	begin_date = paste0(seq(2015, 2095, 20), "0101"), 
	end_date = c(paste0(seq(2034, 2094, 20), "1231"), "21001231"), 
	scenario = "ssp585", 
	model="GFDL-CM4")
### Add URL
ncdf_ssp585 <- ncdf_ssp585 %>% 
	mutate(url = paste0('https://esgdata.gfdl.noaa.gov/thredds/dodsC/gfdl_dataroot4/ScenarioMIP/NOAA-GFDL/GFDL-CM4/ssp585/r1i1p1f1/day/pr/gr1/v20180701/pr_day_GFDL-CM4_ssp585_r1i1p1f1_gr1_', begin_date, '-', end_date, '.nc#fillmismatch'))

ncdf_ssp585

### Merge
ncdf_df <- ncdf_historical %>%
	bind_rows(ncdf_ssp585)
### Check
ncdf_df


###########################################################################
###  Download CMIP6 data
###########################################################################
for (i in seq(1,dim(ncdf_df)[1])){

### Choose a single NCDF file
ncdf_info <- ncdf_df[i,]

### Open the NCDF file
nc_file <- nc_open(ncdf_info$url)

### Extract lat, lon, and time info
lat_list <- ncvar_get(nc_file, "lat")
lon_list <- ncvar_get(nc_file, "lon")
date_list <- ncvar_get(nc_file, "time") + as.Date("1850-01-01")

lat_col <- which(lat_list >= loc$lat[1] & lat_list <= loc$lat[2])
lon_col <- which(lon_list >= loc$lon[1] + 359 & lon_list <= loc$lon[2] + 360)

#lon_list <- lon_list - 360

#prcp <- ncvar_get(nc_file, ncdf_info$short_name, start = c(1,1,175), count = c(-1,-1,1))
#image.plot(lon_list,lat_list, t(prcp))
#maps::map(add = T)


### Recreate the Lat Lon grid
#lat_lon_df <- expand.grid(lat = lat_list, lon=lon_list)
#lat_lon_sf <- st_as_sf(lat_lon_df, coords = c("lon", "lat"),  crs = 4326)

### Extract only the points within a buffer
#buffer_test <- buffer_zone %>%
#	st_transform(crs = modis) %>%
#	st_contains(st_transform(lat_lon_sf, crs = modis), sparse = FALSE)

### Subset to just within the buffer zone
#lat_lon_subset <- lat_lon_sf[buffer_test,] %>%
#	st_transform(crs = modis) ==

# Test plot
#plot(buffer_zone, axes=TRUE, graticule=TRUE)
#plot(st_geometry(us_states), add=TRUE)
#plot(lat_lon_sf[buffer_test,], add=TRUE)

### Find nearest point
#nearest_test <- p1 %>%
#	st_transform(crs = modis) %>%
#	st_nearest_feature(lat_lon_subset)

### Extract lat lon
#nearest_pt <-lat_lon_df[buffer_test,][nearest_test,]
#nearest_pt$a <- p1$a

### Find the grid cell that corresponds
#lat_i <- match(nearest_pt$lat, lat_list)
#lon_i <- match(nearest_pt$lon, lon_list)
#nearest_pt$lat_index <- lat_i
#nearest_pt$lon_index <- lon_i

#j <- 1

### Extract variable data
var_data_j <-  ncvar_get(nc_file, ncdf_info$short_name, 
                         start=c(lon_col,lat_col, 1), 
                          count=c(1,1,-1)) 
### Remove missing values
var_data_j[var_data_j > 1.00e+20] <- NA

### Extract variable
yup <- tibble(site = loc$site[1], date = date_list, 
              variable = as.character(ncdf_info$short_name), 
              value = var_data_j, model = ncdf_info$model, 
              scenario = ncdf_info$scenario)

### Close the nc file
nc_close(nc_file)


if (i == 1) {
	clim_df <- yup
} else {
	clim_df <- rbind(clim_df,yup)
}
}

###########################################################################
###  Do some processing
###########################################################################
### Convert from kg/m2s to mm/day
clim_df <- clim_df %>%
  mutate(precip_mm = value * 24*60*60) %>%
	arrange(date)

###########################################################################
###  Check plot
###########################################################################
ggplot(clim_df, aes(x=date, y=precip_mm, colour=scenario)) + geom_line() + theme_bw() + scale_y_continuous(name="Precip (mm)")

###########################################################################
###  Save results from as an RDS file for analysis
###########################################################################
#gfdl_df <- clim_df
### Save resutls for next step
saveRDS(clim_df, file = paste0(output_path,loc$site[1], "gfdl_results.rds"))


#############################################################
######calculate 3 months average precipitation ##############

#make monthly sum
accum_df <- clim_df %>%
  arrange(date) %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  group_by(year,month, .add = TRUE)%>%
  summarise(mprcp = sum(precip_mm)) %>%
  ungroup()

n_roll <-3
n_roll_min <-2
pred_df <- accum_df %>%
  mutate(roll_mean_3 = rollsumr(x=mprcp, k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(roll_mean_3_notna = rollsumr(x=!is.na(mprcp), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(precip = case_when(roll_mean_3_notna  > n_roll_min ~ roll_mean_3,
                            TRUE ~ NA_real_)
  ) %>%
  mutate(precip = precip/3) %>%
  mutate(units = "mm/month") %>%
  mutate(date = as.Date(paste0(year,'-',month,'-01'))) %>%
  mutate(site = loc$site[1]) %>%
  mutate(variable = as.character("pr")) %>%
  mutate(model = "GFDL-CM4") %>%
  select(date, site, year, variable, model, precip,units, month)
ggplot(pred_df, aes(x=date, y=precip)) + geom_line() + theme_bw() + scale_y_continuous(name="Precip (mm)")

saveRDS(pred_df, file = file.path(write_output_path,paste0(loc$site[1],"pred_df.rds") ))
