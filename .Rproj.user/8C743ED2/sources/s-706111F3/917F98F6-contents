# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: 0.downloadforgam.R
# | DATE: Aug.31.2021
# | CREATED BY:  Kay Sung     
# *--------------------------------------------------------------------------
# | PURPOSE: 
# | Download gridmet, Noaa reanlaysis and NASPA data
# | 
# | 
# *--------------------------------------------------------------------------
require(stars)
require(ncdf4)
require(lubridate)
require(fields)
require(purrr)
require(tidyverse)
require(zoo)
require(mgcv)
require(gratia)
require(dplyr)

select <- dplyr::select
### Path for Data and Output	
data_path <- "../data/"
output_path <- "../output/"

##set up the location
#find the row and column 
# Columbus: Lat: 39.9612 N, Lon: 82.9988 W
# Oklahoma city : 35.4676? N, 97.5164? 
# San Diego: 32.7157? N, 117.1611? W
# San Antonio 29.4814578N,99.073463 W

#loc <- data.frame(site="OKC,OK",
#                  lon=c(-98,-97.5),
#                  lat= c(35.0,35.5))

#loc <- data.frame(site="phoenix,AZ",
#                  lon= c(-112.5,-112.0),
#                  lat=c(33.50,34.00))

#loc <- data.frame(site="SanJuan_CO",
#                 lon= c(-106.5,-106.0),
#                  lat=c(37.50,38.00))
#loc <- data.frame(site = "COL_OH",
#					lon = c(-83,-82.5),
#					lat = c(39.5, 40))


######################################################################
#1. Download Gridmet data
##########################################################################3
inst_download <- function(loc){

ncdf_df <- data.frame(short_name = c("pr"))

ncdf_df <- ncdf_df %>%
  mutate(url = paste0('http://thredds-prod.nkn.uidaho.edu:8080/thredds/dodsC/agg_met_',short_name, 
                      '_1979_CurrentYear_CONUS.nc#fillmismatch'))


### Open the NCDF file
nc_file <- nc_open(ncdf_df$url)
var_prcp <- attributes(nc_file$var)$names

### Extract lat, lon, and time info
lat_list <- ncvar_get(nc_file, "lat")
lon_list <- ncvar_get(nc_file, "lon")
date_list <- ncvar_get(nc_file, "day") + as.Date("1900-01-01")

#find grid cells belongs to one city and average all measures
#lat_col <- which(lat_list > 35 & lat_list < 35.5) #OKC 35.4676? N, 97.5164? W
#lon_col <- which(lon_list > -98 & lon_list < -97.5)
# San antonio 29.4814578,-99.073463

lat_col <- which(lat_list > loc$lat[1] & lat_list < loc$lat[2])
lon_col <- which(lon_list > loc$lon[1] & lon_list < loc$lon[2])

count <- lat_col[length(lat_col)]-lat_col[1]

var_data_j <- ncvar_get(nc_file, as.character(var_prcp), 
                         start=c(lon_col[1], lat_col[1], 1), 
                         count=c(count,count,-1)) 

### Remove missing values
var_data_j[var_data_j ==32767] <- NA

####Gridmet has fine resolution, need to average over all grid cells
#Unit of gridcell: daily accumulated 
a <- colMeans(var_data_j, dim = 2, na.rm = TRUE)

yup <- tibble(site =  loc$site[1], date = date_list, 
                  variable = as.character(var_prcp), value = a, 
              model = "Gridmet")

accum_df <- yup %>%
  arrange(date) %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  group_by(year,month)%>%
  summarise(mprcp = sum(value)) %>%
  ungroup()
  

n_roll <-3
n_roll_min <-2
gridmet_df <- accum_df %>%
  mutate(roll_mean_3 = rollsumr(x=mprcp, k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(roll_mean_3_notna = rollsumr(x=!is.na(mprcp), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(precip = case_when(roll_mean_3_notna  > n_roll_min ~ roll_mean_3,
                                 TRUE ~ NA_real_)
  ) %>%
  mutate(precip = precip/3) %>%
  select(-roll_mean_3_notna,mprcp) %>%
  mutate(units = "mm/month") %>%
  mutate(date = as.Date(paste0(year,'-',month,'-01'))) %>%
  mutate(site = loc$site[1]) %>%
  mutate(variable = as.character(ncdf_df$short_name)) %>%
  mutate(model = "gridmet") %>%
  select(date, site, year, variable, model, precip,units)
	
### Close the nc file
nc_close(nc_file)

#######################################################################
## Choose a CRU file
#####################################################################

nc_cru <- nc_open(filename = "../data/cru_ts4.05.1901.2020.pre.dat.nc")
attributes(nc_cru$var)$names
var_prcp <- attributes(nc_cru$var)$names[[1]]
var_unit <- ncatt_get(nc_cru, attributes(nc_cru$var)$names[[1]], "units")
 
lat_list <- ncvar_get(nc_cru,"lat")
lon_list <- ncvar_get(nc_cru,"lon")

time_scale <- ncatt_get(nc_cru, "time", "units")$value %>%
  str_split_fixed("since",n =2)
date_list <- ncvar_get(nc_cru, "time") + as.Date(time_scale[2])

#Example names as the row names, column 1 as longitude, and column 2 as latitude
lat_col <- which(lat_list > loc$lat[1] & lat_list < loc$lat[2])
lon_col <- which(lon_list > loc$lon[1] & lon_list < loc$lon[2])

var_data_j <- ncvar_get(nc_cru, varid= var_prcp, start = c(lon_col,lat_col,1), 
                        count = c(1,1,-1))

var_data_j[var_data_j < -9.96920996838687e+36] <- NA

yup <- tibble(site = loc$site[1], date = date_list, 
              variable = var_prcp, value = var_data_j)

cru_df <- yup %>%
  arrange(date) %>%
  mutate(year = year(date)) %>%
  mutate(model = "cru") %>%
  mutate(scenario = "observed")
  
n_roll <- 3
n_roll_min <- 2
cru_df <- cru_df %>%
  mutate(roll_mean_3 = rollsumr(x=value, k=n_roll, fill=NA, na.rm=TRUE))  %>%
#  mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
#  mutate(precip = case_when(roll_mean_3_notna > n_roll_min ~ roll_mean_3,
#                              TRUE ~ NA_real_)) %>%
  mutate(precip = roll_mean_3/3) %>%
  mutate(units = "mm") %>% 
  mutate(date = as.Date(floor_date(date, "month"))) %>%
  select(date, site, year, variable,model, precip,units)



################################################################
##########Combine instrumental data 
###############################################################

instrument_df <- rbind(gridmet_df,cru_df)

instrument_df <- instrument_df %>%
  arrange(date) %>%
  drop_na(precip) 
saveRDS(instrument_df,file =paste0(output.p,"/instrument_df.rds"))

p <- ggplot(data = instrument_df %>%filter(month(date) == 1)) + 
  geom_line(aes(x=date, y=precip, group = model, colour = model)) +

  scale_color_manual(values = c("#1A3FD4", "orangered2")) +
  labs(title =  loc$site[1],
       y = "3 months ave. precip.(mm)") +
  theme_classic(base_size = 30)


p
filename <- paste0(output.p,"/Pinst")
ggsave(filename =paste0(filename,".png"), plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename =paste0(filename,".svg"), plot = p, width =12.5, height = 8, dpi = 300)


return()
}

###Discarded

#######################################################################
## Choose a single Reanalysis file
#####################################################################
ncdf_link <- "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/20thC_ReanV3/Monthlies/accumsSI-MO/apcp.mon.mean.nc"

### Read locally
#ncdf_link <- file.path(data_path, "noaa/apcp.mon.mean.nc")

### Open the NCDF file
nc_file <- nc_open(ncdf_link)
var_prcp <- attributes(nc_file$var)$names[[2]]

### Extract lat, lon, and time info
lat_list <- ncvar_get(nc_file, "lat")
lon_list <- ncvar_get(nc_file, "lon")

ncatt_get(nc_file, "time", "units")$value
tm <- ncvar_get(nc_file, "time")
dim(tm)
ncatt_get(nc_file,0,"title")

#date_list <- ncvar_get(nc_file, "time") + ymd_hms("1800-1-1 00:00:0.0")
date_list <- ncvar_get(nc_file, "time")/24 + as.Date("1800-01-01")

lat_col <- which(lat_list >= loc$lat[1] & lat_list <= loc$lat[2])
lon_col <- which(lon_list >= loc$lon[1] + 360 & lon_list <= loc$lon[2] + 360)

count <- lat_col[length(lat_col)] - lat_col[1] + 1

var_data_j <- ncvar_get(nc_file, var_prcp, 
                        start=c(lon_col[1], lat_col[1], 1), 
                        count=c(count,count,-1)) 

### Remove missing values
var_data_j[var_data_j < -9.96920996838687e+36] <- NA

### Extract variable #unit: 3hourly data have to multiply 8 and make monthly
yup <- tibble(site =  loc$site[1], date = date_list, 
              variable = variable, value = var_data_j*8 , 
              model = "Noaa_v3", unit = "mm")
accum_df <- yup %>%
  arrange(date) %>%
  mutate(year = year(date)) %>%
  mutate(month_day = paste0(month(date),"-",day(date))) 

n_roll <- 3
n_roll_min <- 2
reanal_v3 <- accum_df %>%
  mutate(roll_mean_3 = rollsumr(x=value, k=n_roll, fill=NA, na.rm=TRUE))  %>%
  mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(precip = case_when(roll_mean_3_notna > n_roll_min ~ roll_mean_3,
                            TRUE ~ NA_real_)) %>%
  select(-roll_mean_3_notna) %>%
  mutate(precip = precip * 10) %>%
  #filter(month_day == "7-1") %>%
  mutate(units = "mm/month") %>%
  select(date, site, year, variable, model, precip,units, month_day) #date removed
### Close the nc file
reanal_v3


nc_close(nc_file)
