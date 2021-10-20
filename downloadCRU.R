# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: get_ncdf.R
# | DATE: Sep.20.2021
# | CREATED BY:  Kay Sung     
# *--------------------------------------------------------------------------
# | PURPOSE: 
# | Extracting CRU data 
# | 
# | 
# *--------------------------------------------------------------------------
require(zoo)
require(ncdf4)
require(raster)
require(tidyverse)
require(lubridate)
#ncdf_df <- ncdf_df %>%
#  mutate(url = "https://dap.ceda.ac.uk/badc/cru/data/cru_ts/cru_ts_4.05/data/pre/cru_ts4.05.1901.2020.pre.dat.nc.gz")
#url = "https://dap.ceda.ac.uk/badc/cru/data/cru_ts/cru_ts_4.05/data/pre/cru_ts4.05.1901.2020.pre.dat.nc.gz"
#dest.file <- "D:/Multi_decadal_precp/data/cru_ts4.05.1901.2020.pre.dat.nc.gz"
#download.file(url = url,destfile = dest.file)

nc_cru <- nc_open(filename = "../data/cru_ts4.05.1901.2020.pre.dat.nc")
cru_precip <- brick("./data/cru_ts4.05.1901.2020.pre.dat.nc", varname="pre")

var1 <- attributes(nc_cru$var) #variable name (cru: pre)
lat_list <- ncvar_get(nc_cru,"lat")
lon_list <- ncvar_get(nc_cru,"lon")

time_scale <- ncatt_get(nc_cru, "time", "units")$value %>%
  str_split_fixed("since",n =2)
date_list <- ncvar_get(nc_cru, "time") + as.Date(time_scale[2])

#Example names as the row names, column 1 as longitude, and column 2 as latitude
#Lake Oroville 39.6278275,-121.4944494
#loc <- data.frame(site="lakeOroville,CA",lon=c(-121.5,-121.00),lat= c(39.5,40.0))
lat_col <- which(lat_list > loc$lat[1] & lat_list < loc$lat[2])
lon_col <- which(lon_list > loc$lon[1] & lon_list < loc$lon[2])

#samples <- read.csv("samples.csv", header=TRUE, row.names="site", sep=",")
#pre.sites <- data.frame(extract(cru_precip, samples, ncol=2))
var_data_j <- ncvar_get(nc_cru, varid= var1$names[1], start = c(lon_col,lat_col,1), 
                        count = c(1,1,-1))

var_data_j[var_data_j < -9.96920996838687e+36] <- NA

plot(x = date_list, y = var_data_j)

yup <- tibble(site = loc$site[1], date = date_list, 
              variable = as.character(var1$names[1]), value = var_data_j)

cru_df <- yup %>%
  arrange(date) %>%
  mutate(year = year(date)) %>%
  mutate(month_day = paste0(month(date),"-",day(date))) %>%
  mutate(model = "cru") %>%
  mutate(scenario = "observed") %>%
  mutate(units = "mm")
  
n_roll <- 3
n_roll_min <- 2
cru_df <- cru_df %>%
  mutate(roll_mean_3 = rollsumr(x=value, k=n_roll, fill=NA, na.rm=TRUE))  %>%
 # mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
#  mutate(precip = case_when(roll_mean_3_notna > n_roll_min ~ roll_mean_3,
#                            TRUE ~ NA_real_)) %>%
#  select(-roll_mean_3_notna) %>%
  mutate(precip = roll_mean_3) %>%
#  filter(month_day == "7-16") %>%
  mutate(units = "mm") %>%
  select(date, site, year, variable, model, precip,units, month_day)


