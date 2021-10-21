
require(stars)
require(ncdf4)
require(lubridate)
require(fields)
require(purrr)
require(tidyverse)
require(zoo)
require(mgcv)

fn_wm <- "NASPA_WARM_TOTALPRECIP.nc"
short_name <- "pr"
ncdf_df <- data.frame(short_name = c("pr"), 
                      var_name = c("precipitation_amount"))
nc_rm <- nc_open(fn_wm)
print(nc_rm)

var1<-attributes(nc_rm$var)
lat_list <- ncvar_get(nc_rm, "lat")
lon_list <- ncvar_get(nc_rm, "lon")
date_list <- ncvar_get(nc_rm, "time")
tm_orig <- as.Date("0000-08-01")
date_list <- tm_orig %m+% years(date_list)

loc <- data.frame(site="OKC,OK",
                  lon=c(-98,-97.5),
                  lat= c(35.0,35.5))

lat_col <- which(lat_list > loc$lat[1] & lat_list < loc$lat[2])
lon_col <- which(lon_list > loc$lon[1] & lon_list < loc$lon[2])

var_data_j <- ncvar_get(nc_rm, varid= var1$names[1], start = c(1,lat_col,lon_col), 
                        count = c(-1,1,1))
plot(x = date_list, y = var_data_j, type = 'l')

yup <- tibble(site = "Oklahoma City", date = date_list, 
              variable = as.character(ncdf_df$short_name), value = var_data_j)

naspa_df <- yup %>%
  #select(date, precip_mm_day) %>%
  mutate(precip = var_data_j / 3) %>%
  mutate(model = "NASPA") %>%
  mutate(scenario = "reconstruction") %>%
  mutate(units = "mm/month") %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
 # mutate(month_day= paste0(month(date),"-",day(date))) %>%
  select(site, date,year,month, variable, model, precip, units)	

