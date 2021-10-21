

NASPA_COOL_SPI

fn_co <- "NASPA_COOL_SPI.nc"

short_name <- "pr"

ncdf_df <- data.frame(short_name = c("pr"), 
                      var_name = c("precipitation_amount"))
nc_co <- nc_open(fn_co)
print(nc_co)

var1 <- attributes(nc_co$var)
lat_list <- ncvar_get(nc_co, "lat")
lon_list <- ncvar_get(nc_co, "lon")
date_list <- ncvar_get(nc_co, "time")
tm_orig <- as.Date("0000-04-30")
date_list <- tm_orig %m+% years(date_list)

lat_col <- which(lat_list > loc$lat[1] & lat_list < loc$lat[2])
lon_col <- which(lon_list > loc$lon[1] & lon_list < loc$lon[2])

var_data_k <- ncvar_get(nc_co, varid= var1$names[1], start = c(1,lat_col,lon_col), 
                        count = c(-1,1,1))

yuc <- tibble(date = date_list, 
              spi5 = var_data_k)

naspa_spi5 <- yuc %>% 	
  mutate(month = 4) %>%
  mutate(year = year(date)) %>%
  complete(year = seq(0,2016), month = seq(1,12)) %>%
  mutate(date = as.Date(paste0(year, "-", month, "-15"))) %>%
  mutate(date = as.Date(ceiling_date(date, "month")-1)) %>%
  arrange(date)




