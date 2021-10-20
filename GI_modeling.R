# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: get_ncdf.R
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

### Path for Data and Output	
data_path <- "../data"
output_path <- "../output/"

### Set up output folders
#dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)
##set up the location
#find the row and column 
# Columbus: Lat: 39.9612 N, Lon: 82.9988 W
# Oklahoma city : 35.4676° N, 97.5164° 
# San Diego: 32.7157° N, 117.1611° W
# San Antonio 29.4814578N,99.073463 W

loc <- data.frame(site="OKC,OK",
                  lon=c(-98,-97.5),
                  lat= c(35.0,35.5))
####################################################
#Now open Naspa ncdf4 package
# The data is 0.5*0.5 resolution, unit = mm
fn_wm <- "NASPA_WARM_TOTALPRECIP.nc"
fn_co <- "NASPA_COOL_TOTALPRECIP.nc"

short_name <- "pr"

ncdf_df <- data.frame(short_name = c("pr"), 
                      var_name = c("precipitation_amount"))
nc_rm <- nc_open(fn_wm)
nc_co <- nc_open(fn_co)
print(nc_rm)

var1 <- attributes(nc_co$var)
lat_list <- ncvar_get(nc_co, "lat")
lon_list <- ncvar_get(nc_co, "lon")
date_list <- ncvar_get(nc_co, "time")
tm_orig <- as.Date("0000-04-30")
date_list <- tm_orig %m+% years(date_list)

lat_col <- which(lat_list > loc$lat[1] & lat_list < loc$lat[2])
lon_col <- which(lon_list > loc$lon[1] & lon_list < loc$lon[2])

var_data_j <- ncvar_get(nc_rm, varid= var1$names[1], start = c(1,lat_col,lon_col), 
                  count = c(-1,1,1))

var_data_k <- ncvar_get(nc_co, varid= var1$names[1], start = c(1,lat_col,lon_col), 
                        count = c(-1,1,1))


plot(x = date_list, y = var_data_k, type = 'l')

yup <- tibble(site = loc$site[1], date = date_list, 
              variable = as.character(ncdf_df$short_name), value = var_data_j)
yuc <- tibble(date = date_list, 
              value = var_data_k)

naspa_df <- yup %>%
  #select(date, precip_mm_day) %>%
  mutate(precip = var_data_j) %>%
  mutate(model = "NASPA") %>%
  #mutate(scenario = "reconstruction") %>%
  mutate(units = "mm") %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  mutate(month_day= paste0(month(date),"-",day(date))) %>%
  select(site, year,month, variable, model, precip, units)	

naspa_co <- yuc %>%
  #select(date, precip_mm_day) %>%
  mutate(precip = value) %>%
  mutate(model = "NASPA") %>%
  mutate(scenario = "reconstruction") %>%
  mutate(units = "mm") %>%
  mutate(year = year(date)) %>%
  mutate(month_day= paste0(month(date),"-",day(date))) %>%
  select(date,site, year,month_day, variable, model, precip, units)	

naspa_df <- naspa_df %>%
 # bind_rows(naspa_co) %>%
  arrange(date)

write.csv(naspa_df,'oklahomacity_naspa.csv')
saveRDS(naspa_df, file = "oklahomacity_naspa.rds")
######################################################################
#Download Gridmet data
##########################################################################3
ncdf_df <- ncdf_df %>%
  mutate(url = paste0('http://thredds-prod.nkn.uidaho.edu:8080/thredds/dodsC/agg_met_',short_name, 
                      '_1979_CurrentYear_CONUS.nc#fillmismatch'))

### Open the NCDF file
nc_file <- nc_open(ncdf_df$url)

### Extract lat, lon, and time info
lat_list <- ncvar_get(nc_file, "lat")
lon_list <- ncvar_get(nc_file, "lon")
date_list <- ncvar_get(nc_file, "day") + as.Date("1900-01-01")


#find grid cells belongs to one city and average all measures
#lat_col <- which(lat_list > 35 & lat_list < 35.5) #OKC 35.4676° N, 97.5164° W
#lon_col <- which(lon_list > -98 & lon_list < -97.5)
# San antonio 29.4814578,-99.073463

lat_col <- which(lat_list > loc$lat[1] & lat_list < loc$lat[2])
lon_col <- which(lon_list > loc$lon[1] & lon_list < loc$lon[2])

count <- lat_col[length(lat_col)]-lat_col[1]

var_data_j <- ncvar_get(nc_file, as.character(ncdf_df$var_name), 
                         start=c(lon_col[1], lat_col[1], 1), 
                         count=c(count,count,-1)) 


### Remove missing values
var_data_j[var_data_j ==32767] <- NA

a <- colMeans(var_data_j, dim = 2)

yup <- tibble(site =  loc$site[1], date = date_list, 
                  variable = as.character(ncdf_df$short_name), value = a, 
              model = "Gridmet")

n_roll <- 92
n_roll_min <- n_roll - 14

accum_df <- yup %>%
  arrange(date) %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  mutate(jdate = yday(date)) %>%
  mutate(month_day = paste0(month(date),"-",day(date))) 

gridmet_df <- accum_df %>%
  mutate(roll_mean_3 = rollsumr(x=value, k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(precip = case_when(roll_mean_3_notna  > n_roll_min ~ roll_mean_3,
                                 TRUE ~ NA_real_)
  ) %>%
  select(-roll_mean_3_notna) %>%
  filter(day(date) == 1) %>%
  mutate(units = "mm") %>%
  select(date, site, year, month_day, variable, model, precip,units) #date removed
	
### Close the nc file
nc_close(nc_file)


#######################################################################
## Choose a single Reanalysis file
#####################################################################
ncdf_link <- "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/20thC_ReanV3/Monthlies/accumsSI-MO/apcp.mon.mean.nc"

### Read locally
#ncdf_link <- file.path(data_path, "noaa/apcp.mon.mean.nc")

### Open the NCDF file
nc_file <- nc_open(ncdf_link)

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
variable <- as.character("apcp")

var_data_j <- ncvar_get(nc_file, variable, 
                        start=c(lon_col[1], lat_col[1], 1), 
                        count=c(count,count,-1)) 

### Remove missing values
var_data_j[var_data_j < -9.96920996838687e+36] <- NA

### Extract variable #unit: hourly data have to multiply 8
yup <- tibble(site =  loc$site[1], date = date_list, 
              variable = variable, value = var_data_j * 8 , 
              model = "Noaa_v3", unit = "mm")

accum_df <- yup %>%
  arrange(date) %>%
  mutate(year = year(date)) %>%
  mutate(month_day = paste0(month(date),"-",day(date))) 

  n_roll <- 3
  n_roll_min <- 2
reanal_v3 <- accum_df %>%
  mutate(roll_mean_3 = rollsumr(x=value, k=n_roll, fill=NA, na.rm=TRUE))  %>%
 # mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
 # mutate(precip = 30 * case_when(roll_mean_3_notna > n_roll_min ~ roll_mean_3,
 #                                TRUE ~ NA_real_)) %>%
 # select(-roll_mean_3_notna) %>%
  mutate(precip = 30.67 * roll_mean_3) %>%
 #filter(month_day == "7-1") %>%
  mutate(units = "mm") %>%
  select(date, site, year, variable, model, precip,units, month_day) #date removed
### Close the nc file
nc_close(nc_file)

################################################################
##########Combine instrumental data 
#################################################################
instrument_df <- gridmet_df %>%
  bind_rows(cru_df) %>%
  bind_rows(reanal_v3) %>%
  group_by(model) %>%
  arrange(year) %>%
  drop_na(precip)

instrument_df <- instrument_df %>%
  filter(precip >0)

p <- ggplot(instrument_df, aes(x=year, y=precip, colour=model)) + geom_line(size = 1) +
  labs(title =  loc$site[1],
       y = "3 months ave. precip.(mm)") +
  theme_classic(base_size = 30)

p
filename <- paste0(loc$site[1],"precip_inst")
ggsave(filename =paste0(filename,"3.png"), plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename =paste0(filename,"3.svg"), plot = p, width =12.5, height = 8, dpi = 300)

pr_modGI_ins <- gam(precip ~ model + ti(year, bs = "cr", k = 7),
                data=instrument_df, family = Gamma)

precip <- transform(instrument_df, modGI = predict(pr_modGI_ins, type="response"),
                    na.rm = TRUE)
GI_predict_ins <- predict(pr_modGI_ins, se.fit = TRUE, type = "response")

### make a new combined_df with Gridmet as model
min_noaa <- min((instrument_df %>% filter(model == "Noaa_v3"))$year) 
max_noaa <- max((instrument_df %>% filter(model == "Noaa_v3"))$year)
min_gridm <- min((instrument_df %>% filter(model == "Gridmet"))$year) 
max_gridm <- max((instrument_df %>% filter(model == "Gridmet"))$year)
min_cru <- min((instrument_df %>% filter(model == "cru"))$year) 
max_cru <- max((instrument_df %>% filter(model == "cru"))$year)

instrument_predict_df <- data.frame(year = seq(min_noaa,max_noaa), model = "Gridmet", plot_model = "NOAA") %>%
  bind_rows(data.frame(year = seq(min_gridm,max_gridm), model = "Gridmet", plot_model = "Gridmet")) %>%
  bind_rows(data.frame(year = seq(min_cru,max_cru), model = "Gridmet", plot_model = "cru"))


### Make predictions based on this
GI_predict_modified <- predict(pr_modGI_ins, newdata = instrument_predict_df, se.fit = TRUE, type = "response")


instrument_predict_df <- transform(instrument_predict_df,
                         modGI = GI_predict_modified$fit, 
                         modGI_se = GI_predict_modified$se.fit)

p <- ggplot(data=precip, aes(x=year,  group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_ribbon(data = instrument_predict_df,aes( ymin=(modGI-2*modGI_se),
                                      ymax=(modGI+2*modGI_se)), alpha=0.25) +
  geom_line(data = instrument_predict_df, aes(y=modGI)) +
  geom_line(aes(y=precip),size = 1) +
  labs(title = "August1_reanalysis and observation: knots = 7",
       x="year",y="precip(mm)", size = 14)+
  theme_classic(base_size = 30)

p
filename <- paste0("GI_","inst_com",loc$site[1])
ggsave(filename = paste0(filename,"3.png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,"3.svg"), plot = p, width =12.5, height = 8, dpi = 300)

####################################################################
###                  Combined_df
####################################################################
combined_df <- naspa_df %>%
  bind_rows(gridmet_df) %>%
  bind_rows(cru_df) %>%
  bind_rows(reanal_v3) %>%
  group_by(model) %>%
  arrange(date) %>%
  drop_na(precip)

combined_df <- combined_df %>%
  filter(precip >0)

p <- ggplot(combined_df, aes(x=year, y=precip, colour=model)) + geom_line() +
  labs(title = loc$site[1]) + ylab("3 months averaged rainfall(mm)") +
  theme_classic(base_size = 30)
p

filename <- paste0("precip_naspa",loc$site[1])
ggsave(filename = paste0(filename,"3.png"), plot = p,  width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,"3.svg"), plot = p, width =12.5, height = 8, dpi = 300)

########################################################
#########             GI                     ###########
########################################################
pr_modGI <- gam(precip ~ model + ti(year, bs = "cr", k = 80),
                data=combined_df, family = Gamma)

precip <- transform(combined_df, modGI = predict(pr_modGI, type="response"),
                    na.rm = TRUE)
GI_predict <- predict(pr_modGI, se.fit = TRUE, type = "response")

######make a new model based on each model
min_naspa <- min((combined_df %>% filter(model == "NASPA"))$year) 
max_naspa <- max((combined_df %>% filter(model == "NASPA"))$year)

######make a model with seperate interception using GI##########
################################################################

combined_df_ <- data.frame(year = seq(min_naspa,max_naspa), model = "NASPA", plot_model = "NASPA") %>%
  bind_rows(data.frame(year = seq(min_noaa,max_noaa), model = "NOAA", plot_model = "NOAA")) %>%
  bind_rows(data.frame(year = seq(min_gridm,max_gridm), model = "Gridmet", plot_model = "Gridmet"))
### Make predictions based on this
GI_predict <- predict(pr_modGI, newdata = combined_df, se.fit = TRUE, type = "response")

combined_df <- transform(combined_df,
                         modGI = GI_predict$fit, 
                         modGI_se = GI_predict$se.fit)

p <- ggplot(data=precip, aes(x=year,  group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_ribbon(data = combined_df,aes( ymin=(modGI-2*modGI_se),
                                      ymax=(modGI+2*modGI_se)), alpha=0.25) +
  geom_line(data = combined_df, aes(y=modGI), size = 1) +
  geom_line(aes(y=precip),alpha = 0.6) +
  labs(title = "August1_GI model: knots = 80",
       x="year",y="3 months averaged rainfall(mm)") +
  theme_classic(base_size = 30)

p
filename <- paste0(loc$site[1],"GI_k80")
ggsave(filename = paste0(filename,"3.png"), plot = p,  width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,"3.svg"), plot = p,  width =12.5, height = 8, dpi = 300)

#######################################################################

### make a new combined_df with Gridmet as model
calibrated_df <- data.frame(year = seq(min_naspa,max_naspa), model = "Gridmet", plot_model = "NASPA") %>%
    bind_rows(data.frame(year = seq(min_cru,max_cru), model = "Gridmet", plot_model = "cru")) %>% 
    bind_rows(data.frame(year = seq(min_noaa,max_noaa), model = "Gridmet", plot_model = "NOAA")) %>%
    bind_rows(data.frame(year = seq(min_gridm,max_gridm), model = "Gridmet", plot_model = "Gridmet"))

### Make predictions based on this
GI_predict_modified <- predict(pr_modGI, newdata = calibrated_df, se.fit = TRUE, type = "response")

calibrated_df <- transform(calibrated_df,
                    modGI = GI_predict_modified$fit, 
                    modGI_se = GI_predict_modified$se.fit)

p <- ggplot(data=precip, aes(x=year,  group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_ribbon(data = calibrated_df,aes(ymin=(modGI-2*modGI_se),
                  ymax=(modGI+2*modGI_se)), alpha=0.25) +
  geom_line(aes(y=precip),alpha = 0.5) +
  geom_line(data = calibrated_df, aes(y=modGI), size = 1) +
  labs(title = "August1_GI model: knots = 80",
       x="year",y="precip(mm)", size = 14) +
  theme_classic(base_size = 30)

p

filename <- paste0("calibrated",loc$site[1])

ggsave(filename = paste0(filename,"3.svg"), plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,"3.png"), plot = p, width =12.5, height = 8, dpi = 300)

save(calibrated_df, file = paste0(data_path,filename,".RData"))
#######################################################################


