# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: get_ncdf.R
# | DATE: 
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

setwd("C:/Users/kyungmin/Documents/Biascorrection/NASPA")
### Path for Data and Output	
data_path <- "../data"
output_path <- "../output"

### Set up output folders
write_output_path <- output_path
#write_output_path <- file.path(output_path, "flow_model")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

# The data is 0.5*0.5 resolution, unit = mm
fn_wm <- "NASPA_WARM_TOTALPRECIP.nc"

#####################################################
#Now open Naspa ncdf4 package
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

# plot for certain years, timing in US
prcp <- ncvar_get(nc_rm, varid= var1$names[1], start = c(175,1,1), count = c(1,-1,-1))
image.plot(lon_list,lat_list, t(prcp))
maps::map(add = T)

#find the row and colum near the columbus gauge
# Columbus: Lat: 39.9612 N, Lon: 82.9988 W
# Oklahoma city : 35.4676° N, 97.5164° W
# San diego: 32.7157° N, 117.1611° W
lat_col <- which(lat_list > 35 & lat_list < 35.5)
lon_col <- which(lon_list > -98 & lon_list < -97.5)

var_data_j <- ncvar_get(nc_rm, varid= var1$names[1], start = c(1,lat_col,lon_col), 
                  count = c(-1,1,1))
plot(x = date_list, y = var_data_j, type = 'l')

yup <- tibble(site = "Oklahoma City", date = date_list, 
              variable = as.character(ncdf_df$short_name), value = var_data_j)

naspa_df <- yup %>%
  #select(date, precip_mm_day) %>%
  mutate(precip = var_data_j / 92) %>%
  mutate(location = "columbus") %>%
  mutate(model = "NASPA") %>%
  mutate(scenario = "reconstruction") %>%
  mutate(units = "mm") %>%
  mutate(year = year(date)) %>%
  mutate(month_day= paste0(month(date),"-",day(date))) %>%
  select(site, date,year,month_day, variable, model, precip, units)	

######################################################################
#Download Gridmet data
##########################################################################3
ncdf_df <- ncdf_df %>%
  mutate(url = paste0('http://thredds-prod.nkn.uidaho.edu:8080/thredds/dodsC/agg_met_',short_name, 
                      '_1979_CurrentYear_CONUS.nc#fillmismatch'))
#ncdf_info <- ncdf_df[1,]

### Open the NCDF file
nc_file <- nc_open(ncdf_df$url)

### Extract lat, lon, and time info
lat_list <- ncvar_get(nc_file, "lat")
lon_list <- ncvar_get(nc_file, "lon")
date_list <- ncvar_get(nc_file, "day") + as.Date("1900-01-01")

lat_col <- which(lat_list > 35 & lat_list < 35.5)
lon_col <- which(lon_list > -98 & lon_list < -97.5)
#lat_col <- which(lat_list > 39.5 & lat_list < 40)
#lon_col <- which(lon_list > -83 & lon_list < -82.5)

count <- lat_col[length(lat_col)]-lat_col[1]

var_data_j <- ncvar_get(nc_file, as.character(ncdf_df$var_name), 
                         start=c(lon_col[1], lat_col[1], 1), 
                         count=c(count,count,-1)) 

### Remove missing values
var_data_j[var_data_j ==32767] <- NA

a <- colMeans(var_data_j, dim = 2)

yup <- tibble(site = "Oklahoma City", date = date_list, 
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
  mutate(roll_mean_3 = rollmeanr(x=value, k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(precip = case_when(roll_mean_3_notna  > n_roll_min ~ roll_mean_3,
                                 TRUE ~ NA_real_)
  ) %>%
  select(-roll_mean_3_notna) %>%
  filter(month_day == "8-1") %>%
  mutate(units = "mm") %>%
  select(site,date, year, month_day, variable, model, precip,units)
	
### Close the nc file
nc_close(nc_file)

combined_df <-naspa_df %>%
  bind_rows(gridmet_df) %>%
  #bind_rows(forts_df) %>% 
  group_by(model) %>%
  arrange(year) %>%
  drop_na(precip)

p <- ggplot(combined_df, aes(x=year, y=precip, colour=model)) + geom_line() +
  labs(title = "Oklahoma City")
p
filename <- paste0("precip_OK.png")
ggsave(filename = filename, plot = p, path = write_output_path, width =12.5, height = 8, dpi = 300)

#######################################################################
## Choose a single Reanalysis file
#####################################################################
ncdf_link <- "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/20thC_ReanV3/Monthlies/accumsSI-MO/apcp.mon.mean.nc"

### Read locally
ncdf_link <- file.path(data_path, "noaa/apcp.mon.mean.nc")

### Open the NCDF file
nc_file <- nc_open(ncdf_link)

### Extract lat, lon, and time info
lat_list <- ncvar_get(nc_file, "lat")
lon_list <- ncvar_get(nc_file, "lon")

ncatt_get(nc_file, "time", "units")
tm <- ncvar_get(nc_file, "time")
dim(tm)
ncatt_get(nc_file,0,"title")

date_list <- ncvar_get(nc_file, "time") + ymd_hms("1800-1-1 00:00:0.0")

lat_col <- which(lat_list >= 39.5 & lat_list <= 40)
lon_col <- which(lon_list >= -83 + 180 & lon_list <= -82.5 + 180)

count <- lat_col[length(lat_col)] - lat_col[1] + 1

var_data_j <- ncvar_get(nc_file, as.character("apcp"), 
                        start=c(lon_col[1], lat_col[1], 1), 
                        count=c(count,count,-1)) 

var_data_tmp <- ncvar_get(nc_file, as.character("apcp"), 
                        start=c(1, 1, 1000), 
                        count=c(-1,-1,1)) 

image(lon_list, lat_list, var_data_tmp)

### Remove missing values
var_data_j[var_data_j > 1.00e+20] <- NA

### Extract variable
yup <- tibble(site = "columbus", date = date_list, 
              variable = as.character(ncdf_df$short_name), value = var_data_j, 
              model = "Noaa_v3")

accum_df <- yup %>%
  arrange(date) %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  mutate(jdate = yday(date)) 
  
accum_df <- accum_df %>%
  mutate(roll_mean_3 = rollmeanr(x=precip, k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(roll_mean_3_notna = rollsumr(x=!is.na(precip), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(roll_mean_3 = case_when(roll_mean_3_notna > n_roll_min ~ roll_mean_3,
                                 TRUE ~ NA_real_)
  ) %>%
  select(-roll_mean_3_notna, -precip)

### Close the nc file
nc_close(nc_file)

pr_modGS <- gam(precip ~ ti(year, bs = "cr", k = 80),
                  data = combined_df, family= Gamma)
draw(pr_modGS)

precip <- transform(combined_df, modGS = predict(pr_modGS, type="response"),
                    na.rm = TRUE)

p <- ggplot(data=precip, aes(x=modGS, y=precip,color = model)) +
  geom_point(alpha=0.3) +
  scale_x_continuous(limits = c(0,8)) +
  scale_y_continuous(limits = c(0,8)) +
  coord_cartesian() +
  geom_abline() +
  labs(x="Predicted precip (model *GS*)", y= "Observed precip")

filename <- paste0("GS","observed_vs_modeled.png")
ggsave(filename = filename, plot = p, path = write_output_path, width =12.5, height = 8, dpi = 300)

GS_predict <- predict(pr_modGS, se.fit = TRUE, type = "response")

precip <- transform(combined_df,
                    modGS = GS_predict$fit, 
                    modGS_se = GS_predict$se.fit)

p <- ggplot(data=precip, aes(x=year, y=precip, group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_ribbon(aes(ymin=(modGS-2*modGS_se),
                  ymax=(modGS+2*modGS_se)), alpha=0.25) +
  geom_line(aes(y=modGS)) +
  geom_point(alpha = 0.25) +
  labs(title = "August1_GS model: knots = 80",x="year",y="precip(mm)", size = 14)
p
filename <- paste0("GS_","GAM_modeled_k80.png")
ggsave(filename = filename, plot = p, path = write_output_path, width =12.5, height = 8, dpi = 300)


########################################################
#########             GI        ########################
########################################################
pr_modGI <- gam(precip ~ model + ti(year, bs = "cr", k = 80),
                data=combined_df, family = Gamma)

plotpr <- pr_modGI$fitted.values + 
plot(pr_modGI$fitted.values)

precip <- transform(combined_df, modGI = predict(pr_modGI, type="response"),
                    na.rm = TRUE)

p <- ggplot(data=precip, aes(x=modGI, y=precip)) +
  geom_point(alpha=0.1) +
  scale_x_continuous(limits = c(0,8)) +
  scale_y_continuous(limits = c(0,8)) +
  coord_cartesian() +
  geom_abline() +
  labs(title = "Aug.1st GI model: knots = 20",
       x="Predicted precip (model *GI*)", y= "Observed precip")
p

GI_predict <- predict(pr_modGI, se.fit = TRUE, type = "response")

precip <- transform(combined_df,
                    modGI = GI_predict$fit, 
                    modGI_se = GI_predict$se.fit)

p <- ggplot(data=precip, aes(x=year, y=precip, group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_ribbon(aes(ymin=(modGI-2*modGI_se),
                  ymax=(modGI+2*modGI_se)), alpha=0.25) +
  geom_line(aes(y=modGI)) +
  geom_point(alpha = 0.25) +
  labs(title = "August1_GI model: knots = 80",
       x="year",y="precip(mm)", size = 14)
p
filename <- paste0("GI_","GAM_modeled_k80.png")
ggsave(filename = filename, plot = p, path = write_output_path, width =12.5, height = 8, dpi = 300)

##mapping example

mapCDFtemp <- function(lat,lon,tas) #model and perc should be a string
  
{
  titletext <- "title"
  
  expand.grid(lon_list, lat_list) %>%
    
    rename(lon = Var1, lat = Var2) %>%
    
    mutate(lon = ifelse(lon > 180, -(360 - lon), lon),
           
           tas = as.vector(var_data_temp)) %>%  
     ggplot() + 
     geom_point(aes(x = lon, y = lat, color = tas),size = 0.8) + 
     borders("world", colour="black", fill=NA) + 
     scale_color_viridis(na.value="white",name = "precip") + 
     theme(legend.direction="vertical", legend.position="right", 
           legend.key.width=unit(0.4,"cm"), legend.key.heigh=unit(2,"cm")) + 
     coord_quickmap() + 
     ggtitle(titletext) 
  
}

mapCDFtemp(lat = lat_list,lon = lon_list, var_data_tmp)




