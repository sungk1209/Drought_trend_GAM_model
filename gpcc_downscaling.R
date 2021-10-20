# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: gpcc_downscaling.R
# | DATE: Oct.02.2021
# | CREATED BY:  Kay Sung     
# *--------------------------------------------------------------------------
# | PURPOSE: 
# | 
# | 
# | 
# *--------------------------------------------------------------------------
require(lubridate)
require(tidyverse)
require(dplyr)
require(ncdf4)
require(fitdistrplus)
require(zoo)
library(RANN)
select <- dplyr::select

data_path <- "../data"
output_path <- "../output/"

### Read in
#warm_spi <- read_tsv("../data/NASPA_WARM_SPI.txt")
naspa_spi3 <- read_tsv("../data/NASPA_Reconstructed_MJJ_0_2016_SPI.txt", skip = 2)%>% 
	rename(year = Year, spi3=SPI)
	
### Add a month column
naspa_spi3 <- naspa_spi3 %>% 	
  mutate(month = 7) %>%
  complete(year = seq(0,2016), month = seq(1,12)) %>%
  mutate(date = as.Date(paste0(year, "-", month, "-15"))) %>%
  mutate(date = as.Date(ceiling_date(date, "month")-1)) %>%
  arrange(date)

head(naspa_spi3)

### SPI-3 3 months equals 3 months
n_days <- 12
n_roll <- 3

gpcc <- "../data/precip.mon.total.v7.nc"

nc_gpcc <- nc_open(gpcc)
print(nc_gpcc)

#gpcc_df: precipitation: total monthly precipitation
# Columbus: Lat: 39.9612 N, Lon: 82.9988 W
#loc <- data.frame(site="columbus,co",
#                  lon=c(-83,-82.5),
#                  lat= c(39.5,40))
loc <- data.frame(site="OKC,OK",
                  lon=c(-98,-97.5),
                  lat= c(35.0,35.5))

var1 <- attributes(nc_gpcc$var)
lat_list <- ncvar_get(nc_gpcc, "lat")
lon_list <- ncvar_get(nc_gpcc, "lon")
date_list <- ncvar_get(nc_gpcc, "time")
tm_orig <- as.Date("1800-01-01")
date_list <- date_list + tm_orig

lat_col <- which(lat_list > loc$lat[1] & lat_list < loc$lat[2])
lon_col <- which(lon_list > loc$lon[1] + 360 & lon_list < loc$lon[2] + 360)

var_data_j <- ncvar_get(nc_gpcc, varid= var1$names[2], start = c(lon_col,lat_col,1), 
                        count = c(1,1,-1))
plot(x = date_list, y = var_data_j, type = 'l')

yup <- tibble(site = loc$site[1], date = date_list, 
              value = as.numeric(var_data_j))

gpcc_df <- yup %>%
  mutate(model = "GPCC") %>%
  mutate(units = "mm") %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  mutate(month_day= paste0(month(date),"-",day(date))) %>%
  select(date,site, year,month, model, value,units)

n_roll <- 3
n_roll_min <- 2

gpcc_df <- gpcc_df %>%
  mutate(roll_mean_3 = rollmeanr(x=value, k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(precip = case_when(roll_mean_3_notna  > n_roll_min ~ roll_mean_3,
                            TRUE ~ NA_real_)
  ) %>%
  select(-roll_mean_3_notna) %>%
  mutate(date = ceiling_date(date, "month") - 1) %>%
  mutate(units = "mm/day") %>%
  select(date, site, model,precip,units) #date removed

gpcc_df <- gpcc_df %>%
  filter(precip > 0)

gpcc_df <- gpcc_df %>%
  group_by(month(date)) %>%
  mutate(shape = fitdist(precip, "gamma")$estimate[[1]],
            rate = fitdist(precip, "gamma")$estimate[[2]]) %>%
  mutate(prob = pgamma(precip, shape = shape, rate = rate)) %>%
  mutate(spi = qnorm(prob,mean = 0, sd = 1)) %>%
  ungroup()

n_library <- length(gpcc_df$precip)
### Convert library into a dataframe
library_df <- data.frame(i = seq(1,n_library), spi_thisjuly = as.numeric(gpcc_df$spi)) %>%
  
  mutate(spi_nextjuly = dplyr::lead(gpcc_df$spi, n=12,  default = NA))%>% 
  drop_na()

library_df[1:50,]
library_df <- library_df %>% select(-i)

n_neighbors <- 20

for (year_i in c(300:2000)){ 

#### Eventually put this in a loop through years
date_subset <- seq(as.Date(paste0(year_i,"-08-01")), as.Date(paste0(year_i+1,"-08-01")), 
                     by = "month")-1
  
naspa_subset <- naspa_spi3 %>%
	filter(date %in% date_subset)

naspa_points <- data.frame(spi_thisjuly = naspa_subset$spi3[[1]], spi_nextjuly = naspa_subset$spi3[[13]])

### Find the k closest points
#library_df <- library_df %>% select(-i)
closest <- nn2(data= library_df,
               query = naspa_points, k=n_neighbors)

#want to improve: pick only July 
for(k in seq(1,n_neighbors)){
	
	closest_k <- closest$nn.idx[[k]]
	fragment_k <- library_df[seq(closest_k, closest_k+12),]
	
	naspa_subset_k <- naspa_subset %>%
		mutate(iter = paste0("iter_", k)) %>%
		mutate(spi3 = fragment_k$spi_thisjuly)

	if(k == 1){
		naspa_subset_iter <- naspa_subset_k
		
	} else {
		naspa_subset_iter <- naspa_subset_iter %>% bind_rows(naspa_subset_k)
	}
}
  smoothed <- loess(spi3 ~ as.numeric(date),data = naspa_subset_iter)
  predicted_k <- data.frame(date= date_subset,
                            spi3 = predict(smoothed, newdata = naspa_subset$date))
if(year_i == 300){
  predicted_df <- predicted_k
  
  } else {
  predicted_df <- predicted_df %>% bind_rows(predicted_k)
  }

}
naspa_subset_iter_wide <- naspa_subset_iter %>% 
	pivot_wider(names_from = iter, values_from = spi3) 

naspa_subset_iter_wide
temp <- gpcc_df %>% 
  filter(date %in% date_subset) %>%
  rename(spi3 = spi) %>%
  mutate(iter = "GPCC") %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
 select(year,month,spi3,date, iter)

temp2 <- rbind(naspa_subset_iter,temp)

p <-ggplot(temp2, aes(x=date, y=spi3)) +
  geom_point(aes(colour= iter)) + 
  geom_smooth(alpha = 0.5) +
  geom_line(data = temp2 %>% filter(iter == "GPCC"), size = 2) +
  labs(title = paste0(year_i,"years, OKC,OK")) 
p
ggsave(p,filename = paste0(output_path,year_i,"OKC.png"))
p <-ggplot(predicted_df %>% filter(year(date) > 1000 & year(date) <1930), 
           aes(x=date, y=spi3)) +
  geom_line() +
  geom_line(data = gpcc_df %>% filter(year(date) > 1900 & year(date) <1930),
            aes(y = spi), color = "blue" ) +
  labs(title = paste0("1900- 1930,OKC,OK"),
       color = "data") 
p

ggsave(p,filename = paste0(output_path,"1900-1930 Tempa.png"))

p <- ggplot(predicted_df %>% filter(year(date) > 1900 & month(date) == 5),
            aes(x = date, y=spi3)) + 
  geom_line() +
 # geom_point(naspa_spi3 %>% filter(year > 1900 & month == 5), 
#             mapping = aes(y = spi3, color = "red")) +
  geom_line(data = gpcc_df%>%filter(month(date) == 5), aes(y = spi), color = "blue" ) +
  #scale_color_manual(values=c("black", "red", "blue"), 
  #       name="Data",
  #       labels=c("simulated", "naspa", "gpcc")) +
  labs(title ="MAM precip of modeled and naspa,OKC,OK") 
p

ggsave(p,filename = paste0(output_path,"May","OKC2.png"))

p <-ggplot(predicted_df , aes(x=date, y=spi3)) +
  geom_line() +
  geom_line(data = gpcc_df,aes(y = spi), color = "blue" ) +
  labs(title = paste0("OKC,OK"),
       color = "data") 

p

