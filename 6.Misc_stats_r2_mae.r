# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: 5.GI_model_loop.R
# | DATE: Dec.20.2021
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

data_path <- "../data/"
output_path <- "../output/"
select <- dplyr::select

site_list <- list()
site_list[1] <- list(data.frame(site="okc_OK",
                               lon=c(-98,-97.5),
                               lat= c(35.0,35.5)))
site_list[2] <- list(data.frame(site="phx_AZ",
                               lon= c(-112.5,-112.0),
                               lat=c(33.50,34.00)))
site_list[3] <- list(data.frame(site="Den_CO",
                                lon= c(-105.0,-104.5),
                                lat=c(39.50,40.00)))
site_list[4] <- list(data.frame(site="Aber_WA",
                                lon= c(-124.0,-123.5),
                                lat=c(46.50,47.00)))
site_list[5] <- list(data.frame(site="Alb_GA",
                                lon= c(-84.5,-84.0),
                                lat=c(31.50,32.00)))
site_list[6] <- list(data.frame(site="Mrv_OH",
                                lon= c(-83.5,-83.0),
                                lat=c(40.00,40.50)))
site_list[7] <- list(data.frame(site="Mor_MN",
                                lon= c(-96.0,-95.5),
                                lat=c(45.50,46.00)))
site_list[8] <- list(data.frame(site="Nyc_NY",
                                lon= c(-74.5,-74.0),
                                lat=c(40.50,41.00)))
site_list[9] <- list(data.frame(site="Los_CA",
                                lon= c(-118.5,-118.0),
                                lat=c(34.00,34.50)))
site_list[10] <- list(data.frame(site="Wax_TX",
                                lon= c(-97.0,-96.5),
                                lat=c(32.00,32.50)))
site_list[11] <- list(data.frame(site="SJn_CO",
                                 lon= c(-106.5,-106.0),
                                 lat=c(37.50,38.00)))
site_list[12] <- list(data.frame(site="Cwa_NC",
                                 lon= c(-79.0,-78.5),
                                 lat=c(33.50,34.00)))
site_list[13] <- list(data.frame(site="Grd_MT",
                                 lon= c(-111,-110.5),
                                 lat=c(44.5,45.00)))
site_list[14] <- list(data.frame(site="Roe_NM",
                                 lon= c(-109,-108.5),
                                 lat=c(31.0,31.50)))


####################################################
############Read R2 at each location 
###################################################

fn_wm <- "NASPA_WARM_STATS.nc"
fn_co <- "NASPA_COOL_STATS.nc"

short_name <- "CRSQ"

nc_wm <- nc_open(paste0(data_path,fn_wm))
nc_co <- nc_open(paste0(data_path,fn_co))
print(nc_wm)

var1 <- attributes(nc_wm$var)
lat_list <- ncvar_get(nc_wm, "lat")
lon_list <- ncvar_get(nc_wm, "lon")

site_df <- data.frame(
  site = c("okc_OK","phx_AZ","Den_CO","Aber_WA","Alb_GA","Mrv_OH",
            "Mor_MN","Nyc_NY","Los_CA","Wax_TX","SJn_CO","Cwa_NC",
            "Grd_MT","Roe_NM"),
  lon_low = c(-98,-112.5,-105.0,-124.0,-84.5,-83.5,
               -96.0,-74.5,-118.5,-97.0,-106.5,-79.0,
               -111,-109),
  lon_up = c(-97.5,-112.0,-104.5,-123.5,-84.0,-83.0,
              -95.5,-74.0,-118.0,-96.5,-106.0,-78.5,
              -110.5,-108.5),
  lat_low = c(35.0,33.50,39.50,46.50,31.50,40.00,
               45.50,40.50,34.00,32.00,37.50,33.50,
               44.5,31.0),
  lat_up = c(35.5, 34.00,40.00, 47.00,32.00,40.50,
              46.00,41.00,34.50,32.50,38.00,34.00,
              45.00,31.50)
               )

lat_index <- NA
lon_index <- NA
var_cvrq_w <- NA
var_cvrq_c<- NA

for (j in seq(1:dim(site_df)[1])) {
 lat_index[j] = which(lat_list > site_df$lat_low[j] &  lat_list < site_df$lat_up[j])
 lon_index[j] = which(lon_list > site_df$lon_low[j] &  lon_list < site_df$lon_up[j])
 var_cvrq_w[j] <- ncvar_get(nc_wm, varid= var1$names[1], start = c(lat_index[j],lon_index[j]), 
                         count = c(1,1))
 var_cvrq_c[j] <- ncvar_get(nc_co, varid= var1$names[1], start = c(lat_index[j],lon_index[j]), 
                            count = c(1,1))
 }
  
site_df <- site_df %>%
  mutate(lat_index = lat_index) %>%
  mutate(lon_index = lon_index) %>%
  mutate(CRSQ_warm = var_cvrq_w) %>%
  mutate(CRSQ_cool = var_cvrq_c)

write.csv(site_df,paste0(output_path,"cvsq.csv"))

for (j in seq(1:14)) {
  
  loc <- site_list[[j]]
  output.p <- file.path(output_path, loc$site[1])
  
  inst_df <- readRDS(file =paste0(output.p,"/instrument_df.rds"))
  naspa_df <- readRDS(file = paste0(output.p,"/pre_naspa.rds"))
  gpcc_df <- readRDS(file = paste0(output.p,"/gpcc_df.rds"))
  naspa_df <- naspa_df %>%
    group_by(date) %>%
    summarize(precip = mean(precip, na.rm = TRUE))
  
  naspa_df <- naspa_df %>%
    mutate(year = year(date)) %>%
    mutate(month = month(date)) %>%
    mutate(variable = "predicted") %>%
    mutate(model = "naspa")%>%
    mutate(units = "mm/month") %>%
    mutate(site = loc$site[1]) %>%
    select(date,site, year, variable, model , precip, units, month)
  
size <- dim(gpcc_df)[1]

mape_df <- gpcc_df %>%
  select(date,site, precip) %>%
  rename(pr_gp = precip) %>%
  inner_join(naspa_df, by = "date") %>%
  rename(pr_naspa = precip) %>%
  mutate(
    MAE = case_when((pr_gp >= 1) ~ (abs(pr_gp - pr_naspa)/pr_gp),
                      TRUE       ~ (abs(1 - pr_naspa)/1)
    )
  )%>%
  #mutate(MAE = abs((pr_gp - pr_naspa)/pr_gp)) %>%
  mutate(nE = abs((pr_gp - pr_naspa))) %>%
  select(date, site.x,pr_gp, pr_naspa,year, month, MAE,nE) 

mape_df <- mape_df %>%
  group_by(month) %>%
  summarise(SAE = sum(MAE), nSAE = sum(pr_gp), nSE = sum(nE)) %>%
  mutate(MAPE = 100 * SAE/size) %>%
  mutate(nMAE = nSE/nSAE) %>%
  mutate(site = loc$site[1]) %>%
  ungroup()

month_level <- c("7","8","9","10","11","12","1","2","3",
                 "4","5","6")

mape_df <- mape_df %>%
  arrange(factor(month,levels = month_level)) %>%
  rbind(mape_df %>%filter(month == 7)) 
  
mape_df <- mape_df %>% mutate(mon_col = labels)


  
  if (j == 1) {
    stats_2 <- mape_df
  }else {
    stats_2 <- bind_rows(stats_2, mape_df)
  }
}
 stats_ave <- stats_2 %>%
   group_by(month) %>%
   summarise(aveMAE = mean(nMAE), aveMAPE = mean(MAPE),
             var_nMAE = var(nMAE), var_MAPE = var(MAPE),
             medMAE= median(nMAE),medMAPE = median(nMAE)) %>%
   mutate(site = "Average")
 
 stats_ave <- stats_ave %>%
 arrange(factor(month,levels = month_level)) %>%
   rbind(stats_ave %>%filter(month == 7))
 stats_ave <- stats_ave %>%
   mutate(mon_col = labels)
 stats_ave$mon_col <- as.factor(stats_ave$mon_col)
 
 
write.csv(stats_ave,paste0(output_path,"stats_allsites.csv"))
#write.csv(stats_ave,paste0(output_path,"stats_ave.csv"))
stats_2 <- stats_2 %>%
  bind_rows(stats_2, stats_ave)

read.csv(stats,paste0(output_path,"stats_allsites.csv"))

stats_2$month <- as.factor(stats_2$month)
stats_2$mon_col <- as.factor(stats_2$mon_col)

labels <- c("MJJ","JJA",
            "JAS","ASO","SON","OND","NDJ","DJF","JFM","FMA","MAM","AMJ","MJJ")
labels <- c("ASO","SON","OND","NDJ","DJF","JFM","FMA","MAM","AMJ","MJJ","JJA","JAS")

#stats_2$month <- as.factor(stats_2$month)

p <- ggplot(stats_2,aes(x = mon_col, y = nMAE)) +
     geom_line(aes(group = site), color = "grey", size = 0.8) +
     geom_smooth(stats_ave, mapping = aes(x = mon_col, y= medMAE, group = site),
             color = "steelblue", fill = "NA", alpha = 0.5,size = 1.5) +
     labs(x="Month",y="nMAE", size = 14)+
     scale_x_discrete(limits= labels,
                      guide = guide_axis(angle = 45)) +
     scale_y_continuous(limits = c(0, 0.6),breaks = c(0.0,0.2,0.4,0.6))+
     theme_classic(base_size = 20)

p

ggsave(p, file = "nMAE_line_mjj.svg",width = 7, height = 5)
ggplot(combine %>% filter(year(date) > 1901 & month(date) == 6), aes(x = date ,y= precip)) +
  geom_line(aes( color = model))

