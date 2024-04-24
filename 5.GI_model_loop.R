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
output_path <- "../output"
select <- dplyr::select

site_list <- list()
site_list[1] <-list(data.frame(site="Okc_OK",
                               lon=c(-98,-97.5),
                               lat= c(35.0,35.5)))
site_list[2] <-list(data.frame(site="Phx_AZ",
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


#for (j in seq(1, length(site_list))){
for (j in seq(12,14)){
  loc <- site_list[[j]]
  dir.create(file.path(output_path,loc$site[1]), recursive = FALSE)
  output.p <- file.path(output_path, loc$site[1])
  #inst_download(loc)
  #yrs<- naspa_download(loc)
  #knn_result <- knn(loc)
  #GI_model(loc)
  plot_model2(loc)
  }
####################################################################
a <- site_list[[1]]
for (i in seq(1:10)){
  a <- bind_rows(a, site_list[[i+1]])
}

site_df <- a[seq(1,22,2),]

require(sf)
require(maps)
require(spData)
usa = st_as_sf(map('usa', plot = FALSE, fill = TRUE))
usa <- st_transform(usa)
us_states4326 = st_transform(us_states, 4326)

p <- ggplot() + geom_sf(data = us_states4326) + 
  geom_point(data = site_df, aes(x = lon, y = lat),
                           size = 3, shape = 25, fill = "Red") + 
  geom_label(data = site_df, aes(x = lon, y = lat), label= site_df$site,
    nudge_x = 0.25, nudge_y = 1.2) +
  coord_sf(expand = TRUE) +
  theme(
    text = element_text(size = 20),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude")
 
p

ggsave(filename = paste0(output_path,"/map.png"),plot = p, width =10, height = 8, dpi = 300)
ggsave(filename = paste0(output_path,"/map.svg"),plot = p, width =10, height = 8, dpi = 300)

##########################################################################

for (j in seq(1, length(site_list))){
  loc <- site_list[[j]]
  output.p <- file.path(output_path, loc$site[1])

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
 # mutate(site = loc$site[1]) %>%
  select(date, year, variable, model , precip, units, month)

####################################################################
##########     Compare downscaled Naspa and GPCC
####################################################################
size <- dim(gpcc_df)[1]
mape_df <- gpcc_df %>%
  select(date,site, precip) %>%
  rename(pr_gp = precip) %>%
  inner_join(naspa_df, by = "date") %>%
  rename(pr_naspa = precip) %>%
  mutate(MAE = abs((pr_gp - pr_naspa)/pr_gp)) %>%
  mutate(nE = abs((pr_gp - pr_naspa))) %>%
  select(date, site,pr_gp, pr_naspa,year, month, MAE,nE) 

mape_df <- mape_df %>%
  group_by(month) %>%
  summarise(SAE = sum(MAE), nSAE = sum(pr_gp), nSE = sum(nE)) %>%
  mutate(MAPE = 100 * SAE/size) %>%
  mutate(nMAE = nSE/nSAE) %>%
  mutate(site = loc$site[1]) %>%
  ungroup()

if(j == 1){
  nmaes <- mape_df } 
  else{
    nmaes <- rbind(nmaes,mape_df)}
}

nmaes$month <- as.factor(nmaes$month)

p <- ggplot(nmaes, aes(fill=month, y = nMAE, x=site)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = T, option = "viridis") +
  scale_y_continuous(breaks =  c(0, 0.2, 0.4,0.6,0.8,1.0)) + 
  theme_classic(base_size = 17)
p

