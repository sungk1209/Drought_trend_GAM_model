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

site_list <- list()
site_list[1] <-list(data.frame(site="okc_OK",
                               lon=c(-98,-97.5),
                               lat= c(35.0,35.5)))
site_list[2] <-list(data.frame(site="phx_AZ",
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


for (j in seq(7, length(site_list))){
  loc <- site_list[[j]]
  dir.create(file.path(output_path,loc$site[1]), recursive = FALSE)
  output.p <- file.path(output_path, loc$site[1])
  inst_download(loc)
  yrs<- naspa_download(loc)
  knn_result <- knn(loc)
  GI_model(loc)
  plot_model(loc)
  }

