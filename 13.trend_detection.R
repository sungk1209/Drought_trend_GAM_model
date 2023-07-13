
# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: gpcc_downscaling.R
# | DATE: Jul.06.2023
# | CREATED BY: Kay Sung     
# *--------------------------------------------------------------------------
# | PURPOSE: 
# | 
# | 
# | 
# *--------------------------------------------------------------------------

require(lubridate)
require(tidyverse)
require(dplyr)
require(gratia)
require(mgcv)
require(MASS)

select <- dplyr::select

data_path <- "../data/"
output_path <- "../output/"

loc <- site_list[[j]]
output.p <- file.path(output_path, loc$site[1])
gam_fit <- readRDS( file =paste0(output.p,"/model.rds"))

fd <- derivSimulCI(gam_fit, samples = 10000, n = 200)














