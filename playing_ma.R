require(tidyverse)
require(dplyr)
select <- dplyr::select

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

n_days <- 12
### SPI-3 3 months equals 3 months
n_roll <- 3

### Use a 2 month moving average MA(2) with innovations of sqrt(3)/n_roll, which produces N(0,1)
### Technically, the moving average would be on precip, not SPI, but this is a  close approximation
innov_c <- rnorm(n_days, 0, sqrt(n_roll))

### Generate SPI series using an MA(2) model with coef of 1. Need to remove the first 91 values
sim_spi <- arima.sim(list(order = c(0,0,(n_roll - 1)), 
                          ma = rep(1,(n_roll -1))), 
                     n = n_days, innov=innov_c/n_roll)

sim_spi[seq(1,(n_roll-1))] <- NA

#sim_spi <- sim_spi %>%
#	mutate(innov = innov_c/n_roll) # %>%                            
#	mutate(spi = c(sim_spi))

#ggplot(spi_true, aes(x=date, y=spi)) + geom_line() + theme_classic()

year_i <- 2000

#### Brute force it
n_library <- 1000000
### Create a library of Moving Average
library_ts <- arima.sim(n = 1000000, list( ma = rep(1,(n_roll -1))), sd = sqrt(3)/n_roll)

### Check
mean(library_ts)
sd(library_ts)

### Convert library into a dataframe
library_df <- data.frame(i = seq(1,n_library), spi_thisjuly = as.numeric(library_ts)) %>%
	mutate(spi_nextjuly = dplyr::lead(spi_thisjuly, n=12,  default = NA))%>% 
	drop_na()

library_df[1:50,]

#### Eventually put this in a loop through years
year_i <- 2000

date_subset <- seq(as.Date(paste0(year_i,"-08-01")), as.Date(paste0(year_i+1,"-08-01")), 
                   by = "month")-1

naspa_subset <- naspa_spi3 %>%
	filter(date %in% date_subset)

naspa_points <- data.frame(spi_thisjuly = naspa_subset$spi3[[1]], spi_nextjuly = naspa_subset$spi3[[13]])

### Find the k closest points
n_neighbors <- 20

library(RANN)
library_df <- library_df %>% select(spi_thisjuly, spi_nextjuly)
closest <- nn2(data= library_df,
               query = naspa_points, k=n_neighbors)

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

naspa_subset_iter_wide <- naspa_subset_iter %>% 
	pivot_wider(names_from = iter, values_from = spi3)

naspa_subset_iter_wide

ggplot(naspa_subset_iter, aes(x=date, y=spi3)) + geom_line(aes(colour= iter)) + geom_smooth( alpha = 0.5)

