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
output_path <- "../output"
output.p <- file.path(output_path, loc$site[1])

### Read in

gpcc_df <- readRDS(file = paste0(output.p,"/gpcc_df.rds"))

paras<- gpcc_df %>% select(date,shape3, rate3) %>%
  group_by(month(date)) %>%
  summarise(shape = mean(shape3), rate = mean(rate3)) %>%
  rename(month =`month(date)`)

n_library <- length(gpcc_df$precip)
### Convert library into a dataframe
library_df <- data.frame(i = seq(1,n_library), spi_thisjuly = as.numeric(gpcc_df$spi3)) %>%
  mutate(spi_nextmay = dplyr::lead(gpcc_df$spi5, n=9,  default = NA))%>% 
  mutate(spi_nextjuly = dplyr::lead(gpcc_df$spi3, n=12,  default = NA))%>% 
  drop_na()

library_df[1:50,]
library_df <- library_df %>% select(-i)

n_neighbors <- 10

for (year_i in c(begin.y:end.y)){ 

#### Eventually put this in a loop through years
date_subset <- seq(as.Date(paste0(year_i,"-08-01")),
                   as.Date(paste0(year_i+1,"-08-01")), 
                     by = "month")-1
  
naspa_subset <- naspa_spi3 %>%
	filter(date %in% date_subset)

naspa_points <- data.frame(spi_thisjuly = naspa_subset$spi3[[1]],
                           spi_nextapr = naspa_subset$spi5[[10]],
                           spi_nextjuly = naspa_subset$spi3[[13]])

#if (sum(is.na(naspa_points))>0) {
#  next
#}
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
  #smoothed <- loess(spi3 ~ as.numeric(date),data = naspa_subset_iter)
  #predicted_k <- data.frame(date= date_subset,
  #                         spi3 = predict(smoothed, newdata = naspa_subset$date))
if(year_i == begin.y){
  predicted_df <- naspa_subset_iter
  
  } else {
  predicted_df <- predicted_df %>% bind_rows(naspa_subset_iter)
  }

}
#naspa_subset_iter_wide <- naspa_subset_iter %>% 
#	pivot_wider(names_from = iter, values_from = spi3) 

#naspa_subset_iter_wide
predicted_df <- predicted_df %>%
  mutate(month = month(date)) %>%
  right_join(paras, by= "month") %>%
  mutate(prob = pnorm(spi3)) %>%
  mutate(precip = qgamma(prob, shape = shape, rate = rate)) 
predicted_df$iter <- as.factor(predicted_df$iter)

saveRDS(predicted_df, file = paste0(output.p,"/pre_naspa.rds"))

###########################################################
###                     plotting                      #####
###########################################################
temp <- gpcc_df %>% 
  filter(date %in% date_subset) %>%
  mutate(iter = "GPCC") %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
 select(year,month,spi3,date, iter,spi5)

temp2 <- rbind(naspa_subset_iter,temp)

p <-ggplot(temp2, aes(x=date, y=spi3)) +
  geom_point(aes(colour= iter)) + 
    geom_smooth(alpha = 0.5) +
  geom_line(data = temp2 %>% filter(iter == "GPCC"), size = 2) +
  labs(title = paste0(year_i,"years with, OKC,OK")) +
  coord_cartesian(ylim = c(-2.5,2.5)) +
  theme_classic(base_size = 18)

p
ggsave(p,filename = paste0(output_path,year_i,loc$site[1],"1930_5OKC.png"))

p <-ggplot(predicted_df %>% filter(year(date) > 1900 & year(date) < 1930), 
           aes(x=date, y=spi3)) +
  geom_line(aes(colour = "black")) + 
  geom_line(data = gpcc_df %>% filter(year(date) > 1900 & year(date) < 1930),
            aes(y = spi3), color = "blue" ) +
 
  labs(title = paste0("1900 - 1930",loc$site),
       color = "data") 
p

ggsave(p,filename = paste0(output_path,loc$site[1],".png"),width = 6, height = 4)

###########################################################
###                  plotting all NNs and neighbors   #####
###########################################################
tyear <- seq(1930,1950)
p <-ggplot(predicted_df %>% filter(year(date) %in% tyear),
           aes(x=date)) +
  geom_line(aes(y = precip, group = iter, color = iter, alpha = 0.6), size = 0.8) +
  geom_line(data = gpcc_df %>% filter(year(date) %in% tyear),
            aes(x = date, y = precip, color = "GPCC") ) +
  scale_color_manual(values = c("orangered2",rep("grey76",10))) +
  labs(title = paste0("3months ave.mean and all NNs"),
       y = "3-months Prcp.(mm/month)",
       color = "data") +
  theme_classic(base_size = 20)

p <- p + stat_summary(aes(y = precip),fun = "mean", colour = "green4", geom = "line")
p

ggsave(p,filename = paste0(output.p,"/",tyear[1],"ds_gpcc.png"),width = 10, height = 6 )
ggsave(p,filename = paste0(output.p,"/",tyear[1],"ds_gpcc.svg"),width = 10, height = 6 )

tyear <- seq(1900,2020)
tmonth <- 10

p <-ggplot(predicted_df %>% filter(year(date) %in% tyear & month(date) == tmonth),
           aes(x=date)) +
  geom_line(aes(y = precip, group = iter, color = iter, alpha = 0.6), size = 0.8) +
  geom_line(data = gpcc_df %>% filter(year(date) %in% tyear & month(date) == tmonth),
            aes(x = date, y = precip, color = "GPCC") ) +
  scale_color_manual(values = c("orangered2",rep("grey76",10))) +
  labs(title = paste0("3months ave.mean and all NNs"),
       y = "3-months Prcp.(mm/month)") +
  theme_classic(base_size = 12)

p <- p + stat_summary(aes(y = precip),fun = "mean", colour = "green4", geom = "line")
p

ggsave(p,filename = paste0(output.p,"/",tmonth,"month_ds.png"),width = 6, height = 4)
ggsave(p,filename = paste0(output.p,"/",tmonth,"month_ds.svg"),width = 6, height = 4)


