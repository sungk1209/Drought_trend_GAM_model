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

### Read in
knn <- function(loc){
  
gpcc_df <- readRDS(file = paste0(output.p,"/gpcc_df.rds"))
naspa_spi3<- readRDS(file = paste0(output.p,"/naspa_spi3.rds"))

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

#library_df[1:50,]
library_df <- library_df %>% select(-i)

n_neighbors <- 10
begin.y <- yrs[[1]][1]
end.y <- yrs[[2]][1] -1
predicted_df <- data.frame(year = NA)

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

if (sum(is.na(naspa_points))>0) {
  next
}
### Find the k closest points

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

if(is.na(predicted_df$year[1]) == TRUE){

  predicted_df <- naspa_subset_iter
  
  } else {
  predicted_df <- predicted_df %>% bind_rows(naspa_subset_iter)
  }

}

predicted_df <- predicted_df %>%
  right_join(paras, by= "month") %>%
  mutate(prob = pnorm(spi3)) %>%
  mutate(precip = qgamma(prob, shape = shape, rate = rate)) 
predicted_df$iter <- as.factor(predicted_df$iter)

saveRDS(predicted_df, file = paste0(output.p,"/pre_naspa.rds"))

return()

}
###########################################################
###                     plotting                      #####
###########################################################
predicted_df <- readRDS(paste0(output.p,"downscaled_preip.rds"))

temp <- gpcc_df %>% 
  mutate(iter = "GPCC") %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
 select(year,month,spi3,date, iter,spi5)
temp2<- predicted_df %>%
  select(year,month,spi3,date, iter,spi5)
naspa_temp <- naspa_spi3 %>%
 mutate(iter = "naspa")

temp3 <- rbind(temp, temp2, naspa_temp)

p <-ggplot(temp3 %>%filter(year >1901 & month == 7), aes(x=date, y=spi3)) +
  geom_line(aes(group = iter, colour= "grey")) + 
  geom_line(data = temp3 %>% filter(iter == "GPCC") %>%
              filter(year >1901 & month == 7), size =1) +
  geom_line(data = temp3 %>% filter(iter == "naspa") %>%
              filter(year >1901 & month == 7), size =1,aes(color =  "red")) +
  labs(title = paste0("")) +
  coord_cartesian(ylim = c(-2.5,2.5)) +
  theme_classic(base_size = 18)

p
ggsave(p,filename = paste0(output_path,loc$site[1],"ds_spi3.png"))

p <-ggplot(predicted_df %>% filter(year(date) > 2000 & year(date) < 2010), 
           aes(x=date, y=spi3)) +
  geom_line(aes(group = iter, color = 'grey')) + 
  geom_line(data = gpcc_df %>% filter(year(date) > 2000 & year(date) < 2010),
            aes(y = spi3), color = "blue" ) +
 
  labs(title = paste0("1900 - 1930",loc$site),
       color = "data") 
p

ggsave(p,filename = paste0(output_path,loc$site[1],".png"),width = 6, height = 4)

###########################################################
###                  plotting all NNs and neighbors   #####
###########################################################
years<- c(1500, 1950,2010)

for (i in c(2:3)) {

fyear <- years[i]
fdate<- seq(as.Date(paste0(fyear,"-08-01")),
                          as.Date(paste0(fyear+1,"-08-01")), 
                          by = "month")-1
fd_lab<- seq(as.Date(paste0(fyear,"-08-01")),
            as.Date(paste0(fyear+1,"-08-01")), 
            by = "2 months")-1


p <-ggplot(predicted_df %>% filter(date %in% fdate),
           aes(x=date, y= precip)) +
  geom_line(aes(group = iter, color = "10NNs")) +
  geom_line(data = gpcc_df %>% filter(date %in% fdate),
         aes(x = date, y = precip, color = "GPCC"),size = 1.2, linetype = "dashed") +
  stat_summary(aes(y = precip), fun = "mean", size = 1.2, geom = "line", na.rm = TRUE) +
  labs(title = paste0(loc$site[1],fyear),
       y = "3-months Prcp.(mm/Month)",
       X = "Month",
       color = "data") +
 # scale_x_date(date_breaks = "2 months", date_minor_breaks = "months", date_labels = "%b") +
  scale_x_date(breaks = fd_lab, date_labels = "%b") +
  #scale_y_continuous(limits = c(0,400)) +
  scale_color_discrete() + 
  theme_classic(base_size = 18) 
p
ggsave(p,filename = paste0(output.p,"/",fyear,"_ds.png"),width = 6, height = 4 )
ggsave(p,filename = paste0(output.p,"/",fyear,"_ds.svg"),width = 6, height = 4 )

}

p <-ggplot(predicted_df %>% filter(year(date) %in% fyear),
           aes(x=date)) +
  geom_line(aes(y = precip, group = iter, color = iter, alpha = 0.6), size = 0.8) +
  geom_line(data = gpcc_df %>% filter(year(date) %in% fyear),
            aes(x = date, y = precip, color = "GPCC") ) +
  scale_color_manual(values = c("orangered2",rep("grey76",10))) +
  labs(title = paste0("3months ave.mean and all NNs"),
       y = "3-months Prcp.(mm/month)") +
  theme_classic(base_size = 12)

p <- p + stat_summary(aes(y = precip),fun = "mean", colour = "green4", geom = "line")
p



p <- ggplot(predicted_df %>% 
              filter(year(date) > 1902 & year(date) < 2020 & month == 4),
            aes(x = date, y= precip)) + 
  geom_line(aes(group = iter), color = "skyblue3") +
  #geom_point(naspa_df %>%
  #            filter(year(date) > 900 & year(date) < 1000), 
  #           mapping = aes(y = precip, color = "naspa")) +
  geom_line(data = gpcc_df %>%filter(year(date) > 1902, month(date) == 4), 
            aes(y = precip,color = "gpcc"),size = 0.8) +
  labs(title ="MJJ precip, naspa and predicted",
       y = "3months ave.prcp(mm/m)") +
  scale_color_brewer(palette = "Set1") +
  theme_classic(base_size = 18)
p <- p + stat_summary(fun = "mean", colour = "black", size = 1.5, geom = "line")
p
ggsave(p,filename = paste0(output_path,loc$site[1],"10.png"),width = 6, height = 4 )


==