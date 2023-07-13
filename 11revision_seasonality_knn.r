# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: gpcc_downscaling.R
# | DATE: Jul.02.2023
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

data_path <- "../data/"
output_path <- "../output/"

### Read in

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
  
  n_neighbors <- 15
  #begin.y <- yrs[[1]][1]
  #end.y <- yrs[[2]][1] -1
  begin.y <- min(naspa_spi3$year)
  end.y <- max(naspa_spi3$year) -1

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
    mutate(month = month(date)) %>%
    right_join(paras, by= "month") %>%
    mutate(prob = pnorm(spi3)) %>%
    mutate(precip = qgamma(prob, shape = shape, rate = rate)) 
  
  predicted_df$iter <- as.factor(predicted_df$iter)


#saveRDS(predicted_df, file = paste0(output_path,loc$site[1],"15_nn.rds"))


###########################################################
###                     plotting                      #####
###########################################################
 
  
   title <- c("NDJ", "DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND")
  
  for (j in c(1:14)){

  loc <- site_list[[j]]
  
  output.p <- file.path(output_path, loc$site[1])
 
#Load neighbors data with GPCC  
  gpcc_df <- readRDS(file = paste0(output.p,"/gpcc_df.rds"))  
  gpcc_df <- gpcc_df %>% mutate(month = month(date))
  
  predicted_df <- readRDS(file = paste0(output.p,"/pre_naspa.rds"))
  dsnaspa_df <- predicted_df %>%
    group_by(date) %>%
    summarize(precip = mean(precip, na.rm = TRUE))

  dsnaspa_df <- dsnaspa_df %>%
    mutate(model = "naspa") %>%
    filter(year(date) > 1901)

  compare_df <- dsnaspa_df %>% rbind(data.frame(date= gpcc_df$date, precip = gpcc_df$precip, model = "GPCC"))

  compare_df <- pivot_wider(compare_df,names_from = model, values_from = precip )
  compare_df <- compare_df %>% mutate(month = month(date))

###########Plot for long-time period for one season#################
  
  predicted_df$label <- factor(predicted_df$month, labels = title)
  gpcc_df$label <- factor(gpcc_df$month, labels = title)
  compare_df$label <- factor(compare_df$month, labels = title)

#  
    p <- ggplot(predicted_df %>% 
                  filter(year(date) > 1980 & year(date) < 2013 & month %in% tmonth ),
                aes(x = date, y= precip)) + 
      geom_line(aes(group = iter, color = "Neighbors")) +
      geom_line(data = gpcc_df %>% filter(year(date) > 1980 & year(date) < 2013 & month(date) %in% tmonth), 
                aes(y = precip,color = "GPCC"),size = 1.0) +
      facet_wrap(~label, nrow = 1) +
      labs(y = "Prcp.(mm/Month)",x = "Year", 
           color = "Dataset") +
      scale_color_manual(values = c("Neighbors" = "skyblue3", "GPCC" = "red", "Average"= "black")) +
      theme_bw(base_size = 18) +
      theme(legend.position = "none")
    
    p <- p + stat_summary(fun = "mean", aes(colour = "Average"), size = 1.0, geom = "line")
    
    
    p
    
  ggsave(p,filename = paste0(output.p,"/knn_4mons",".png"), width = 12, height = 4 )
   
####################create gpcc VS dsnaspa plot
  
  
  p <- ggplot(compare_df %>% 
                filter( month(date) %in% tmonth ),
              aes(x = naspa, y= GPCC)) + 
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~label, nrow = 1, scales = "free") +
    coord_cartesian(xlim = c(0, NA),
                    ylim = c(0, NA)) +
    labs(#title =paste0(title[tmonth]," ", loc$site[1]),
         y = "Obs. (GPCC)", x = "Downscaled (NASPA)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "none")
  p
  
  ggsave(p,filename = paste0(output.p,"/obs_scatter_4mons.png"), width = 12, height = 4)
    
######################Plot for comparing seasonality
  
  summary_data <- compare_df %>%
    group_by(month(date)) %>%
    drop_na()
  
  summary_data <- summary_data %>%
    summarize(avg_naspa = mean(naspa),
              lower_naspa = mean(naspa) - 1.96 * sd(naspa),
              upper_naspa = mean(naspa) + 1.96 * sd(naspa),
              
              avg_gpcc = mean(GPCC),
              lower_gpcc = mean(GPCC) - 1.96 * sd(GPCC),
              upper_gpcc = mean(GPCC) + 1.96 * sd(GPCC))
  
  # Create a scatter plot with average line and shaded confidence interval ribbon
  
  p <- ggplot(compare_df) +
    geom_point(aes(x = month(date), y = GPCC, color = "GPCC"), alpha = 0.3) +
   
    geom_ribbon(data = summary_data,
                aes(x = `month(date)`, ymin = lower_gpcc, ymax = upper_gpcc, fill = "GPCC"),
                alpha = 0.3) +
    geom_point(aes(x = month(date), y = naspa, color = "NASPA"), alpha = 0.3) +
    geom_ribbon(data = summary_data,
                aes(x =`month(date)`, ymin = lower_naspa, ymax = upper_naspa, fill = "NASPA"),
                alpha = 0.3) +
    geom_line(data = summary_data, aes(x = `month(date)`, y = avg_gpcc, color = "GPCC"),
              size = 1) +
    geom_line(data = summary_data, aes(x = `month(date)`, y = avg_naspa, color= "NASPA"),
              size = 1) +
  #  scale_color_manual(values = c("GPCC" = "#00AFBB", "NASPA"= "#FC4E07")) +
  #  scale_fill_manual(values = c("GPCC" = "#00AFBB", "NASPA"= "#FC4E07")) +
    scale_x_continuous(breaks = c(1:12), labels = title,name = "Months") +
    scale_y_continuous(name = "Precipitation (mm/Month)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  p
  
  ggsave(p,filename = paste0(output.p,"/seasonality",".png"), width = 3.75, height = 3.3)
  
  }
  