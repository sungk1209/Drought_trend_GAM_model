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

title <- c("NDJ", "DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND")

for (j in c(1:14)){
  
  loc <- site_list[[j]]
  
  output.p <- file.path(output_path, loc$site[1])
  
  #Load neighbors data with GPCC  
  gpcc_df <- readRDS(file = paste0(output.p,"/gpcc_df.rds"))  
  gpcc_df <- gpcc_df %>% mutate(month = month(date))
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
  
  
  n_neighbors <- i
 
  begin.y <- min(naspa_spi3$year)
  end.y <- max(naspa_spi3$year) - 2

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

  dsnaspa_df <- predicted_df %>%
    group_by(date) %>%
    summarize(precip = mean(precip, na.rm = TRUE))
  
  dsnaspa_df <- dsnaspa_df %>%
    mutate(year = year(date)) %>%
    mutate(month = month(date)) %>%
    mutate(variable = "predicted") %>%
    mutate(model = "naspa")%>%
    mutate(units = "mm/month") %>%
    # mutate(site = "OKC,OK") %>%
    select(date,year, variable, model , precip, units, month)
  
  
#saveRDS(predicted_df, file = paste0(output_path,loc$site[1],"15_nn.rds"))

###########Plot for long-time period for one season#################
  
  size <- dim(gpcc_df)[1]
  mape_df <- gpcc_df %>%
    select(date,site, precip) %>%
    rename(pr_gp = precip) %>%
    inner_join(dsnaspa_df, by = "date") %>%
    rename(pr_naspa = precip) %>%
    mutate(MAE = abs((pr_gp - pr_naspa)/pr_gp)) %>%
    mutate(nE = abs((pr_gp - pr_naspa))) %>%
    select(date, pr_gp, pr_naspa,year, month, MAE,nE) 
  mape_df <- mape_df %>%
    group_by(month) %>%
    summarise(nSAE = sum(pr_gp), nSE = sum(nE)) %>%
    mutate(nMAE = nSE/nSAE) %>%
    mutate(iterate = i) %>%
    mutate(loc = loc$site[1]) %>%
    ungroup()
  
  if (i == 4){
    nmae_df <- mape_df
  }else{nmae_df <- rbind(nmae_df, mape_df)}
  
  }
  
  if(j == 1){nmae_all <- nmae_df}else{
  nmae_all <- rbind(nmae_all, nmae_df)}
  
  saveRDS(nmae_all, file = paste0(output_path,"/nmae_all.rds"))

plabel <- c("NDJ",  "JFM", "MAM",  "MJJ",  "JAS", "SON")

  p <- ggplot(data=nmae_all ,
              aes(x=month, y = nMAE, group = iterate, color = iterate)) +
    facet_wrap(~loc, scale = "free", ncol = 3) +
    geom_line() +
    geom_line(data= nmae_all %>% filter(iterate == 10), color = "red", linetype = "dashed") +
    labs(x="Month",y="nMAE", color = "Num of \nNeighbors") +
    scale_y_continuous(limits = c(0,1.1)) +
    scale_x_continuous(breaks = seq(1,12,2), labels = plabel) +
    scale_color_viridis() +
    theme_classic()

  p

  ggsave(p, file = paste0(output_path,"nMAE.png"), height = 12, width = 8)
  
  nmae_all$label <- factor(nmae_all$month, labels = title)
  
  p <- ggplot(data=nmae_all,
              aes(x=iterate, y = nMAE, group= loc, color = loc)) + #, group = iterate, color = iterate)) +
    facet_wrap(~label, scale = "free", ncol = 3) +
    geom_line() +
    geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
    #  geom_line(data= nmae_all %>% filter(iterate == 10), color = "red", linetype = "dashed") +
      labs(title = paste0("nMAE"),
           x="Num of Neighbors",y="nMAE", color = "Location") +
    scale_y_continuous(limits = c(0,1.2)) +
    scale_x_continuous(breaks = seq(0,30,5)) +
    theme_classic()
  
  p
  
  ggsave(p, file = paste0(output_path,"nMAE_iterate.png"))
  

  # ##################################################################################### 
  
  
  predicted_df$label <- factor(predicted_df$month, labels = title)
  gpcc_df$label <- factor(gpcc_df$month, labels = title)
  compare_df$label <- factor(compare_df$month, labels = title)
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
   

  