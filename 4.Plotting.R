# *-----------------------------------------------------------------------# | PROGRAM NAME: 
# | FILE NAME: 4.Plotting.R
# | DATE: NOV.19.2021
# | CREATED BY:  Kay Sung     
# *--------------------------------------------------------------------------
# | PURPOSE: 
# | 
# | 
# | 
# *--------------------------------------------------------------------------
require(lubridate)
require(tidyverse)
require(mgcv)
require(dplyr)
require(fitdistrplus)
require(viridis)

select <- dplyr::select

data_path <- "../data"
output_path <- "../output/"

#################################################################
#######                         Plotting               ##########
#################################################################

plotm <- 1

p <- ggplot(data= modeled_df %>% filter(month %in% plotm), 
            aes(x = year, y = modGI, group = plot_model)) +
  geom_line(data= prcp_pos %>%filter(month %in% plotm), 
            aes(y = precip, group=model,colour = model), size = 0.7, alpha = 0.45) +
  geom_ribbon(data = modeled_df%>%filter(month %in% plotm), 
              aes(ymin= lowp, ymax= upp,fill = plot_model), alpha = 0.5) +
  geom_ribbon(data = modeled_df%>%filter(month %in% plotm),
              aes(ymin= lower,ymax= upper), alpha = 0.25, fill = "grey") +
  geom_line(aes(colour = plot_model), size = 1.3) +
  scale_fill_manual(values = colours) +
  scale_color_manual(values = colours) +
  labs(title = paste(plotm,loc$site[1]),
       y="3-Month ave.prcp(mm/month)",
       color ='datasets')+
  theme_classic(base_size = 20)  
p

ggsave(filename = paste0(output.p,"/",plotm,"LT.png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(output.p,"/",plotm,"LT.svg"),plot = p, width =12.5, height = 8, dpi = 300)

tyear <- seq(1900,2020)
plotm <- 10
p <- ggplot(data= modeled_df %>% filter(year %in% tyear, month %in% plotm), 
            aes(x = year, y = modGI, group = interaction(month,plot_model))) +
  geom_line(color = "black", size = 1) +
  geom_line(data= prcp_pos %>%filter(year %in% tyear,month %in% plotm), 
            aes(y = precip, group=model,colour = model)) +
  scale_color_manual(values = colours) +
  geom_ribbon(data = modeled_df%>%filter(year %in% tyear,month %in% plotm), 
              aes(ymin= lowp, ymax= upp), alpha = 0.25, fill = "blue") +
  geom_ribbon(data = modeled_df%>%filter(year %in% tyear, month %in% plotm),
              aes(ymin= lower,ymax= upper), alpha = 0.25, fill = "grey") +
  #facet_wrap(~ month, ncol = 2) +
  #scale_colour_viridis() +
  labs(title = paste(plotm,loc$site[1]),
       y="3-Month ave.prcp(mm/month)",
       color ='datasets')+
  theme_classic(base_size = 20)  
p

ggsave(filename = paste0(output.p,"/",plotm,"20LT.png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(output.p,"/",plotm,"20LT.svg"),plot = p, width =12.5, height = 8, dpi = 300)

tyear <- seq(1,2020,20)
p <- ggplot(data= modeled_df %>% filter (year %in% tyear), 
            aes(x = month, y = modGI, 
                group = interaction(year,plot_model), color = year)) +
	geom_line(size = 0.7) +
	scale_colour_viridis() +
	labs(title = loc$site[1],
       y="3-Month ave.prcp(mm/month)", size = 20)+
	scale_x_discrete(name ="Month",limits=seq(1,12)) +
	theme_classic(base_size = 20)
p

ggsave(filename = paste0(output.p, "combinedm.png"),plot = p, width =18, height = 12, dpi = 300)


### Loop through each site
for (j in seq(1, length(site_list))){

short_name_j <- site_list[j]

cat(j)
cat(short_name_j)

### Quick plot of Initial Values
month_breaks <- c(yday(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month")), 365)
month_labels <- c(as.character(month(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month"), label=TRUE)), "Jan")

### Find all files within folder
#write_output_path <- "/run/media/jhstagge/Seagate Backup Plus Drive/spibayes_paper/output/gauge_tensor"
site_folder <- file.path(write_output_path, short_name_j)

### Quick plot of Initial Values
month_breaks <- c(yday(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month")), 365)
month_labels <- c(as.character(month(seq(as.Date("1900-01-01"), as.Date("1900-12-31"), by = "1 month"), label=TRUE)), "Jan")

### Find all files within folder
#write_output_path <- "/run/media/jhstagge/Seagate Backup Plus Drive/spibayes_paper/output/gauge_tensor"
site_folder <- file.path(write_output_path, short_name_j)

estimate_nc <- nc_open(file.path(site_folder, paste0(short_name_j, "_estimate.nc")))

date_vec <- ncvar_get(estimate_nc, "time") + as.Date("1900-01-01")
year_vec <- ncvar_get(estimate_nc, "year") 
jdate_vec <- ncvar_get(estimate_nc, "jdate") 

mean_mat <- ncvar_get(estimate_nc, "mean") 
shape_mat <- ncvar_get(estimate_nc, "shape") 
scale_mat <- ncvar_get(estimate_nc, "scale") 

dist_df <- data.frame(year = year_vec, jdate = jdate_vec, shape = apply(est_shape, 1, mean), scale = apply(est_scale, 1, mean))

dist_df <- bind_cols(gridmet_pred_df,GAM_predict)


year_seq <-  c(seq(0, 2020, 200),2020)
month_seq <- seq(1,12)

plot_dist <- tibble(year = NA, month = NA,  ymin = NA, lower = NA, middle = NA, upper = NA, ymax = NA)

for (k in seq(1, length(year_seq))) {
              for(i in seq(1,12)){

              dist_temp <- dist_df %>%
                            filter(year == year_seq[k] & month %in%  month_seq[i])
              
              if(dim(dist_temp)[1] == 0) {next}
            
              dist_temp_full <- tibble(year = year_seq[k], 
                             month = month_seq[i], 
                            # date = x_date[i],
                             ymin = qgamma(0.025, shape = dist_temp$est_shape[1], scale = dist_temp$est_scale[1]),
                             lower =qgamma(0.25, shape = dist_temp$est_shape[1], scale = dist_temp$est_scale[1]),
                             middle = qgamma(0.5, shape = dist_temp$est_shape[1], scale = dist_temp$est_scale[1]),
                             upper = qgamma(0.75, shape = dist_temp$est_shape[1], scale = dist_temp$est_scale[1]),
                             ymax = qgamma(0.975, shape = dist_temp$est_shape[1], scale = dist_temp$est_scale[1]))
              
              plot_dist <- plot_dist %>%
                             bind_rows(dist_temp_full)

              rm(dist_temp_full)
              rm(dist_temp)
              }
}

plot_dist <- plot_dist %>%
              mutate(plot_month = month) %>%
              mutate(plot_month = factor(month, levels = unique(month)))
plot_dist <- plot_dist %>%
  drop_na()

p <- ggplot(plot_dist, aes(x=plot_month, lower=lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax, fill = factor(year))) %>%
              + geom_boxplot(stat = "identity", colour = "grey30", alpha = 0.9, size = 0.35) %>%
              + scale_fill_viridis(name = "Year", option = "inferno", discrete = TRUE) %>%
              + scale_y_continuous(name = "3-month prcp (mm/month)") %>%
              + scale_x_discrete(name = "Month") %>%
              + theme_classic(11, base_size = 20)


### Save Figure
ggsave(paste0(output.p,"/plot.png"), p, width =12.5, height = 4.5, dpi = 300)
ggsave(paste0(output.p,"/plot.svg"), p, width =12.5, height = 4.5)

}

