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
require(viridis)

select <- dplyr::select

data_path <- "../data"
output_path <- "../output/"

#################################################################
#######                         Plotting               ##########
#################################################################
plot_model <- function(loc) {

prcp_pos<- readRDS(file = paste0(output.p,"/prcp_df.rds"))
modeled_df<- readRDS(file = paste0(output.p,"/modeled_df.rds"))
plot_ms <- c(1,4,7,10)
colours = c("naspa" = "green4", "cru" = "#1A3FD4", "gridmet" = "orangered3") 
for (i in 1:length(plot_ms)) {
  
  plotm <- plot_ms[i]
 
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


ggsave(filename = paste0(output.p,"/",plotm,"_LT2.png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(output.p,"/",plotm,"_LT2.svg"),plot = p, width =12.5, height = 8, dpi = 300)
}
return()
}

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

return()
}


# tyear <- seq(1,2020,20)
# p <- ggplot(data= modeled_df %>% filter (year %in% tyear),
#             aes(x = month, y = modGI,
#                 group = interaction(year,plot_model), color = year)) +
# 	geom_line(size = 0.7) +
# 	scale_colour_viridis() +
# 	labs(title = loc$site[1],
#        y="3-Month ave.prcp(mm/month)", size = 20)+
# 	scale_x_discrete(name ="Month",limits=seq(1,12)) +
# 	theme_classic(base_size = 20)
# 
# 
# ggsave(filename = paste0(output.p, "combinedm.png"),plot = p, width =18, height = 12, dpi = 300)

# 
# tyear <- seq(1900,2020)
# for (i in 1:length(plot_ms)) {
#   
#   plotm <- plot_ms[i]
# p <- ggplot(data= modeled_df %>% filter(year %in% tyear, month %in% plotm), 
#             aes(x = year, y = modGI, group = interaction(month,plot_model))) +
#   geom_line(color = "black", size = 1.2) +
#   geom_line(data= prcp_pos %>%filter(year %in% tyear,month %in% plotm), 
#             aes(y = precip, group=model,colour = model)) +
#   scale_color_manual(values = colours) +
#   geom_ribbon(data = modeled_df%>%filter(year %in% tyear,month %in% plotm), 
#               aes(ymin= lowp, ymax= upp), alpha = 0.25, fill = "blue") +
#   geom_ribbon(data = modeled_df%>%filter(year %in% tyear, month %in% plotm),
#               aes(ymin= lower,ymax= upper), alpha = 0.25, fill = "grey") +
#   #facet_wrap(~ month, ncol = 2) +
#   #scale_colour_viridis() +
#   labs(title = paste(plotm,loc$site[1]),
#        y="3-Month ave.prcp(mm/month)",
#        color ='datasets')+
#   theme_classic(base_size = 20)  
# p 

#ggsave(filename = paste0(output.p,"/",plotm,"20LT.png"),plot = p, width =12.5, height = 8, dpi = 300)
#ggsave(filename = paste0(output.p,"/",plotm,"20LT.svg"),plot = p, width =12.5, height = 8, dpi = 300)

#}
# p <- ggplot(data= modeled_df %>% filter(year > 1800 & month %in% plotm), 
            #             aes(x = year, y = modGI, group = plot_model)) +
            #   geom_line(data= prcp_pos %>%filter(year > 1800 & month %in% plotm), 
            #             aes(y = precip, group=model,colour = model), size = 0.7, alpha = 0.45) +
            #   geom_ribbon(data = modeled_df%>%filter(year > 1800 & month %in% plotm), 
            #               aes(ymin= lowp, ymax= upp,fill = plot_model), alpha = 0.5) +
            #   geom_ribbon(data = modeled_df%>%filter(year > 1800 &month %in% plotm),
            #               aes(ymin= lower,ymax= upper), alpha = 0.25, fill = "grey") +
            #   geom_line(aes(colour = plot_model), size = 1.3) +
            #   scale_fill_manual(values = colours) +
            #   scale_color_manual(values = colours) +
            #   labs(title = paste(plotm,loc$site[1]),
            #        y="3-Month ave.prcp(mm/month)",
            #        color ='datasets')+
            #   theme_classic(base_size = 20)  
            # 
            # ggsave(filename = paste0(output.p,"/",plotm,"18002.png"),plot = p, width =12.5, height = 8, dpi = 300)
            # ggsave(filename = paste0(output.p,"/",plotm,"18002.svg"),plot = p, width =12.5, height = 8, dpi = 300)
            
            