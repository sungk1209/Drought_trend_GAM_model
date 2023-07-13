
# *-----------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: 12.revision_data_bias_1900.R
# | DATE: Jul.04.2023
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

for (j in c(1:14)){
  
  loc <- site_list[[j]]
  
  
  output.p <- file.path(output_path, loc$site[1])
  
gam_fit <- readRDS( file =paste0(output.p,"/model.rds"))

prcp_pos<- readRDS(file = paste0(output.p,"/prcp_df.rds"))
### make a new combined_df with Gridmet as model

min_gridm <- min((prcp_pos %>% filter(model == "gridmet"))$year) 
max_gridm <- max((prcp_pos %>% filter(model == "gridmet"))$year)
min_cru <- min((prcp_pos %>% filter(model == "cru"))$year) 
max_cru <- max((prcp_pos %>% filter(model == "cru"))$year)
min_naspa <- min((prcp_pos %>% filter(model == "naspa"))$year) 
max_naspa <- max((prcp_pos %>% filter(model == "naspa"))$year)
#min_gcm <- min((prcp_pos %>% filter(model == "GFDL-CM4"))$year) 
#max_gcm <- max((prcp_pos %>% filter(model == "GFDL-CM4"))$year)

pred_df <- data.frame(expand.grid(year = seq(min_naspa,max_naspa),month = seq(1,12)), model = "naspa", plot_model = "naspa") %>%
  bind_rows(data.frame(expand.grid(year = seq(min_gridm,max_gridm),month = seq(1,12)),model = "gridmet", plot_model = "gridmet")) %>%
  bind_rows(data.frame(expand.grid(year = seq(min_cru,max_cru),month = seq(1,12)), model = "cru", plot_model = "cru")) 

GAM_predict <- predict(gam_fit, newdata = pred_df, se.fit = TRUE, type = "response")

GAM_predict  <- GAM_predict %>%
  data.frame() %>%
  as_tibble() %>%
  rename(est_mean = 1) %>%
  rename(est_shape = 2) 

GAM_predict <- GAM_predict %>%
  mutate(est_shape = 1/exp(est_shape)) %>%
  mutate(est_scale = est_mean/est_shape) %>%
  mutate(est_sd = (est_shape * est_scale^2)^0.5)

modeled_df <- transform(pred_df,
                        modGI = GAM_predict$est_mean,
                        lower = qgamma(0.025, shape = GAM_predict$est_shape,
                                       scale = GAM_predict$est_scale),
                        upper = qgamma(0.975, shape = GAM_predict$est_shape,
                                       scale = GAM_predict$est_scale),
                        # lowcrit = GAM_predict$est_mean + crit * GAM_predict$se.fit.1,
                        #  upcrit = GAM_predict$est_mean - crit * GAM_predict$se.fit.1,
                        lowp = GAM_predict$est_mean -  2* GAM_predict$se.fit.1,
                        upp = GAM_predict$est_mean +  2* GAM_predict$se.fit.1,
                        est_sd = GAM_predict$est_sd)

modeled_df <- modeled_df %>%
  mutate(date = as.Date(paste0(year,"-",month,"-01")))

title <- c("NDJ", "DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND")

modeled_df$label <- factor(modeled_df$month, labels = title)
prcp_pos$label <- factor(prcp_pos$month, labels = title)

plotm <- c(1,4,7,10)
colours = c("naspa" = "green4", "cru" = "#1A3FD4", "gridmet" = "orangered3") 
  
  p <- ggplot(data= modeled_df %>% filter(year > 1880 & month %in% plotm), 
              aes(x = year, y = modGI, group = plot_model)) +
    geom_line(data= prcp_pos %>%filter(year > 1880 &month %in% plotm), 
              aes(y = precip, group=model,colour = model), size = 0.7, alpha = 0.20) +
   # geom_ribbon(data = modeled_df%>%filter(year > 1880 & month %in% plotm), 
  #              aes(ymin= lowp, ymax= upp,fill = plot_model), alpha = 0.5) +
    geom_ribbon(data = modeled_df%>%filter(year > 1880 & month %in% plotm),
                aes(ymin= lower,ymax= upper, fill = model), alpha = 0.20) +
    geom_line(aes(colour = plot_model), size = 1) +
    scale_fill_manual(values = colours) +
    scale_color_manual(values = colours) +
    facet_wrap(~label, nrow = 1, scales = "free") +
    labs(y="3-Month ave.prcp(mm/month)",
         color ='datasets', fill = 'datasets')+
    theme_classic(base_size = 20)  +
    theme(legend.position = "none")
  
  p

  ggsave(p, filename = paste0(output.p,"/_sep_1900.png"), width =12, height = 4, dpi = 300)
}
