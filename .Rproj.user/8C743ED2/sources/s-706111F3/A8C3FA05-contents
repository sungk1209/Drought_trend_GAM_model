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
supp.labs <- c("NDJ","FMA","MJJ","ASO")
names(supp.labs) <- plot_ms

#for (i in 1:length(plot_ms)) {
  
  plotm <- plot_ms
 
p <- ggplot(data= modeled_df %>% filter(year > 1400 & month %in% plotm), 
            aes(x = year, y = modGI, group = plot_model)) +
   geom_ribbon(data = modeled_df%>% filter(year > 1400 & month %in% plotm),
              aes(ymin= lowp, ymax = upp), alpha = 0.3) +
  geom_ribbon(data = modeled_df%>% filter(year > 1400 & month %in% plotm),
              aes(ymin= lower,ymax= upper), alpha = 0.25, fill = "grey") +
  geom_line(data= prcp_pos%>% filter(year > 1400 & month %in% plotm) , 
            aes(y = precip, group = model, colour = model), size = 0.7, alpha = 0.6) +
  geom_line(aes(color = "Modeled Trend"), color = "black",size = 1.4, alpha = 0.7) +
  scale_fill_manual() +
  scale_colour_manual(values = colours) +
  scale_x_continuous(expand= c(0,0), limits = c(1400, NA)) +
  scale_y_continuous(expand= c(0,0), limits = c(0, NA)) +
  facet_wrap(~month, ncol = 1,
             labeller = labeller(month = supp.labs)) +
  labs(#title = paste(plotm,loc$site[1]),
       y="3-Month ave.prcp(mm/Month)",
       color ='datasets')+
  theme(text = element_text(colour = 'black',size = 17),
        strip.text = element_blank(),
        #strip.text = element_text(colour = 'black',size =15),
        axis.ticks.length=unit(0.3, "lines"),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank())
  
p

ggsave(filename = paste0(output.p,"/LT.png"),plot = p, width =6, height = 10, dpi = 300)
ggsave(filename = paste0(output.p,"/LT.svg"),plot = p, width =6, height = 10, dpi = 300)

#ggsave(filename = paste0(output.p,"/",plotm,"_LT2.png"),plot = p, width =10, height = 8, dpi = 300)
#ggsave(filename = paste0(output.p,"/",plotm,"_LT2.svg"),plot = p, width =10, height = 8, dpi = 300)
}
return()

####################################################################################
#########Plot seasonally varied long-term change             #######################
####################################################################################

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


p### Save Figure
ggsave(paste0(output.p,"/plot.png"), p, width =12.5, height = 4.5, dpi = 300)
ggsave(paste0(output.p,"/plot.svg"), p, width =12.5, height = 4.5)

return()
}

##############################################################################
################plot for all region
#############################################################################
plot_model2 <- function(loc) {
  
  prcp_pos<- readRDS(file = paste0(output.p,"/prcp_df.rds"))
  modeled_df<- readRDS(file = paste0(output.p,"/modeled_df.rds"))
  plotm <- c(1,4,7,10)
  supp.labs <- c("NDJ","FMA","MJJ","ASO")
  names(supp.labs) <- plotm

#plotm <- plot_ms

modeled_df$month <- as.factor(modeled_df$month)  
p <- ggplot(data= modeled_df %>% filter(month %in% plotm), 
            aes(x = year, y = modGI)) +
  geom_line(aes(group = month, color = month), size = 0.6) +
  #geom_ribbon(data = modeled_df %>% filter(month %in% plotm),
  #            aes(ymin = lower, ymax = upper, group = month, color = month, fill = month), 
  #            outline.type = "full", linetype = 2, alpha = 0.1, size = 0.6) +
  geom_line(data = modeled_df %>% filter(month %in% plotm),
              aes(y = lower, group = month, color = month), 
              linetype = 2, size = 0.9) +
  #scale_colour_viridis_d() +
  scale_colour_manual(values = c("#0836C1","#3E88EF","#d62828","#f77f00")) +
  scale_x_continuous(expand= c(0,0), limits = c(0, NA)) +
  #labs(title = loc$site[1],
  #    color ='Month')+
  theme(text = element_text(colour = 'black',size = 13),
        strip.text = element_text(colour = 'black',size = 13),
        axis.ticks.length=unit(0.3, "lines"),
        panel.border = element_rect(fill = NA, color = "black"),
      #  legend.position="none",
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank())

p
ggsave(paste0(output.p,"/LTspi_low.png"), p, width = 3, height = 2, dpi = 300)
ggsave(paste0(output.p,"/LTspi_low.svg"), p, width = 3, height = 2)
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

#ggsave(filename = paste0(output.p,"/",plotm,"20LT.png"),plot = p, width =12.5, height = 8, dpi = 300)
#ggsave(filename = paste0(output.p,"/",plotm,"20LT.svg"),plot = p, width =12.5, height = 8, dpi = 300)

#}
plotm <- plot_ms

p <- ggplot(data= modeled_df %>% filter(year > 1875 & month %in% plotm), 
            aes(x = year, y = modGI, group = model)) +
  geom_line(aes(color = model), size = 1.2) +
  geom_line(data= prcp_pos %>% filter(year > 1875 & month %in% plotm), 
            aes(y = precip, group = model, colour = model), size = 0.5, alpha = 0.7) +
  geom_ribbon(data = modeled_df %>% filter(year > 1875 & month %in% plotm), 
              aes(ymin= lowp, ymax = upp, group = model,fill = model), alpha = 0.25) +
  geom_ribbon(data = modeled_df %>% filter(year > 1875 & month %in% plotm),
              aes(ymin= lower,ymax= upper), alpha = 0.25, fill = "grey") +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(expand= c(0,0), limits = c(1875, NA)) +
  facet_wrap(~month, ncol = 2,
             labeller = labeller(month = supp.labs)) +
  labs(#title = paste(plotm,loc$site[1]),
    y="3-Month ave.prcp(mm/Month)",
    color ='datasets')+
  theme(text = element_text(colour = 'black',size = 20),
        strip.text = element_text(colour = 'black',size =15),
        axis.ticks.length=unit(0.3, "lines"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank())

p
ggsave(paste0(output.p,"/1875LT.png"), p, width =4, height = 3)
ggsave(paste0(output.p,"/1875LT.svg"), p, width = 4, height = 3)
