# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: ModelingGAM.R
# | DATE: Oct.20.2021
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
require(ncdf4)
require(fitdistrplus)
require(zoo)
library(RANN)
require(viridis)

select <- dplyr::select

data_path <- "../data"
output_path <- "../output/"

instrument_df <- readRDS(file =paste0(output_path,"instrument_df.rds"))
predicted_df <- readRDS(file = paste0(output_path,"predictedpreip.rds"))

inst_df<- instrument_df %>%
  mutate(month = month(date))

predict_df <- predicted_df %>%
  group_by(date) %>%
  summarize(precip = mean(precip, na.rm = TRUE))

predict_df <- predict_df %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  mutate(variable = "predicted") %>%
  mutate(model = "naspa")%>%
  mutate(units = "mm/month") %>%
 # mutate(site = "OKC,OK") %>%
  select(date,site, year, variable, model , precip, units, month)

prcp_df <- as.data.frame(rbind(inst_df,predict_df))
prcp_df$model <- as.factor(prcp_df$model)
prcp_df <- prcp_df %>%
  mutate(decimal_year = decimal_date(date)) %>%
  arrange(date)

p <- ggplot(data=prcp_df %>%filter(year > 1900 & month == 9), 
            aes(x=date, y = precip, group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_line(size = 1) +
  labs(title = "JAS prcp",
       x="year",y="precip(mm)", size = 20) +
  theme_classic(base_size = 20)

p
ggsave(filename = paste0(output_path,loc$site[1],"JASprcp.png"),plot = p, width =12.5, height = 8, dpi = 300)

##############################################################
####                      Modeling                       #####
##############################################################

nknot <- (max(prcp_df$year) - min(prcp_df$year) + 1)/30
nknot <- round(nknot) + 1

#
#first model : never run
#pr_mod <- gam(precip ~ model + ti(year, bs = "cr", k = 15) + ti(month, bs = "cc", k = 4) +
#                ti(year,month, bs = c("cc","cr"), k = c(15,4))),
#  ~ model + ti(year, bs = "cr", k = 15) + ti(month, bs = "cc", k = 4) +
#  ti(year,month, bs = c("cc","cr"), k = c(15,4)))
#                data=prcp_df, family = Gamma)

#second model
prcp_pos <- prcp_df %>%
  filter(precip >0)

start_time <- Sys.time()
gam_fit <- gam(list(
  precip ~ model + s(month, by = model, bs = "cc", k = 6) + 
    te(year, month, bs = c("cr","cc"), k = c(20,6)), 
  ~model + s(month, by = model, bs = "cc", k = 6) + 
    te(year, month, bs = c("cr","cc"), k = c(20,6))),
  family=gammals, link=list("identity","log"), data = prcp_pos)
end_time <- Sys.time()
end_time - start_time

plot(gam_fit,all.terms=TRUE,scheme=2, pages=1)
summary(gam_fit)

### Predict results with the new data
gam_predict <- predict(gam_fit, type="response")

gam_predict  <- gam_predict  %>%
  data.frame() %>%
  as_tibble() %>%
  rename(est_mean = 1) %>%
  rename(est_shape = 2) %>%
  bind_cols(prcp_pos) %>%
  select(date, month,  model, est_mean, est_shape) 

seasonal_df <- gam_predict %>%
  mutate(est_shape = 1/exp(est_shape)) %>%
  mutate(est_scale = est_mean/est_shape) %>%
  mutate(est_sd = (est_shape * est_scale^2)^0.5)

### Plot seasonal pattern for three months accumulated precipitation
tmonth <- 4
p <- ggplot(filter(seasonal_df, month(date) == tmonth), aes(x=date, y= est_mean * 3, colour=model)) %>%
  + geom_line(size = 1) %>%
 # + geom_point(data=seasonal_df %>% filter(model=="naspa" & month == 7)) %>%
  + scale_x_date(name = "year") %>% 
  + scale_y_continuous(name= "3-Month Precip Estimated Mean (mm)") %>%
  + scale_colour_brewer(name="", type = "qual", palette = "Set2") %>%
  + labs("JFM 3-month precip(mm)") %>%
  + theme_classic(base_size = 25) %>%
  + theme(legend.position="bottom")
p
ggsave(filename = paste0(output_path,loc$site[1],tmonth,"temodel.png"),plot = p, width =12.5, height = 8, dpi = 300)
#ggsave(filename = paste0(filename,"3.svg"), plot = p, width =12.5, height = 8, dpi = 300)

### make a new combined_df with Gridmet as model

min_gridm <- min((prcp_pos %>% filter(model == "gridmet"))$year) 
max_gridm <- max((prcp_pos %>% filter(model == "gridmet"))$year)
min_cru <- min((prcp_pos %>% filter(model == "cru"))$year) 
max_cru <- max((prcp_pos %>% filter(model == "cru"))$year)
min_naspa <- min((prcp_pos %>% filter(model == "naspa"))$year) 
max_naspa <- max((prcp_pos %>% filter(model == "naspa"))$year)

gridmet_pred_df <- data.frame(expand.grid(year = seq(min_naspa,max_naspa),month = seq(1,12)), model = "gridmet", plot_model = "naspa") %>%
  bind_rows(data.frame(expand.grid(year = seq(min_gridm,max_gridm),month = seq(1,12)),model = "gridmet", plot_model = "gridmet")) %>%
  bind_rows(data.frame(expand.grid(year = seq(min_cru,max_cru),month = seq(1,12)), model = "gridmet", plot_model = "cru"))

pred_df <- data.frame(expand.grid(year = seq(min_naspa,max_naspa),month = seq(1,12)), model = "naspa", plot_model = "naspa") %>%
  bind_rows(data.frame(expand.grid(year = seq(min_gridm,max_gridm),month = seq(1,12)),model = "gridmet", plot_model = "gridmet")) %>%
  bind_rows(data.frame(expand.grid(year = seq(min_cru,max_cru),month = seq(1,12)), model = "cru", plot_model = "cru"))

### Make predictions based on this
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
                        lwr = GAM_predict$est_mean - 2*GAM_predict$se.fit.1,
                        upr = GAM_predict$est_mean + 2*GAM_predict$se.fit.1)
                       
modeled_df <- modeled_df %>%
  mutate(date = as.Date(paste0(year,"-",month,"-01")))

#modeled_df <- transform(pred_df,
#                        modGI = GAM_predict$fit,
#                        modGI_se = GAM_predict$se.fit)

#modeled_df <- modeled_df %>%
#  mutate(date = as.Date(paste0(year,"-",month,"-01")))

#Plot annual data
tmonth <- seq(1:12)
p <- ggplot(data=prcp_df%>%filter(month==7), aes(x=year, group=model,colour = model)) +
  geom_line(aes(y=precip), alpha = 0.5) +
  geom_ribbon(data = modeled_df %>% filter(month(date)==7), aes(ymin=lwr,
                                     ymax=upr, fill = model), alpha=0.8) +
  geom_line(data = modeled_df%>%filter(month(date)==7), aes(y=modGI, colour = model)) +
  scale_colour_brewer(palette = "Set2", aesthetics = "colour", name = "model") +
 # scale_colour_OkabeIto(aesthetics = "colour", name = "model") +
  scale_fill_brewer(name = "model") +
  labs(title = paste0(tmonth," modeled"),
       x="year",y="3-Month ave.precip(mm)", size = 10)+
  theme_classic(base_size = 10)

p
filename <- paste0("lterm",loc$site[1])
ggsave(filename = paste0(output_path,filename,".png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,".svg"), plot = p, width =12.5, height = 8, dpi = 300)

#Plot monthly data
tyear <- seq(2000,2010)
p <- ggplot(data=prcp_df %>%filter(year %in%  tyear), aes(x= date, group=model,colour = model)) +
    geom_ribbon(data = modeled_df %>% filter(year %in%  tyear), aes(ymin=(modGI-2*modGI_se.1),
                                                                 ymax=(modGI+2*modGI_se.1)), alpha=0.25) +
  geom_line(data = modeled_df%>% filter(year %in%  tyear), aes(y=modGI.1)) +
  geom_line(aes(y=precip), alpha = 0.8) +
  #facet_wrap(~year) +
  labs(title = paste(tyear,"year modeled"),
       x="date",y="3-Month ave.precip(mm)", size = 14)+
  scale_x_date(name = "date") +
  theme_classic(base_size = 20)

p
filename <- paste0(tyear[1], "allseason","modeled",loc$site[1])
ggsave(filename = paste0(filename,".png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,".svg"), plot = p, width =12.5, height = 8, dpi = 300)


p <- ggplot(modeled_df, aes(x=year, y=month)) %>%
  + geom_raster(aes(fill=modGI.1)) %>%
  + scale_fill_viridis(name = "Distr\nMean" ) %>%
  + theme_bw() %>%
  + scale_x_continuous(name = "Year", expand = c(0, 0)) %>%
  + scale_y_continuous(name = "Month", expand = c(0, 0)) 
  #+ coord_cartesian(xlim=c(1,365), ylim = c(1900,1999))

### Plot
p	
filename <- paste0("raster",loc$site[1])
ggsave(filename = paste0(filename,".png"),plot = p, width =12.5, height = 8, dpi = 300)





