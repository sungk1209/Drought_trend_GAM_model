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
library(ggthemes)
library(cowplot)
require(colorblindr)

select <- dplyr::select

data_path <- "../data"
output_path <- "../output/"

inst_df <- readRDS(file =paste0(output_path,loc$site[1],"instrument_df.rds"))
naspa_df <- readRDS(file = paste0(output_path,loc$site[1],"predictedpreip.rds"))
pred_df<- readRDS(file = paste0(output_path,loc$site[1],"pred_df.rds"))

naspa_df <- naspa_df %>%
  group_by(date) %>%
  summarize(precip = mean(precip, na.rm = TRUE))

inst_df <- inst_df %>%
  mutate(month = month(date))

naspa_df <- naspa_df %>%
  mutate(year = year(date)) %>%
  mutate(month = month(date)) %>%
  mutate(variable = "predicted") %>%
  mutate(model = "naspa")%>%
  mutate(units = "mm/month") %>%
  mutate(site = "COL-OH") %>%
  select(date,site, year, variable, model , precip, units, month)

prcp_df <- as.data.frame(rbind(inst_df,naspa_df,pred_df)) %>%
  arrange(date)
prcp_df$model <- as.factor(prcp_df$model)

#second model
prcp_pos <- prcp_df %>%
  filter(precip >0)

p <- ggplot(data=prcp_df %>% filter(year >1900), 
            aes(x=date, y = precip, group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_line() +
  labs(title = "prcp",
       x="year",y="precip(mm)", size = 20) +
  theme_classic(base_size = 20)

p
ggsave(filename = paste0(output_path,loc$site[1],"prcp.png"),plot = p, width =12.5, height = 8, dpi = 300)

##############################################################
####                      Modeling                       #####
##############################################################
nknot <- (max(prcp_df$year) - min(prcp_df$year) + 1)/30
nknot <- round(nknot) + 1

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

#gam_fit <- gam(precip ~ model + s(month, by = model, bs = "cc", k = 6) + 
#       family=gamma, data = prcp_pos)


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
p <- ggplot(filter(seasonal_df, month(date) == tmonth), aes(x=date, y= est_mean, colour=model)) %>%
  + geom_line(size = 1) %>%
  + geom_point(data=seasonal_df %>% filter(model=="naspa" & month == tmonth)) %>%
  + scale_x_date(name = "year") %>% 
  + scale_y_continuous(name= "3-Month Precip Estimated Mean (mm)") %>%
  #+ scale_colour_brewer(name="", type = "qual", palette = "Set2") %>%
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
min_gcm <- min((prcp_pos %>% filter(model == "GFDL-CM4"))$year) 
max_gcm <- max((prcp_pos %>% filter(model == "GFDL-CM4"))$year)

gridmet_pred_df <- data.frame(expand.grid(year = seq(min_naspa,max_naspa),month = seq(1,12)), model = "gridmet", plot_model = "naspa") %>%
  bind_rows(data.frame(expand.grid(year = seq(min_gridm,max_gridm),month = seq(1,12)),model = "gridmet", plot_model = "gridmet")) %>%
  bind_rows(data.frame(expand.grid(year = seq(min_cru,max_cru),month = seq(1,12)), model = "gridmet", plot_model = "cru")) %>%
  bind_rows(data.frame(expand.grid(year = seq(min_gcm,max_gcm),month = seq(1,12)), model = "gridmet", plot_model = "GFDL-CM4"))

pred_df <- data.frame(expand.grid(year = seq(min_naspa,max_naspa),month = seq(1,12)), model = "naspa", plot_model = "naspa") %>%
  bind_rows(data.frame(expand.grid(year = seq(min_gridm,max_gridm),month = seq(1,12)),model = "gridmet", plot_model = "gridmet")) %>%
  bind_rows(data.frame(expand.grid(year = seq(min_cru,max_cru),month = seq(1,12)), model = "cru", plot_model = "cru")) %>%
  bind_rows(data.frame(expand.grid(year = seq(min_gcm,max_gcm),month = seq(1,12)), model = "GFDL-CM4", plot_model = "GFDL-CM4"))

### Make predictions based on this
GAM_predict <- predict(gam_fit, newdata = Gridmet_pred_df, se.fit = TRUE, type = "response")

GAM_predict  <- GAM_predict %>%
  data.frame() %>%
  as_tibble() %>%
  rename(est_mean = 1) %>%
  rename(est_shape = 2) 

GAM_predict <- GAM_predict %>%
  mutate(est_shape = 1/exp(est_shape)) %>%
  mutate(est_scale = est_mean/est_shape) %>%
  mutate(est_sd = (est_shape * est_scale^2)^0.5)

###############################################################
########calculate simultaneous 95% confidence interval########
##############################################################
#https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

set.seed(42)
N <- 10000

rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}
Vb <- vcov(gam_fit) #Returns the variance-covariance matrix of the main parameters of a fitted model object
BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
Cg <- predict(gam_fit, gridmet_pred_df, type = "lpmatrix")
simDev <- Cg %*% t(BUdiff)

absDev <- abs(sweep(simDev, 1, GAM_predict$se.fit.1, FUN = "/"))
masd <- apply(absDev, 2L, max)
crit <- quantile(masd, prob = 0.95, type = 8)

modeled_df <- transform(gridmet_pred_df,
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

#Plot annual data
basis_func_palette = colorblind_pal()(8)
basis_func_palette = basis_func_palette[-c(1,5)]

tmonth <- 7
p <- ggplot(data=prcp_pos %>% filter(year > 200 & month == tmonth), aes(x=year, group=model,colour = model)) +
  geom_line(aes(y=precip)) +
  geom_ribbon(data = modeled_df%>%filter(year > 200 & month == tmonth), aes(ymin= lowp,
                                     ymax= upp), alpha = 0.25, fill = "black") +
 # geom_ribbon(data = modeled_df%>%filter(year > 1980 & month == 4), aes(ymin= lowcrit, ymax= upcrit), fill = "black", 
#              alpha = 0.5) +
  geom_ribbon(data = modeled_df%>%filter(year > 200 & month == tmonth), aes(ymin= lower,
                                                          ymax= upper), alpha = 0.25, fill = "black") +
   geom_line(data = modeled_df%>%filter(year > 200 & month == tmonth), aes(y=modGI)) +
   # facet_wrap(~month) +
  labs(title = paste0("MJJ"," modeled"),
       subtitle = "",
       x="year",y="3-Month ave.precip(mm)", size = 10)+
    theme_classic(base_size = 20)

p
filename <- paste0(tmonth,"lterm",loc$site[1])
ggsave(filename = paste0(output_path,filename,"1.png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,".svg"), plot = p, width =12.5, height = 8, dpi = 300)

#Plot monthly data
tyear <- seq(2000,2020)
p <- ggplot(data=prcp_pos %>%filter(year %in%  tyear), aes(x= date, group=model,colour = model)) +
    #geom_ribbon(data = modeled_df %>% filter(year %in%  tyear), 
    #            aes(ymin=lower,ymax=upper), alpha=0.25) +
  #geom_line(data = modeled_df%>% filter(year %in%  tyear), aes(y=modGI, group = model, colour = model)) +
  geom_line(aes(y=precip), alpha = 0.8) +
  #facet_wrap(~year) +
 # scale_colour_brewer(palette = "Set1", aesthetics = "colour", name = "Data set") +
  labs(title = paste(tyear,"year modeled"),
       x="date",y="3-Month ave.precip(mm)", size = 14)+
  scale_x_date(name = "date") +
  theme_classic(base_size = 20)

p
filename <- paste0(tyear[1], "allseason","modeled",loc$site[1])
ggsave(filename = paste0(filename,".png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,".svg"), plot = p, width =12.5, height = 8, dpi = 300)


basis_func_palette = colorblind_pal()(8)
basis_func_palette = basis_func_palette[-c(1,5)]

p <- ggplot(modeled_df, aes(x= month, y = year)) %>%
  + geom_raster(aes(fill=modGI)) %>%
  + scale_fill_viridis(name = "Mean" ) %>%
  + theme_bw() %>%
  + scale_x_continuous(name = "Year", expand = c(0, 0)) %>%
  + scale_y_continuous(name = "Month", expand = c(0, 0)) 
  #+ coord_cartesian(xlim=c(1,365), ylim = c(1900,1999))

### Plot
p	
filename <- paste0("raster",loc$site[1])
ggsave(filename = paste0(filename,".png"),plot = p, width =12.5, height = 8, dpi = 300)
####

tyear <- seq(500,2020,100)

p <- ggplot(data= modeled_df %>% filter ( year %in% tyear), 
            aes(x = month, y = modGI, group = interaction(year,plot_model), color = year)) +
  geom_line(se = FALSE) +
  #geom_ribbon(data = modeled_df %>% filter(year %in%  tyear), 
  #            aes(ymin=lower,ymax=upper), alpha=0.25) +
  #geom_line(data = modeled_df%>% filter(year %in%  tyear), aes(y=modGI, group = model, colour = model)) +
  #facet_wrap(~year) +
  scale_colour_viridis() +
  labs(title = paste("Average precipitation"),
      y="3-Month ave.precip(mm)", size = 14)+
    theme_classic(base_size = 20)

p


filename <- paste0(tyear[1], "allseason","modeled",loc$site[1])
ggsave(filename = paste0(filename,".png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,".svg"), plot = p, width =12.5, height = 8, dpi = 300)





