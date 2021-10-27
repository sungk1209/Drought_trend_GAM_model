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
select <- dplyr::select

data_path <- "../data"
output_path <- "../output/"

inst_df <- readRDS(file =paste0(output_path,"instrument_df.rds"))
predict_df <- readRDS(file = paste0(output_path,"predictedpreip.rds"))

inst_df<- instrument_df %>%
  mutate(month = month(date)) %>%
   select(-month_day)
predict_df <- predict_df %>%
  mutate(year = year(date)) %>%
  mutate(variable = "predicted") %>%
  mutate(model = "naspa")%>%
  mutate(units = "mm/month") %>%
  mutate(site = "OKC,OK") %>%
  select(date,site, year, variable, model , precip, units, month)

prcp_df <- as.data.frame(rbind(inst_df,predict_df))
prcp_df$model <- as.factor(prcp_df$model)
prcp_df <- prcp_df %>%
  mutate(decimal_year = decimal_date(date)) %>%
  arrange(date)

p <- ggplot(data=prcp_df %>%filter(month == 7, year > 1850), 
            aes(x=date, y = precip, group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_line() +
  labs(title = "prcp",
       x="year",y="precip(mm)", size = 14) +
  theme_classic(base_size = 10)

p

##############################################################
####                      Modeling                       #####
##############################################################

nknot <- (max(prcp_df$year) - min(prcp_df$year) + 1)/30
nknot <- round(nknot) + 1

#first model
#pr_mod <- gam(precip ~ model + ti(year,month, bs = c("cc","cr"), k = c(80,4)) +
#                  ti(year,month, bs = c("cc","cr"), k = c(80,4)),
#                data=prcp_df, family = Gamma)

#second model

prcp_nozero <- prcp_df %>%
  filter(precip > 0)
start_time <- Sys.time()
gam_fit <- gam(list(
  precip ~ model + ti(year, month, bs = c("cc","cr"), k = c(10,4)), 
                     ~ model + te(year, month, bs = c("cc","cr"), k = c(10,4))),
                family=gammals, link=list("identity","log"), data = prcp_nozero)
end_time <- Sys.time()
end_time - start_time

#third model
comb_gam <- gam(list(
  precip~s(year,bs="cc", k=15) + s(month, bs="cr", k=4) + model, 
  ~s(year,bs="cc", k=15) + s(month, bs="cr", k=4) + model),
  family=gammals, 
  link=list("identity","log"), 
  data = prcp_df, 
  method = "REML")
end_time <- Sys.time()
end_time - start_time

plot(comb_gam,all.terms=TRUE,scheme=2, pages=1)
summary(comb_gam)

seasonal_data <- data.frame(date=seq(as.Date("1700-01-01"), as.Date("1700-12-31"), by= "1 month")) %>%
  mutate(month = month(date)) %>%
  mutate(year = year(date)) %>%
  mutate(model = "grimet")


### Predict results with the new data
gam_prdict <- predict(gam_fit, type="response")

gam_prdict  <- gam_prdict  %>%
  data.frame() %>%
  as_tibble() %>%
  rename(est_mean = 1) %>%
  rename(est_shape = 2) %>%
  bind_cols(prcp_df) %>%
  select(date, month, decimal_year, model, est_mean, est_shape) 


seasonal_df <- gam_prdict %>%
  mutate(est_shape = 1/exp(est_shape)) %>%
  mutate(est_scale = est_mean/est_shape) %>%
  mutate(est_sd = (est_shape * est_scale^2) ^ 0.5)


### Plot seasonal pattern for three months accumulated precipitation
p <- ggplot(filter(seasonal_df, month(date) == 3), aes(x=date, y= est_mean * 3, colour=model)) %>%
  + geom_line() %>%
 # + geom_point(data=seasonal_df %>% filter(model=="naspa" & month == 7)) %>%
  + scale_x_date(name = "year") %>% 
  + scale_y_continuous(name= "3-Month Precip Estimated Mean (mm)") %>%
  + scale_colour_brewer(name="", type = "qual", palette = "Set2") %>%
  + labs("JFM 3-month precip(mm)") %>%
  + theme_classic(base_size = 15) %>%
  + theme(legend.position="bottom")
p
ggsave(filename = paste0(output_path,"temodel.png"),plot = p, width =12.5, height = 8, dpi = 300)
#ggsave(filename = paste0(filename,"3.svg"), plot = p, width =12.5, height = 8, dpi = 300)

precip <- transform(combined_df, modGI = predict(pr_modGI, type="response"),
                    na.rm = TRUE)
GI_predict <- predict(pr_modGI, se.fit = TRUE, type = "response")

pr_modGI_ins <- gam(precip ~ model + ti(year, bs = "cr", k = 7),
                    data=instrument_df, family = Gamma)

precip <- transform(instrument_df, modGI = predict(pr_modGI_ins, type="response"),
                    na.rm = TRUE)
GI_predict_ins <- predict(pr_modGI_ins, se.fit = TRUE, type = "response")

### make a new combined_df with Gridmet as model
min_noaa <- min((prcp_df %>% filter(model == "Noaa_v3"))$year) 
max_noaa <- max((prcp_df %>% filter(model == "Noaa_v3"))$year)
min_gridm <- min((prcp_df %>% filter(model == "Gridmet"))$year) 
max_gridm <- max((prcp_df %>% filter(model == "Gridmet"))$year)
min_cru <- min((prcp_df %>% filter(model == "cru"))$year) 
max_cru <- max((prcp_df %>% filter(model == "cru"))$year)
min_naspa <- min((prcp_df %>% filter(model == "naspa"))$year) 
max_naspa <- max((prcp_df %>% filter(model == "naspa"))$year)

gridmet_pred_df <- data.frame(expand.grid(year = seq(min_noaa,max_noaa), month = seq(1,12)), model = "Gridmet", plot_model = "NOAA") %>%
  bind_rows(data.frame(expand.grid(year = seq(min_gridm,max_gridm),month = seq(1,12)),model = "Gridmet", plot_model = "Gridmet")) %>%
  bind_rows(data.frame(expand.grid(year = seq(min_cru,max_cru),month = seq(1,12)), model = "Gridmet", plot_model = "cru")) %>%
  bind_rows(data.frame(expand.grid(year = seq(min_naspa,max_naspa),month = seq(1,12)), model = "Gridmet", plot_model = "naspa"))


### Make predictions based on this
GAM_predict <- predict(gam_fit, newdata = gridmet_pred_df, se.fit = TRUE, type = "response")

modeled_df <- transform(gridmet_pred_df,
                        modGI = GAM_predict$fit,
                        modGI_se = GAM_predict$se.fit)
modeled_df <- modeled_df %>%
  mutate(date = as.Date(paste0(year,"-",month,"-01")))
#Plot annual data
tmonth <- 12
p <- ggplot(data=prcp_df %>%filter(month == tmonth), aes(x=year, group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_ribbon(data = modeled_df %>% filter(month == tmonth), aes(ymin=(modGI.1-2*modGI_se.1),
                                                ymax=(modGI.1+2*modGI_se.1)), alpha=0.25) +
  geom_line(data = modeled_df%>% filter(month == tmonth), aes(y=modGI.1)) +
  geom_line(aes(y=precip), alpha = 0.8) +
 labs(title = "DJF modeled",
       x="year",y="3-Month ave.precip(mm)", size = 14)+
  theme_classic(base_size = 30)

p
filename <- paste0(tmonth, "te_","modeled",loc$site[1])
ggsave(filename = paste0(filename,".png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,".svg"), plot = p, width =12.5, height = 8, dpi = 300)

#Plot monthly data
tyear <- seq(1970,1970)
p <- ggplot(data=prcp_df %>%filter(year %in%  tyear), aes(x= date, group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_ribbon(data = modeled_df %>% filter(year %in%  tyear), aes(ymin=(modGI.1-2*modGI_se.1),
                                                                 ymax=(modGI.1+2*modGI_se.1)), alpha=0.25) +
  geom_line(data = modeled_df%>% filter(year %in%  tyear), aes(y=modGI.1)) +
  geom_line(aes(y=precip), alpha = 0.8) +
  labs(title = paste(tyear,"year modeled"),
       x="date",y="3-Month ave.precip(mm)", size = 14)+
  scale_x_date(name = "date") +
  theme_classic(base_size = 20)

p
filename <- paste0(tyear, "GI_","modeled",loc$site[1])
ggsave(filename = paste0(filename,".png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,".svg"), plot = p, width =12.5, height = 8, dpi = 300)













########################################################
######### Misc.Old code ModGI                ###########
########################################################
pr_modGI <- gam(precip ~ model + ti(year, bs = "cr", k = 80),
                data=combined_df, family = Gamma)

precip <- transform(combined_df, modGI = predict(pr_modGI, type="response"),
                    na.rm = TRUE)
GI_predict <- predict(pr_modGI, se.fit = TRUE, type = "response")

######make a new model based on each model
min_naspa <- min((combined_df %>% filter(model == "NASPA"))$year) 
max_naspa <- max((combined_df %>% filter(model == "NASPA"))$year)

######make a model with seperate interception using GI##########
################################################################

combined_df_ <- data.frame(year = seq(min_naspa,max_naspa), model = "NASPA", plot_model = "NASPA") %>%
  bind_rows(data.frame(year = seq(min_noaa,max_noaa), model = "NOAA", plot_model = "NOAA")) %>%
  bind_rows(data.frame(year = seq(min_gridm,max_gridm), model = "Gridmet", plot_model = "Gridmet"))
### Make predictions based on this
GI_predict <- predict(pr_modGI, newdata = combined_df, se.fit = TRUE, type = "response")

combined_df <- transform(combined_df,
                         modGI = GI_predict$fit, 
                         modGI_se = GI_predict$se.fit)

p <- ggplot(data=precip, aes(x=year,  group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_ribbon(data = combined_df,aes( ymin=(modGI-2*modGI_se),
                                      ymax=(modGI+2*modGI_se)), alpha=0.25) +
  geom_line(data = combined_df, aes(y=modGI), size = 1) +
  geom_line(aes(y=precip),alpha = 0.6) +
  labs(title = "August1_GI model: knots = 80",
       x="year",y="3 months averaged rainfall(mm)") +
  theme_classic(base_size = 30)

p
filename <- paste0(loc$site[1],"GI_k80")
ggsave(filename = paste0(filename,"3.png"), plot = p,  width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,"3.svg"), plot = p,  width =12.5, height = 8, dpi = 300)

#######################################################################

### make a new combined_df with Gridmet as model
calibrated_df <- data.frame(year = seq(min_naspa,max_naspa), model = "Gridmet", plot_model = "NASPA") %>%
  bind_rows(data.frame(year = seq(min_cru,max_cru), model = "Gridmet", plot_model = "cru")) %>% 
  bind_rows(data.frame(year = seq(min_noaa,max_noaa), model = "Gridmet", plot_model = "NOAA")) %>%
  bind_rows(data.frame(year = seq(min_gridm,max_gridm), model = "Gridmet", plot_model = "Gridmet"))

### Make predictions based on this
GI_predict_modified <- predict(pr_modGI, newdata = calibrated_df, se.fit = TRUE, type = "response")

calibrated_df <- transform(calibrated_df,
                           modGI = GI_predict_modified$fit, 
                           modGI_se = GI_predict_modified$se.fit)

p <- ggplot(data=precip, aes(x=year,  group=model,colour = model)) +
  #facet_wrap(~model) +
  geom_ribbon(data = calibrated_df,aes(ymin=(modGI-2*modGI_se),
                                       ymax=(modGI+2*modGI_se)), alpha=0.25) +
  geom_line(aes(y=precip),alpha = 0.5) +
  geom_line(data = calibrated_df, aes(y=modGI), size = 1) +
  labs(title = "August1_GI model: knots = 80",
       x="year",y="precip(mm)", size = 14) +
  theme_classic(base_size = 30)

p

filename <- paste0("calibrated",loc$site[1])

ggsave(filename = paste0(filename,"3.svg"), plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0(filename,"3.png"), plot = p, width =12.5, height = 8, dpi = 300)

save(calibrated_df, file = paste0(data_path,filename,".RData"))
#######################################################################


