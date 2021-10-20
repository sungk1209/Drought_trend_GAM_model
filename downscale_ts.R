# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: get_ncdf.R
# | DATE: Sep.26.2021
# | CREATED BY:  Kay Sung     
# *--------------------------------------------------------------------------
# | PURPOSE: 
# | Download gridmet, Noaa reanlaysis and NASPA data
# | 
# | 
# *--------------------------------------------------------------------------
require(stars)
require(ncdf4)
require(lubridate)
require(fields)
require(purrr)
require(tidyverse)
require(zoo)
require(mgcv)
require(gratia)
require(dplyr)
require(leaps)
require(MASS)

### Path for Data and Output	
data_path <- "./data"
output_path <- "./output/"

########work with one location: 

#loc <- data.frame(site="lakeOroville,CA",lon=c(-121.5,-121.00),lat= c(39.5,40.0))
loc <- data.frame(site="OKC,OK",
                  lon=c(-98,-97.5),
                  lat= c(35.0,35.5))

filename <- paste0("calibrated",loc$site[1])
load(paste0(data_path,filename,".RData"))
ts_df <- calibrated_df %>%
  mutate(month_day = "7-31") %>%
  select(-model,-modGI_se)

#Think about AR with trend / real precipitation
#function: prcip_yi,mk = precip_y(i-1),mk + precip_yi,m(k-1) + precip_yi,m(k+1)
acf(cru_df$precip[3:length(cru_df$precip)])
pacf(cru_df$precip[3:length(cru_df$precip)])

acf(gridmet_df$precip[5:length(gridmet_df$precip)])
pacf(gridmet_df$precip[100:length(gridmet_df$precip)])

acf(reanal_v3$precip[3:length(reanal_v3$precip)])
pacf(reanal_v3$precip[3:length(reanal_v3$precip)])

monthly_gm <- gridmet_df %>%
  group_by(month_day)

data_temp <- monthly_gm$data[[8]]$precip[2:43]

fit <- lm(data[[8]]$precip ~ 
            data[[7]]$precip + 
            data[[6]]$precip + 
            data[[5]]$precip +
            c(data_temp,NA), na.action=na.exclude,
          data= monthly_gm)

step <- stepAIC(fit, direction="both")
step$anova # display results

leaps <- regsubsets(data[[8]]$precip ~ 
            data[[7]]$precip + 
            data[[6]]$precip + 
            data[[5]]$precip +
            c(data_temp,3),
          data= monthly_gm)

pr_naspa <- combined_df %>% filter(model == "NASPA") %>% select(precip)

fit_gam <- gam(list(precip
                  ~ te(year,month, bs = c("cc","cr")) + model,
                  ~ te(year,month, bs = c("cc","cr")) + model), 
             data= instrument_df, 
             family=gammals,
             select=TRUE,
             method="REML")

comb_gam <- gam(list(roll_mean_3~te(decimal_year, jdate, bs = c("cr","cc"), k = c(0,0) ) + model, ~te(decimal_year, jdate, bs = c("cr","cc"), k = c(0,0)) + model),family=gammals, link=list("identity","log"), data = yup_data, method = "REML")


fit_gam <- gam(list(precip
                    ~ s(year, bs = c("cc"),k = 100) ,
                    ~ s(year, bs = c("cc"), k = 100)), 
               data= naspa_df%>%filter(year> 1500), 
               family=gammals,
               select=TRUE,
               method="REML")

summary(fit_gam)
plot(fit_gam)

month_test <- 7
combined_df <- expand.grid(year = seq(1836,2000),month = (1:12), model = c("cru","Gridmet","Noaa_v3"))

  ### Make predictions based on this
predict <- predict(fit_gam, newdata = combined_df, se.fit = TRUE)

combined_df <- transform(combined_df,
                         modGI = predict$fit[,1]+predict$fit[,2], 
                         modGI_se = predict$se.fit[,1]+predict$se.fit[,2])

p <- ggplot(data=combined_df, aes(x=year)) +
  #facet_wrap(~model) +
  geom_ribbon(data = combined_df,aes(ymin=(modGI-2*modGI_se),
                                    ymax=(modGI+2*modGI_se)), alpha=0.25) +
  geom_line(aes(y=modGI), size = 1) +
  geom_line(data = instrument_df, aes(y=precip, group = model, color = model),alpha = 0.6) +
 # geom_line(data = naspa_df, aes(y=precip),alpha = 0.6) +
  labs(title = "Apr-May-Jun sum precip prediction", y = "precip(mm)") +
  theme_classic(base_size = 30)
p

filename <- "OKC,OK,apr-may-jun"
ggsave(filename = paste0("../output/",filename,".png"),plot = p, width =12.5, height = 8, dpi = 300)
ggsave(filename = paste0("../output/",filename,".svg"), plot = p, width =12.5, height = 8, dpi = 300)


