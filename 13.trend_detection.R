
# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: gpcc_downscaling.R
# | DATE: Jul.06.2023
# | CREATED BY: Kay Sung     
# *--------------------------------------------------------------------------
# | PURPOSE: 
# | 
# | 
# | 
# *--------------------------------------------------------------------------

require(lubridate)
require(tidyverse)
require(dplyr)
require(gratia)
require(mgcv)
require(MASS) 
require(gridExtra)

select <- dplyr::select

data_path <- "../data/"
output_path <- "../output/"

#call list file in 5.GI_model_loop.R 
j <- 1
loc <- site_list[[j]]

output.p <- file.path(output_path, loc$site[1])
gam_fit <- readRDS(file =paste0(output.p,"/model.rds"))
prcp_df <- readRDS(file =paste0(output.p,"/prcp_df.rds"))

prcp_pos <- prcp_df %>%
  filter(precip >0)


min_gridm <- min((prcp_pos %>% filter(model == "gridmet"))$year) 
max_gridm <- max((prcp_pos %>% filter(model == "gridmet"))$year)
min_cru <- min((prcp_pos %>% filter(model == "cru"))$year) 
max_cru <- max((prcp_pos %>% filter(model == "cru"))$year)
min_naspa <- min((prcp_pos %>% filter(model == "naspa"))$year) 
max_naspa <- max((prcp_pos %>% filter(model == "naspa"))$year)

tmonth <- 7

gridmet_pred_df <- data.frame(expand.grid(year = seq(min_naspa,max_naspa),month = tmonth), model = "gridmet", plot_model = "naspa") %>%
  bind_rows(data.frame(expand.grid(year = seq(min_gridm,max_gridm),month = tmonth),model = "gridmet", plot_model = "gridmet")) %>%
  bind_rows(data.frame(expand.grid(year = seq(min_cru,max_cru),month = tmonth), model = "gridmet", plot_model = "cru"))

### Make predictions based on this
GAM_predict <- predict(gam_fit, newdata = gridmet_pred_df, se.fit = TRUE, type = "response")

GAM_predict  <- GAM_predict %>%
  data.frame() %>%
  as_tibble() %>%
  rename(est_mean = 1) %>%
  rename(est_shape = 2) 

GAM_predict <- GAM_predict %>%
  mutate(est_shape = 1/exp(est_shape)) %>%
  mutate(est_scale = est_mean/est_shape) %>%
  mutate(est_sd = (est_shape * est_scale^2)^0.5)

modeled_df <- transform(gridmet_pred_df,
                        meanGI = GAM_predict$est_mean,
						shapeGI = GAM_predict$est_shape,
						scaleGI = GAM_predict$est_scale)
                        
						
modeled_df <- modeled_df %>%
  mutate(date = as.Date(paste0(year,"-",month,"-01"))) %>%
  arrange(date) %>%
  group_by(date) %>%
  summarize(meanGI = mean(meanGI), shapeGI = mean(shapeGI), scaleGI = mean(scaleGI))
  

derive_df <- modeled_df %>%
	mutate(year = year(date)) %>%
	mutate(slope = meanGI - lag(meanGI)) %>%
	mutate(draw = 10000)
	
	draw_total <- data.frame()	
	
  for( i in c(1:100)) {
  
  
  draw_df <- draw_fromgammals(model = gam_fit, n = 1, new_data = gridmet_pred_df)
  draw_df <- draw_df %>%
	mutate(year = year(date)) %>%
	mutate(slope = meanGI - lag(meanGI)) %>%
	mutate(draw = i)
	
   draw_total <- bind_rows(draw_total, draw_df)
  
  }
  
draw_total <- bind_rows(derive_df, draw_total)  

p1 <- ggplot(draw_total, aes(x = date, y = slope, group = draw, color = "draws")) + 
	geom_line() +
	#geom_line(aes(color = ifelse(slope < 0, "Below Zero", "Above Zero")), size = 1) +
	
	geom_line(data = draw_total %>% filter(draw == 10000), aes(x = date, y= slope, color = "Modeled")) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_bw()
	
p2 <- ggplot(draw_total, aes(x = date, y = meanGI, group = draw, color = "draws")) + 
	geom_line() +
	geom_line(data = draw_total %>% filter(draw == 10000), aes(x = date, y= meanGI, color = "Modeled")) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_bw()
	
combined_plot <- grid.arrange(p1, p2, nrow = 2, top = "Oklahoma city 3 month prcp. in MJJ (top: slope, bottom: prcp)")

# Display the combined plot
print(combined_plot)
	
	

disp_to_shape <- function(disp){
  1/ exp(-7 + log(1 + exp(disp)))
}

### Function to draw from gammals


rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n), m, n))
}

draw_fromgammals <- function(model, n, new_data){
  
  a <- nrow(new_data)
  
  Cg <- predict(model, new_data, type = "lpmatrix")
  Vb <- vcov(model)
  sims <- rmvn(n, mu = coef(model), sig = Vb)
  
  coef_dims <- dim(sims)
  mean_dims <- 1:(coef_dims[2]/2)
  shape_dims <- (max(mean_dims)+1):coef_dims[2]
  fits_mean <- Cg[,mean_dims] %*% t(sims)[mean_dims,]
  fits_shape <- Cg[,shape_dims] %*% t(sims)[shape_dims,]
  
  pred_values <- data.frame(draw = seq(1,coef_dims[1]), mean_link = c(fits_mean), shape_link = c(fits_shape))
  pred_values <- pred_values %>%
    mutate(meanGI = exp(mean_link)) %>%
    mutate(shapeGI = disp_to_shape(shape_link)) %>%
    mutate(scaleGI = meanGI/ shapeGI)
  
  pred_values <- bind_cols(new_data, pred_values) %>%
	mutate(date = as.Date(paste0(year,"-",month,"-01")))  %>%
	group_by(date) %>%
    summarize(meanGI = mean(meanGI), shapeGI = mean(shapeGI), scaleGI = mean(scaleGI)) %>%
  	select(date,meanGI,shapeGI,scaleGI)
}
