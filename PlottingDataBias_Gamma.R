
# *--------------------------------------------------------------------------
# | PROGRAM NAME: 
# | FILE NAME: ModelingGAM.R
# | DATE:May.20.2022
# | CREATED BY:  Kay Sung     
# *--------------------------------------------------------------------------
# | PURPOSE: 
# | 
# | 
# | 
# *--------------------------------------------------------------------------


output_path <- "../output/"

loc <- "Monte Vista,CO"
dir.create(file.path(output_path,loc), recursive = FALSE)
output.p <- file.path(output_path, loc)
write_figures_path <- file.path(output.p, "figure")                  

inst_df <- readRDS(file =paste0(output.p,"/instrument_df.rds"))
naspa_df <- readRDS(file = paste0(output.p,"/pre_naspa.rds"))
gpcc_df <- readRDS(file = paste0(output.p,"/gpcc_df.rds"))
#pred_df<- readRDS(file = paste0(output.p,"/pmip_3m.rds"))

inst_df <- inst_df %>%
  mutate(month = month(date))

# >>>>>>> Stashed changes

naspa_df <- naspa_df %>%
  group_by(date) %>%
  summarize(precip = mean(precip, na.rm = TRUE))

naspa_df <- naspa_df %>%
  mutate(year = year(date)) %>%
  mutate(variable = "pr") %>%
  mutate(month = month(date)) %>%
  mutate(model = "naspa")%>%
  mutate(units = "mm/month") %>%
  mutate(site = loc) %>%
  select(date,site, year,variable, model, precip, units, month)

prcp_df <- data.frame(rbind(inst_df,naspa_df)) %>%
  arrange(date)
prcp_df$model <- as.factor(prcp_df$model)

#second model
prcp_pos <- prcp_df %>%
  filter(precip >0)

prcp_pos$model <- as.factor(prcp_pos$model)

#dist_df <- prcp_pos %>%
#         filter(month == i) %>%
#         group_by(model) %>%
#         mutate(shape = coef(fitdist(precip, "gamma"))[1], rate = coef(fitdist(precip, "gamma"))[2])

pr<- data.frame(prcp = c(0,70))
plotm <- c("NDJ","FMA","MJJ","ASO")
plot_month <- c(1,4,7,10)

  for(i in c(1:4)){
    
    
    dist_df <- prcp_pos %>%
      filter(month == plot_month[i]) %>%
      group_by(model) %>%
      summarise(shape = coef(fitdist(precip, "gamma"))[1], rate= coef(fitdist(precip, "gamma"))[2]) %>%
      mutate(month = plot_month[i])
               
    p <- ggplot(data=pr,aes(x=prcp)) +
      stat_function(fun=dgamma, args=list(shape= dist_df$shape[1], rate = dist_df$rate[1]), 
                    size = 1.2, aes(colour = "CRU"))+
      stat_function(fun=dgamma, args=list(shape= dist_df$shape[2], rate= dist_df$rate[2]), 
                    size = 1.2, aes(colour = "Gridmet"))+
      stat_function(fun=dgamma, args=list(shape= dist_df$shape[3], rate= dist_df$rate[3]), 
                    size = 1.2, aes(colour = "NASPA"))+
      labs(x = "3-month averaged prcp.(mm)", y = "f(x)", 
           color='Datasets',
           title = paste0("Gamma Distribution ",
                          "in ",plotm[i]," at Monte Vista,CO")) + 
      theme(plot.title = element_text(hjust = 0.5), 
            axis.title.x = element_text(face="bold", colour="blue", size = 12),
            axis.title.y = element_text(face="bold", colour="blue", size = 12),
            legend.title = element_text(face="bold", size = 10),
            legend.position = "right") +
      theme_classic(base_size = 20)
    p
    
    ggsave(filename = paste0(write_figures_path,"/Gamma","_",plot_month[i],".png"),plot = p, width =12.5, height = 8, dpi = 300)
    ggsave(filename = paste0(write_figures_path,"/Gamma","_",plot_month[i],".svg"),plot = p, width =12.5, height = 8, dpi = 300)
    
      
  }












