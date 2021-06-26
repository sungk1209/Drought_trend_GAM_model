setwd("C:/Users/kyungmin/Documents/Biascorrection/NASPA")
### Path for Data and Output	
data_path <- "../data"
output_path <- "../output"

### Set up output folders
write_output_path <- output_path
#write_output_path <- file.path(output_path, "flow_model")
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

# The data is 0.5*0.5 resolution, unit = mm
fn_wv <- "OH_Westerville.data.csv"

accum_df <- read.csv(fn_wv)
accum_df <- accum_df %>%
  mutate(variable = "PRCP.99") %>%
  mutate(value = as.numeric(as.character(PRCP.99))) %>%
  select(YEAR, MO,DAY,PRCP.99, variable, value)
  
accum_df <- accum_df %>%
  mutate(roll_mean_3 = rollmeanr(x=value, k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(roll_mean_3_notna = rollsumr(x=!is.na(value), k=n_roll, fill=NA, na.rm=TRUE)) %>%
  mutate(precip_inch = case_when(roll_mean_3_notna> n_roll_min ~ roll_mean_3,
                            TRUE ~ NA_real_)) %>%
  select(-roll_mean_3_notna) %>%
  mutate(precip = precip_inch * 25.4)%>%
  filter( MO  == 8 & DAY == 1) 

forts_df <- accum_df %>%
  mutate(units = "mm") %>%
  mutate(site = "weterville") %>%
  mutate(model = "forts") %>%
  mutate(month_day = "8-1") %>%
  rename(year = YEAR) %>%
  select(site, year, month_day, variable, model, precip,units)
  

  
  