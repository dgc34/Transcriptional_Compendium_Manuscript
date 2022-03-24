library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(patchwork)

data_list <- readRDS('~/Documents/darpa-manuscript-data/data_list_final/data_list_v1.RDS')
source('~/Dropbox/plot_palette.R')

# goal of this script: 
# panel C: create hist of sample sizes for bacterial and viral datasets
# panel D: plot number of time series datasets that are viral and bacterial
# panel E: plot number of each platform in bact & vir

# there's a chance we might need this for other categories, so maybe let's try to write functions?

### HISTOGRAMS OF SAMPLE SIZE
findStudyWithClass <- function(MIobj, pos_class, neg_class, N_threshold = 4){
  pos_ix <- sum(MIobj$pheno$Class %in% pos_class)
  neg_ix <- sum(MIobj$pheno$Class %in% neg_class)
  
  if(sum(pos_ix) >= N_threshold & sum(neg_ix) >= N_threshold){
    return(T)
  } else {
    return(F)
  }
}

bac_ix <- sapply(data_list, findStudyWithClass, pos_class = 'Bacteria', neg_class = c('Healthy', 'Other.Infectious', 'Other.NonInfectious', 'Convalescent'))
sapply(data_list[bac_ix], function(x){table(x$pheno$Class)})
vir_ix <- sapply(data_list, findStudyWithClass, pos_class = 'Virus', neg_class = c('Healthy', 'Other.Infectious', 'Other.NonInfectious', 'Convalescent'))
sapply(data_list[vir_ix], function(x){table(x$pheno$Class)})

getNClass <- function(MIobj, input_classes){
  return(sum(MIobj$pheno$Class %in% input_classes))
}

bac_N <- sapply(data_list[bac_ix], getNClass, input_classes = c('Bacteria', 'Healthy', 'Other.Infectious', 'Other.NonInfectious'))
vir_N <- sapply(data_list[vir_ix], getNClass, input_classes = c('Virus', 'Healthy', 'Other.Infectious', 'Other.NonInfectious'))

hist_plot_dat <- data.frame('Type' = c(rep('Virus', length(vir_N)),
                      rep('Bacteria', length(bac_N))),
           'N' = c(vir_N, bac_N))
ggplot(hist_plot_dat, aes(x = N, fill = Type)) + geom_histogram()
p2c <- ggplot(hist_plot_dat, aes(x = Type, y = N, fill = Type)) + 
  geom_boxplot() + 
  stat_compare_means(method = 'wilcox.test') + 
  theme_clean

### TIME SERIES DATSET COMPARISON
findTS <- function(MIobj){
  if(sum(names(MIobj$pheno) == 'Time.Point') > 0){
    return(T)
  } else{
    return(F)
  }
}
time_series_ix <- sapply(data_list, findTS)

ts_list <- data_list[time_series_ix]
ts_bac_ix <- sapply(ts_list, findStudyWithClass, pos_class = 'Bacteria', neg_class = c('Healthy', 'Other.NonInfectious', 'Convalescent'))
ts_vir_ix <- sapply(ts_list, findStudyWithClass, pos_class = 'Virus', neg_class = c('Healthy', 'Other.NonInfectious', 'Convalescent'))
ts_non_ix <- sapply(ts_list, findStudyWithClass, pos_class = 'Other.Infectious', neg_class = c('Healthy', 'Other.NonInfectious', 'Convalescent'))
sum(ts_bac_ix)
sum(ts_vir_ix)
sum(ts_non_ix)


sapply(ts_list[ts_non_ix], function(x){table(x$pheno$Pathogen)})
# GSE28405_GPL2700 is viral
# GSE68310_GPL10558 is viral
# GSE7000_GPL4857 is bacterial and other.infectious
# GSE35859_GPL15241 is other.infectious
sapply(ts_list[ts_vir_ix], function(x){table(x$pheno$Pathogen)})
# these are all correct
sapply(ts_list[ts_bac_ix], function(x){table(x$pheno$Pathogen)})
N_other <- 2 # 7000 is common with bacteria (manually checked)
N_vir <- sum(ts_vir_ix) # all look correct
N_bac <- sum(ts_bac_ix) # 2 are in common with the virus list (40012 and 20346)

ts_N_df <- data.frame('Type' = c('Virus', 'Bacteria', 'Other Infectious'), 
                      'Time Series Studies' = c(N_vir, N_bac, N_other))
p2d <- ts_N_df %>%
  filter(Type %in% c('Bacteria', 'Virus'))  %>% 
  ggplot(aes(x = Type, y = Time.Series.Studies)) + 
  geom_bar(stat = 'identity') + 
  theme_clean + 
  labs(x = 'Infection Type', y = 'Number of Time Series')

unused_ix <- !(ts_bac_ix|ts_vir_ix|ts_non_ix)
names(which(unused_ix))
#data_list$GSE35859_GPL15241$pheno$Class %>% table()
#data_list$GSE35859_GPL15241$pheno$Pathogen%>% table()

# snapshot studies are the unused from each category
vir_accs <- names(data_list[vir_ix])
bac_accs <- names(data_list[bac_ix])
if(mean(names(ts_list[ts_vir_ix]) %in% vir_accs) != 1){
  stop('not all viral ts in full ts list')
}
N_vir_snapshot <- sum(!vir_accs %in% names(ts_list[ts_vir_ix]))
N_bac_snapshot <- sum(!bac_accs %in% names(ts_list[ts_bac_ix]))

ts_N_df <- data.frame(Infection = c('Virus', 'Virus', 'Bacteria', 'Bacteria'), 
                      Type = c('Time Series', 'Snapshot', 'Time Series', 'Snapshot'),
                      N = c(N_vir, N_vir_snapshot, N_bac, N_bac_snapshot))
p2d <- ts_N_df %>%
  ggplot(aes(x = Type, y = N, fill = Infection)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_clean + 
  labs(x = 'Study Design', y = 'Number of Studies')

### Platforms by pathogen type
parsePlatform <- function(str){
  str_extract(pattern = 'GPL\\d+', string = str) %>%
    return()
}
getPlatformsFromList <- function(input_list, type){
  input_list %>%
    names() %>%
    parsePlatform() %>%
    table() %>%
    sort(decreasing = T) %>%
    as.data.frame() %>%
    mutate('Type' = type) %>%
    return()
}
vir_platforms <- getPlatformsFromList(data_list[vir_ix], 'Virus')
bac_platforms <- getPlatformsFromList(data_list[bac_ix], 'Bacteria')

platform_df <- bind_rows(vir_platforms, bac_platforms)
names(platform_df)[1] <- 'Platform'

annotation_df <- openxlsx::read.xlsx('~/Documents/darpa-echo/Input/platform_specific/platform_types.xlsx')
annotation_df <- platform_df %>%
  dplyr::select(Platform) %>% 
  left_join(annotation_df) %>%
  distinct()
appending_platforms <- 
data.frame('Platform' = c('GPL6254',   'GPL14604', 'GPL15615', 'GPL17586', 'GPL17692', 'GPL21947', 'GPL23126', 'GPL8300', 'GPL9392', 'GPL13287', 'GPL16951', 'GPL20844', 'GPL6254'),
       'Manufacturer' = c('other',     'affy',     'other',    'affy',     'affy',     'other',    'affy',     'affy',    'affy',    'other',    'other',    'other',    'other'))
annotation_df <- openxlsx::read.xlsx('~/Documents/darpa-echo/Input/platform_specific/platform_types.xlsx') %>%
  filter(!is.na(Manufacturer)) %>% 
  bind_rows(appending_platforms)

platform_df <- platform_df %>%
  left_join(annotation_df)
platform_df %>% filter(is.na())

manufacturer_df <- platform_df %>%
  mutate(Manufacturer = ifelse(!Manufacturer %in% c('ilmn', 'affy'), 'other', Manufacturer)) %>% 
  group_by(Manufacturer, Type) %>%
  summarize('Freq' = sum(Freq)) 

p2e <- ggplot(manufacturer_df, aes(x = Type, y = Freq, fill = Manufacturer)) + 
  geom_bar(stat = 'identity') + 
  theme_clean + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(y = 'n. of studies')
p2e
p2c + p2d + p2e
setwd('~/Dropbox/compendium_manuscript/figures/ggplot objects/')
saveRDS(object = p2c, 'p2c.RDS')
saveRDS(object = p2d, 'p2d.RDS')
saveRDS(object = p2e, 'p2e.RDS')
