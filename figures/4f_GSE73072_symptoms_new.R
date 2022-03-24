# load time series studies
# test all signatures, get geomMeans for all studies (don't split symptomatic and asymptomatic before z-scoring)
# group by subject, Class
# mutate max - Healthy
# compare between two groups for each study and each signature
library(dplyr)
library(ggplot2)
library(MetaIntegrator)
library(scales)
library(gridExtra)
library(tidyr)
library(patchwork)
library(ggpubr)
source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))
source('~/Dropbox/time_series_fxn_manuscript.R')
source('~/Dropbox/plot_palette.R')

signatures_file <- '~/Dropbox/compendium_manuscript/tables/supp_1.xlsx'
signatures <- read.xlsx(signatures_file) %>%
  filter(Type == 'Virus') 

accessions_of_interest <- 'GSE73072'
input_file <- '~/Dropbox/GSE73072_v1.RDS'
if(!file.exists(input_file)){
  data_list <- readRDS('~/Documents/darpa-manuscript-data/data_list_final/data_list_v1.RDS')
  keep_ix <- grepl(pattern = paste(accessions_of_interest, collapse = '|'), x = names(data_list))
  data_list <- data_list[keep_ix]
  saveRDS(object = data_list, file = input_file)
} else {
  data_list <- readRDS(input_file)
}

names(data_list)

# generate scores for each signature and each subject
data_list <- data_list %>%
  lapply(function(MIobj){
    MIobj$class <- as.numeric(MIobj$pheno$Class == 'Virus')
    names(MIobj$class) <- rownames(MIobj$pheno)
    if(checkDataObject(MIobj, 'Dataset')){
      return(MIobj)
    } else {
      stop()
    }
  })

full_df <- list()
for(i in 1:nrow(signatures)){
  current_sig <- sig2Meta(signatures[i,])
  score_df <- lapply(data_list, function(MIobj){
    df <- MIobj$pheno %>%
      dplyr::select(SubjectID, Time.Point, Symptoms, Class) %>%
      mutate('Score' = calculateScoreRobust(datasetObject = MIobj, filterObject = current_sig, method = 'geomMean', zScore = F)) %>%
      mutate('Study' = MIobj$formattedName, 'Signature' = current_sig$filterDescription)
    return(df)
  })
  full_df[[i]] <- score_df %>% bind_rows()
}
plot_dat <- bind_rows(full_df) %>%
  mutate(Signature = factor(Signature, levels = paste0('V', 1:15), ordered = T))

## verify that asymptomatic individuals still demonstrate an increase in viral signature scores
## compared to baseline (requires z-scoring for facet grid. could facet wrap, but relative changes are preserved either way)
plot_dat_scaled <- plot_dat %>%
  group_by(Signature, Study) %>%
  mutate(Score = scale(Score)) %>%
  mutate(N = n())



#saveRDS(object = p4f, file = '~/Dropbox/compendium_manuscript/figures/ggplot objects/4f_time_series.RDS')
### AUCs asymptomatic vs symptomatic
# take the maximum value over time for Virus and Healthy
summarized_dat <- plot_dat %>%
  group_by(SubjectID, Study, Signature, Symptoms, Class) %>%
  summarize('Max' = max(Score)) %>%
  ungroup() %>%
  group_by(SubjectID, Study, Signature, Symptoms) %>%
  mutate(N = n()) %>%
  filter(N == 2)

N_dat <- plot_dat %>%
  filter(Signature == 'V1') %>%
  dplyr::select(SubjectID, Study, Signature, Symptoms) %>%
  distinct() %>%
  group_by(Study, Signature, Symptoms) %>%
  summarize('N' = n()) %>%
  group_by(Study) %>%
  summarize('N' = min(N))

stat_studies <- N_dat %>%
  filter(N >= 4) %>%
  pull(Study)

auc_dat <- summarized_dat %>%
  filter(Study %in% stat_studies) %>%
  group_by(Symptoms, Study, Signature) %>%
  summarize('AUC' = auroc(bool = Class == 'Virus', score = Max))

study_lab <- 'DEE5'

AUC_threshold = 0.7
auc_plot_dat <- auc_dat %>%
  mutate(Symptoms = ifelse(Symptoms, 'symptomatic', 'asymptomatic')) %>%
  mutate(Symptoms = factor(Symptoms, levels = c('symptomatic', 'asymptomatic'), ordered = T)) %>% 
  mutate(Study = ifelse(grepl(pattern = 'DEE3', x = Study), 'H1N1 cohort 1', Study)) %>% 
  mutate(Study = ifelse(grepl(pattern = 'DEE4', x = Study), 'H1N1 cohort 2', Study)) %>% 
  mutate(Study = ifelse(grepl(pattern = 'HRV', x = Study), 'hRV cohort 3', Study))
symp_colors <- c('gray 80', 'white')
names(symp_colors) <- c('symptomatic', 'asymptomatic')
p4i <- ggplot(auc_plot_dat %>%
                filter(grepl(pattern = 'cohort 3', x = Study)), aes(x = Symptoms, y = AUC)) + 
  #geom_line(aes(group = Signature),color = 'gray80') + 
  geom_boxplot(outlier.shape = NA) + 
  #geom_point(aes(shape = Study), position = position_jitter(width = 0.1, height = 0, seed = 0208)) +
  #geom_point(aes(shape = Study)) +
  geom_point() + 
  #ylim(c(0.5, 1)) +
  #geom_line(aes(group = Signature),color = 'gray50', position = position_jitter(width = 0.1, height = 0, seed = 0208)) +
  
  scale_fill_manual(values = symp_colors) + 
  theme_clean + 
  theme(legend.title = element_blank(), legend.position = 'bottom')  + 
  labs(x = '') + 
  #facet_grid(. ~ Study, scales = 'free', space = 'free') + 
  scale_y_continuous(breaks = seq(to = 1, from = 0.5, by = 0.1), limits = c(0.45,1.05)) + 
  geom_hline(yintercept = 0.5, color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = AUC_threshold, color = 'red', linetype = 'dashed') + 
  theme(strip.placement = 'none') + 
  stat_compare_means(method = 'wilcox.test', paired = T, position = position_nudge(y = -0.025))
p4i
saveRDS(object = p4i, file = '~/Dropbox/ggplot_objects/fig4_020922/fig4i.RDS')

supp_p4i <- ggplot(auc_plot_dat %>%
                filter(!grepl(pattern = 'cohort 3', x = Study)), aes(x = Symptoms, y = AUC)) + 
  geom_boxplot() + 
  #geom_point(aes(shape = Study), position = position_jitter(width = 0.1, height = 0, seed = 0208)) +
  #geom_point(aes(shape = Study)) +
  geom_point() + 
  #ylim(c(0.5, 1)) +
  #geom_line(aes(group = Signature),color = 'gray50', position = position_jitter(width = 0.1, height = 0, seed = 0208)) +
  geom_line(aes(group = Signature),color = 'gray80', alpha = 0.5) + 
  scale_fill_manual(values = symp_colors) + 
  theme_clean + 
  theme(legend.title = element_blank(), legend.position = 'none')  + 
  labs(x = '') + 
  facet_grid(. ~ Study, scales = 'free', space = 'free') + 
  scale_y_continuous(breaks = seq(to = 1, from = 0.5, by = 0.1), limits = c(0.45,1.05)) + 
  #geom_hline(yintercept = 0.5, color = 'black', linetype = 'dashed') +
  #geom_hline(yintercept = AUC_threshold, color = 'red', linetype = 'dashed') + 
  theme(strip.background = element_rect(fill = 'white')) + 
  stat_compare_means(method = 'wilcox.test', paired = F, position = position_nudge(y = -0.01))
supp_p4i
pdf('~/Dropbox/compendium_manuscript/figures/supplements/supp6.pdf', height = 3.5, width = 5)
supp_p4i
dev.off()
auc_plot_dat %>%
  filter(grepl(pattern = 'cohort 1', x = Study), AUC > 0.8, Symptoms == 'asymptomatic')

auc_plot_dat %>%
  filter(grepl(pattern = 'cohort 3', x = Study), AUC > 0.8, Symptoms == 'asymptomatic')
