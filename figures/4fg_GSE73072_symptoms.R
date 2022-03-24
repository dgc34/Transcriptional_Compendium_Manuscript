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

plotTrajectory <- function(plot_dat, fname, h = 24, w = 15){
  grout <- ggplot(plot_dat, aes(x = Time.Point, y = Score, group = SubjectID, color = Symptoms)) + 
    geom_point(alpha = 0.5) + 
    geom_line(alpha = 0.5) + 
    geom_smooth(method = 'loess', aes(group = Symptoms)) + 
    facet_grid(Signature  ~ Study, scales= 'free', space = 'free') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text.y.right = element_text(angle = 0))
  #ggsave(filename = fname, plot = grout, height = h, width = w) 
  return(grout)
}

plotTrajectory(plot_dat_scaled, '~/Dropbox/gse73072_asymptomatic_trends.png')
for(study in unique(plot_dat_scaled$Study)){
  if(study == "GSE73072_GPL14604_HRV DUKE"){
    p4f <- plotTrajectory(plot_dat = plot_dat_scaled %>% 
                            filter(grepl(pattern = study, x = Study), Signature == 'V1'), 
                          fname = paste0('~/Dropbox/', study, '.png'), h = 6, w = 6)
  }
}

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

auc_plot_dat <- auc_dat %>%
  mutate(Symptoms = ifelse(Symptoms, 'symptomatic', 'asymptomatic')) %>%
  mutate(Study = ifelse(grepl(pattern = 'DEE3', x = Study), 'H1N1 cohort 1', Study)) %>% 
  mutate(Study = ifelse(grepl(pattern = 'DEE4', x = Study), 'H1N1 cohort 2', Study)) %>% 
  mutate(Study = ifelse(grepl(pattern = 'HRV', x = Study), 'hRV cohort 3', Study))
symp_colors <- c('gray 80', 'white')
names(symp_colors) <- c('symptomatic', 'asymptomatic')
p4i <- ggplot(auc_plot_dat, aes(x = Symptoms, y = AUC)) + 
  geom_boxplot() + 
  geom_point(aes(shape = Study), position = position_jitter(width = 0.1, height = 0, seed = 0208)) + 
  #ylim(c(0.5, 1)) +
  scale_fill_manual(values = symp_colors) + 
  theme_clean + 
  theme(legend.title = element_blank(), legend.position = 'bottom')  + 
  labs(x = '') + 
  scale_y_continuous(breaks = seq(to = 1, from = 0.5, by = 0.1), limits = c(0.45,1.05)) + 
  geom_hline(yintercept = 0.5, color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = AUC_threshold, color = 'red', linetype = 'dashed')
p4i
for(study_lab in unique(auc_dat$Study)){
  min_aucs <- auc_dat %>% 
    filter(grepl(pattern = study_lab, x = Study)) %>% 
    ungroup() %>% 
    group_by(Symptoms) %>% 
    summarize('med' = min(AUC)) %>% 
    pull(med)
  
  p1a <- ggplot(auc_dat %>% 
                  filter(grepl(pattern = study_lab, x = Study)), 
                aes(x = Symptoms, y = AUC)) +
    geom_violin() + 
    #geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(seed = 208, width = 0.15, height = 0)) + 
    #geom_point() + 
    #facet_grid(. ~ Signature, scales = 'free', space = 'free') + 
    theme_clean + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_y_continuous(breaks = seq(to = 1, from = 0, by = 0.1), limits = c(0,1.01)) + 
    geom_hline(yintercept = 0.5, color = 'red', linetype = 'dashed') + 
    geom_hline(yintercept = min_aucs[1], color = 'blue', linetype = 'dashed') + 
    geom_hline(yintercept = min_aucs[2], color = 'gray20', linetype = 'dashed') + 
    stat_compare_means(method = 'wilcox', label = 'p.signif')
  #p1a
  
  fc_dat <- summarized_dat %>%
    filter(Study %in% stat_studies) %>%
    filter(grepl(pattern = study_lab, x = Study)) %>%
    group_by(Signature, Study, SubjectID, Symptoms) %>%
    summarize('FC' = Max[Class == 'Virus'] - Max[Class == 'Healthy'])
  
  label_dat <- fc_dat %>%
    group_by(Study, Signature) %>%
    summarize('p.value' = wilcox.test(x = FC[Symptoms == T], y = FC[Symptoms == F], alternative = 'greater')$p.value) %>%
    mutate('pval' = p.value) %>%
    mutate(p.value = paste0('p = ', round(p.value, 3)), Symptoms = T)
  
  max_vals <- fc_dat %>%
    filter(Symptoms == T) %>%
    group_by(Study, Signature) %>%
    summarize(FC = max(FC))
  
  label_dat <- label_dat %>%
    left_join(max_vals) %>%
    mutate('sig_lab' = ifelse(pval < 0.1, '', '')) %>%
    mutate('sig_lab' = ifelse(pval < 0.05, '*', sig_lab)) %>%
    mutate('sig_lab' = ifelse(pval < 0.01, '**', sig_lab)) %>%
    mutate(FC = 1.02 * FC)
  
  p1b <- ggplot(fc_dat, aes(x = Symptoms, y = FC)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    facet_wrap(. ~ Signature, scales = 'free', nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_text(data = label_dat, aes(label = sig_lab), nudge_y = 0, nudge_x = -0.5, size = 10) + 
    labs(y = 'Signature Score Difference\n Peak vs. Pre-infection')
  #p1b
  
  
  pout <- p1a + p1b + plot_layout(widths = c(1,6))
  plot(pout)
  if(grepl(pattern = 'HRV DUKE', x = study_lab)){
    fig_out_1a <- p1a 
    fig_out_1b <- p1b
    saveRDS(object = fig_out_1a, file = '~/Dropbox/compendium_manuscript/figures/ggplot objects/4f.RDS')
    saveRDS(object = fig_out_1b, file = '~/Dropbox/compendium_manuscript/figures/ggplot objects/4g.RDS')
  }
  ggsave(pout, filename = paste0('~/Dropbox/compendium_manuscript/figures/extra/', study_lab, '.png'), width = 14, height = 5)
}

fc_dat <- summarized_dat %>%
  group_by(Signature, Study, SubjectID, Symptoms) %>%
  summarize('FC' = Max[Class == 'Virus'] - Max[Class == 'Healthy'])

label_dat <- fc_dat %>%
  group_by(Study, Signature) %>%
  summarize('p.value' = wilcox.test(x = FC[Symptoms == T], y = FC[Symptoms == F], alternative = 'greater')$p.value) %>%
  mutate('pval' = p.value) %>%
  mutate(p.value = paste0('p = ', round(p.value, 3)), Symptoms = T)

max_vals <- fc_dat %>%
  filter(Symptoms == T) %>%
  group_by(Study, Signature) %>%
  summarize(FC = max(FC))

label_dat <- label_dat %>%
  left_join(max_vals) %>%
  mutate('sig_lab' = ifelse(pval < 0.1, '', '')) %>%
  mutate('sig_lab' = ifelse(pval < 0.05, '*', sig_lab)) %>%
  mutate('sig_lab' = ifelse(pval < 0.01, '**', sig_lab)) %>%
  mutate(FC = 1.02 * FC)

grid_dat <- label_dat %>% 
  dplyr::select(Study, Signature, sig_lab) %>%
  mutate(Signature = factor(Signature, levels = paste0('V', 1:15), ordered = T)) %>%
  pivot_wider(id_cols = Study, names_from = Signature, values_from = sig_lab) %>%
  filter(Study %in% stat_studies)
#arrange(Signature)
grid_dat <- grid_dat %>%
  dplyr::select(Study, contains('V')) %>%
  mutate(Study = as.character(Study)) %>%
  mutate(Study = gsub(Study, pattern = 'GSE73072_GPL14604_', replacement = '')) %>%
  mutate(Study = gsub(Study, pattern = 'DEE', replacement = 'Cohort ')) %>%
  mutate(Study = gsub(Study, pattern = 'DUKE', replacement = 'Cohort 6')) %>%
  mutate(Study = gsub(Study, pattern = 'UVA', replacement = 'Cohort 7')) %>%
  mutate(Study = factor(Study, levels = c('RSV Cohort 1', 'H3N2 Cohort 2', 'H1N1 Cohort 3', 'H1N1 Cohort 4', 'H3N2 Cohort 5', 'HRV Cohort 6', 'HRV Cohort 7'), ordered = T)) %>%
  arrange(Study)
tt = ttheme_default(colhead=list(fg_params=list(rot=0)))
p2 <- tableGrob(grid_dat, theme = tt)
saveRDS(object = p2, file = '~/Dropbox/compendium_manuscript/figures/ggplot objects/4h.RDS')
mat <- c('
         ABBBBBB
         ABBBBBB
         ABBBBBB
         CCCCCCC
         CCCCCCC
         ')
pout <- fig_out_1a + fig_out_1b + p2 + plot_layout(design = mat)
ggsave(pout, filename = '~/Dropbox/compendium_manuscript/figures/4ghi.png', height = 7, width = 10, dpi = 300)
