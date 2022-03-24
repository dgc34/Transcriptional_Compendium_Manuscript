# the goal of this figure is to evaluate all signatures in all appropriate datasets
library(dplyr)
library(ggplot2)
library(openxlsx)
library(ggpubr)
library(stringr)

source('~/Dropbox/plot_palette.R')

sig_order <- lapply(c('VxB', 'V', 'B'), function(x){
  paste0(x, 1:15)}) %>%
  unlist()

gen_df <- readRDS('~/Dropbox/big_eval_heatmap_list_generalizability_021322_healthy_sbj_ts.RDS') %>%
  bind_rows() %>%
  #filter(Type == 'Discovery') %>%
  mutate(Comparison = str_extract(string = Signature, pattern = 'VxB|V|B')) %>%
  mutate(Comparison = factor(Comparison, levels = c('VxB', 'V', 'B'), labels = c('VxB', 'Virus', 'Bacteria'), ordered = T)) %>%
  mutate(Signature = factor(Signature, levels = sig_order, ordered = T)) %>%
  filter(N > 4) %>%
  filter(Comparison != 'VxB', !TS)

ann1 <- read.xlsx('~/Dropbox/bacterial_acute_annotation.xlsx') %>%
  dplyr::select(Accession, Time)
ann2 <- read.xlsx('~/Dropbox/viral_timing_annotations.xlsx') %>%
  dplyr::select(Accession, Time)
ann3 <- read.xlsx('~/Dropbox/unannotated_time.xlsx') %>%
  dplyr::select(Accession, Time)
ann <- bind_rows(ann1, ann2, ann3) %>%
  distinct() #%>%
  # group_by(Accession) %>%
  # mutate(N = n()) %>%
  # filter(N > 1)
study_plot_dat <- gen_df %>%
  left_join(ann) #%>%
  #filter(!grepl(pattern = 'GSE36238|GSE13699|GSE74027', x = Accession))#%>%
study_plot_dat %>%
  filter(is.na(Time)) %>%
  dplyr::select(Study, Comparison) %>%
  distinct()

ggplot(study_plot_dat %>%
         filter(Comparison == 'Virus', Time != 'Acute & Chronic'), 
       aes(x = Time, y = Scores)) + 
  geom_boxplot() + 
  facet_grid(. ~ Signature) +
  stat_compare_means()

fig4e_plot_dat <- study_plot_dat %>%
  group_by(Signature, Time, Comparison) %>% # group by study so each observation is a signature
  filter(Comparison == 'Virus', Time != 'Acute & Chronic', Signature != 'V9') %>% # remove V9 because there are no chronic datasets with this control group
  summarize(Scores = median(Scores, na.rm = T)) %>% # find the median signature scores
  mutate(Time = factor(Time, levels = c('Acute', 'Chronic', 'Acute & Chronic'), ordered = T))

p4h <- ggplot(fig4e_plot_dat,
             aes(x = Time, y = Scores)) +
  theme_clean + 
  geom_boxplot() +
  #facet_grid(. ~ Signature, space = 'free', scales = 'free') +
  #facet_grid(. ~ Comparison, space = 'free', scales = 'free') + 
  stat_compare_means(method = 'wilcox.test', position = position_nudge(x = 0.2, y = -0.02)) +
  theme(axis.text.x = element_text(angle= 90, hjust = 1)) +
  geom_hline(yintercept = 0.5, color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = AUC_threshold, color = 'red', linetype = 'dashed') + 
  scale_y_continuous(breaks = seq(to = 1, from = 0.5, by = 0.1), limits = c(0.45,1.05)) + 
  labs(y = 'Average Dataset AUC', x = 'Infection Timing')
p4h
saveRDS(object = p4h, file = '~/Dropbox/ggplot_objects/fig4_020922/fig4h.RDS')
ggsave(filename = '~/Dropbox/compendium_manuscript/figures/4e_acuteVSchronic.png', plot = p4, height = 6, width = 3, dpi = 300)


fig4e_plot_dat %>% group_by(Time, Comparison) %>% summarize('med' = median(Scores))
