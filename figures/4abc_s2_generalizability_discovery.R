library(dplyr)
library(ggplot2)
library(stringr)
library(openxlsx)
library(ggridges)
library(AnnotationDbi)
library(org.Hs.eg.db)

source('~/Dropbox/plot_palette.R')
FLIP_V10_V11 <- F

sig_order <- lapply(c('VxB', 'V', 'B'), function(x){
  paste0(x, 1:15)}) %>%
  unlist()

input_file <- '~/Dropbox/big_eval_heatmap_list_generalizability_021322_healthy_sbj_ts.RDS'
gen_df <- readRDS(input_file) %>%
  bind_rows() %>%
  #filter(Type == 'Discovery') %>%
  filter(!grepl(pattern = 'I', x = Signature)) %>% 
  mutate(Comparison = str_extract(string = Signature, pattern = 'VxB|V|B')) %>%
  mutate(Comparison = factor(Comparison, levels = c('VxB', 'V', 'B'), labels = c('VxB', 'Virus', 'Bacteria'), ordered = T)) %>%
  mutate(Signature = factor(x = Signature, levels = sig_order, ordered = T))

contrast <- str_extract(pattern = 'healthy|noninf|orig', string = input_file)

gen_df %>% filter(Signature == 'VxB1', Accession == 'GSE6269')

ggplot(gen_df %>% filter(Type != 'Discovery', grepl(pattern = '^B', x = Signature), N > 4), aes(x = Signature, y = Scores)) + 
  geom_boxplot() + geom_hline(yintercept = 0.7) 

if(FLIP_V10_V11){
  table(gen_df$Signature)
  gen_df$Signature[gen_df$Signature == 'V11'] <- 'V12'
  gen_df$Signature[gen_df$Signature == 'V10'] <- 'V11'
  gen_df$Signature[gen_df$Signature == 'V12'] <- 'V10'
  table(gen_df$Signature)
}

# verify that sweeney AKA V10 is IFI27, JUP, LAX1
gen_df %>% filter(Signature == 'V10') %>% pull(pos_genes) %>% unique() %>% strsplit(split = ' ') %>% unlist() %>% mapIds(x = org.Hs.eg.db, keytype = 'ENTREZID', column = 'SYMBOL')

modified_disc <- read.xlsx('~/Dropbox/compendium_manuscript/tables/supp_1.xlsx') %>%
  filter(grepl(pattern = 'GSE73072', x = Discovery.Accessions))

for (sig in 1:nrow(modified_disc)){
  signature_row <- modified_disc[sig, ]
  sdys <- str_split(signature_row$Discovery.Accessions, pattern = ';') %>% unlist()
  gen_df <- gen_df %>%
    mutate(Type = ifelse((Signature == signature_row$Signature.Label & Study %in% sdys), 'Discovery', Type))
}

gen_df <- gen_df %>%
  filter(N >= 5) %>%
  mutate(Signature = str_replace(pattern = 'VxB', replacement = 'V/B', string = Signature)) %>% 
  mutate(Comparison = str_replace(pattern = 'VxB', replacement = 'V/B', string = Comparison)) 



disc_df <- gen_df %>%
  filter(Type == 'Discovery')

nondisc_df <- gen_df %>%
  filter(Type %in% c('Validation', 'New'))

outlier.studies <- nondisc_df %>%
  filter(N > 4) %>% 
  group_by(Comparison, Accession) %>% 
  summarize('med' = median(Scores)) %>% 
  filter(med < 0.5) 

annotateOutliers <- function(input_df, outlier.studies){
  outlier.studies <- outlier.studies %>%
    mutate(Pt.Color = Accession)
  output_df <- input_df %>%
    left_join(outlier.studies) %>%
    mutate('Pt.Color' = ifelse(is.na(Pt.Color), '', Pt.Color)) %>%
    mutate(Pt.Color = factor(Pt.Color)) 
  return(output_df)
}

nondisc_df2 <- annotateOutliers(nondisc_df, outlier.studies)
disc_df2 <- annotateOutliers(disc_df, outlier.studies)

plotGenDf <- function(input_df, text_flag = T, AUC_threshold = 0.6, fname){
  color_dat <- input_df %>%
    group_by(Signature) %>%
    summarize('med' = median(Scores)) %>% # find the fraction of values below this point
    arrange(med) %>%
    mutate('Specificity' = factor(ifelse(med <= AUC_threshold, 'not robust', 'robust'))) %>% # label signatures with more than 25% of their values above the threshold red
    dplyr::select(-med)
  color_vals <- c("#FFFFFF", "gray50")#, "#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#F0E442")
  shape_vals <- c(8, 3, 4, 16, 17)
  names(shape_vals) <- unique(outlier.studies$Accession)
  names(color_vals) <- c('robust', 'not robust')
  #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  #color_vals <- c('black', 'red3')
  
  
  text_df <- input_df %>%
    group_by(Comparison, Signature) %>%
    summarize('sig_med' = median(Scores)) %>%
    ungroup() %>%
    group_by(Comparison) %>%
    summarize('comp_med' = median(sig_med), 'N_sigs' = n()) %>%
    mutate(Signature = N_sigs/1.6) %>%
    mutate(Scores = 0.3) %>%
    mutate(plot_label = paste0('median AUC: ', round(comp_med, 2)))
  
  input_df <- input_df %>%
    left_join(color_dat)
  
  text_geom <- theme()
  if(text_flag){
    text_geom <- geom_text(data = text_df, aes(label = plot_label, color = NULL)) 
  } 
  
  pout <- ggplot(input_df %>% 
                   mutate(Signature = factor(Signature, levels = gsub(pattern = 'x', replacement = '/', x = sig_order), ordered = T)) %>%
                   mutate(Comparison = factor(Comparison, levels = c('V/B', 'Virus', 'Bacteria'), ordered = T)), 
                 aes(x = Signature, y = Scores)) + 
    geom_boxplot(outlier.shape = NA, aes(fill = Specificity, alpha = 0.5)) + 
    geom_point(data = input_df %>% filter(Pt.Color != ''), aes(shape = Pt.Color)) + 
    #text_geom +
    geom_hline(yintercept = c(AUC_threshold), color = 'red', linetype = 'dashed') + 
    geom_hline(yintercept = c(0.5), color = 'black', linetype = 'dashed') + 
    facet_grid(. ~ Comparison, scales = 'free', space = 'free') + 
    scale_shape_manual(values = shape_vals, name = 'Outlier Points') + 
    scale_fill_manual(values = color_vals[1:2], name = 'Performance') + 
    scale_alpha(guide = 'none') +
    theme_clean + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = 'right',
          panel.spacing = unit(2, "lines")) + 
    scale_y_continuous(breaks = seq(to = 1, from = 0, by = 0.1), limits = c(0,1)) + 
    labs(y = 'AUC') 
  
  #ggsave(fname, width = 8, height = 5, dpi = 300)
  return(list(input_df, pout))
}


s2 <- plotGenDf(disc_df2, text_flag = T, AUC_threshold = 0.7, fname = paste0('~/Dropbox/compendium_manuscript/figures/s2_generalizability_discovery_', contrast, '.png'))
pabc <- plotGenDf(nondisc_df2, text_flag = F, AUC_threshold = 0.7, fname = paste0('~/Dropbox/compendium_manuscript/figures/4abc_generalizability_', contrast, '.png'))
#saveRDS(object = pabc[[2]], file = '~/Dropbox/compendium_manuscript/figures/ggplot objects/4abc.RDS')
pabc[[2]]

color_vals <- c('#F2BE84', '#C4E8F0', '#8FD0BA', 'gray50')
names(color_vals) <- c('B', 'V', 'V/B', 'not robust')
sig_order <- c(
  paste0('V', 1:11), paste0('B', 1:7), paste0('V/B', 1:7)
) %>% rev()

AUC_threshold <- 0.7
plot_dat <- pabc[[1]] %>% 
  mutate(Comparison = ifelse(Comparison == 'Bacteria', 'B', Comparison)) %>% 
  mutate(Comparison = ifelse(Comparison == 'Virus', 'V', Comparison)) %>%
  mutate(Comparison = factor(Comparison, levels = c("V", "B", "V/B"), ordered = T)) %>%
  mutate(Signature = factor(Signature, levels = sig_order, ordered = T)) %>%
  mutate(Specificity = ifelse(Specificity == 'robust', as.character(Comparison), 'not robust'))
shape_vals <- c(2, 4, 1, 16, 5)
names(shape_vals) <- outlier.studies$Accession

point_color_vals <- c(rep('black', length(outlier.studies$Accession)),  NA)
names(point_color_vals) <- c(outlier.studies$Accession, '')

op3 <- list()
for (comp in unique(plot_dat$Comparison)){
  op3[[comp]] <- ggplot(plot_dat %>% filter(Comparison == comp), 
                        aes(x = Scores, y = Signature, fill = Specificity, color = Specificity)) +
    
    #geom_density_ridges(stat = 'binline', binwidth = 0.05) +
    #geom_density_ridges(quantile_lines = TRUE, quantiles = 2) +
    geom_density_ridges(alpha = 0.9) +
    #geom_point(data = plot_dat %>% filter(Pt.Color != '', Comparison == comp), aes(shape = Pt.Color)) + 
    #geom_vline(xintercept = AUC_threshold, color = 'red', linetype = 'dashed') + 
    #geom_vline(xintercept = 0.5, color = 'black', linetype = 'dashed') + 
    #scale_shape_manual(values = shape_vals, name = 'Outlier Points') + 
    #facet_grid(Comparison ~ ., space = 'free', scales = 'free') +
    scale_fill_manual(values = color_vals) +
    scale_color_manual(values = color_vals) +
    #xlim(c(0,1)) + 
    labs(x = 'AUROC', y = '') +
    theme_clean +
    theme(legend.position = 'none') + 
    scale_x_continuous(breaks = seq(to = 1, from = 0, by = 0.2), limits = c(0,1))
}
op3[[3]]
op1 <- ggplot(plot_dat, aes(x = Scores, y = Signature, fill = Specificity)) +
  #geom_vline(xintercept = AUC_threshold, color = 'red', linetype = 'dashed') + 
  #geom_vline(xintercept = 0.5, color = 'black', linetype = 'dashed') + 
  #geom_density_ridges(stat = 'binline', binwidth = 0.05) +
  geom_density_ridges() +
  facet_grid(Comparison ~ ., space = 'free', scales = 'free') +
  scale_fill_manual(values = color_vals) +
  xlim(c(0,1)) + 
  labs(x = 'AUC') +
  theme_clean +
  theme(legend.position = 'none')
if(contrast == 'noninf'){
  op1 <- ggplot(plot_dat %>%
                  filter(Comparison != 'V/B'), 
                aes(x = Scores, y = Signature, fill = Specificity, color = Specificity)) +
    #geom_vline(xintercept = AUC_threshold, color = 'red', linetype = 'dashed') + 
    #geom_vline(xintercept = 0.5, color = 'black', linetype = 'dashed') + 
    #geom_density_ridges(stat = 'binline', binwidth = 0.05) +
    geom_density_ridges(alpha = 0.95) +
    facet_grid(Comparison ~ ., space = 'free', scales = 'free') +
    scale_fill_manual(values = color_vals) +
    scale_color_manual(values = color_vals) +
    xlim(c(0,1)) + 
    labs(x = 'AUROC (validation)') +
    theme_clean +
    theme(legend.position = 'none')
}
#saveRDS(object = op3, file = paste0('~/Dropbox/ggplot_objects/fig4_020922/fig4abc_0217_', contrast, '.RDS'))

s_table_4 <- gen_df %>%
  dplyr::select(Study, Signature, Comparison, Type, Scores, TS, N_neg, N_pos) %>%
  dplyr::rename(AUC = Scores, N_samples_positive = N_pos, N_samples_negative = N_neg, Time.Series = TS)
head(s_table_4)
#write.csv(x = s_table_4, file = paste0('~/Dropbox/compendium_manuscript/tables/supp4_', contrast, '.csv'), quote = F, row.names = F)


#shape_vals <- c(2, 4, 1)
#names(shape_vals) <- outlier.studies$Accession
#color_vals <- c('#F2BE84', '#C4E8F0', '#8FD0BA', 'gray50')
#names(color_vals) <- c('B', 'V', 'V/B', 'not robust')
#AUC_threshold = 0.7
#op3_disc <- list()
disc_df <- s2[[1]] %>% 
  mutate(Comparison = str_extract(string = Signature, pattern = 'V/B|^V|^B')) %>%
  mutate(Comparison = factor(Comparison, levels = c('V', 'B', 'V/B'), ordered = T)) 
extra_rows <- disc_df %>%
  filter(Signature == 'V1') %>%
  mutate(Signature = 'V/B1', Comparison = 'V/B')
disc_df <- disc_df %>%
  bind_rows(extra_rows)
op3_disc <- list()
for (comp in unique(disc_df$Comparison)){
  op3_disc[[comp]] <- ggplot(disc_df %>% 
                               filter(Comparison == comp), 
                        aes(x = Scores, y = Signature, fill = Comparison)) +
    
    #geom_density_ridges(stat = 'binline', binwidth = 0.05) +
    #geom_density_ridges(quantile_lines = TRUE, quantiles = 2) +
    #geom_boxplot() + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(pch = 21, size = 5) + 
    #geom_density_ridges() +
    #geom_boxplot() + 
    #geom_boxplot(outlier.shape = NA) + 
    
    #geom_point(data = s2[[1]] %>% filter(Pt.Color != '', Comparison == comp), aes(shape = Pt.Color)) + 
    geom_vline(xintercept = AUC_threshold, color = 'red', linetype = 'dashed') + 
    geom_vline(xintercept = 0.5, color = 'black', linetype = 'dashed') + 
    #scale_shape_manual(values = shape_vals, name = 'Outlier Points') + 
    #facet_grid(Comparison ~ ., space = 'free', scales = 'free') +
    scale_fill_manual(values = color_vals) +
    #scale_color_manual(values = color_vals) + 
    xlim(c(0,1)) + 
    labs(x = 'AUC', y = '') +
    theme_clean +
    theme(legend.position = 'none')
}
op3_disc[[3]]

pout <- ggplot(disc_df %>% 
                 mutate(Signature = factor(Signature, levels = gsub(sig_order, pattern = 'x', replacement = '/'), ordered = T)) %>%
                 mutate(Comparison = factor(Comparison, levels = c('V', 'B', 'V/B'), ordered = T)), 
       aes(x = Scores, y = Signature, fill = Comparison)) +
  
  geom_boxplot(outlier.shape = NA) + 
  geom_point(pch = 21, size = 3) + 
  
  geom_vline(xintercept = AUC_threshold, color = 'red', linetype = 'dashed') + 
  geom_vline(xintercept = 0.5, color = 'black', linetype = 'dashed') + 
  
  facet_grid(Comparison ~ ., space = 'free', scales = 'free') +
  scale_fill_manual(values = color_vals) +
  xlim(c(0,1)) + 
  labs(x = 'AUROC (discoery)', y = '') +
  theme_clean +
  theme(legend.position = 'none', strip.background = element_rect(fill = 'white'))
pout
saveRDS(pout, file = paste0('~/Dropbox/ggplot_objects/3b_discovery_AUROCs_', contrast, '.RDS'))
