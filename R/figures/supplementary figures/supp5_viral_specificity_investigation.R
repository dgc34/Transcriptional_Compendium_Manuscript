library(ggplot2)
library(dplyr)
library(openxlsx)
library(ggthemes)
library(stringr)
library(scales)
library(ggh4x)
library(ggpubr)
library(ggradar)
library(ggridges)
library(grid)

FLIP_V10_V11 <- F

source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))
source('~/Dropbox/plot_palette.R')

filterGeneFraction <- function(input_df, frac_threshold){
  input_df %>% 
    mutate(frac_pos = ifelse(is.na(frac_pos), 1, frac_pos)) %>%
    mutate(frac_neg = ifelse(is.na(frac_neg), 1, frac_pos)) %>% 
    filter(frac_pos >= frac_threshold, frac_neg >= frac_threshold) %>%
    return()
}

getSignificantDf <- function(df2){
  df2 %>% 
    ungroup() %>%
    filter(N >= 4, !Signature %in% sigs_to_omit) %>% #, !is.na(P.Value)) %>%
    group_by(Signature, Pathogen) %>%
    #summarize("N_sig" = sum(P.Value <= 0.05, na.rm = T), "N" = sum(!is.na(P.Value)), 'median' = median(Scores, na.rm = T)) %>%
    summarize('median' = median(Scores, na.rm = T)) %>%
    mutate('Fraction' = N_sig/N) %>%
    ungroup() %>%
    #mutate('Specific' = ifelse(median <= 0.6, 1, 0)) %>%
    mutate('Specific' = ifelse(Fraction < 0.5, 1, 0)) %>%
    mutate(Specific = as.logical(Specific)) %>%
    mutate(Specific.Logical = as.logical(Specific)) %>%
    mutate(Specific = ifelse(Specific, 'Specific', 'Non-Specific')) %>%
    mutate(Specific = factor(Specific, levels = c('Specific', 'Non-Specific'), ordered = T)) %>%
    mutate('Slabel' = ifelse(Specific == 'Specific', '', '.'))  %>%
    # mutate(Specific = ifelse(Fraction == 0, 1, 0)) %>%
    # mutate(Specific.Logical = as.logical(Specific)) %>%
    # mutate(Specific = ifelse(Specific, 'Specific', 'Non-Specific')) %>%
    # mutate(Specific = factor(Specific, levels = c('Specific', 'Non-Specific'), ordered = T)) %>%
    return()
}

getOrder <- function(input_data, col_name){
  col_name <- enquo(col_name)
  order <- input_data %>%
    group_by(!!col_name) %>%
    summarize('S' = sum(Specific.Logical)) %>%
    arrange(desc(S)) %>%
    pull(!!col_name)
  return(order)
}

### virus
# load phenodata labeling
phenoData <- read.xlsx(paste0(Sys.getenv('DARPA'), 'Input/classLabels/phenoData.xlsx'))
# load signatures
#signatures_file <- paste0(Sys.getenv("DARPA"), '/Input/signatures/bacterial_viral_healthy.xlsx')
signatures_file <- '~/Dropbox/compendium_manuscript/tables/supp_1.xlsx'
signatures <- read.xlsx(signatures_file)

#plot_file <- '~/Documents/darpa-manuscript-data/figures/bact_vir_specificity_heatmap_6a_v2.RDS'
plot_file <- '~/Dropbox/big_eval_heatmap_list_specificity_021322_healthy_sbj_ts.RDS'
out <- readRDS(plot_file)

plot_dat <- bind_rows(out) %>%
  filter(Type != 'Null') %>% # remove permuted labels
  mutate('Type' = factor(Type, levels = c('Discovery', 'Validation', 'New', 'Null'), ordered = T)) %>%
  filter(!grepl(pattern = 'NA|DARPA|Chawla', x = Signature)) %>% # remove the single gene testing variants and DARPA signatures
  mutate(Comparison = str_extract(pattern = 'VxB|V|B', string = Signature)) # identify signature comparison

if(FLIP_V10_V11){
  table(plot_dat$Signature)
  plot_dat$Signature[plot_dat$Signature == 'V11'] <- 'V12'
  plot_dat$Signature[plot_dat$Signature == 'V10'] <- 'V11'
  plot_dat$Signature[plot_dat$Signature == 'V12'] <- 'V10'
  table(plot_dat$Signature)
  
}
plot_dat %>% 
  filter(Signature == 'V10') %>% 
  pull(pos_genes) %>% 
  unique() %>% 
  strsplit(split = ' ') %>% 
  unlist() %>% 
  mapIds(x = org.Hs.eg.db, keytype = 'ENTREZID', column = 'SYMBOL')

table(plot_dat$Comparison)
table(plot_dat$Signature)

N_threshold = 5
colors <- createColorScale()

input_tmp_plot_dat <- plot_dat %>%
  filter(N >= N_threshold, Type != 'Null') %>%
  mutate(Signature = as.character(Signature))

order_dat <- input_tmp_plot_dat %>%
  group_by(Study) %>%
  summarize('N' = n()) %>%
  arrange(N)

signature_labels <- rev(c(paste0('VxB', 1:10), paste0('V', 1:12), paste0('B', 1:10)))

tmp_plot_dat <- input_tmp_plot_dat %>%
  mutate('Signature' = factor(Signature, levels = signature_labels, ordered = T)) %>%
  mutate('Study' = factor(Study, levels = order_dat$Study), ordered = T)

study_plot_dat <- tmp_plot_dat %>%
  mutate(Train_Label = ifelse(Type == 'Discovery', '○', '')) 

color_vals <- c("#F2BE84", "#C4E8F0", "#8FD0BA", "white")
names(color_vals) <- c('B', 'V', 'V/B', 'non-specific') # swap bacterial and viral for this case 
specificity_threshold = 0.6
specificity_plot_dat <- study_plot_dat %>%
  filter(!Signature %in% c('B1', 'B7', 'V3')) %>% 
  filter(!grepl(pattern = 'VxB', x = Signature)) %>%
  group_by(Signature) %>%
  mutate('med' = median(Scores)) %>% 
  mutate(Color = as.character(Comparison)) %>%
  mutate(Color = ifelse(med <= specificity_threshold, Color, 'non-specific')) %>% 
  #mutate(Signature = factor(Signature, levels = sig_order, ordered = T)) %>%
  mutate(Comparison = factor(Comparison, levels = c('V', 'B'), ordered = T))

p3_option1 <- ggplot(specificity_plot_dat, 
                     aes(x = Signature, y = Scores, fill = Color)) + 
  geom_boxplot() +
  #geom_point(position = position_jitter(seed = 1013, width = 0.1, height = 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15)) + 
  #labs(y = 'median AUC', title = 'Bacterial Signatures\nEvaluted in Viral Datasets') + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.03)) +
  labs(y = 'AUC') + 
  scale_fill_manual(values = color_vals) + 
  facet_grid(. ~ Comparison, scales = 'free', space = 'free') + 
  geom_hline(yintercept = specificity_threshold, linetype = 'dashed', color = 'blue') + 
  geom_hline(yintercept = c(0.5), linetype = 'dashed', color = 'black') + 
  theme_clean + 
  theme(legend.position = 'none')
#geom_text(data = fraction_labs, aes(label = round(Fraction, 2), y = 1.03), size = 4)


p3_option2 <- ggplot(specificity_plot_dat, aes(x = Scores, y = Signature, fill = Color, color = Comparison)) + 
  geom_density_ridges() + 
  #geom_vline(xintercept = specificity_threshold, linetype = 'dashed', color = 'blue') + 
  #geom_vline(xintercept = c(0.5), linetype = 'dashed', color = 'black') + 
  scale_fill_manual(values = color_vals) + 
  scale_color_manual(values = color_vals) + 
  scale_x_continuous(breaks = seq(to = 1, from = 0, by = 0.1), limits = c(0,1)) + 
  labs(x = 'AUC') + 
  facet_grid(Comparison ~ ., scales = 'free', space = 'free') +  
  #xlim(c(0,1)) + 
  theme_clean + 
  theme(legend.position = 'none', strip.background = element_rect(fill = 'white'), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

vir_annotation_df <- read.xlsx('~/Dropbox/annotated_pathogen_df_unique.xlsx')
specificity_plot_dat2 <- specificity_plot_dat %>%
  left_join(vir_annotation_df) %>%
  dplyr::select(Study, Signature, Envelope, Genome, N, Scores) %>%
  filter(grepl(pattern = 'B', x = Signature))
specificity_plot_dat2$Envelope %>% table()
specificity_plot_dat2 %>% filter(Envelope == '')
ggplot(specificity_plot_dat2 %>% filter(Envelope %in% c('enveloped', 'nonenveloped')), aes(x = Signature, y = Scores, fill = Envelope)) + 
  geom_boxplot()

color_vals <- c('#f2be84', 'white')
names(color_vals) <- c('specific', 'non-specific')
table(specificity_plot_dat2$Genome, specificity_plot_dat2$Signature)
pa <- ggplot(specificity_plot_dat2 %>% 
         group_by(Signature) %>% 
         mutate(Signature = factor(Signature, levels = c('B1', 'B2', 'B3', 'B4', 'B5', 'B6'), ordered = T)) %>% 
         mutate(color = ifelse(median(Scores) <= 0.6, 'specific', 'non-specific')), aes(x = Signature, y = Scores, fill = color)) + 
  geom_boxplot() + 
  scale_fill_manual(values = color_vals)  +
  facet_grid(. ~ Signature, space = 'free', scales = 'free') + 
  theme_clean + 
  geom_hline(yintercept = 0.6, linetype = 'dashed', color = 'gray50') + 
  labs(y = '') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  theme(legend.position = 'none')

pb <- ggplot(specificity_plot_dat2 %>% 
         mutate(Envelope = ifelse(!Envelope %in% c('enveloped', 'nonenveloped'), 'other', Envelope)) %>% 
         mutate(Envelope = factor(Envelope, levels = c('enveloped', 'nonenveloped', 'other'), ordered = T)) %>% 
           filter(!grepl(pattern = 'other', x = Envelope)) %>% 
         group_by(Signature, Envelope) %>% 
         mutate(color = ifelse(median(Scores) <= 0.6, 'specific', 'non-specific')) %>% 
         mutate(Signature = factor(Signature, levels = c('B1', 'B2', 'B3', 'B4', 'B5', 'B6'), ordered = T)), 
         aes(x = Envelope, y = Scores, fill = color)) + 
  geom_boxplot() + 
  scale_fill_manual(values = color_vals)  +
  facet_grid(. ~ Signature, space = 'free', scales = 'free') + 
  theme_clean + 
  geom_hline(yintercept = 0.6, linetype = 'dashed', color = 'gray50') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(y = 'AUROC (bacterial specificity)')+ 
  theme(legend.position = 'none')

pc <- ggplot(specificity_plot_dat2 %>% 
         mutate(Genome = ifelse(!Genome %in% c('dsDNA', 'ssRNA'), 'other', Genome)) %>% 
         mutate(Genome = factor(Genome, levels = c('dsDNA', 'ssRNA', 'other'), ordered = T)) %>%
           filter(!grepl(pattern = 'other', x = Genome)) %>% 
         group_by(Signature, Genome) %>% 
         mutate(color = ifelse(median(Scores) <= 0.6, 'specific', 'non-specific')) %>% 
         mutate(Signature = factor(Signature, levels = c('B2', 'B3', 'B4', 'B5', 'B6'), ordered = T)), 
         aes(x = Genome, y = Scores, fill = color)) + 
  geom_boxplot() + 
  scale_fill_manual(values = color_vals)  +
  facet_grid(. ~ Signature, space = 'free', scales = 'free') + 
  theme_clean + 
  geom_hline(yintercept = 0.6, linetype = 'dashed', color = 'gray50') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(y = '') + 
  theme(legend.position = 'none') 

pdf('~/Dropbox/compendium_manuscript/figures/supplements/supp7.pdf', height = 6, width = 7)
pb / pc
dev.off()


tmp <- specificity_plot_dat2 %>% 
  group_by(Signature, Genome) %>% 
  mutate(color = ifelse(median(Scores) <= 0.6, 'specific', 'non-specific')) %>% 
  mutate(Signature = factor(Signature, levels = c('B1', 'B2', 'B3', 'B4', 'B5', 'B6'), ordered = T)) %>% 
  filter(Genome %in% c('dsDNA', 'ssRNA'))
table(tmp$Signature, tmp$Genome)
tmp %>% filter(Signature == 'B3')

for(sig in unique(tmp$Signature)){
  print(sig)
  print(wilcox.test(x = tmp$Scores[tmp$Signature == sig & tmp$Genome == 'ssRNA'], 
                    tmp$Scores[tmp$Signature == sig & tmp$Genome == 'dsDNA'], alternative = 'greater')$p.value)
}


## statistical test: would we expect to see a median shift this large in a sampling of that size?
set.seed(0608)
pval_vect <- rep(0, length(unique(tmp$Signature)))
names(pval_vect) <- unique(tmp$Signature)
for(sig in unique(tmp$Signature)){
  obs_val <- tmp %>%
    filter(Signature == sig, Genome == 'dsDNA')
  sampling_pool <- tmp %>%
    filter(Signature == sig)
  perm <- sapply(1:1000, function(ind){return(median(sample(x = sampling_pool$Scores, size = nrow(obs_val), replace = F)))})
  pval <- mean(median(obs_val$Scores) >= perm)
  pval_vect[[sig]] <- pval
}
pval_vect

obs_vals <- rep(0, length(unique(tmp$Signature)))
names(obs_vals) <- unique(tmp$Signature)
for(sig in unique(tmp$Signature)){
  obs_vals[sig] <- tmp %>% filter(Signature == sig, Genome == 'dsDNA') %>% pull(Scores) %>% median()
}
obs_vals
plot(obs_vals, pval_vect)

text_df <- data.frame('Signature' = names(pval_vect), 'P.Value' = paste0('p = ', as.numeric(pval_vect))) %>%
  mutate('Genome' = 'dsDNA')
pc <- ggplot(specificity_plot_dat2 %>% 
               mutate(Genome = ifelse(!Genome %in% c('dsDNA', 'ssRNA'), 'other', Genome)) %>% 
               mutate(Genome = factor(Genome, levels = c('dsDNA', 'ssRNA', 'other'), ordered = T)) %>% 
               group_by(Signature, Genome) %>% 
               mutate(color = ifelse(median(Scores) <= 0.6, 'specific', 'non-specific')) %>% 
               mutate(Signature = factor(Signature, levels = c('B2', 'B3', 'B4', 'B5', 'B6'), ordered = T)), 
             aes(x = Genome, y = Scores, fill = color)) + 
  geom_boxplot() + 
  scale_fill_manual(values = color_vals)  +
  facet_grid(. ~ Signature, space = 'free', scales = 'free') + 
  theme_clean + 
  geom_hline(yintercept = 0.6, linetype = 'dashed', color = 'gray50') + 
  geom_text(data = text_df, aes(label = P.Value, y = 0.05, fill = NULL), position = position_nudge(x = 1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(y = '') + 
  theme(legend.position = 'none') 

pa / pb / pc


pval_vect <- rep(0, length(unique(tmp$Signature)))
names(pval_vect) <- unique(tmp$Signature)
for(sig in unique(tmp$Signature)){
  obs_val <- tmp %>%
    filter(Signature == sig, Genome == 'dsDNA') 
  if(nrow(obs_val) >= 4){
    sampling_pool <- tmp %>%
      filter(Signature == sig, Genome != 'dsDNA')
    pval_vect[[sig]] <- wilcox.test(x = obs_val$Scores, sampling_pool$Scores, alternative = 'less')$p.value
  } else {
    pval_vect[[sig]] <- NA
  }
}
pval_vect

text_df <- data.frame('Signature' = names(pval_vect), 'P.Value' = paste0('p = ', round(as.numeric(pval_vect), 3))) %>%
  mutate('Genome' = 'dsDNA')
pc <- ggplot(specificity_plot_dat2 %>% 
               mutate(Genome = ifelse(!Genome %in% c('dsDNA', 'ssRNA'), 'other', Genome)) %>% 
               mutate(Genome = factor(Genome, levels = c('dsDNA', 'ssRNA', 'other'), ordered = T)) %>% 
               group_by(Signature, Genome) %>% 
               mutate(color = ifelse(median(Scores) <= 0.6, 'specific', 'non-specific')) %>% 
               mutate(Signature = factor(Signature, levels = c('B2', 'B3', 'B4', 'B5', 'B6'), ordered = T)), 
             aes(x = Genome, y = Scores, fill = color)) + 
  geom_boxplot() + 
  scale_fill_manual(values = color_vals)  +
  facet_grid(. ~ Signature, space = 'free', scales = 'free') + 
  theme_clean + 
  geom_hline(yintercept = 0.6, linetype = 'dashed', color = 'gray50') + 
  geom_text(data = text_df, aes(label = P.Value, y = 0.05, fill = NULL), position = position_nudge(x = 1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(y = '') + 
  theme(legend.position = 'none') 

pa / pb 
pa / pc
