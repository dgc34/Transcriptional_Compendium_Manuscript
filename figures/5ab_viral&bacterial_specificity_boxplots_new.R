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

# order_dat <- input_tmp_plot_dat %>%
#   group_by(Study) %>%
#   summarize('N_star' = sum(P.Value < 0.05)) %>%
#   arrange(N_star)
# 
# signature_labels <- c(paste0('VxB', 1:10), paste0('V', 1:12), paste0('B', 1:10))
# 
# tmp_plot_dat <- input_tmp_plot_dat %>%
#   mutate('Signature' = factor(Signature, levels = signature_labels, ordered = T)) %>%
#   mutate('Study' = factor(Study, levels = order_dat$Study), ordered = T)
# 
# study_plot_dat <- tmp_plot_dat %>%
#   mutate(Train_Label = ifelse(Type == 'Discovery', '○', '')) %>%
#   mutate(Significance_Label = ifelse(P.Value <= 0.05, '*', '')) 
# 
# study_order <- study_plot_dat %>%
#   filter(Comparison == 'V') %>%
#   group_by(Study) %>%
#   summarize('AUC' = mean(Scores, na.rm = T)) %>%
#   ungroup() %>%
#   arrange(desc(AUC)) %>%
#   mutate(Study = factor(Study, levels = Study, ordered = T)) %>% 
#   pull(Study)
# 
# study_plot_dat <- study_plot_dat %>%
#   mutate(Study = factor(Study, levels = study_order, ordered = T))
# 
# sig_order <- study_plot_dat %>%
#   group_by(Signature) %>%
#   summarize(Scores = median(Scores)) %>%
#   arrange(desc(Scores)) %>%
#   pull(Signature)
sig_order <- rev(c(paste0('V', 1:11), paste0('B', 1:7),  paste0('V/B', 1:7)))
study_plot_dat <- input_tmp_plot_dat %>%
  mutate(Signature2 = factor(Signature, levels = sig_order, ordered = T)) %>%
  filter(!grepl(pattern = 'I', x = Signature)) 


# fraction_labs <- study_plot_dat %>%
#   group_by(Signature) %>%
#   summarize(Fraction = mean(P.Value <= 0.05))
# max_locs <- study_plot_dat %>%
#   group_by(Signature) %>%
#   summarize(Scores = max(Scores))
# fraction_labs <- fraction_labs %>%
#   left_join(max_locs) %>%
#   mutate(Scores = ifelse(Scores > 0.95, Scores-0.1, Scores))

color_vals <- c("#79d2f0", "#ffa31f", "#1fe58c", "gray50")
names(color_vals) <- c('V', 'B', 'V/B', 'non-specific') # swap bacterial and viral for this case 
specificity_threshold = 0.6
specificity_plot_dat <- study_plot_dat %>%
  filter(!Signature %in% c('B1', 'V3')) %>% 
  filter(!grepl(pattern = 'V/B', x = Signature)) %>%
  group_by(Signature) %>%
  mutate('med' = median(Scores)) %>% 
  mutate(Color = as.character(Comparison)) %>%
  mutate(Color = ifelse(med <= specificity_threshold, Color, 'non-specific')) %>% 
  mutate(Signature = factor(Signature, levels = sig_order, ordered = T)) %>%
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


p3_option2 <- ggplot(specificity_plot_dat, aes(x = Scores, y = Signature, fill = Color, color = Color)) + 
  geom_density_ridges(quantile_lines = F, alpha = 0.9) + 
  #geom_vline(xintercept = specificity_threshold, linetype = 'dashed', color = 'blue') + 
  #geom_vline(xintercept = c(0.5), linetype = 'dashed', color = 'black') + 
  scale_fill_manual(values = color_vals) + 
  scale_color_manual(values = color_vals) + 
  labs(x = 'AUC') + 
  facet_grid(Comparison ~ ., scales = 'free', space = 'free', switch = 'y') +  
  #xlim(c(0,1)) + 
  theme_clean + 
  theme(legend.position = 'none', strip.background = element_rect(fill = 'white')) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2), limits = c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p3_option1
p3_option2

# #color_vals2 <- c("#F2BE84", "#C4E8F0", "#8FD0BA", "white")
# #names(color_vals2) <- c('B', 'V', 'V/B', 'specific') 
# g <- ggplot_gtable(ggplot_build(p3_option2))
# strip_both <- which(grepl('strip-', g$layout$name))
# fills <- c("#C4E8F0", "#F2BE84")
# k <- 1
# for (i in strip_both) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid.draw(g)
# p5ab <- g

#ggsave(p3, filename = paste0('~/Dropbox/compendium_manuscript/figures/s5_specificity_boxplots.png'), height = 6, width = 4.5)

study_plot_dat %>% 
  group_by(Signature) %>%
  summarize('median' = median(Scores)) %>%
  pull(median) %>%
  median()
### bacteria
source_node <- 'Bacteria'
plot_dat <- bind_rows(out)
if(FLIP_V10_V11){
  table(plot_dat$Signature)
  plot_dat$Signature[plot_dat$Signature == 'V11'] <- 'V12'
  plot_dat$Signature[plot_dat$Signature == 'V10'] <- 'V11'
  plot_dat$Signature[plot_dat$Signature == 'V12'] <- 'V10'
  table(plot_dat$Signature)
}


plot_dat <- plot_dat %>%
  filter(!grepl(pattern = 'NA|DARPA|Chawla', x = Signature)) %>% # remove the single gene testing variants and DARPA signature
  mutate(Comparison = str_extract(pattern = 'VxB|V|B', string = Signature)) 
table(plot_dat$Comparison)


N_threshold = 5
colors <- createColorScale()

# facet by pathogen and cluster rows/columns
gram_stains <- read.xlsx('~/Dropbox/bact_median_aucs_annotated.xlsx') %>%
  mutate(Include = as.logical(Include)) %>%
  dplyr::select(Gram.stain, Study, Include, Age.Group, Tissue) %>%
  mutate(Gram.stain = gsub(pattern = 'acid-fast', replacement = 'Acid-fast', x = Gram.stain)) %>%
  mutate(Gram.stain = gsub(pattern = 'Pos', replacement = 'Gram positive', x = Gram.stain)) %>%
  mutate(Gram.stain = gsub(pattern = 'Neg', replacement = 'Gram negative', x = Gram.stain)) %>%
  mutate(Gram.stain = gsub(pattern = 'Mixed', replacement = 'Multiple gram types', x = Gram.stain)) %>%
  mutate(Gram.stain = gsub(pattern = 'not reported', replacement = 'Unreported', x = Gram.stain)) %>%
  mutate(Gram.stain = factor(Gram.stain, levels = c('Acid-fast', 'Gram positive', 'Multiple gram types', 'Gram negative', 'Unreported'), ordered = T))
# gse77528 is parasitic
# gse73464 is non-infectious
accessions_to_remove <- c('GSE77528', 'GSE73464')
plot_dat <- plot_dat %>%
  filter(!grepl(pattern = paste(accessions_to_remove, collapse = '|'), x = Study))

input_tmp_plot_dat <- plot_dat %>%
  filter(N >= N_threshold, Type != 'Null', Comparison == 'V') %>%
  mutate(Signature = as.character(Signature))

# order_dat <- input_tmp_plot_dat %>%
#   group_by(Study) %>%
#   summarize('N_star' = sum(P.Value < 0.05)) %>%
#   arrange(N_star)

signature_labels <- unique(input_tmp_plot_dat$Signature)

# tmp_plot_dat <- input_tmp_plot_dat %>%
#   mutate('Signature' = factor(Signature, levels = signature_labels, ordered = T)) %>%
#   mutate('Study' = factor(Study, levels = order_dat$Study), ordered = T)
# 
study_plot_dat <- input_tmp_plot_dat %>%
   mutate(Train_Label = ifelse(Type == 'Discovery', '○', '')) %>%
#   #mutate(Significance_Label = ifelse(P.Value <= 0.05, '*', '')) %>%
   left_join(gram_stains) %>%
   mutate('Signature' = factor(Signature, levels = signature_labels, ordered = T)) #%>%
   #mutate('Study' = factor(Study, levels = order_dat$Study), ordered = T)
# 
# getSignificantDf <- function(df2){
#   df2 %>% 
#     ungroup() %>%
#     filter(N >= 4, !Signature %in% sigs_to_omit, !is.na(P.Value)) %>%
#     group_by(Signature, Gram.stain) %>%
#     summarize("N_sig" = sum(P.Value <= 0.05, na.rm = T), "N" = sum(!is.na(P.Value)), 'median' = median(Scores, na.rm = T)) %>%
#     mutate('Fraction' = N_sig/N) %>%
#     ungroup() %>%
#     #mutate('Specific' = ifelse(median <= 0.6, 1, 0)) %>%
#     mutate('Specific' = ifelse(Fraction < 0.5, 1, 0)) %>%
#     mutate(Specific = as.logical(Specific)) %>%
#     mutate(Specific.Logical = as.logical(Specific)) %>%
#     mutate(Specific = ifelse(Specific, 'Specific', 'Non-Specific')) %>%
#     mutate(Specific = factor(Specific, levels = c('Specific', 'Non-Specific'), ordered = T)) %>%
#     mutate('Slabel' = ifelse(Specific == 'Specific', '', '.'))  %>%
#     # mutate(Specific = ifelse(Fraction == 0, 1, 0)) %>%
#     # mutate(Specific.Logical = as.logical(Specific)) %>%
#     # mutate(Specific = ifelse(Specific, 'Specific', 'Non-Specific')) %>%
#     # mutate(Specific = factor(Specific, levels = c('Specific', 'Non-Specific'), ordered = T)) %>%
#     return()
# }
# sigs_to_omit <- ''
# tmp <- study_plot_dat %>%
#   filterGeneFraction(0.5) %>%
#   getSignificantDf()
# ggplot(tmp %>% mutate(Signature = factor(Signature, levels = paste0('V', 1:11), ordered = T)), 
#        aes(x = Signature, y = Gram.stain, fill = Specific)) + 
#   #geom_tile(alpha = 0.1, fill = 'white', color = 'black') +
#   geom_tile(alpha = 0.1) +
#   #geom_tile() +
#   geom_text(aes(label = Slabel), size = 15, position = position_nudge(y = 0.5)) + 
#   geom_vline(xintercept = (0:length(unique(tmp$Signature)))+0.5, color = 'gray95') + 
#   geom_hline(yintercept = (0:length(unique(tmp$Gram.stain)))+0.5) + 
#   #scale_fill_manual(values = color_bar) + 
#   #scale_fill_gradient2(low = scales::muted('blue'), high = scales::muted('red'), mid = 'white', midpoint = 0.5) + 
#   theme_minimal() + 
#   theme(panel.grid = element_blank()) + 
#   #facet_grid(. ~ Comparison, scales = 'free', space = 'free') + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# ggsave(filename = '~/Dropbox/compendium_manuscript/figures/sX_viralSpecificity_fractional.png', height = 2.5, width = 5, dpi = 300)

gram_plot_dat <- study_plot_dat %>%
  group_by(Gram.stain, Signature) %>%
  #summarize(Fraction = mean(P.Value <= 0.05), N = n(), AUC = median(Scores)) 
  summarize(N = n(), AUC = median(Scores))
# sig_order <- gram_plot_dat %>%
#   ungroup() %>%
#   group_by(Signature) %>%
#   summarize(AUC = median(AUC)) %>%
#   arrange(desc(AUC)) %>%
#   pull(Signature)

head(gram_plot_dat)
# radar_plot_dat <- gram_plot_dat %>%
#   filter(!Signature %in% c('V3','B1')) %>% 
#   mutate(Signature = factor(Signature, levels = paste0('V', 11:1), ordered = T)) %>% 
#   mutate(Gram.stain = as.character(Gram.stain)) %>% 
#   filter(Gram.stain %in% c('Acid-fast', 'Gram positive', 'Gram negative')) %>% 
#   dplyr::select(Signature, Gram.stain, AUC) %>% 
#   tidyr::pivot_wider(names_from = Signature, values_from = AUC)  #%>%
#   #dplyr::select(c('Gram.stain', paste0('V', c(1:2,4:11))))
# p5c <- ggradar(radar_plot_dat %>% dplyr::select(-V2), 
#         background.circle.colour = 'white', 
#         grid.mid = 0.6, label.gridline.mid = F,
#         grid.min = 0, label.gridline.min = F,
#         grid.max = 1, label.gridline.max = F,
#         group.colours = c('#c47115', '#F2BE84', '#5a340a'),
#         gridline.min.colour = 'black',
#         gridline.max.colour = 'black',
#         axis.line.colour = 'black',
#         gridline.mid.colour = 'blue') + 
#   theme(legend.position = 'bottom')
shape_vals <- c(22, 21, 24)
names(shape_vals) <- c('Acid-fast', 'Gram positive', 'Gram negative')
color_vals <- c('#79d2f0', 'gray50')
names(color_vals) <- c('NXR', 'XR')
p5c <- ggplot(gram_plot_dat %>%
                filter(!Signature %in% c('V3')) %>% 
                mutate(Signature = factor(Signature, levels = paste0('V', 11:1), ordered = T)) %>% 
                filter(grepl(pattern = 'positive|negative|fast', x = Gram.stain)) %>%
                mutate(XR = ifelse(AUC <= specificity_threshold, 'NXR', 'XR')), 
              aes(x = Signature, y = AUC, shape = Gram.stain, group = Gram.stain, fill = XR, color = XR)) + 
  geom_line(size = 2, color = 'gray80') + 
  geom_point(size = 4) +
  scale_y_continuous(limits = c(0,1), breaks = seq(to = 1, from = 0, by = 0.2)) + 
  theme_clean + 
  scale_shape_manual(values = shape_vals, name = 'bacterial class') + 
  scale_fill_manual(values = color_vals, name = '') + 
  scale_color_manual(values = color_vals, name = '') + 
  theme(legend.position = 'top') + 
  coord_flip() + 
  guides(fill = 'none')
p5c
gram_plot_dat %>% filter(grepl(pattern = 'Acid', x = Gram.stain)) %>% filter(AUC > 0.6) %>% pull(AUC) %>% summary()
p5c2 <- p5c
p5c2 <- p5c2 + theme(legend.text = element_text(size = 10), legend.position = 'bottom')
#grid.draw(p5ab) 
p5c_legend <- get_legend(p5c2)
p5c <- p5c + theme(legend.position = 'none')
p5c

