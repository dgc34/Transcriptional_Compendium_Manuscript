library(dplyr)
library(ggplot2)
library(stringr)
library(ggpattern)

source('~/Dropbox/plot_palette.R')
getSignificantDf <- function(df2){
  df2 %>% 
    ungroup() %>%
    filter(N >= 4, !Signature %in% sigs_to_omit, !is.na(P.Value)) %>%
    group_by(Signature, Pathogen) %>%
    summarize("N_sig" = sum(P.Value <= 0.05, na.rm = T), "N" = sum(!is.na(P.Value)), 'median' = median(Scores, na.rm = T)) %>%
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

df2 <- readRDS('~/Documents/darpa-manuscript-data/parasite_eval_heatmap.RDS')
sigs_to_omit <- c('V3', 'B1')

sig_dat <- getSignificantDf(df2)

specificity_threshold <- 0.6
med_plot_sig <- sig_dat %>%
  #mutate(Specific = ifelse(median > specificity_threshold, 'Non-Specific', 'Specific')) %>% 
  dplyr::select(Signature, Specific) %>%
  filter(Specific == 'Non-Specific') %>%
  distinct()

leishmania_sdy <- 'GSE77528_GPL10558'
df2 %>% filter(Study == 'GSE77528_GPL10558') %>%
  group_by(Sig.Type) %>% summarize('med' = median(Scores)) 

VxB_order <- df2 %>%
  filter(Sig.Type == 'VxB') %>%
  mutate(Signature = gsub(pattern = 'x', replacement = '/', x = Signature)) %>% 
  group_by(Signature) %>%
  summarize('med' = median(Scores)) %>%
  arrange(desc(med)) %>%
  pull(Signature) 

#sig_order_para <- c(as.character(sig_order), as.character(VxB_order))
sig_order_para <- sig_order <- c(paste0('V', 1:11), paste0('V/B', 1:7), paste0('B', 1:7))
plot_dat_para <- df2 %>%
  dplyr::filter(Sig.Type != 'Influenza') %>%
  dplyr::filter(!Signature %in% sigs_to_omit) %>%
  mutate(Sig.Type = gsub(Sig.Type, pattern = 'x', replacement = '/')) %>%
  mutate(Sig.Type = gsub(Sig.Type, pattern = 'Virus', replacement = 'V')) %>%
  mutate(Sig.Type = gsub(Sig.Type, pattern = 'Bacteria', replacement = 'B')) %>%
  mutate(Sig.Type = factor(Sig.Type, levels = c('V', 'B', 'V/B'))) %>% 
  group_by(Signature) %>%
  mutate(med = median(Scores)) %>% 
  mutate(color = ifelse(med > specificity_threshold, 'not specific', as.character(Sig.Type))) %>%
  ungroup() %>%
  mutate(Signature = gsub(pattern = 'x', replacement = '/', x = Signature)) %>% 
  mutate(Signature = factor(Signature, levels = sig_order_para, ordered = T))
outlier_df <- plot_dat_para %>% filter(Study == 'GSE77528_GPL10558')


color_vals <- c("#79d2f0", "#ffa31f", "#1fe58c", "gray50")
names(color_vals) <- c('V', 'B', 'V/B', 'not specific')
color_vals_para <- color_vals
#color_vals_para[c('V', 'B', 'V/B')] <- c(rep('gold1', 3))
f5b <- ggplot(plot_dat_para %>% mutate('Exposure' = 'Parasite'), 
       aes(x = Signature, y = Scores, fill = color)) + 
  geom_boxplot(outlier.shape = NA) +
  #geom_boxplot() +
  #geom_density_ridges() + 
  #geom_point(data = outlier_df, color = 'black', pch = 4) + 
  theme_clean + 
  #geom_point(data = df %>% 
  #             filter(grepl(pattern = 'GSE1124', x = Study), Sig.Type == 'VxB'), 
  #           aes(color = 'red')) + 
  #geom_line(aes(group = Study), alpha = 0.3) + 
   facet_grid(Exposure ~ Sig.Type, space = 'free', scales = 'free', switch = 'y') + 
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = color_vals_para) + 
   #geom_hline(yintercept = c(0.6), linetype = 'dashed', color = 'blue') +
   #geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'black') +
   #labs(x = 'AUROC') + 
  theme(strip.background = element_rect(fill = 'white')) + 
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2), limits = c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = 'none')
f5b  

p3_option2 + f5b + plot_layout(design = c('AB
                                           AB
                                           AB
                                           AB
                                           AB
                                           AB
                                           AB
                                           #B'))
# pathogens <- lapply(data_list, function(x){return(data.frame('Study' = x$formattedName, 
#                                                              'Pathogen' = setdiff(unique(x$pheno$Pathogen), c('Healthy', 'Unknown'))
# ))}) %>%
#   bind_rows() %>%
#   group_by(Study) %>%
#   summarize(Pathogen = paste0(Pathogen, collapse = ';')) %>%
#   mutate(Pathogen = ifelse(grepl(pattern = 'Plasmodium', x = Pathogen), 'Plasmodium', Pathogen))
# 
# df2 <- df %>%
#   left_join(pathogens)

color_vals <- c("#79d2f0", "#ffa31f", "#1fe58c", "gray50")
names(color_vals) <- c('V', 'B', 'V/B', 'not specific')
ggplot(data.frame('type' = c('B', 'V', 'V/B'), 'value' = c(1, 1, 1)), aes(x = type, y = value, color = type)) + 
  geom_point() + 
  scale_color_manual(values = color_vals)

#sig_order <- getOrder(sig_dat, Signature)
#sig_order <- rev(c(paste0('VxB', 1:12), paste0('V', 1:12), paste0('B', 1:12)))
path_order <- getOrder(sig_dat, Pathogen)

sig_dat_plot <- sig_dat %>%
  mutate('Signature' = factor(Signature, levels = rev(sig_order), ordered = T)) %>%
  mutate('Pathogen' = factor(Pathogen, levels = path_order, ordered = T)) %>%
  mutate(Comparison = str_extract(pattern = 'I|VxB|V|B', string = Signature)) %>%
  mutate(Comparison = factor(Comparison, levels = c('VxB', 'V', 'B'), ordered = T))

f5c <- ggplot(sig_dat_plot %>%
                filter(Comparison %in% c('VxB', 'V', 'B')) %>%
                mutate(Comparison = as.character(Comparison)) %>%
                mutate(Comparison = gsub(pattern = 'VxB', replacement = 'V/B', x = Comparison)) %>%
                mutate(Comparison = factor(Comparison, levels = c('V', 'B', 'V/B'), ordered = T)), 
              aes(x = Signature, y = Pathogen)) + 
  #geom_tile(alpha = 0.1, fill = 'white', color = 'black') +
  geom_tile(alpha = 0.1, fill = 'white') +
  #geom_tile() +
  geom_text(aes(label = Slabel), size = 15, position = position_nudge(y = 0.65)) + 
  geom_vline(xintercept = (0:length(unique(sig_dat$Signature)))+0.5, color = 'gray95') + 
  geom_hline(yintercept = (0:length(unique(sig_dat$Pathogen)))+0.5) + 
  #scale_fill_manual(values = color_bar) + 
  #scale_fill_gradient2(low = scales::muted('blue'), high = scales::muted('red'), mid = 'white', midpoint = 0.5) + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) + 
  facet_grid(. ~ Comparison, scales = 'free', space = 'free') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = 'Parasite')
f5c

df2 %>% filter(grepl(pattern = 'GSE5418', x = Study), Sig.Type == 'VxB')

out <- readRDS('~/Documents/darpa-manuscript-data/figures/aging_obesity_heatmap_012522.RDS') 

plot_dat <- bind_rows(out) %>%
  filter(!grepl(pattern = '^I', x = Signature))
plot_dat[grepl(pattern = 'GSE65219', x = plot_dat$Study) & grepl(pattern = 'Sampson', x = plot_dat$Signature), 'Type'] <- 'Discovery'

plot_dat <- plot_dat %>%
  mutate('Type' = factor(Type, levels = c('Discovery', 'Validation', 'New', 'Null'), ordered = T)) %>%
  mutate('Comparison' = str_extract(pattern = 'VxB|V|B', string = Signature))  %>%
  arrange(Comparison) %>%
  mutate(Signature = factor(Signature, levels = unique(Signature), ordered = T))

table(plot_dat$Comparison)

aging_list <- readRDS('~/Documents/darpa-manuscript-data/data_list_final/aging_list.RDS')
obesity_list <- readRDS('~/Documents/darpa-manuscript-data/data_list_final/obesity_list.RDS')

exposure_df <- lapply(c(aging_list, obesity_list), function(x){
  data.frame('Study' = x$formattedName, 'Exposure' = setdiff(unique(x$pheno$Standardized.Exposure), 'Healthy')) %>%
    return()
}) %>%
  bind_rows()

#write.xlsx('~/Dropbox/compendium_manuscript/tables/supp_table5.xlsx', x = exposure_df)

plot_dat <- plot_dat %>%
  left_join(exposure_df)

N_threshold = 0
#colors <- createColorScale()

# facet by pathogen and cluster rows/columns
plot_list <- list()
head(plot_dat)
table(plot_dat$Exposure.Class)
exposures <- unique(plot_dat$Comparison)

pattern_vals <- c('none', 'stripe')
names(pattern_vals) <- c('Aging', 'Obesity')
f5efg_list <- list()
plot_dat <- plot_dat %>% 
  mutate(Comparison = gsub(pattern = 'VxB', replacement = 'V/B', x = Comparison)) %>%
  mutate(Comparison = factor(Comparison, levels = c('V', 'B', 'V/B'), ordered = T))


#color_vals_noninf <- color_vals_para[2:4]
#names(color_vals_noninf) <- c('Aging', 'Obesity', 'not specific')
f5efg <- ggplot(plot_dat %>%
                  filter(!Signature %in% sigs_to_omit) %>%
                  filter(!is.na(Comparison)) %>% 
                  #filter(Exposure == 'Aging') %>% 
                  mutate(Exposure = factor(Exposure, levels = c('Obesity', 'Aging'), ordered = T)) %>% 
                  group_by(Exposure, Signature) %>%
                  mutate('med' = median(Scores)) %>% 
                  mutate('specific' = ifelse(med <= specificity_threshold, as.character(Exposure), 'not specific')),
                aes(x = Signature, y = Scores, fill = Comparison)) + 
  geom_boxplot() + 
  #scale_pattern_manual(values = pattern_vals) + 
  facet_grid(Exposure ~ Comparison, scales = 'free', space = 'free') +
  #facet_grid(. ~ Comparison, scales = 'free', space = 'free') +
  #ylim(c(0,1)) + 
  scale_fill_manual(values = color_vals) +
  geom_hline(yintercept = c(0.6), color = 'blue', linetype = 'dashed') +
  geom_hline(yintercept = c(0.5), color = 'black', linetype = 'dashed') + 
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = '', y = '') + 
  theme(legend.position = 'none', strip.background = element_rect(fill = 'white'))
f5efg

sig_order <- c(paste0('V', 1:11), paste0('V/B', 1:7), paste0('B', 1:7))
f5efg_list <- list()
limits_list <- list('Aging' = c(0,1), 'Obesity' = c(0, 1))
for (exp in c('Obesity', 'Aging')){
  print(exp)
  f5efg_list[[exp]] <- ggplot(plot_dat %>%
                                 filter(!Signature %in% c('V3', 'B1')) %>%
                                 mutate(Signature = gsub(Signature, pattern = 'x', replacement = '/')) %>% 
                                 mutate(Signature = factor(Signature, levels = sig_order, order = T)) %>% 
                                 group_by(Exposure, Signature) %>%
                                 mutate('med' = median(Scores)) %>% 
                                 mutate('specific' = ifelse(med <= specificity_threshold, as.character(Comparison), 'not specific')) %>% 
                                 filter(Exposure == exp), 
         aes(x = Signature, y = Scores, fill = specific)) + 
    #geom_boxplot_pattern(color = 'black', pattern_fill = 'white', pattern_angle = 45, pattern_density = 1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) + 
    geom_boxplot() + 
    #scale_pattern_manual(values = pattern_vals) + 
    #facet_grid(Exposure ~ Comparison, scales = 'free', space = 'free') +
    facet_grid(Exposure ~ Comparison, scales = 'free', space = 'free', switch="y") +
    #ylim(c(0,1)) + 
    #geom_hline(yintercept = c(0.6), color = 'blue', linetype = 'dashed') +
    #geom_hline(yintercept = c(0.5), color = 'black', linetype = 'dashed') + 
    theme_clean + 
    theme(strip.background = element_rect(fill = 'white')) + 
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2), limits = limits_list[[exp]]) + 
    scale_fill_manual(values = color_vals) +
    labs(x = '', y = '') + 
    theme(legend.position = 'none')
}
pout <- f5efg_list[[1]] / f5efg_list[[2]]
#pout
#ggsave(pout, filename = '~/Dropbox/5efg.png', width = 8.5, height = 3, dpi = 300)

f5_obesity <- f5efg_list[[1]]
f5_aging <-   f5efg_list[[2]]

radar_plot <- p5c2

design = c(
                                          'AAB
                                           AAB
                                           AAB
                                           AAB
                                           AAB
                                           AAB
                                           AAB
                                           CCB
                                           CCB
                                           CC#
                                           DDD
                                           DDD
                                           EEE
                                           EEE
                                           EEE
                                           EEE')

design <- c('AABB
             AABB
             AABB
             AABB
             AABB
             AABB
             AABB
             CCCC
             CCCC
             DDDD
             DDDD
             EEEE
             EEEE')
f5b2 <- f5b + theme(axis.text.x = element_blank()) + labs(x = '')
f5_obesity2 <- f5_obesity + theme(axis.text.x = element_blank(), strip.background.x = element_blank(), strip.text.x = element_blank()) + labs(x = '')
f5_aging2 <- f5_aging + theme(strip.background.x = element_blank(), strip.text.x = element_blank())
pout <- p3_option2 + plot_spacer() + inset_element(wrap_elements(radar_plot), align_to = 'panel', left = -0.2, 
                                 bottom = -0.7, 
                                 right = 1.2, 
                                 top = 1.2, on_top = F) + 
  f5b2 + f5_obesity2 + f5_aging2 + plot_layout(design = design) + plot_annotation(tag_level = 'A')

pout <- p3_option2 + plot_spacer() + inset_element(wrap_elements(plot_spacer()), align_to = 'panel', left = -0.2, 
                                                   bottom = -0.7, 
                                                   right = 1.2, 
                                                   top = 1.2, on_top = F) + 
  f5b2 + f5_obesity2 + f5_aging2 + plot_layout(design = design) + plot_annotation(tag_level = 'A')

pout <- p3_option2 + p5c2 + 
  f5b2 + f5_obesity2 + f5_aging2 + plot_layout(design = design) + plot_annotation(tag_level = 'A')

pout
pdf(file = '~/Dropbox/fig5_assembled_sans_radar.pdf', height = 11, width = 8.5); plot(pout); dev.off()
ggsave(pout, filename = '~/Dropbox/compendium_manuscript/figures/5_assembled.png', width = 8.5, height = 11, dpi = 300)
 # 
# plot_dat %>%
#   filter(!Signature %in% c('V3', 'B1', 'B5', 'B7')) %>%
#   group_by(Signature, Comparison, Exposure) %>%
#   summarize('min.p' = min(P.Value)) %>%
#   arrange(Exposure) %>%
#   write.xlsx(file = '~/Dropbox/aging_obesity_significance.xlsx')

aodf <- plot_dat %>%
  filterGeneFraction(frac_threshold = 0.5) %>% 
  dplyr::rename(Pathogen = Exposure) %>% 
  getSignificantDf() %>%
  mutate(Comparison = str_extract(pattern = 'I|VxB|V|B', string = Signature)) %>%
  mutate(Comparison = factor(Comparison, levels = c('VxB', 'V', 'B'), ordered = T))

# pd <- filterGeneFraction(plot_dat, frac_threshold = 0.5) 
# pd %>%
#   filter(Signature == 'V10') %>%
#   dplyr::select(Study, Signature, frac_pos, frac_neg, Exposure, Scores, P.Value)

f5d <- ggplot(aodf %>%
          mutate(Comparison = gsub(pattern = 'VxB', replacement = 'V/B', x = as.character(Comparison))) %>%
            mutate(Comparison = factor(Comparison, levels = c('V', 'B', 'V/B'), ordered = T)) %>%
            filter(Pathogen == 'Obesity'), 
          aes(x = Signature, y = Pathogen)) + 
  #geom_tile(alpha = 0.1, fill = 'white', color = 'black') +
  geom_tile(alpha = 0.1, fill = 'white') +
  #geom_tile() +
  geom_text(aes(label = Slabel), size = 15, position = position_nudge(y = 0.15)) + 
  #geom_vline(xintercept = (0:length(unique(sig_dat$Signature)))+0.5, color = 'gray95') + 
  #geom_hline(yintercept = (0:length(unique(sig_dat$Pathogen)))+0.5) + 
  #scale_fill_manual(values = color_bar) + 
  #scale_fill_gradient2(low = scales::muted('blue'), high = scales::muted('red'), mid = 'white', midpoint = 0.5) + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) + 
  facet_grid(. ~ Comparison, scales = 'free', space = 'free') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
f5d

plot_list <- list(f5c, f5d)
for (i in 1:length(plot_list)){
  current_plot <- plot_list[[i]]
  current_plot <- egg::set_panel_size(current_plot, height=unit((length(unique(current_plot$data$Pathogen))/2), "cm"),
                                      width=unit((length(unique(current_plot$data$Signature))/1.5), "cm") )
  plot_list[[i]] <- current_plot
}
cowplot::plot_grid(f5c, f5d, nrow=2)


combined_plot_dat <- f5c$data %>%
  mutate(Signature = as.character(Signature)) %>%
  bind_rows(f5d$data %>% mutate(Signature = as.character(Signature))) %>%
  mutate('Infectious' = ifelse(Pathogen %in% c('Obesity', 'Aging'), F, T)) %>%
  mutate(Signature = factor(Signature, levels = rev(sig_order), ordered = T))

ggplot(combined_plot_dat %>%
         mutate(Infectious = !Infectious) %>%
         filter(Comparison %in% c('VxB', 'V', 'B')), aes(x = Signature, y = Pathogen, fill = Specific)) + 
  #geom_tile(alpha = 0.1, fill = 'white', color = 'black') +
  geom_tile(alpha = 0.1) +
  #geom_tile() +
  geom_text(aes(label = Slabel), size = 15, position = position_nudge(y = 0.75)) + 
  geom_vline(xintercept = (0:length(unique(sig_dat$Signature)))+0.5, color = 'gray95') + 
  geom_hline(yintercept = (0:length(unique(sig_dat$Pathogen)))+0.5) + 
  #scale_fill_manual(values = color_bar) + 
  #scale_fill_gradient2(low = scales::muted('blue'), high = scales::muted('red'), mid = 'white', midpoint = 0.5) + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) + 
  facet_grid(Infectious ~ Comparison, scales = 'free', space = 'free') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = '')
ggsave(filename = '~/Dropbox/compendium_manuscript/figures/5cd_parasite_aging_obesity_specificity.png', width = 8, height = 2.4, dpi = 300)

f5c_out <- f5c + theme(axis.text.x = element_blank()) + labs(x = '', y = '')
f5d_out <- f5d + theme(axis.text.x = element_blank()) + labs(x = '', y = '')
f5efg_out <- f5efg + labs(title = 'Aging') + theme(strip.background =element_rect(fill="white"))

lay <- rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,3,3),
             c(3,3,3,3),
             c(4,4,4,4),
             c(5,5,5,5),
             c(5,5,5,5),
             c(5,5,5,5),
             c(5,5,5,5))
pout <- grid.arrange(p5ab, p5c, f5c_out, f5d_out, f5efg_out, layout_matrix = lay)

ggsave(plot = pout, filename = '~/Dropbox/5_assembled.png', dpi = 300, height = 11, width = 8.5)

mat <- c(
  'AAA###
   AAA###
   AAABBB
   AAABBB
   AAABBB
   AAABBB
   AAABBB
   AAABBB
   AAA###
   AAA##F
   AAA###
   AAA###
   CCCCCC
   CCCCCC
   CCCCCC
   DDDDDD
   EEEEEE
   EEEEEE
   EEEEEE
   EEEEEE
   EEEEEE
   EEEEEE
   EEEEEE'
  )

mat <- c('AAA
          AAA
          AAA
          BBB
          BBB
          CCC
          DDD
          DDD
          DDD')
 
top_row <- wrap_elements(p5ab) + plot_spacer() + inset_element(wrap_elements(p5c), align_to = 'panel', left = -0.2, 
                                                               bottom = -0.5, 
                                                               right = 1.2, 
                                                               top = 1.5, on_top = F)
#top_row
pout <- top_row + f5c_out + f5d_out + f5efg_out + wrap_elements(p5c_legend) + plot_layout(design = mat)
pout
ggsave(filename = '~/Dropbox/5_assembled.png', plot = pout, height = 11, width = 8.5)
