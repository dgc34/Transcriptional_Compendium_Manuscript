library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(ggpubr)

source('~/Dropbox/plot_palette.R')

output_dir <- '~/Dropbox/figS1/'
corr_files <- list.files(output_dir, pattern = 'RDS', full.names = T)
corr_list <- lapply(corr_files, readRDS) 
null_ix <- sapply(corr_list, function(x){is.null(unlist(x))})
corr_list <- corr_list[!null_ix]
corr_list <- corr_list %>%
  lapply(FUN = function(x){x$Signature <- as.character(x$Signature); return(x)})

geom_mean_scores <- readRDS('~/Dropbox/big_eval_heatmap_list_generalizability_100621sdy_ts.RDS') %>%
  bind_rows() %>%
  #dplyr::select(Study, Signature, Type, Scores) %>%
  mutate('Method' = 'zScore_geomMean')

plot_dat <- bind_rows(corr_list) %>%
  bind_rows(geom_mean_scores) %>%
  filter(N >= 15) %>%
  dplyr::select(Study, Signature, Comparison, Method, Scores, Type) %>%
  pivot_wider(names_from = 'Method', values_from = 'Scores') %>%
  mutate(Type = factor(Type)) %>%
  mutate(Null = Type == 'Null') %>% 
  filter(!is.na(logit_), !is.na(zScore_geomMean)) %>%
  mutate(Comparison = str_extract(pattern = 'VxB|V|B', string = Signature))

label_dat <- plot_dat %>%
  group_by(Comparison, Null) %>%
  summarize('corr' = round(cor.test(x = zScore_geomMean, y = logit_)$estimate,3),
            'pval' = round(cor.test(x = zScore_geomMean, y = logit_)$p.value,3))


comparisons <- c("VxB", "V", "B")
p_list <- list()
for (i in seq_along(comparisons)){
  comparison <- comparisons[i]
  
  p_list[[i]] <- ggplot(plot_dat %>% 
                          filter(Comparison %in% comparison) %>%
                          pivot_longer(cols = matches(c('zScore_geomMean', 'logit_')), names_to = 'Method', values_to = 'AUC')
                          , aes(x = Signature, y = AUC, color = Method)) + 
    geom_boxplot() +
    stat_compare_means(method = 't.test', paired = T, label = "p.signif") + 
    # geom_line(data = data.frame(x = c(0,1), 
    #                             y = c(0,1)), 
    #           aes(x = x, 
    #               y = y, 
    #               color = 'y=x')) + 
    # labs(title = comparison, y = 'logistic regression AUC', x = 'geometric mean zscore AUC') + 
    # geom_text(data = label_dat %>% 
    #             filter(Comparison %in% comparison), 
    #           aes(x = 0.5, y = 0.2, color = NULL,
    #               label = paste0('r = ', corr,'\np.value = ', pval))) + 
    #facet_grid(Null ~ Signature, scales = 'free', space = 'free') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylim(c(0,1))
  
  #if(i != length(comparisons)){
  if(i != 1){
    p_list[[i]] <- p_list[[i]] + 
      theme(legend.position = 'none')
  }
}


layout <- "
AAA###
BBBBBB
CCC###
"
p4 <- p_list[[1]] + p_list[[2]] + p_list[[3]] + plot_layout(design = layout)
ggsave(p4, filename = '~/Dropbox/fig_s4_061321.png', width = 18, height = 12, units = 'in')

plot_dat <- plot_dat %>% mutate(Comparison = gsub(pattern = 'x', replacement = '/', x = Comparison))
label_dat <- label_dat %>% mutate(Comparison = gsub(pattern = 'x', replacement = '/', x = Comparison))
comparisons <- c("V/B", "V", "B")
p_list <- list()
for (i in seq_along(comparisons)){
  comparison <- comparisons[i]
  
  p_list[[i]] <- ggplot(plot_dat %>% 
                          filter(Comparison %in% comparison
                          ), aes(x = zScore_geomMean, y = logit_)) + 
    geom_point() +
    geom_line(data = data.frame(x = c(0,1), 
                                y = c(0,1)), 
              aes(x = x, 
                  y = y),
              color = 'black', linetype = 'dashed') + 
    labs(title = comparison, y = 'logistic regression AUROC', x = 'geometric mean AUROC') + 
    geom_text(data = label_dat %>% 
                filter(Comparison %in% comparison), 
              aes(x = 0.5, y = 0.2, color = NULL,
                  label = paste0('r = ', corr,'\np.value = ', pval))) + 
    #facet_grid(Null ~ Signature, scales = 'free', space = 'free') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlim(c(0,1)) + ylim(c(0,1)) + 
    theme_clean
  
  #if(i != length(comparisons)){
  if(i != 1){
    p_list[[i]] <- p_list[[i]] + 
      theme(legend.position = 'none')
  }
}


layout <- "BCA"
p4 <- p_list[[1]] + p_list[[2]] + p_list[[3]] + plot_layout(design = layout)
pdf('~/Dropbox/compendium_manuscript/figures/supplements/supp1_ABC.pdf', width = 8.5, height = 3.5)
p4
dev.off()
ggsave(p4, filename = '~/Dropbox/compendium_manuscript/figures/s1_logit_vs_geomMean.png', width = 18, height = 12, units = 'in')

label_dat %>% 
  filter(!Null) %>% 
  group_by(Comparison) %>% 
  summarize(corr = median(corr), pval = median(pval))
plot_dat %>%
  filter(!Null) %>%
  group_by(Signature, Comparison) %>%
  summarize(N = n()) %>%
  arrange(desc(Comparison))

corr_values <- plot_dat %>%
  #filter(grepl(pattern = 'Sweeney', x = Signature)) %>%
  mutate(Null = ifelse(Type == 'Null', 'Permuted', 'Real Labels')) %>%
  group_by(Signature, Null) %>%
  filter(!(is.na(zScore_geomMean) | is.na(logit_))) %>%
  mutate('N' = n()) %>%
  group_by(Comparison, Signature, Null, N) %>%
  summarize('corr' = cor(x = zScore_geomMean, y = logit_), 'pval' = cor.test(x = zScore_geomMean, y = logit_)$p.value)

corr_values
colors <- createColorScale()
ggplot(corr_values %>% filter(Comparison %in% c("Virus vs. Bacteria", "Virus vs. Healthy", "Bacteria vs. Healthy")), aes(y = corr, x = Null)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point() + 
  #scale_color_manual(values = colors) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(. ~ Comparison) +
  ylim(-1,1)

corr_values %>% filter(Null == 'Real Labels') %>% group_by(Comparison) %>% summarize('median' = median(corr))



corr_values_filtered <- plot_dat %>%
  #filter(grepl(pattern = 'Sweeney', x = Signature)) %>%
  #filter(N >= 15) %>% # this filtering was already done
  mutate(Null = ifelse(Type == 'Null', 'Permuted', 'Real Labels')) %>%
  filter(Null == 'Real Labels') %>%
  group_by(Signature, Null) %>%
  filter(!(is.na(zScore_geomMean) | is.na(logit_))) %>%
  mutate('N' = n()) %>%
  ungroup() %>%
  group_by(Comparison, Signature, N) %>%
  summarize('corr' = cor(x = zScore_geomMean, y = logit_), 'pval' = cor.test(x = zScore_geomMean, y = logit_)$p.value)

for (comp in unique(corr_values_filtered$Comparison)){
  corr_values_filtered %>% filter(Comparison == comp) %>% print()
}


