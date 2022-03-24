library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(ggpubr)

source("~/Dropbox/plot_palette.R")
output_dir <- '~/Dropbox/figS1/'
corr_files <- list.files(output_dir, full.names = T, pattern = 'RDS')
corr_list <- lapply(corr_files, readRDS) 
null_ix <- sapply(corr_list, function(x){is.null(unlist(x))})
corr_list <- corr_list[!null_ix]
corr_list <- corr_list %>%
  lapply(FUN = function(x){x$Signature <- as.character(x$Signature); return(x)})

color_vals <- c('#F2BE84', '#C4E8F0', '#8FD0BA')
names(color_vals) <- c('V', 'B', 'V/B')

geom_mean_scores <- readRDS('~/Dropbox/big_eval_heatmap_list_generalizability_100621sdy_ts.RDS') %>%
  bind_rows() %>%
  #dplyr::select(Study, Signature, Type, Scores) %>%
  mutate('Method' = 'zScore_geomMean')

full_df <- bind_rows(corr_list) %>%
  bind_rows(geom_mean_scores) %>%
  filter(N >= 15) %>%
  dplyr::select(Study, Signature, Comparison, Method, Scores, Type) %>%
  pivot_wider(names_from = 'Method', values_from = 'Scores') %>%
  mutate(Type = factor(Type)) %>%
  mutate(Null = Type == 'Null') %>% 
  filter(!is.na(logit_), !is.na(zScore_geomMean)) %>%
  mutate(Comparison = str_extract(pattern = 'VxB|V|B', string = Signature)) %>%
  mutate(Comparison = gsub(Comparison, pattern = 'x', replacement = '/')) %>%
  mutate(Signature = gsub(Signature, pattern = 'x', replacement = '/')) %>%
  mutate(Signature = factor(Signature, levels = c(paste0('V', 1:11), paste0('B', 1:7), paste0('V/B', 1:6)), ordered = T))


label_dat <- full_df %>%
  group_by(Signature, Comparison, Null) %>%
  summarize('corr' = round(cor.test(x = zScore_geomMean, y = logit_)$estimate,3),
            'pval' = round(cor.test(x = zScore_geomMean, y = logit_)$p.value,3))

plot_list <- list()
for(comparison in unique(full_df$Comparison)){
  plot_list[[comparison]] <- ggplot(full_df %>% filter(Comparison == comparison), 
                                    aes(y = logit_, x = zScore_geomMean, color = Comparison)) + 
    geom_point() + 
    theme_clean + 
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray50') +
    facet_grid(. ~ Signature) + 
    labs(x = 'geometric mean', y = 'logistic regression') + 
    xlim(c(0,1)) + 
    ylim(c(0,1)) + 
    scale_color_manual(values = color_vals) + 
    geom_text(data = label_dat %>% 
                filter(Comparison %in% comparison), 
              aes(x = 0.5, y = 0.2, color = NULL,
                  label = paste0('r = ', corr,'\np.value = ', pval)), size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)) + 
    theme(axis.text.y = element_text(size = 5)) + 
    theme(text = element_text(size = 7)) + 
    labs(x = 'AUROC\n(geometric mean score)', y = 'AUROC\n(logistic regression score)')
  
}
plot_list[['V']] <- plot_list[['V']] + theme(legend.position = 'none') + labs(y = '')
plot_list[['V/B']] <-  plot_list[['V/B']] + theme(legend.position = 'none') + labs(y = '')

pout <- 
plot_list[['V']] + plot_list[['B']] + plot_list[['V/B']] + plot_layout(design = c('AAAAAAAAAAA
                                                                                  BBBBBBBB###
                                                                                  CCCCCCC####'))
ggsave(pout, filename = '~/Dropbox/geomMean_vs_logit_supp.png', height = 5.5, width = 12)
pdf('~/Dropbox/compendium_manuscript/figures/supplements/supp1_full.pdf', height = 5, width = 8.5)
pout
dev.off()

summarized_df <- full_df %>%
  group_by(Signature, Comparison) %>%
  summarize('logit' = mean(logit_), 'geom' = mean(zScore_geomMean))
plot_list <- list()
for(comparison in unique(full_df$Comparison)){
  plot_list[[comparison]] <- ggplot(summarized_df %>% filter(Comparison == comparison), 
                                    aes(y = logit, x = geom)) + 
    geom_point() + 
    theme_clean + 
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray50') +
    facet_grid(. ~ Signature) + 
    xlim(c(0,1)) + 
    ylim(c(0,1)) + 
    #scale_color_manual(values = color_vals) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = 'AUROC (geometric mean score', y = 'AUROC (logistic regression score)')
  
}

pout <- 
  plot_list[['V']] + plot_list[['B']] + plot_list[['V/B']] + plot_layout(design = c('AAAAAAAAAAA
                                                                                  BBBBBBBB###
                                                                                  CCCCCCC####'))
ggsave(pout, filename = '~/Dropbox/geomMean_vs_logit_summ.png', height = 5.5, width = 12)

ggplot(summarized_df, aes(y = logit, x = geom, color = Comparison)) + 
  geom_point() + 
  theme_clean + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray50') +
  xlim(c(0,1)) + 
  ylim(c(0,1)) + 
  scale_color_manual(values = color_vals) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

rlm_eqn <- function(df){
  m <- rlm(logit ~ geom, df, psi = psi.bisquare);
  #eq <- substitute(~~italic(r)^2~"="~r2, 
  #                 list(r2 = format(summary(m)$r.squared, digits = 3)))
  #as.character(as.expression(eq));
  return(m)
}

obj <- rlm_eqn(summarized_df)
obj

lm_eqn <- function(df){
  m <- lm(logit ~ geom, df);
  # eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
  #                  list(a = format(unname(coef(m)[1]), digits = 2),
  #                       b = format(unname(coef(m)[2]), digits = 2),
  #                       r2 = format(summary(m)$r.squared, digits = 3)))
  eq <- substitute(~~italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

option1 <- ggplot(summarized_df, aes(y = logit, x = geom)) + 
  geom_point(size = 2) + 
  theme_clean + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray50') +
  theme(legend.position = 'none') + 
  #geom_text(x = 0.8, y = 0.65, color = 'black', label = lm_eqn(summarized_df), parse = TRUE) + 
  #geom_smooth(aes(fill = NULL, color = NULL), method = 'lm', color = 'black')  +
  xlim(c(0,1)) + 
  ylim(c(0,1)) + 
  scale_fill_manual(values = color_vals) + 
  labs(x = 'AUROC (geometric mean score)', y = 'AUROC (logistic regression score)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs()
option1
pdf('~/Dropbox/compendium_manuscript/figures/supplements/supp1_condensed.pdf', height = 3.5, width = 3.5); option1; dev.off()

option2 <- ggplot(summarized_df, aes(y = logit, x = geom, fill = Comparison)) + 
  geom_point(pch = 21, size = 2) + 
  theme_clean + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray50') +
  geom_text(x = 0.8, y = 0.65, color = 'black', label = lm_eqn(summarized_df), parse = TRUE) + 
  geom_smooth(aes(fill = NULL, color = NULL), method = 'lm', color = 'black')  +
  xlim(c(0,1)) + 
  ylim(c(0,1)) + 
  scale_fill_manual(values = color_vals) + 
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

option3 <- ggplot(summarized_df, aes(y = logit, x = geom, fill = Comparison)) + 
  geom_point(pch = 21, size = 2) + 
  theme_clean + 
  #geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray50') +
  geom_text(x = 0.8, y = 0.65, color = 'black', label = lm_eqn(summarized_df), parse = TRUE) + 
  geom_smooth(aes(fill = NULL, color = NULL), method = 'lm', color = 'black')  +
  xlim(c(0,1)) + 
  ylim(c(0,1)) + 
  scale_fill_manual(values = color_vals) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = 'AUROC (geometric mean score)', y = 'AUROC (logistic regression score') + 
  theme(legend.position = 'none')

option4 <- ggplot(summarized_df, aes(y = logit, x = geom, fill = Comparison)) + 
  geom_point(pch = 21, size = 2) + 
  theme_clean + 
  #geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray50') +
  #geom_text(x = 0.8, y = 0.65, color = 'black', label = rlm_eqn(summarized_df), parse = TRUE) + 
  #geom_smooth(aes(fill = NULL, color = NULL), method = 'lm', color = 'black')  +
  geom_abline(slope = 0.5320738, intercept = 0.4373441) +
  xlim(c(0,1)) + 
  ylim(c(0,1)) + 
  scale_fill_manual(values = color_vals) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = 'AUROC (geometric mean score)', y = 'AUROC (logistic regression score', title = 'robust regression') + 
  theme(legend.position = 'none')

pout2 <- option1 + option2 + option3 + option4 + plot_layout(nrow = 2)
ggsave(plot = pout2, filename = '~/Dropbox/geomMean_vs_logit_4options.png', height = 8, width = 8)

plot_dat <- full_df %>%
  pivot_longer(names_to = 'method', values_to = 'score', cols = c(logit_, zScore_geomMean))

lm_eqn <- function(df){
  m <- lm(logit ~ geomMean, df);
  # eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
  #                  list(a = format(unname(coef(m)[1]), digits = 2),
  #                       b = format(unname(coef(m)[2]), digits = 2),
  #                       r2 = format(summary(m)$r.squared, digits = 3)))
  eq <- substitute(~~italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}




summarized_df <- plot_dat %>%
  #filter(Type == 'Discovery') %>% 
  group_by(Signature, method, Comparison) %>%
  summarize('med_score' = median(score))
sum_plot_dat <- summarized_df %>%
  pivot_wider(names_from = method, values_from = med_score) %>%
  dplyr::rename('logit' = logit_, 'geomMean' = zScore_geomMean) %>%
  mutate(Comparison = gsub(pattern = 'x', replacement = '/', x = Comparison)) 
option1 <- ggplot(sum_plot_dat, aes(x = logit, y = geomMean, fill = Comparison)) + 
  theme_clean +
  labs(x = 'AUROC (logistic regression)', y = 'AUROC (geometric mean score)') + 
  scale_fill_manual(values = color_vals) + 
  #xlim(c(0.5, 1)) + 
  #ylim(c(0.5,1)) + 
  geom_smooth(aes(fill = NULL, color = NULL), method = 'lm', color = 'black')  +
  geom_point(pch = 21, size = 4) + 
  coord_flip() + 
  geom_text(x = 0.8, y = 0.65, color = 'black', label = lm_eqn(sum_plot_dat), parse = TRUE) + 
  theme(legend.position = 'none', )  + 
  scale_x_continuous(limits = c(0.5, 1), breaks = seq(to = 1, from = 0.5, by = 0.1), minor_breaks = seq(to = 1, from = 0.5, by = 0.05)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(to = 1, from = 0.5, by = 0.1), minor_breaks = seq(to = 1, from = 0.5, by = 0.05))
option1
saveRDS(object = option1, file = '~/Dropbox/ggplot_objects/3c_logit_vs_geommean.RDS')


option2 <- ggplot(sum_plot_dat, aes(x = logit, y = geomMean, color = Comparison)) + 
  geom_point() + 
  theme_clean +
  labs(x = 'Logistic Regression AUC', y = 'Geometric Mean AUC') + 
  scale_color_manual(values = color_vals) + 
  #xlim(c(0.5, 1)) + 
  #ylim(c(0.5,1)) + 
  scale_x_continuous(limits = c(0.5, 1), breaks = seq(to = 1, from = 0.5, by = 0.1)) +
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(to = 1, from = 0.5, by = 0.1)) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') + 
  labs(title = 'option 2') +
  coord_flip()
option2

option1 + option2
saveRDS('3c_opt1.RDS', object = option1)
saveRDS('3c_opt2.RDS', object = option2)



pout <- ggplot(disc_df2, aes(x = Signature, y = Scores)) + 
  geom_boxplot() + 
  #geom_point(data = input_df %>% filter(Pt.Color != ''), aes(color = Pt.Color)) + 
  #text_geom +
  #geom_hline(yintercept = c(AUC_threshold), color = 'red', linetype = 'dashed') + 
  #geom_hline(yintercept = c(0.5), color = 'black', linetype = 'dashed') + 
  facet_grid(. ~ Comparison, scales = 'free', space = 'free') + 
  scale_alpha(guide = 'none') +
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'right',
        panel.spacing = unit(2, "lines")) + 
  scale_y_continuous(breaks = seq(to = 1, from = 0.5, by = 0.1), limits = c(0.5,1)) + 
  labs(y = 'AUC') 

pout + option1 + plot_layout(widths = c(4,1.5))
saveRDS(object = pout, file = '3d.RDS')
saveRDS(object = option1, file = '3c.RDS')
