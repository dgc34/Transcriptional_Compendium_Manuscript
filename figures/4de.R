source('~/Dropbox/plot_palette.R')
library(dplyr)

AUC_threshold <- 0.7

sig_order <- rev(c(paste0('V', 1:15), paste0('B', 1:10)))

plot_dat_hiv <- readRDS('~/Dropbox/hiv_eval_heatmap_list_generalizability_021422_healthy_sbj_ts.RDS') %>%
  bind_rows() %>%
  mutate(Signature = factor(Signature, levels = sig_order, ordered = T)) %>% 
  group_by(Signature) %>% 
  mutate(med = median(Scores)) %>%
  ungroup() %>% 
  mutate('Specificity' = factor(ifelse(med <= AUC_threshold, 'not robust', 'V'))) %>% # label signatures with more than 25% of their values above the threshold red
  #dplyr::select(-med) %>%
  filter(grepl(pattern = 'V', x = Signature))

plot_dat_hiv %>%
  filter(Signature != 'V3') %>%
  group_by(Signature) %>%
  summarize('med' = median(Scores)) %>%
  pull(med) %>%
  summary()

plot_dat_bps <- readRDS('~/Dropbox/bp_eval_heatmap_list_generalizability_021422_healthy_sbj_ts.RDS') %>%
  bind_rows() %>% 
  mutate(Signature = factor(Signature, levels = sig_order, ordered = T)) %>% 
  group_by(Signature) %>%
  mutate(med = median(Scores)) %>%
  ungroup() %>% 
  mutate('Specificity' = factor(ifelse(med <= AUC_threshold, 'not robust', 'B'))) %>% # label signatures with more than 25% of their values above the threshold red
  dplyr::select(-med)
color_vals <- c('#F2BE84', '#C4E8F0', '#8FD0BA', 'gray50')
names(color_vals) <- c('B', 'V', 'V/B', 'not robust')#"#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#F0E442")
#names(color_vals) <- c('robust', 'not robust')

noninf_bps_df <- readRDS('~/Dropbox/bp_eval_heatmap_list_generalizability_021422_noninf_sbj_ts.RDS') %>%
  bind_rows() %>% 
  mutate(Signature = factor(Signature, levels = sig_order, ordered = T)) %>% 
  group_by(Signature) %>%
  mutate(med = median(Scores)) %>%
  ungroup() %>% 
  mutate('Specificity' = factor(ifelse(med <= AUC_threshold, 'not robust', 'B'))) %>% # label signatures with more than 25% of their values above the threshold red
  dplyr::select(-med)

plot_dat_bps %>%
  filter(!Signature %in% c('B1')) %>%
  group_by(Signature) %>%
  summarize('med' = median(Scores)) #%>%
  pull(med) %>%
  summary()

genPlot <- function(df, title_str){
  plot_grob <- df %>%
    ggplot(aes(y = Signature, x = Scores, fill = Specificity)) + 
    geom_boxplot() +
    theme_clean + 
    xlim(0, 1) + 
    scale_fill_manual(values = color_vals) +
    #geom_vline(xintercept = 0.5, linetype = 'dashed', color = 'black') +
    #geom_vline(xintercept = AUC_threshold, linetype = 'dashed', color = 'red') + 
    theme(legend.position = 'none') + 
    labs(x = 'AUC', title = title_str)
  
}

#dlims <- layer_scales(op3[[2]])$y$range_c$range[2]

p4d <- genPlot(plot_dat_hiv, 'HIV datasets') #+ scale_y_discrete(expand= expansion(add = c(0.1,0.2)))
p4d <- p4d + scale_y_discrete(expand = expansion(mult = c(0.07, .14)))
p4e <- genPlot(plot_dat_bps, 'B. pseudomallei datasets') + 
  scale_y_discrete(expand = expansion(mult = c(0.1, .2)))
p4f <- readRDS('~/Dropbox/ggplot_objects/fig4_020922/fig4f.RDS')
#p4f <- p4f + labs(title = 'Mycoplasma in Cohort', x = '') + 
#  theme(axis.text.y = element_text(size = 16, hjust = -0.1))

supp <- genPlot(noninf_bps_df, 'B. pseudomallei non-infectious controls')
pdf('~/Dropbox/compendium_manuscript/figures/supplements/supp3.pdf', height = 5.5, width = 4)
supp
dev.off()

p4d/p4e/p4f

saveRDS('~/Dropbox/ggplot_objects/fig4_020922/fig4d.RDS', object = p4d)
saveRDS('~/Dropbox/ggplot_objects/fig4_020922/fig4e.RDS', object = p4e)
# 
# design = c('AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             BBE
#             BBE
#             BBE
#             BBE
#             BBE
#             BBE
#             BBE
#             CCF
#             CCF
#             CCF
#             CCF
#             CCF
#             CCF
#             GHI
#             GHI
#             GHI
#             GHI
#             GHI
#             GHI
#             GHI
#             GHI')
# op3[[2]] <- op3[[2]] + labs(x = '', title = 'Independent validation datasets')
# op3[[3]] <- op3[[3]] + labs(x = '')
# p4h <- p4h + labs(y = 'Average AUC')
# p4d <- p4d + labs(y = '')
# p4e <- p4e + labs(y = '')
# p4f <- p4f + theme(axis.text.y = element_text(size = 18, hjust = -0.1))
# p4f2 <- p4f + 
#   geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'black') + 
#   geom_hline(yintercept = AUC_threshold, linetype = 'dashed', color = 'red')
# 
# fig_components <- list(op3[[2]], op3[[3]], op3[[1]], p4d, p4e, p4f2, p4g, p4h, p4i)
# saveRDS(object = fig_components, file = '~/Dropbox/compendium_manuscript/figures/ggplot objects/4_list.RDS')
# pout <- op3[[2]] + op3[[3]] + op3[[1]] + p4d+ p4e + p4f2 + p4g + p4h + p4i + plot_layout(design = design) + plot_annotation(tag_levels = 'A')
# pout
# saveRDS(object = pout, file = '~/Dropbox/ggplot_objects/4_assembled.RDS')
# ggsave(filename = '~/Dropbox/compendium_manuscript/figures/4_assembled.png', plot = pout, dpi = 300, height = 11, width = 8.5)
# 
# design = c('AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             AAD
#             BBE
#             BBE
#             BBE
#             BBE
#             BBE
#             BBE
#             BBE
#             CCF
#             CCF
#             CCF
#             CCF
#             CCF
#             CCF
#             GHI
#             GHI
#             GHI
#             GHI
#             GHI
#             GHI
#             GHI
#             GHI')
# design = rbind(
#   c(1, 1, 4),
#   c(1, 1, 4),
#   c(1, 1, 4),
#   c(1, 1, 4),
#   c(1, 1, 4),
#   c(1, 1, 4),
#   c(1, 1, 4),
#   c(1, 1, 4),
#   c(1, 1, 4),
#   c(1, 1, 4),
#   c(2, 2, 5), 
#   c(2, 2, 5),
#   c(2, 2, 5), 
#   c(2, 2, 5),
#   c(2, 2, 5), 
#   c(2, 2, 5),
#   c(2, 2, 5),
#   c(3, 3, 6),
#   c(3, 3, 6),
#   c(3, 3, 6),
#   c(3, 3, 6),
#   c(3, 3, 6),
#   c(3, 3, 6),
#   c(7, 8, 9),
#   c(7, 8, 9),
#   c(7, 8, 9),
#   c(7, 8, 9),
#   c(7, 8, 9),
#   c(7, 8, 9),
#   c(7, 8, 9),
#   c(7, 8, 9))
#          
# fig_components_gg <- fig_components
# fig_components_gg[[2]] <- fig_components_gg[[2]] + labs(title = ' ')
# fig_components_gg[[3]] <- fig_components_gg[[3]] + labs(title = ' ')
# fig_components_gg[[8]] <- fig_components_gg[[8]] + theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0))
# pout <- grid.arrange(grobs = fig_components_gg, layout_matrix = design)
# saveRDS(object = pout, file = '~/Dropbox/ggplot_objects/4_assembled_gridarrange.RDS')
# # op41 <- op3[[1]] + coord_flip()
# # op42 <- op3[[2]] + coord_flip()
# # op43 <- op3[[3]] + coord_flip()
# # 
# # p4d2 <- p4d + coord_flip()
# # p4e2 <- p4e + coord_flip()
# # p4f3 <- p4f2 + coord_flip()
# # 
# # 
# # d_vert <- c('
# #             AAAAAAAAAAABBBBBBBCCCCCC
# #             AAAAAAAAAAABBBBBBBCCCCCC
# #             AAAAAAAAAAABBBBBBBCCCCCC
# #             AAAAAAAAAAABBBBBBBCCCCCC
# #             AAAAAAAAAAABBBBBBBCCCCCC
# #             DDDDDDDDDDDEEEEEEEFFFFFF
# #             DDDDDDDDDDDEEEEEEEFFFFFF
# #             GGGGGGGGHHHHHHHHIIIIIIII
# #             GGGGGGGGHHHHHHHHIIIIIIII
# #             ')
# # pout2 <- op42 + op43 + op41 + p4d2 + p4e2 + p4f3 + p4g + p4h + p4i + plot_layout(design = d_vert) + plot_annotation(tag_levels = 'A')
# # ggsave(filename = '~/Dropbox/compendium_manuscript/figures/4_assembled_vert.png', plot = pout2, dpi = 300, height = 11, width = 8.5)
