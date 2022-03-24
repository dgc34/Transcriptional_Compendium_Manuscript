

male <- readRDS('~/Dropbox/male_eval_heatmap_list_generalizability_021422_healthy_sbj_ts.RDS')
female <- readRDS('~/Dropbox/female_eval_heatmap_list_generalizability_021422_healthy_sbj_ts.RDS')
# plot male vs female along with correlations?

male_df <- bind_rows(male) %>%
  #mutate('Imputed.Sex' = 'male') %>%
  dplyr::select(Study, Signature, Type, Scores, N) %>%
  dplyr::rename('Male.Scores' = Scores, 'Male.N' = N)
female_df <- bind_rows(female) %>%
  #mutate('Imputed.Sex' = 'female') %>%
  dplyr::select(Study, Signature, Type, Scores, N) %>%
  dplyr::rename('Female.Scores' = Scores, 'Female.N' = N)


sig_order <- c(paste0('VxB', 1:7), paste0('V', 1:12), paste0('B', 1:9))
plot_df <- male_df %>%
  left_join(female_df) %>%
  mutate('Comparison' = str_extract(pattern = 'VxB|V|B', string = Signature)) %>%
  mutate(Signature = factor(Signature, levels = sig_order, ordered = T)) 

mini_df <- plot_df %>%
  group_by(Comparison, Signature) %>%
  summarize('Male.AUC' = median(Male.Scores, na.rm = T), 'Female.AUC' = median(Female.Scores, na.rm = T)) %>%
  ungroup() %>%
  mutate(Comparison = gsub(pattern = 'VxB', replacement = 'V/B', x = Comparison))

lm_eqn <- function(df){
  m <- lm(Female.AUC ~ Male.AUC, df);
  # eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
  #                  list(a = format(unname(coef(m)[1]), digits = 2),
  #                       b = format(unname(coef(m)[2]), digits = 2),
  #                       r2 = format(summary(m)$r.squared, digits = 3)))
  eq <- substitute(~~italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


p4g <- ggplot(mini_df, aes(x = Male.AUC, y = Female.AUC)) + 
  geom_text(x = 0.65, y = 0.85, color = 'black', label = lm_eqn(mini_df), parse = TRUE) + 
  geom_point(pch = 21, fill = 'gray50') + 
  xlim(c(0.5,1)) + 
  ylim(c(0.5, 1)) + 
  scale_fill_manual(values = color_vals, name = '') + 
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray50') + 
  labs(x = 'Male median AUC', y = 'Female median AUC') +
  theme_clean +
  theme(legend.position = 'none')
p4g
saveRDS(object = p4g, '~/Dropbox/ggplot_objects/fig4_020922/fig4g.RDS')
