library(ggplot2)
library(dplyr)
library(scales)

initial_signature <- readRDS('~/Dropbox/flu_vs_virus.RDS')
df <- read.csv('~/Dropbox/random_flu_vs_virus_121121.csv')

df <- bind_rows(df, scoreSignature(initial_signature, sum(length(initial_signature$posGeneNames), length(initial_signature$negGeneNames))) %>% dplyr::select(-c(posGenes, negGenes)))
# evaluate full signature and append with label

max_ix <- df %>%
  #mutate('diff' = flu_test - other_train) %>%
  mutate('diff' = flu_test - other_test) %>% 
  pull() %>%
  which.max()

pt_color = 'steelblue2'
df <- df %>%
  #mutate('diff' = flu_test - other_train) %>%
  mutate('diff' = flu_test - other_test) %>% 
  #mutate(col_pt = ifelse(diff != max(diff), 'black', 'darkorange2')) %>%
  mutate('col_pt' = ifelse(k != max(k), 'black', pt_color)) %>%
  #mutate(lab_pt = ifelse(diff != max(diff), '', 'Maximum Difference (Specificity)\n 5 genes')) %>%
  mutate('lab_pt' = ifelse(k != max(k), '', paste0('Full Signature\n', df[nrow(df),'k'],' genes')))
  #mutate(col_pt = ifelse(k != max(k), '', 'steelblue3')) %>%
  #mutate(lab_pt = ifelse(k != max(k), '', 'Full Signature\n152 genes')) 

## merge other_test and other_train based on the number of samples in each
other_train_list <- c(readRDS('~/Documents/darpa-manuscript-data/flu_data_v3/other_disc.RDS'), readRDS('~/Documents/darpa-manuscript-data/flu_data_v3/other_dev.RDS'))
other_test_list <- readRDS('~/Documents/darpa-manuscript-data/flu_data_v3/other_test.RDS')
N1 <- sum(sapply(other_train_list, function(x){nrow(x$pheno)}))
N2 <- sum(sapply(other_test_list, function(x){nrow(x$pheno)}))

# df2 <- df %>%
#   ungroup() %>% 
#   mutate('other' = (N1 * other_train + N2 * other_test) /(N1 + N2))

p_out <- ggplot(df, aes(x = flu_test, y = other_test))  +
#p_out <- ggplot(df2, aes(x = flu_test, y = other))  +
  theme_bw() + 
  geom_point(color = df$col_pt, alpha = 0.8) + 
  #geom_abline(intercept = 0, slope = 1) + 
  labs(x = 'Avg. AUC Influenza Datasets', y = 'Avg. AUC Non-Influenza Datasets') +
  scale_x_continuous(limits = c(0.5, 1.01), expand = c(0,0)) +
  scale_y_continuous(limits = c(0.5, 1.01), expand = c(0,0)) +
  #geom_smooth(method = 'lm') + 
  geom_line(data = data.frame(x = seq(to = 2, from = -1, by = 0.01), y = seq(to = 2, from = -1, by = 0.01)), aes(x = x, y = y), color = 'red', linetype = 'dashed') + 
  geom_label(data = df %>% 
               filter(df$col_pt != 'black'), 
             aes(x = flu_test, y = other_test, label = lab_pt), color = pt_color, position = position_nudge(x = -0.035, y = 0.07)) + 
  theme(axis.text = element_text(size = 12))
p_out
saveRDS(object = p_out, file = '~/Dropbox/compendium_manuscript/figures/ggplot objects/6d_flu_vs_virus_scatter.RDS')
ggsave(filename = '~/Dropbox/compendium_manuscript/figures/6c_flu_sig_cloud_vs_other.png', plot = p_out, height = 5, width = 5, units = 'in')

hist(df$diff)

### investigation of points that are strong non-flu predictors but poor flu-predictors
df %>%
  filter(other_train > 0.75, flu_test < 0.76)

density_plots <- df %>%
  reshape2::melt(id.vars = c('posGenes', 'negGenes', 'k', 'col_pt', 'lab_pt'), variable.name = 'type', value.name = 'auc') %>%
  filter(type %in% c('flu_train', 'flu_test', 'other_train', 'other_test'))
p_out <- ggplot(density_plots, aes(x = auc, color = type))  +
  geom_density() + 
  xlim(c(0,1))
p_out
