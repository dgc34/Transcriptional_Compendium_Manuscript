library(ggpubr)
library(tidyr)
source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))
source('~/Dropbox/plot_palette.R')

signatures_file <- '~/Dropbox/compendium_manuscript/tables/supp_1.xlsx'
signatures <- read.xlsx(signatures_file)

MIobj <- readRDS('~/Dropbox/GSE103119.RDS')
table(MIobj$pheno$Pathogen)
table(MIobj$pheno$Class)
p <- MIobj$pheno
p <- p %>%
  mutate('Gram' = ifelse(grepl(pattern = 'Mycoplasma', x = Pathogen), 'Negative', 'Positive')) %>%
  mutate('Gram' = ifelse(`bacterial organism:ch1` %in% c('N/A', 'NONE'), NA, Gram)) %>%
  mutate('Gram' = ifelse(grepl(pattern = '+', x = `bacterial organism:ch1`, fixed = T), 'Mixed', Gram))
table(p$Class, p$Gram)
table(p$`bacterial organism:ch1`, p$Gram)
rownames(p) <- rownames(MIobj$pheno)
MIobj$pheno <- p
checkDataObject(MIobj, 'Dataset')

repairIndex <- function(index){
  ix <- is.na(index)
  index[ix] <- F
  return(index)
}

data_list <- list()
index_pos <- repairIndex((p$Gram == 'Positive' & p$Class %in% c('Bacteria')) | p$Class %in% c('Virus'))
index_neg <- repairIndex((p$Gram == 'Negative' & p$Class %in% c('Bacteria')) | p$Class %in% c('Virus'))
index_both <- repairIndex((p$Gram %in% c('Positive', 'Negative') & p$Class == 'Bacteria') | p$Class %in% c('Virus'))

MIobj_positive <- filterMIobj(MIobj, index_pos)
MIobj_negative <- filterMIobj(MIobj, index_neg)
MIobj_both <- filterMIobj(MIobj, index_both)

nrow(MIobj_positive$pheno)
nrow(MIobj_negative$pheno)
nrow(MIobj_both$pheno)

table(MIobj_positive$pheno$Pathogen, MIobj_positive$pheno$Class)
table(MIobj_negative$pheno$Pathogen, MIobj_negative$pheno$Class)
table(MIobj_both$pheno$Pathogen, MIobj_both$pheno$Class)

# barplot of subjects for positive and negative gram stains
signatures <- signatures %>%
  filter(Type == 'VxB')
nrow(signatures)

filterObjs <- list()
for(i in 1:nrow(signatures)){
  filterObjs[[i]] <- sig2Meta(signatures[i,])
}

score_df <- sapply(filterObjs, FUN = calculateScore, datasetObject = MIobj_both) %>%
  cbind()
signature_names <- sapply(filterObjs, function(x){return(x$filterDescription)}) 
colnames(score_df) <- signature_names
score_df <- score_df %>%
  as.data.frame()

plot_dat <- MIobj_both$pheno %>%
  dplyr::select(geo_accession, Pathogen, Class, Gram) %>%
  bind_cols(score_df) %>%
  reshape2::melt(id.vars = c('geo_accession', 'Pathogen', 'Class', 'Gram'), value.name = 'Score', variable.name = 'Signature') %>%
  mutate(Gram = ifelse(is.na(Gram), '', Gram)) %>%
  mutate(Infection = ifelse(Class == 'Virus', 'Virus', NA)) %>%
  mutate(Infection = ifelse(Gram == 'Positive', 'Gram positive bacteria', Infection)) %>%
  mutate(Infection = ifelse(Gram == 'Negative', 'Mycoplasma', Infection)) %>%
  mutate(Infection = factor(Infection, levels = c('Gram positive bacteria', 'Virus', 'Mycoplasma'), ordered = T))

signature_labels <- unique(plot_dat$Signature)

plot_dat2 <- plot_dat %>%
  mutate('Signature' = factor(Signature, level = signature_labels, ordered = T)) %>%
  mutate('xl' = factor(Infection, levels = c('Gram positive bacteria',
                                               'Mycoplasma', 
                                               'Virus'), 
                         labels = c('+', '-', 'vir'), 
                         ordered = T)) %>%
  group_by(Signature) %>%
  summarize('+' = auroc(bool = (xl == 'vir'), score = Score), '-' = auroc(bool = (xl == 'vir')[xl != '-'], score = Score[xl != '-'])) %>%
  pivot_longer(names_to = 'Type', values_to = 'AUC', cols = !Signature)

p4d <- ggplot(plot_dat2, aes(x = Type, y = AUC)) + 
  #geom_bar(stat = 'identity', position = position_dodge()) + 
  geom_boxplot() + 
  geom_point() + 
  theme_clean +
  geom_line(aes(group = Signature), color = 'gray70') +
  # theme(axis.text.x = element_text(angle = 0, hjust = 0.5), 
  #       text = element_text(size = 14)) + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1), 
  #      text = element_text(size = 14)) + 
  labs(y = 'VxB Signature AUC', title = 'GSE103119', x = '') +
  scale_y_continuous(breaks = seq(to = 1, from = 0, by = 0.1), limits = c(0,1)) + 
  geom_hline(yintercept = 0.5, color = 'red', linetype = 'dashed') + 
  stat_compare_means(paired = T, method = 'wilcox.test', position = position_nudge(x = 0.18, y = -0.08))
p4d

p4f <- ggplot(plot_dat2 %>% mutate(Comparison = 'V/B'), aes(x = Type, y = AUC, fill = Comparison)) + 
  geom_boxplot() + 
  geom_point() + 
  #geom_vline(xintercept = 0.5, linetype = 'dashed', color = 'black') + 
  scale_fill_manual(values = color_vals) +
  stat_compare_means(paired = T, method = 'wilcox.test', position = position_nudge(y = -0.85, x = 0.4)) + 
  ylim(c(0,1)) + 
  coord_flip() + 
  labs(title = 'Mycoplasma', x = 'Included in Cohort') + 
  theme(axis.text.y = element_text(size = 20)) + 
  theme_clean + 
  theme(legend.position = 'none')
p4f
saveRDS(object = p4f, file = '~/Dropbox/ggplot_objects/fig4_020922/fig4f.RDS')
saveRDS(object = p4d, file = '~/Dropbox/compendium_manuscript/figures/ggplot objects/4d.RDS')
ggsave(filename = '~/Dropbox/compendium_manuscript/figures/4d_mycoplasmaRemoval.png', height = 6, width = 3, dpi = 300)
colors <- createColorScale()
names(colors) <- unique(plot_dat$Infection)

# class_names <- names(MIobj_positive$class)
# MIobj_positive$class <- ifelse(MIobj_positive$pheno$Class == 'Virus', 1, 0)
# names(MIobj_positive$class) <- class_names
# lab_df <- data.frame('AUC' = sapply(filterObjs, 
#                                     calculateAUROC, 
#                                     MIobj = MIobj_positive, 
#                                     method = 'geomMean'),
#                      Signature = sapply(filterObjs, function(x){return(x$filterDescription %>% 
#                                                                          sapply(FUN = function(x){strsplit(x, split = '\n')[[1]][1]}) %>%
#                                                                          gsub(pattern = 'Viral_vs_Non.Viral', replacement = '', fixed = T) %>%
#                                                                          trimws())})) %>% 
#   mutate('Signature' = factor(Signature, level = signature_labels, ordered = T))
# 
# lab_df %>% summarize('median' = median(AUC))
# 
# p_out <- ggplot(plot_dat, aes(x = Infection, y = Score)) + 
#   geom_boxplot() + 
#   #geom_point() + 
#   facet_grid(. ~ Signature) + 
#   labs(y = 'Signature Score') +
#   stat_compare_means(method = 'wilcox.test', 
#                      comparisons = list(c('Virus', 
#                                           'Gram positive bacteria'),
#                                         c('Mycoplasma', 
#                                           'Gram positive bacteria'),
#                                         c('Mycoplasma', 
#                                           'Virus')),
#                      label = 'p.signif') +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 20)) +
#   theme(panel.grid.major.x = element_blank()) +
#   theme(panel.grid.minor.y = element_blank()) +
#   #theme(axis.text.x = element_blank()) +
#   scale_color_manual(values = colors, name = 'Infection Type')  + 
#   scale_x_discrete(labels = c('+ gram','virus', 'mycoplasma')) + 
#   labs(x = 'Infection Type') + 
#   geom_text(data = lab_df, aes(x = 2, y = -2.6, label = paste0('AUC = ', round(AUC, 2))))
# p_out
# ggsave(plot = p_out, filename = '~/Dropbox/3b_gse103119.png', height = 7.5, width = 9)
#   
#   #coord_flip()
# # auc plot for positive and negative gram stains
