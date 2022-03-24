# new plan
# 1. select top k signatures per bin
# 2. signature size vs. flu_test
# 3. type I IFN -log(p.value) for everything in the blue region vs. flu_test
# 4. any single gene frequencies
library(enrichR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(nsga2R)
library(stringr)

source('~/Dropbox/plot_palette.R')
color_vals <- c('gray80', 'black')
names(color_vals) <- c('', 'pareto front')
color_vals2 <- c('gray80', 'black')
names(color_vals2) <- c('excluded', 'included')
rank_threshold <- 20

#p1 <- readRDS('~/Dropbox/ggplot_objects/6c_flu_vs_healthy_scatter.RDS')
p2 <- readRDS('~/Dropbox/ggplot_objects/6d_flu_vs_virus_scatter.RDS')

#df1 <- p1$data
df2 <- p2$data

d2 = df2 %>%
  dplyr::select(flu_test, other_test) %>%
  mutate(x = 1-flu_test, y = other_test) %>%
  dplyr::select(x,y)
D2 = d2[order(d2$x,d2$y,decreasing=FALSE),]
front2 = D2[which(!duplicated(cummin(D2$y))),]
dim(front2)

front2_ix <- rownames(front2)

df2 <- df2 %>% 
  mutate(id = 1:nrow(df2)) %>% 
  mutate(color = ifelse(id %in% front2_ix, 'pareto front', ''))

# find points within range of the trend line fit to the pareto front
df2 <- df2 %>%
  filter(flu_test >= 0.65)
pf <- df2 %>% 
  filter(color == 'pareto front')

ggplot(df2, aes(x = flu_test, y = other_test, color = color)) + 
  geom_point() + 
  geom_smooth(data = pf, method = 'loess')

mod <- loess(other_test ~ flu_test, data = pf)

df5 <- df2 %>%
  mutate(pareto = ifelse(color == 'pareto front', 'pareto front', '')) %>% 
  dplyr::select(id, flu_test, other_test, posGenes, negGenes, pareto, k) %>% 
  mutate('pred' = predict(mod, newdata = df2)) %>%
  mutate('resid' = other_test-pred) %>%
  filter(flu_test >= 0.7) %>%
  mutate(bin = cut(flu_test, seq(to = 0.95, from = 0.7, by = 0.05), right = FALSE)) %>%
  #mutate(bin2 = cut(flu_test, seq(to = 0.95, from = 0.7, by = 0.05), right = FALSE)) %>%
  group_by(bin) %>%
  arrange(resid) %>%
  mutate(ranking = 1:n()) %>% 
  ungroup() %>%
  mutate(color = ifelse(ranking <= rank_threshold, 'included', 'excluded')) 


type_i_ifn <- read.table('~/Dropbox/0060337.tsv', header = T)$gene

# df5 <- df2 %>%
#   mutate(pareto = ifelse(color == 'pareto front', 'pareto front', '')) %>% 
#   dplyr::select(id, flu_test, other_test, posGenes, negGenes, pareto, k) %>% 
#   mutate('pred' = predict(mod, newdata = df2)) %>%
#   mutate('resid' = other_test-pred) %>%
#   filter(flu_test >= 0.7) %>%
#   mutate(bin = cut(flu_test, seq(to = 0.95, from = 0.7, by = 0.025), right = FALSE)) %>%
#   group_by(bin) %>%
#   arrange(resid) %>%
#   mutate(ranking = 1:n()) %>% 
#   ungroup() %>%
#   mutate(color = ifelse(pareto == 'pareto front', 'included', 'excluded')) 
fill_vals <- c('white', 'gray50')
names(fill_vals) <- c('pareto front', '')

cor(df5$flu_test, df5$other_test)
p1 <- ggplot(df5[sample(nrow(df5)/20, x = 1:nrow(df5)),], aes(x = flu_test, y = other_test)) + 
  #stat_density_2d(data = df5 %>% filter(color != 'included'),color = 'black') +
  geom_point(data = df5, alpha = 0.01, size = 1) + 
  geom_smooth(method = 'lm', se = F, color = 'black') + 
  geom_text(aes(x = 0.9, y = 0.6, label = 'r = 0.69')) + 
  theme_clean + 
  scale_color_manual(values = color_vals2) + 
  scale_fill_manual(values = fill_vals) + 
  theme(legend.position = 'none') + 
  labs(x = 'mean AUROC in influenza studies', y = 'mean AUC in non-influenza studies')

pdf('~/Dropbox/compendium_manuscript/figures/supplements/supp10.pdf', height = 3, width = 3.5)
p1
dev.off()

cor(df5$flu_test[df5$color == 'included'], df5$k[df5$color == 'included'])
p2 <- ggplot(df5 %>% 
         #filter(pareto == 'pareto front'), #%>%
         filter(color == 'included'),
         aes(x = flu_test, y = k)) + 
  geom_point(pch = 21, size = 3, aes(fill = pareto, color = NULL)) + 
  geom_smooth(method = 'lm', se=F, aes(color = color)) + 
  geom_text(data = data.frame(label = 'r = 0.49', flu_test = 0.75, k = 50), aes(label = label, color = NULL)) + 
  ylim(0,62) +  
  theme_clean + 
  scale_color_manual(values = color_vals2) + 
  scale_fill_manual(values = fill_vals) + 
  labs(y = 'signature size', x = 'mean AUC in influenza studies') + 
  theme(legend.position = 'none')

cor(df5$other_test[df5$color == 'included'], df5$k[df5$color == 'included'])^2
p2_other <- ggplot(df5 %>% 
               #filter(pareto == 'pareto front'), #%>%
               filter(color == 'included'),
             aes(x = other_test, y = k)) + 
  geom_point(pch = 21, size = 3, aes(fill = pareto, color = NULL)) + 
  geom_smooth(method = 'lm', se=F, aes(color = color)) + 
  geom_text(data = data.frame(label = 'r = 0.52', other_test = 0.5, k = 50), aes(label = label, color = NULL)) + 
  ylim(0,62) +  
  theme_clean + 
  scale_color_manual(values = color_vals2) + 
  scale_fill_manual(values = fill_vals) + 
  labs(y = 'signature size', x = 'mean AUC in non-influenza studies') + 
  theme(legend.position = 'none')

pdf('~/Dropbox/compendium_manuscript/figures/supplements/supp11.pdf', width = 3.5, height = 3)
p2_other
dev.off()


getGeneVect <- function(gene_df_row){
  if(nrow(gene_df_row) != 1){
    stop('only pass in one row')
  }
  gene_df_row %>%
    mutate(genes = paste0(posGenes, ';', negGenes)) %>%
    mutate(genes = posGenes) %>%
    pull(genes) %>%
    strsplit(split = ';') %>%
    unlist() %>%
    mapIds(x = org.Hs.eg.db, keytype = 'ENTREZID', column = 'SYMBOL') %>%
    unique() %>%
    return()
}

input_df <- df5 %>% 
  filter(color == 'included')
input_genes <- list()
for (i in 1:nrow(input_df)){
  input_genes[[i]] <- getGeneVect(input_df[i,])
}

num_ifn <- sapply(input_genes, function(x){return(sum(x %in% type_i_ifn))})
num_genes <- sapply(input_genes, length)
plot_dat <- input_df %>%
  mutate('num_overlap' = num_ifn, 'num_genes' = num_genes) %>%
  mutate('fraction_overlap' = num_overlap/num_genes) %>%
  mutate('num_non0' = num_overlap > 0)

plot_dat2 <- plot_dat %>%
  group_by(bin) %>%
  summarize('num_non0' = mean(num_non0))
p3 <- ggplot(plot_dat2, aes(x = bin, y = num_non0, group = 1)) +
  #geom_bar(stat = 'identity', color = 'black', fill = 'black') +
  geom_line() +
  geom_point() +
  theme_bw() +
  #ylim(c(0,.75)) +
  theme_clean +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = 'fraction of signatures\nwith type I ifn genes', x = 'mean AUC in influenza studies')
p3

# what does specificity look like across the range of bins?
# check aging first
source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))
sig_input_df <- input_df %>% 
  dplyr::rename(Positive.Genes = posGenes, Negative.Genes = negGenes) %>%
  mutate('Author' = 'perm', 'Positive.Class' = 'influenza', 'Negative.Class' = 'healthy', 'Discovery.Accessions' = '', 'Validation.Accessions' = '')
sig_list <- list()
for (i in 1:nrow(sig_input_df)){
  current_sig <- list()
  current_sig$posGeneNames <- as.character(parseSets(sig_input_df$Positive.Genes[i]))
  current_sig$negGeneNames <- as.character(parseSets(sig_input_df$Negative.Genes[i]))
  current_sig$filterDescription <- i
  sig_list[[i]] <- current_sig
}

aging_list <- readRDS(file = '~/Documents/darpa-manuscript-data/data_list_final/aging_list.RDS')

sig_input_df <- sig_input_df %>%
  mutate(sig_id = 1:nrow(sig_input_df))
output_list <- list()
for(i in 1:length(sig_list)){
  output_list[[i]] <- plotDatSignature(signature = sig_list[[i]], data_list = aging_list) %>%
    mutate(sig_id = i)
}

aging_plot_dat <- bind_rows(output_list) %>%
  left_join(sig_input_df %>% dplyr::select(k, sig_id, bin))

aging_plot_dat2 <- aging_plot_dat %>%
  group_by(sig_id, bin) %>%
  summarize('med' = median(Scores))

ggplot(aging_plot_dat2, aes(x = bin, y = med)) + 
  geom_boxplot() + 
  theme_clean + 
  #ylim(c(0,1)) + 
  geom_hline(yintercept = 0.6, color = 'gray50', linetype = 'dashed') + 
  labs('title' = 'Influenza Signatures, Aging Data')
