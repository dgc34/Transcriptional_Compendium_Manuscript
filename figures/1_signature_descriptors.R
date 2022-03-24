library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(enrichR)
library(ggwordcloud)
library(Biobase)
library(org.Hs.eg.db)
library(patchwork)
library(stringr)
library(gridExtra)
library(grid)

source('~/Dropbox/plot_palette.R')
# first figure: describe the signature curation
# options to explore:
## histogram of signature sizes
## word clouds with most frequent genes
## histogram of top enrichment terms
### may split each of these based on the class of signature (V, B, VxB)

# read signatures
sigs <- read.xlsx('~/Dropbox/compendium_manuscript/tables/supp_1.xlsx') %>%
  mutate(genes = paste0(Positive.Genes,';', Negative.Genes))
head(sigs$genes)

# split into character vector, remove white space, dashes, and blanks
sig_list <- lapply(sigs$genes, strsplit, split = ';') %>%
  lapply(function(x){as.character(setdiff(sapply(x, trimws), c('—', '')))})

# map onto latest symbols
sig_list <- sig_list %>%
  lapply(mapIds, x = org.Hs.eg.db, keytype = 'SYMBOL', column = 'SYMBOL')
names(sig_list) <- sigs$Signature.Label

# count the number of each sig
p1b <- ggplot(sigs %>% 
         filter(Type != 'Influenza') %>%
         mutate(Type = factor(Type, levels = c('VxB', 'Virus', 'Bacteria'), ordered = T)), 
       aes(x = Type)) + 
  geom_bar(stat = 'count') +
  theme_clean + 
  scale_y_continuous(breaks = 0:11) + 
  #coord_flip()+ 
  labs(y = 'Count', x = 'Signature Type') 

## signature sizes
sizes <- sapply(sig_list, length)
length(sizes)
nrow(sigs)
plot_dat <- sigs %>% 
  mutate('Size' = sizes) #%>%
  #filter(Type != 'Influenza') #%>%
  #filter(Size < 200)
# option 1: 
p1c <- ggplot(plot_dat, aes(x = Size)) + 
  geom_histogram(breaks = seq(to = 400, from = 0, by = 10)) + 
  #scale_x_continuous(breaks = c(seq(to = 80, from = 0, by = 20), 100, 200, 300, 400)) + 
  scale_x_continuous(breaks = c(seq(to = 400, from = 0, by = 50))) + 
  scale_y_continuous(breaks = 0:8) + 
  labs(x = 'Genes per Signature', y = 'Count') +
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p1c
# option 2:
p2 <- ggplot(plot_dat, aes(x = log2(Size))) + 
  geom_histogram() + 
  labs(x = 'log2(Genes per Signature)', y = 'Count') +
  theme_clean
# option 3:
p3 <- ggplot(plot_dat, aes(x = Size)) + 
  geom_histogram() + 
  labs(x = 'Genes per Signature', y = 'Count') +
  facet_grid(Type ~ .) +
  theme_clean
# option 4:
p4 <- ggplot(plot_dat, aes(x = log2(Size))) + 
  geom_histogram() + 
  labs(x = 'log2(Genes per Signature)', y = 'Count') +
  facet_grid(Type ~ .) +
  theme_clean
#p1 + p2
#p3 | p4

## histogram of top genes
gene_count_df <- sigs %>%
  filter(Type != 'Influenza') %>% 
  separate_rows(genes, sep = ';') %>%
  mutate(gene = trimws(genes)) %>% 
  filter(!gene %in% c('—', '', 'NA')) %>%
  dplyr::select(Signature.Label, Type, gene) %>% 
  mutate(gene = mapIds(x = org.Hs.eg.db, 
                        keys = gene, 
                        keytype = 'SYMBOL', 
                        column = 'SYMBOL')) 
top_genes <- gene_count_df %>%
  group_by(gene) %>%
  summarize('count' = n()) %>%
  arrange(desc(count)) %>%
  head(25) %>%
  pull(gene)
# gene_count_df %>%
#   group_by(Type, gene) %>%
#   summarize('count' = n()) %>%
#   arrange(desc(count))

ggplot(gene_count_df %>% 
         filter(gene %in% top_genes) %>%
         mutate(gene = factor(gene, levels = top_genes, ordered = T))
       , aes(x = gene)) + 
  geom_bar(stat = 'count') + 
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  labs(x = 'Most Frequent Signature Genes', y = 'Count')

# easier to do this as a plot loop than a facet
# what do we want? top positive genes and top negative genes per signature type

getTopGenes <- function(input_df, N = 5){
  gene_count_df <- input_df %>%
    #filter(Type != 'Influenza') %>% 
    group_by(Signature.Label) %>% 
    separate_rows(Positive.Genes, sep = ';') %>%
    mutate(gene = trimws(Positive.Genes)) %>% 
    filter(!gene %in% c('—', '', 'NA')) %>%
    dplyr::select(Signature.Label, Type, gene) %>% 
    mutate(gene = mapIds(x = org.Hs.eg.db, 
                         keys = gene, 
                         keytype = 'SYMBOL', 
                         column = 'SYMBOL')) %>%
    distinct()
  top_genes <- gene_count_df %>%
    group_by(gene) %>%
    summarize('count' = n()) %>%
    arrange(desc(count)) %>%
    top_n(N, count) %>%
    arrange(desc(count)) %>%
    pull(gene) 
  return(gene_count_df %>% 
           filter(gene %in% top_genes) %>% 
           group_by(gene) %>%
           summarize('fraction' = n()/length(unique(gene_count_df$Signature.Label))))
}
N <- 3

# positive and negative mixed is confusing: let's separate them
input_df1 <- sigs %>% filter(Type == 'Virus')
input_df2 <- sigs %>% filter(Type == 'Bacteria')
input_df3 <- sigs %>% filter(Type == 'VxB')
input_df4 <- sigs %>% filter(Type == 'Influenza')
top_vir <- getTopGenes(input_df1, N = 11)
top_bac <- getTopGenes(input_df2)
top_vxb <- getTopGenes(input_df3)
top_flu <- getTopGenes(input_df4)

p1 <- ggplot(top_vir %>% arrange(desc(fraction)) %>% mutate(gene = factor(gene, levels = gene, ordered = T)), aes(x = gene, y = fraction)) + 
  geom_bar(stat = 'identity') +
  ylim(c(0,1)) + 
  labs(x = 'Gene', y = 'Fraction of Signatures', title = 'Virus') + 
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2 <- ggplot(top_bac %>% arrange(desc(fraction)) %>% mutate(gene = factor(gene, levels = gene, ordered = T)), aes(x = gene, y = fraction)) + 
  geom_bar(stat = 'identity') +
  ylim(c(0,1)) +  
  labs(x = 'Gene', y = 'Fraction of Signatures', title = 'Bacteria') + 
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p3 <- ggplot(top_vxb %>% arrange(desc(fraction)) %>% mutate(gene = factor(gene, levels = gene, ordered = T)), aes(x = gene, y = fraction)) + 
  geom_bar(stat = 'identity') +
  ylim(c(0,1)) + 
  labs(x = 'Gene', y = 'Fraction of Signatures', title = 'VxB') + 
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p4 <- ggplot(top_flu %>% arrange(desc(fraction)) %>% mutate(gene = factor(gene, levels = gene, ordered = T)), aes(x = gene, y = fraction)) + 
  geom_bar(stat = 'identity') +
  ylim(c(0,1)) + 
  labs(x = 'Gene', y = 'Fraction of Signatures', title = 'Flu') + 
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


p1 + p2 + p3 + p4

ggplot(gene_count_df %>% 
         filter(gene %in% top_genes) %>% 
         mutate(gene = factor(gene, levels = top_genes, ordered = T))
       , aes(x = gene)) + 
  geom_bar(stat = 'count') + 
  theme_clean + 
  facet_grid(. ~ Type, scales = 'free', space = 'free') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = 'Most Frequent Signature Genes', y = 'Count')

## two sided bar plot of counts
### let's demo with bacteria first bc fewest signatures
plot_list <- list()
for (type in c('Virus', 'Bacteria', 'VxB', 'Influenza')){
  pos_gene_counts <- sigs %>%
    filter(Type == type) %>%
    separate_rows(Positive.Genes, sep = ';') %>%
    mutate(gene = trimws(Positive.Genes)) %>% 
    filter(!gene %in% c('—', '', 'NA')) %>%
    dplyr::select(Signature.Label, Type, gene) %>% 
    mutate(gene = mapIds(x = org.Hs.eg.db, 
                         keys = gene, 
                         keytype = 'SYMBOL', 
                         column = 'SYMBOL')) %>%
    distinct() %>%
    group_by(gene) %>%
    count() %>%
    arrange(desc(n)) %>%
    mutate('direction' = 'positive')
  
  neg_gene_counts <- sigs %>%
    filter(Type == type) %>%
    separate_rows(Negative.Genes, sep = ';') %>%
    mutate(gene = trimws(Negative.Genes)) %>% 
    filter(!gene %in% c('—', '', 'NA')) %>%
    dplyr::select(Signature.Label, Type, gene) %>% 
    mutate(gene = mapIds(x = org.Hs.eg.db, 
                         keys = gene, 
                         keytype = 'SYMBOL', 
                         column = 'SYMBOL')) %>%
    distinct() %>%
    group_by(gene) %>%
    count() %>%
    arrange(desc(n)) %>%
    mutate('direction' = 'negative')
  df <- bind_rows(pos_gene_counts, neg_gene_counts) %>%
    mutate('Type' = type)
  plot_list[[type]] <- df
}

all_counts <- bind_rows(plot_list) %>%
  filter(n > 1) %>%
  #mutate(n = ifelse(direction == 'negative', -1*n, n)) %>%
  mutate(gene = ifelse(Type == 'Virus', paste0(' ', gene), gene))

gene_order <- all_counts %>%
  group_by(gene) %>%
  summarize(n = sum(n)) %>% 
  arrange(desc(n))

color_vals <- c('gray20', 'gray80')
names(color_vals) <- c('Positive', 'Negative')
ggplot(all_counts %>% 
         mutate(gene = factor(gene, levels = gene_order$gene, ordered = T)) %>%
         #mutate('Direction' = ifelse(n > 0, 'Positive', 'Negative')) %>%
         mutate('Direction' = ifelse(direction == 'positive', 'Positive', 'Negative')) %>%
         mutate(Direction = factor(Direction, levels = c('Positive', 'Negative'), ordered = T)), 
       aes(x = gene, y = n, fill = Direction)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(. ~ Type, scales = 'free', space = 'free') + 
  scale_fill_manual('Gene Direction', values = color_vals) + 
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = 'Recurring Genes', y = 'Number of Signatures') + 
  scale_y_continuous(limits = c(-2.5,6.5), breaks = -2:6, labels = c(2,1,0:6))

p1dx <- ggplot(all_counts %>% 
         mutate(Type = factor(Type, levels = c('VxB', 'Virus', 'Bacteria'), ordered = T)) %>% 
         mutate(gene = factor(gene, levels = gene_order$gene, ordered = T)) %>%
         #mutate('Direction' = ifelse(n > 0, 'Positive', 'Negative')) %>%
         mutate('Direction' = ifelse(direction == 'positive', 'Positive', 'Negative')) %>%
         mutate(Direction = factor(Direction, levels = c('Positive', 'Negative'), ordered = T)), 
       aes(x = gene, y = n, fill = Direction)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(. ~ Type, scales = 'free', space = 'free') + 
  scale_fill_manual('Gene Direction', values = color_vals) + 
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = 'Recurring Genes', y = 'Number of Signatures') + 
  scale_y_continuous(breaks = 0:6) 

gene_order <- all_counts %>%
  group_by(gene) %>%
  summarize(n = sum(n)) %>% 
  arrange(n)
p1dy <- ggplot(all_counts %>% 
                mutate(Type = factor(Type, levels = c('VxB', 'Virus', 'Bacteria'), ordered = T)) %>% 
                mutate(gene = factor(gene, levels = gene_order$gene, ordered = T)) %>%
                #mutate('Direction' = ifelse(n > 0, 'Positive', 'Negative')) %>%
                mutate('Direction' = ifelse(direction == 'positive', 'Positive', 'Negative')) %>%
                mutate(Direction = factor(Direction, levels = c('Positive', 'Negative'), ordered = T)), 
              aes(x = n, y = gene, fill = Direction)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(Type ~ ., scales = 'free', space = 'free') + 
  scale_fill_manual('Gene Direction', values = color_vals) + 
  theme_clean + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(y = 'Recurring Genes', x = 'Number of Signatures') + 
  scale_x_continuous(breaks = 0:6) 


## gene frequencies as word clouds for each class (not enough variation in frequency for this tbh)

## histogram of top enrichment terms per category
minimum_enrichment_size <- 7
ix <- sizes >= minimum_enrichment_size
table(ix)
sizes

sigs <- sigs %>%
  mutate('enrichment' = ix)
table(sigs$Type, sigs$enrichment)

enrichment_list <- sig_list[ix]
names(enrichment_list) <- sigs$Signature.Label[ix]

current_sig <- 1

enrichment_df_list <- list()
for (current_sig in seq_along(enrichment_list)){
  current_gene_list <- enrichment_list[[current_sig]]
  enriched_term_list <- enrichr(genes = current_gene_list, databases = c('GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021', 'KEGG_2021_Human'))
  df <- enriched_term_list %>%
    bind_rows() %>%
    mutate('Signature.Label' = names(enrichment_list)[current_sig]) 
  enrichment_df_list[[current_sig]] <- df
}

enrichment_highlights <- enrichment_df_list %>%
  bind_rows() %>%
  filter(Adjusted.P.value <= 0.01) %>% 
  mutate(Type = str_extract(pattern = 'VxB|V|B|I', string = Signature.Label)) %>% 
  #filter(Type != 'I') %>% # remove influenza signatures
  group_by(Type) %>%
  mutate('Total.Sigs' = length(unique(Signature.Label))) %>%
  ungroup() %>% 
  group_by(Type, Term, Total.Sigs) %>%
  count()

#table(enrichment_highlights$Type, enrichment_highlights$Signature.Label)

hist(enrichment_highlights %>% filter(Type == 'VxB') %>% pull(n)) 
hist(enrichment_highlights %>% filter(Type == 'V') %>% pull(n)) 
hist(enrichment_highlights %>% filter(Type == 'B') %>% pull(n)) 
hist(enrichment_highlights %>% filter(Type == 'I') %>% pull(n))

# maybe pull out the top enrichment terms per type?

output_table <- enrichment_highlights %>% 
  group_by(Type) %>%
  filter(n >= (max(n)-0), grepl(pattern = 'GO:', x = Term)) %>%
  mutate('Fraction of Signatures Enriched' = paste0(n, '/', Total.Sigs)) %>%
  dplyr::select(Type, Term, `Fraction of Signatures Enriched`) %>%
  mutate(Type = gsub(pattern = '^B$', replacement = 'Bacteria', x = Type)) %>% 
  mutate(Type = gsub(pattern = '^V$', replacement = 'Virus', x = Type)) %>% 
  mutate(Type = factor(Type, levels = c('VxB', 'Virus', 'Bacteria'), ordered = T)) %>%
  arrange(Type) %>%
  dplyr::rename(`Signature Type` = Type, `GO Term` = Term)


### aggregated gene enrichments
bact_genes <- sig_list[str_detect(pattern = '^B\\d+$', string = names(sig_list))] %>%
  unlist() %>%
  unique()
vir_genes <- sig_list[str_detect(pattern = '^V\\d+$', string = names(sig_list))] %>%
  unlist() %>%
  unique()
vxb_genes <- sig_list[str_detect(pattern = '^VxB\\d+$', string = names(sig_list))] %>%
  unlist() %>%
  unique()

procEnr <- function(input_list){
  input_list %>%
    bind_rows() %>%
    arrange(Adjusted.P.value) %>%
    return()
}

GO_vir <- enrichr(genes = vir_genes, databases = c('GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021')) %>% procEnr()
kegg_vir <- enrichr(genes = vir_genes, databases = c('KEGG_2021_Human')) %>% procEnr()

GO_bac <- enrichr(genes = bact_genes, databases = c('GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021')) %>% procEnr()
kegg_bac <- enrichr(genes = bact_genes, databases = c('KEGG_2021_Human')) %>% procEnr()

GO_vxb <- enrichr(genes = vxb_genes, databases = c('GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021')) %>% procEnr()
kegg_vxb <- enrichr(genes = vxb_genes, databases = c('KEGG_2021_Human')) %>% procEnr()

output_xlsx <- list('GO vir' = GO_vir,
                    'KEGG_vir' = kegg_vir,
                    'GO_bac' = GO_bac,
                    'KEGG_bac' = kegg_bac,
                    'GO_vxb' = GO_vxb,
                    'KEGG_vxb' = kegg_vxb
                    )
write.xlsx(output_xlsx, "~/Dropbox/aggregated_enrichment_terms.xlsx")

tt = ttheme_default(colhead=list(fg_params=list(rot=0)))
p1e <- tableGrob(output_table, theme = tt)

layout <- c('AAAB
             AAAC
             DDDD
             EEE#')

p0 <- ggplot() + 
  theme_clean + 
  theme(panel.border = element_blank())
pout <- p0 + p1b + p1c + p1dx + p1e + plot_layout(design = layout) + plot_annotation(tag_levels = 'A')
ggsave(pout, filename = '~/Dropbox/compendium_manuscript/figures/1_signature_curation.png', height = 18, width = (8.5/11*18), dpi = 300)


#### signed enrichment 
# read signatures
signed_sigs <- read.xlsx('~/Dropbox/compendium_manuscript/tables/supp_1.xlsx') 

# split into character vector, remove white space, dashes, and blanks
pos_sig_list <- lapply(sigs$Positive.Genes, strsplit, split = ';') %>%
  lapply(function(x){as.character(setdiff(sapply(x, trimws), c('—', '')))})
neg_sig_list <- lapply(sigs$Negative.Genes, strsplit, split = ';') %>%
  lapply(function(x){as.character(setdiff(sapply(x, trimws), c('—', '')))})

# map onto latest symbols
pos_sig_list <- pos_sig_list %>%
  lapply(mapIds, x = org.Hs.eg.db, keytype = 'SYMBOL', column = 'SYMBOL')
names(pos_sig_list) <- paste0(sigs$Signature.Label, '_pos')

neg_sig_list <- neg_sig_list %>%
  lapply(FUN = function(keys){
    tryCatch(mapIds(keys = keys, x = org.Hs.eg.db, keytype = 'SYMBOL', column = 'SYMBOL'),
             error = function(e){
               print(e)
               return(NULL)
               }) %>%
      return()})
names(neg_sig_list) <- paste0(sigs$Signature.Label, '_neg')

pos_vxb <- pos_sig_list[grepl(pattern = 'VxB', x = names(pos_sig_list))] %>%
  unlist() %>% unique()
neg_vxb <- neg_sig_list[grepl(pattern = 'VxB', x = names(neg_sig_list))] %>%
  unlist() %>% unique()

pos_v <- pos_sig_list[grepl(pattern = '^V\\d', x = names(pos_sig_list))] %>%
  unlist() %>% unique()
neg_v <- neg_sig_list[grepl(pattern = 'V\\d', x = names(neg_sig_list))] %>%
  unlist() %>% unique()

pos_b <- pos_sig_list[grepl(pattern = '^B', x = names(pos_sig_list))] %>%
  unlist() %>% unique()
neg_b <- neg_sig_list[grepl(pattern = '^B', x = names(neg_sig_list))] %>%
  unlist() %>% unique()


pos_flu <- pos_sig_list[grepl(pattern = '^I', x = names(pos_sig_list))] %>%
  unlist() %>% unique()
neg_flu <- neg_sig_list[grepl(pattern = '^I', x = names(neg_sig_list))] %>%
  unlist() %>% unique()

enrichment_list <- list(pos_vxb, neg_vxb, pos_v, neg_v, pos_b, neg_b, pos_flu, neg_flu)
names(enrichment_list) <- c('pos_vxb', 'neg_vxb', 'pos_v', 'neg_v', 'pos_b', 'neg_b', 'pos_flu', 'neg_flu')


enrichment_df_list <- list()
for (current_sig in seq_along(enrichment_list)){
  current_gene_list <- enrichment_list[[current_sig]]
  enriched_term_list <- enrichr(genes = current_gene_list, databases = c('GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021'))
  df <- enriched_term_list %>%
    bind_rows() %>%
    mutate('Signature.Label' = names(enrichment_list)[current_sig]) 
  enrichment_df_list[[current_sig]] <- df
}
length(enrichment_df_list)
names(enrichment_df_list) <- names(enrichment_list)

enrichment_df_list <- lapply(enrichment_df_list, procEnr)
enrichment_df_list[[1]]

output_xlsx <- list('GO vir pos' = enrichment_df_list$pos_v,
                    'GO vir neg' = enrichment_df_list$neg_v,
                    'GO bac pos' = enrichment_df_list$pos_b,
                    'GO bac neg' = enrichment_df_list$neg_b,
                    'GO vxb pos' = enrichment_df_list$pos_vxb,
                    'GO vxb neg' = enrichment_df_list$neg_vxb,
                    'GO flu pos' = enrichment_df_list$pos_flu,
                    'GO flu neg' = enrichment_df_list$neg_flu 
)
write.xlsx(output_xlsx, "~/Dropbox/aggregated_signed_enrichment_terms_with_flu.xlsx")

### find mutually exclusive terms
# filter to q <= 0.05
enrichment_df_list_filtered <- enrichment_df_list %>%
  lapply(function(df){df %>% filter(Adjusted.P.value <= 0.05)})

# get all terms per group 
vxb_pos_terms <- enrichment_df_list_filtered$pos_vxb %>% pull(Term)
v_pos_terms <- enrichment_df_list_filtered$pos_v %>% pull(Term)
b_pos_terms <- enrichment_df_list_filtered$pos_b %>% pull(Term)

vxb_neg_terms <- enrichment_df_list_filtered$neg_vxb %>% pull(Term) # no significant terms
v_neg_terms <- enrichment_df_list_filtered$neg_v %>% pull(Term)
b_neg_terms <- enrichment_df_list_filtered$neg_b %>% pull(Term) # no significant terms

# make a table of filtered terms
output_xlsx <- list('GO vir pos' = enrichment_df_list$pos_v %>% filter(Term %in% v_pos_terms),
                    'GO vir neg' = enrichment_df_list$neg_v %>% filter(Term %in% v_neg_terms),
                    'GO bac pos' = enrichment_df_list$pos_b %>% filter(Term %in% b_pos_terms),
                    'GO bac neg' = enrichment_df_list$neg_b %>% filter(Term %in% b_neg_terms),
                    'GO vxb pos' = enrichment_df_list$pos_vxb %>% filter(Term %in% vxb_pos_terms),
                    'GO vxb neg' = enrichment_df_list$neg_vxb %>% filter(Term %in% vxb_neg_terms) 
)
write.xlsx(output_xlsx, "~/Dropbox/aggregated_signed_enrichment_terms_significant.xlsx")


# get mutually exclusive positive terms
vxb_pos_terms_ME <- setdiff(vxb_pos_terms, c(v_pos_terms, b_pos_terms))
v_pos_terms_ME <- setdiff(v_pos_terms, c(vxb_pos_terms, b_pos_terms))
b_pos_terms_ME <- setdiff(b_pos_terms, c(vxb_pos_terms, v_pos_terms))

vxb_neg_terms_ME <- setdiff(vxb_neg_terms, c(v_neg_terms, b_neg_terms))
v_neg_terms_ME <- setdiff(v_neg_terms, c(vxb_neg_terms, b_neg_terms))
b_neg_terms_ME <- setdiff(b_neg_terms, c(vxb_neg_terms, v_neg_terms))

# make a table of the exclusive terms
output_xlsx <- list('GO vir pos' = enrichment_df_list$pos_v %>% filter(Term %in% v_pos_terms_ME),
                    'GO vir neg' = enrichment_df_list$neg_v %>% filter(Term %in% v_neg_terms_ME),
                    'GO bac pos' = enrichment_df_list$pos_b %>% filter(Term %in% b_pos_terms_ME),
                    'GO bac neg' = enrichment_df_list$neg_b %>% filter(Term %in% b_neg_terms_ME),
                    'GO vxb pos' = enrichment_df_list$pos_vxb %>% filter(Term %in% vxb_pos_terms_ME),
                    'GO vxb neg' = enrichment_df_list$neg_vxb %>% filter(Term %in% vxb_neg_terms_ME) 
)
write.xlsx(output_xlsx, "~/Dropbox/aggregated_signed_enrichment_terms_mutually_exclusive.xlsx")
