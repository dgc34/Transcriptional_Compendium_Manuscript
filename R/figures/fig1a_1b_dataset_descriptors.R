library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MetaIntegrator)
library(RColorBrewer)
library(networkD3)
library(tibble)

source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))
source('~/Dropbox/plot_palette.R')

setwd('~/Dropbox/compendium_manuscript/figures/')

### Setting up colors 
cols <- brewer.pal(8, 'Set3')
plot(1:length(cols), rep(1, length(cols)), pch = 16, col = cols, cex = 5)
color_df <- data.frame('Class' = c('Healthy', 
                                   'Bacteria', 
                                   'Virus', 
                                   'Other.Infectious', 'Other Infectious', 
                                   'Other.NonInfectious', 'Other Non-Infectious', 'Other NonInfectious', 'Non-Infectious', 
                                   'Multiple', 'Multiple Classes', 'Bact+Vir'), 
                       'Color' = c(cols[7], 
                                   cols[5], 
                                   cols[4], 
                                   cols[1], cols[1], 
                                   cols[6], cols[6], cols[6], cols[6], 
                                   cols[3], cols[3], cols[3]))
class_colors <- as.character(color_df$Color)
names(class_colors) <- color_df$Class

ggplot(color_df, aes(x = Class, fill = Class)) + 
  scale_fill_manual(values = class_colors) + 
  geom_bar() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90))


### Panel A: 21 signatures in 3 types
sig_nums <- read.xlsx('~/Dropbox/compendium_manuscript/tables/supp_1.xlsx') %>%
  group_by(Type) %>%
  summarize('N' = n()) %>%
  mutate(Type = gsub(pattern = 'VxB', replacement = 'Virus vs. Bacteria\n(VxB)', x = Type)) %>%
  filter(Type != 'Influenza') %>%
  arrange(desc(N)) %>%
  mutate(Type = factor(Type, levels = Type, ordered = T))
p1a <- ggplot(sig_nums, aes(x = Type, y = N)) + 
  geom_bar(stat = 'identity') + 
  theme_bw() +
  theme_clean + 
  labs(x = '', y = 'Number of Signatures') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 20)) +
  scale_y_continuous(breaks = 0:11)
p1a 
ggsave(filename = '1a_signature_type_barplot.png', height = 8, width = 3.5, dpi = 300)
### Panel C: Number of subjects in each class
# get pheno data of all subjects
full_pheno <- read.csv('~/Dropbox/full_pheno_dlv1.csv') 

noninf_dl <- readRDS('~/Documents/darpa-manuscript-data/noninfectious_data_list.RDS')
noninf_df <- lapply(noninf_dl, function(x){
  df <- x$pheno %>%
    dplyr::select(geo_accession, Class, Standardized.Exposure) %>%
    mutate('Study' = x$formattedName) %>%
    dplyr::rename(geo = geo_accession, Pathogen = Standardized.Exposure) 
  return(df)
}) %>% bind_rows() %>% filter(Pathogen %in% c('Aging', 'Obesity'))

full_pheno <- bind_rows(full_pheno, noninf_df) %>%
  mutate(Class = gsub(Class, pattern = 'Other.NonInfectious', replacement = 'Non-Infectious')) %>%
  mutate(Class = gsub(Class, pattern = 'Other.Infectious', replacement = 'Other Infectious'))

# count classes 
class_orders <- c('Virus', 'Bacteria', 'Other Infectious', 'Multiple', 'Non-Infectious', 'Healthy')
plot_dat <- table(full_pheno$Class) %>%
  reshape2::melt('Class', value.name = 'Count') %>%
  mutate('Class' = gsub(pattern = 'Convalescent', replacement = 'Healthy', x = Class)) %>%
  mutate('Class' = gsub(pattern = 'Bact\\+Vir', replacement = 'Multiple', x = Class)) %>%
  mutate('Class' = gsub(pattern = '\\.', replacement = ' ', x = Class)) %>%
  mutate('Class' = gsub(pattern = 'Other NonInfectious', replacement = 'Non-Infectious', x = Class)) %>%
  filter(!grepl(pattern = '+', fixed = T, x = Class)) %>%
  group_by(Class) %>%
  summarize(Count = sum(Count)) %>%
  ungroup() %>%
  mutate(Class = factor(Class, levels = class_orders, ordered = T))

ggplot(plot_dat, aes(x = Class, y = Count, fill = Class)) + 
  theme_bw() + 
  theme_clean +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(hjust = 1, angle = 90), 
        text = element_text(size = 20),
        legend.position = 'none') + 
  scale_fill_manual(values = class_colors) + 
  scale_y_continuous(breaks = seq(to = 5500, from = 0, by = 1000)) +
  labs(x = element_blank(), y = 'Number of Subjects')
ggsave(filename = '1c_subjectN_barplot_exposures.png', height = 8, width = 7, dpi = 300)

### Panel D: sankey diagram of pathogens / classes

paths_to_plot <- full_pheno %>%
  filter(Class %in% c('Virus', 'Bacteria', 'Other Infectious', 'Non-Infectious')) %>%
  separate_rows(Pathogen, sep = ';') %>%
  filter(!Pathogen == '') %>%
  mutate(Pathogen = trimws(Pathogen)) %>%
   mutate(Pathogen = ifelse(str_detect(pattern = 'Acinetobacter', string = Pathogen), 'Acinetobacter', Pathogen)) %>%
   mutate(Pathogen = ifelse(str_detect(pattern = 'Enterobacter', string = Pathogen), 'Enterobacter', Pathogen)) %>%
  mutate(Pathogen = ifelse(str_detect(pattern = 'Enterococcus', string = Pathogen), 'Enterococcus', Pathogen)) %>%
  # mutate(Pathogen = ifelse(str_detect(pattern = 'Viridans', string = Pathogen), 'Viridans', Pathogen)) %>%
  # mutate(Pathogen = ifelse(str_detect(pattern = 'Streptococcus', string = Pathogen), 'Streptococcus', Pathogen)) %>%
  # mutate(Pathogen = ifelse(str_detect(pattern = 'Serratia', string = Pathogen), 'Serratia', Pathogen)) %>%
   mutate(Pathogen = ifelse(str_detect(pattern = 'Salmonella', string = Pathogen), 'Salmonella', Pathogen)) %>%
   mutate(Pathogen = ifelse(str_detect(pattern = 'Neisseria', string = Pathogen), 'Neisseria', Pathogen)) %>%
  # mutate(Pathogen = ifelse(str_detect(pattern = 'Chlamydia', string = Pathogen), 'Chlamydia', Pathogen)) %>%
  mutate(Pathogen = ifelse(str_detect(pattern = 'Plasmodium', string = Pathogen), 'Plasmodium', Pathogen)) %>%
  mutate(Pathogen = ifelse(str_detect(pattern = 'Staphylococcus', string = Pathogen), 'Staphylococcus', Pathogen)) %>%
  mutate(Pathogen = ifelse(str_detect(pattern = 'Streptococcus', string = Pathogen), 'Streptococcus', Pathogen)) %>%
  mutate(Pathogen = ifelse(str_detect(pattern = 'Klebsiella', string = Pathogen), 'Klebsiella', Pathogen)) %>%
  mutate(Pathogen = ifelse(str_detect(pattern = 'Aeromonas', string = Pathogen), 'Aeromonas', Pathogen)) %>%
  mutate(Pathogen = ifelse(str_detect(pattern = 'Aeromonas$|Klebsiella$|Acinetobacter$|Enterobacter$|Enterococcus$|Salmonella$|Neisseria$|Plasmodium$|Streptococcus$|Staphylococcus$', string = Pathogen), paste0(Pathogen, ' spp.'), Pathogen))

paths_to_plot2 <- paths_to_plot %>%
  filter(Class == 'Bacteria')

table(paths_to_plot2$Pathogen, paths_to_plot2$Class)

# input to sankey: 
## links data frame: source, target, value IDsource, IDTarget
### source and IDsource identify the left node
### target and IDtarget identify the right node
### value is the number of elements to move left to right

total = length(unique(paths_to_plot$Study))
unique_paths <- paths_to_plot %>%
  group_by(Pathogen) %>%
  filter(Pathogen != 'Unknown') %>%
  mutate('value' = n()) %>% # get the number of subjects
  ungroup() %>%
  mutate(Pathogen = ifelse(value <= 10, paste0('Other ', Class), Pathogen)) %>% # filter to pathogens with 5 or more samples
  group_by(Pathogen, Class) %>%
  summarize('value' = length(unique(Study)))  # count the number of unique studies with each label

unique_classes <- data.frame('Class' = c('Virus', 'Bacteria', 'Other Infectious', 'Non-Infectious'),
                             'Type' = c(rep('Infectious', 3), 'Non-Infectious'))
unique_classes <- unique_paths %>%
  left_join(unique_classes) %>%
  group_by(Class, Type) %>%
  summarize('value' = sum(value)) 
  

unique_types <- unique_classes %>%
  group_by(Type) %>%
  summarize('value' = sum(value)) %>%
  mutate('source' = 'Total') %>%
  dplyr::select(Type,source, value)

col_names <- c('target', 'source', 'value')
names(unique_paths) <- col_names
names(unique_classes) <- col_names
names(unique_types) <- col_names

link_df <- bind_rows(list(unique_paths, unique_classes, unique_types)) %>%
  arrange(desc(value)) %>%
   mutate_at(vars(matches('target', 'source')), function(x){
     gsub(pattern = 'Other Other', replacement = 'Other ', fixed = T, x = x)
     }) %>%
  filter(target != source)

term_ID_df <- data.frame('source' = unique(c(link_df$source, link_df$target))) %>%
  mutate('ID' = 1:length(unique(source))-1)

link_df <- link_df %>%
  left_join((term_ID_df %>% dplyr::rename('IDsource' = ID))) %>%
  left_join((term_ID_df %>% dplyr::rename('target' = source, 'IDtarget' = ID))) %>%
  dplyr::relocate(source, target) 

nodes %>% filter(grepl(pattern = '^Virus$|immunodef|tuberculosis', x = name))
link_df2 <- link_df %>%
  filter(!(IDsource == 2 & IDtarget == 7))

nodes <- term_ID_df %>%
  dplyr::rename('name' = source)

# Make the Network
p <- sankeyNetwork(Links = link_df2 %>% arrange(desc(IDtarget)), Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", fontFamily = 'Arial', fontSize = 12,
                   sinksRight=F)

p

  
# sankey labels
level_1 <- c('Infectious', 'Non-Infectious')
level_2 <- c('Bacteria', 'Virus', 'Parasite')
level_3 <- unique(paths_to_plot$Pathogen)

# sankey diagrams at study level — take all standardized labels from annotated data frame
# manually annotate links, regenerate diagram
input_file <- '~/Dropbox/pathogen_table_annotated.xlsx'
sheet_ix <- grep(pattern = 'pathogen_table_annotated', x = getSheetNames(input_file))
std_exp_df <- read.xlsx(input_file, sheet = sheet_ix)
path_replacements <- read.xlsx('~/Dropbox/ncbi_taxonomy_names.xlsx') %>%
  dplyr::rename(Standardized.Exposure = Standardized.Pathogen)

study_names <- unique(full_pheno$Study) %>% parseAccession()
filtered_df <- std_exp_df %>%
  filter(Accession %in% study_names) %>%
  mutate(Standardized.Exposure = gsub(x = Standardized.Exposure, pattern = '|', replacement = ';', fixed = T)) %>%
  separate_rows(Standardized.Exposure, sep = ';') %>%
  mutate(Standardized.Exposure = trimws(Standardized.Exposure)) %>%
  dplyr::select(Accession, Platform, Exposure, Standardized.Exposure) %>%
  distinct()

filtered_df <- filtered_df %>%
  left_join(path_replacements) %>%
  dplyr::select(-Standardized.Exposure) %>%
  dplyr::rename(Standardized.Exposure = NCBI.Taxonomy)

# level 3 numbers
terms_to_omit <- c('Bacteria', 'Virus', 'Sepsis', 'Pneumonia', 'Upper Respiratory Infection', 'Meningitis')
terms <- filtered_df %>%
  filter(!Standardized.Exposure %in% terms_to_omit) %>%
  pull(Standardized.Exposure) %>%
  table() %>%
  sort(decreasing = T) 

single_pathogens <- terms#[terms > 3]



# level 2 numbers
class_label_df <- std_exp_df %>%
  filter(Accession %in% study_names) %>%
  mutate(Class = gsub(Class, pattern = '+', replacement = ';', fixed = T)) %>%
  separate_rows(Class, sep = ';') %>%
  mutate(Class = trimws(Class)) %>% 
  filter(Class %in% level_2) %>%
  dplyr::select(Accession, Platform, Class, Exposure, Standardized.Exposure) %>%
  distinct() 
class_labels <- class_label_df %>%
  pull(Class) %>%
  table() %>%
  sort(decreasing = T)

# level 1 numbers
inf_terms <- c('Bact', 'Vir', 'Par')
noninf_terms <- c('Non')
index_label_df <- std_exp_df %>%
  filter(Accession %in% study_names) %>%
  mutate(Class = gsub(Class, pattern = '+', replacement = ';', fixed = T)) %>%
  separate_rows(Class, sep = ';') %>%
  mutate(Class = trimws(Class)) %>% 
  filter(grepl(pattern = paste(c(inf_terms, noninf_terms), collapse = '|'), x = Class)) %>%
  dplyr::select(Accession, Platform, Class, Exposure, Standardized.Exposure) %>%
  distinct() 
index_labels <- index_label_df %>%
  mutate(Index = ifelse(grepl(pattern = paste(inf_terms, collapse = '|'), x = Class), 'Infectious', 'Non-Infectious')) %>%
  pull(Index) %>%
  table() %>%
  sort(decreasing = T)

# level 0 numbers
N0 <- std_exp_df %>%
  dplyr::select(Accession, Platform) %>%
  filter(Accession %in% study_names) %>%
  nrow()
names(N0) <- 'Total'

#sankey numbers
level_0N <- N0
level_1N <- index_labels
level_2N <- class_labels
level_3N <- single_pathogens

all_terms <- c(names(c(level_0N, level_1N, level_2N, level_3N)))
input_matrix <- matrix(data = 0, nrow = length(all_terms), ncol = length(all_terms))
rownames(input_matrix) <- colnames(input_matrix) <- all_terms

input_node <- 'Total'
for (i in 1:length(level_1N)){
  print(i)
  current_output_node <- level_1N[i]
  #print(current_output_node)
  input_matrix[input_node, names(current_output_node)] <- as.numeric(current_output_node)
}

input_node <- 'Infectious'
for (i in 1:length(level_2N)){
  print(i)
  current_output_node <- level_2N[i]
  #print(current_output_node)
  input_matrix[input_node, names(current_output_node)] <- as.numeric(current_output_node)
}

level_3N
viruses <- grep(pattern = paste(c('^Influenza', 'Virus', 'Hep',
                                  'Hepatitis', 'Measles'), collapse = '|'), x = names(level_3N), ignore.case = T, value = T)
bacteria <- grep(pattern = paste(c('bact', 'staph', 'strep', 'coli', 'influenzae', 'salmonella', 'burkholderia', 'mycoplasma'), collapse = '|'), x = names(level_3N), ignore.case = T, value = T)
parasites <- grep(pattern = paste(c('Leish', 'Malaria|Plasmodium', 'Hookworm'), collapse= '|'), x = names(level_3N), ignore.case = T, value = T)
                  
input_node <- 'Virus'
for (i in 1:length(viruses)){
  print(i)
  current_output_node <- level_3N[names(level_3N) %in% viruses][i]
  #print(current_output_node)
  input_matrix[input_node, names(current_output_node)] <- as.numeric(current_output_node)
}

input_node <- 'Bacteria'
for (i in 1:length(bacteria)){
  print(i)
  current_output_node <- level_3N[names(level_3N) %in% bacteria][i]
  #print(current_output_node)
  input_matrix[input_node, names(current_output_node)] <- as.numeric(current_output_node)
}

input_node <- 'Parasite'
for (i in 1:length(parasites)){
  print(i)
  current_output_node <- level_3N[names(level_3N) %in% parasites][i]
  #print(current_output_node)
  input_matrix[input_node, names(current_output_node)] <- as.numeric(current_output_node)
}

colSums(input_matrix)
links <- input_matrix %>% 
  as.data.frame() %>% 
  rownames_to_column(var="source") %>% 
  gather(key="target", value="value", -1) %>%
  filter(value != 0)
nodes <- data.frame(
  name=c(as.character(links$source), as.character(links$target)) %>% 
    unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1


