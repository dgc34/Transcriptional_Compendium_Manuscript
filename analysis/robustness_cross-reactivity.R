# the goal of this figure is to evaluate all signatures in all appropriate datasets
library(dplyr)
library(MetaIntegrator)
library(ggplot2)
library(openxlsx)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(aveytoolkit)
library(singscore)
library(pROC)
library(ggh4x)
library(scales)

# initialization
data_dir <-  '~/Documents/darpa-manuscript-data/'
source(paste0(Sys.getenv("DARPA"), 'helperFunctions.R'))
source(paste0(Sys.getenv('DARPA'), 'Analysis/scoreStudyTimeSeriesFxn_2021.R'))
source('~/Dropbox/FUN_heatmap.R')

# read in data / phenoData / time series annotations
time_series_df <- read.xlsx(paste0(Sys.getenv('DARPA'), 'Input/time_series.xlsx')) 
signatures_file <- paste0(Sys.getenv("DARPA"), '/Input/signatures/bacterial_viral_healthy.xlsx')
signatures <- read.xlsx(signatures_file)
data_list <- readRDS(paste0(data_dir, 'data_list.RDS'))

# need to identify time series studies for secparate processing
time_series_df <- time_series_df %>%
	mutate(Series = as.logical(Series)) %>%
	filter(Series)

# settings for figure generation
metric <- 'AUROC'
method <- 'zScore'
specificMethod = 'geomMean'

# generate plotDat table
out <- list()
for (signature_ID in 1:nrow(signatures)){
  
  if(!signatures[signature_ID,'Positive.Class'] %in% c("Bacteria", "Virus", "Non-Infectious", "Resp.Virus", "Resp.Bact")){
  	next
  } else {
  	signatures[signature_ID,'Positive.Class'] <- gsub(pattern = 'Resp.Virus', replacement = 'Virus', x = signatures[signature_ID,'Positive.Class'])
  	signatures[signature_ID,'Positive.Class'] <- gsub(pattern = 'Resp.Bact', replacement = 'Bacteria', x = signatures[signature_ID,'Positive.Class'])
  }

  tmp <- lapply(X = data_list, function(x){
    parseClass(x, signature_row = signatures[signature_ID,])
  })

  # remove studies that don't have multiple classes
  ix <- sapply(X = tmp, filterClass)
  tmp <- tmp[ix]
  if(length(tmp) == 0){
    next
  }
  # remove studies that don't have all signature genes
  sig <- sig2Meta(signatures[signature_ID,])
  ix <- sapply(X = tmp, filterGenes, signature = sig, threshold = 0.5)
  tmp <- tmp[ix]
  if(length(tmp) == 0){
    next
  }

  time_ix <- sapply(names(tmp), parseAccession) %in% time_series_df$Accession
  time_tmp <- tmp[time_ix]
  tmp <- tmp[!time_ix]

  # evaluate non-time-series data
  sapply(tmp, function(x){
    print(x$formattedName);
    return(checkDataObject(x, objectType = 'Dataset'))
  })

  current_sig <- sig2Meta(signatures[signature_ID,])
  df <- plotDatSignature(signature = current_sig, 
  		data_list = tmp, 
  		performance_metric = metric, 
  		method = method, 
  		specific_method = specificMethod, pval_calc = T, parallel = T)

  # evaluate time-series data
  ts <- summarizeTimeSeriesPlotDat(signature = current_sig, data_list = time_tmp, pval_calc = T)
  ts <- ts %>%
  	mutate('TS' = T)
  df <- df %>%
  	mutate('TS' = F)
  df <- bind_rows(df, ts)

  out[[signature_ID]] <- df
}
plot_dat <- bind_rows(out)
saveRDS('~/Documents/darpa-manuscript-data/figures/big_eval_heatmap_list.RDS', object = out)
#plot_dat <- plot_dat %>%
#  mutate('P.Value' = ifelse(is.na(P.Value), pval, P.Value)) %>%
#  dplyr::select(-pval)

plot_dat <- plot_dat %>%
  mutate('Type' = factor(Type, levels = c('Discovery', 'Validation', 'New', 'Null'), ordered = T)) %>%
  mutate('Comparison' = sapply(X = Signature, FUN = function(x){
    strsplit(x, split = '\n')[[1]][2]
  }))  %>%
  arrange(Comparison) %>%
  mutate(Signature = factor(Signature, levels = unique(Signature), ordered = T))

plot_dat <- plot_dat %>%
	filter(!grepl(pattern = 'NA|DARPA|Chawla', x = Signature)) %>% # remove the single gene testing variants and DARPA signatures
	mutate('Study.Design' = ifelse(TS, 'Time Series', 'Acute Exposure')) %>%
	mutate('Study.Design' = factor(Study.Design)) %>%
	mutate('Exposure.Class' = sapply(Signature, str_extract, pattern = '\n.*vs')) %>%
	mutate(Exposure.Class = gsub(x = Exposure.Class, pattern = '\\n| |vs', replacement = '')) %>%
	mutate('Author' = str_extract(pattern = '^.* [:digit:]+', string = Signature)) %>%
	#mutate(Author = str_replace(pattern = '', replacement = '', string = Author)) %>%
	mutate(Author = trimws(Author)) %>%
  mutate(Comparison = str_extract(pattern = '\\n.*', string = Signature)) %>%
  mutate(Comparison = str_remove(pattern = '\\n', string = Comparison))
table(plot_dat$Comparison)

N_threshold = 5
colors <- createColorScale()

# facet by pathogen and cluster rows/columns
head(plot_dat)
table(plot_dat$Exposure.Class)
exposures <- unique(plot_dat$Exposure.Class) %>% grep(pattern = 'Non', invert = T, value = T)
for(exposure in exposures){
  input_tmp_plot_dat <- plot_dat %>%
    filter(N >= N_threshold, Study.Design == 'Acute Exposure', Type != 'Null', Exposure.Class == exposure) %>%
    mutate(Signature = as.character(Signature))
  
  # apply clustering on studies
  data_mat <- input_tmp_plot_dat %>%
    pivot_wider(names_from = 'Signature', values_from = 'Scores', id_cols = 'Study', values_fill = -10)
  
  dist <- dist(as.matrix(data_mat))
  orders <- hclust(dist)$order
  clust_dat <- data.frame('Study' = data_mat$Study, 'clust' = orders)
  
  # arrange studies by number of blank rows and clustering
  order_dat <- input_tmp_plot_dat %>%
    complete(Signature, nesting(Study)) %>%
    group_by(Study) %>%
    mutate('blank_row_sum' = sum(is.na(Scores))) %>%
    ungroup() %>%
    dplyr::select(Study, blank_row_sum) %>%
    distinct() %>%
    left_join(clust_dat) %>% 
    arrange(desc(blank_row_sum), clust)
  
  tmp_plot_dat <- input_tmp_plot_dat %>%
    mutate('Study' = factor(Study, levels = order_dat$Study, ordered = T))
  
  # apply clustering on signatures
  data_mat <- input_tmp_plot_dat %>%
    pivot_wider(names_from = 'Study', values_from = 'Scores', id_cols = 'Signature', values_fill = -10)
  dist <- dist(as.matrix(data_mat))
  orders <- hclust(dist)$order
  clust_dat <- data.frame('Signature' = data_mat$Signature, 'clust' = orders)
  order_dat <- input_tmp_plot_dat %>%
    complete(Study, nesting(Signature)) %>%
    group_by(Signature) %>%
    mutate('blank_row_sum' = sum(is.na(Scores))) %>%
    ungroup() %>%
    dplyr::select(Signature, blank_row_sum) %>%
    distinct() %>%
    left_join(clust_dat) %>% 
    arrange(blank_row_sum, desc(clust))
  
  tmp_plot_dat <- tmp_plot_dat %>%
    mutate('Signature' = factor(Signature, levels = order_dat$Signature, ordered = T))
  
  study_plot_dat <- tmp_plot_dat %>%
    mutate(Train_Label = ifelse(Type == 'Discovery', '○', '')) %>%
    mutate(Significance_Label = ifelse(P.Value <= 0.05, '*', '')) 
  
  
  p3 <- ggplot(study_plot_dat, aes(x = Signature, y = Study, fill = Scores)) + 
    geom_tile() + 
    scale_fill_gradient2(high = muted('red'), low = muted('blue'), mid = 'white', midpoint = 0.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    facet_nested(. ~ Exposure.Class, scales = 'free', space = 'free', labeller = label_value) + 
    geom_text(aes(label = Significance_Label), size = 6, vjust = 0.77) +
    geom_text(aes(label = Train_Label), size = 9, vjust = 0.37) 
  p3
  filename <- paste0('~/Dropbox/fig4_heatmap_', exposure, '.png')
  width_num <- length(unique(tmp_plot_dat$Signature))
  ggsave(p3 , filename = filename, dpi = 300, width = (4+0.5*width_num), height = 16)
}



# facet by comparison
## pro: non-sparse plots
## con: some signatures are in their own columns
head(plot_dat)
for(comparison in comparisons){
  tmp_plot_dat <- plot_dat %>%
    filter(N >= N_threshold, Study.Design == 'Acute Exposure', Type != 'Null', Comparison == comparison)
  p2 <- ggplot(tmp_plot_dat, aes(x = Signature, y = Study, fill = Scores)) + 
    geom_tile() + 
    scale_fill_gradient2(high = muted('red'), low = muted('blue'), mid = 'white', midpoint = 0.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    facet_nested(. ~ Exposure.Class, scales = 'free', space = 'free', labeller = label_value)
  
  filename <- paste0('~/Dropbox/fig4_heatmap_', comparison, '.png')
  width_num <- length(unique(tmp_plot_dat$Signature))
  ggsave(p2, filename = filename, dpi = 300, width = (4+0.5*width_num), height = 16)
}


# new version for each comparison, split into blocks based on the number of studies that 
# have evaluations calculated: 3 blocks (all 12, vir*, vir vs bact)
# sort rows by number of stars

# why are studies missing evaluations? e.g. GSE6269_GPL570
plot_dat %>% filter(grepl(pattern = 'GSE6269_GPL570', x = Study)) # because different comparisons, different exposures per dataset

# why are some viral vs. healthy studies missing evaluations in nearly all signatures?
plot_dat %>% filter(grepl(pattern = 'GSE81926|30310', x = Study), Type != 'Null')
# need to go back to debug plot_dat generation

# generate plotDat table
out <- list()
for (signature_ID in 1:nrow(signatures)){
  #signature_ID = 6
  if(!signatures[signature_ID,'Positive.Class'] %in% c("Bacteria", "Virus", "Non-Infectious", "Resp.Virus", "Resp.Bact")){
    next
  } else {
    signatures[signature_ID,'Positive.Class'] <- gsub(pattern = 'Resp.Virus', replacement = 'Virus', x = signatures[signature_ID,'Positive.Class'])
    signatures[signature_ID,'Positive.Class'] <- gsub(pattern = 'Resp.Bact', replacement = 'Bacteria', x = signatures[signature_ID,'Positive.Class'])
  }
  
  if(signatures[signature_ID,'Positive.Class'] == 'Virus' & signatures[signature_ID,'Negative.Class'] == 'Healthy'){
    print(signature_ID)
  }
  
  tmp <- data_list[c('GSE30310_GPL9392', 'GSE81926_GPL21947')]
  tmp <- lapply(X = tmp, function(x){
    parseClass(x, signature_row = signatures[signature_ID,])
  })
  
  # remove studies that don't have multiple classes
  ix <- sapply(X = tmp, filterClass)
  tmp <- tmp[ix]
  if(length(tmp) == 0){
    next
  }
  # remove studies that don't have all signature genes
  sig <- sig2Meta(signatures[signature_ID,])
  ix <- sapply(X = tmp, filterGenes, signature = sig, threshold = 0.5)
  
  if(sum(ix) == 0){
    sig_genes <- c(sig$posGeneNames, sig$negGeneNames)
    print(sum(tmp[[1]]$keys %in% sig_genes))
  }
  
  tmp <- tmp[ix]
  if(length(tmp) == 0){
    next
  }
  
  time_ix <- sapply(names(tmp), parseAccession) %in% time_series_df$Accession
  time_tmp <- tmp[time_ix]
  tmp <- tmp[!time_ix]
  
  # evaluate non-time-series data
  sapply(tmp, function(x){
    print(x$formattedName);
    return(checkDataObject(x, objectType = 'Dataset'))
  })
  
  current_sig <- sig2Meta(signatures[signature_ID,])
  df <- plotDatSignature(signature = current_sig, 
                         data_list = tmp, 
                         performance_metric = metric, 
                         method = method, 
                         specific_method = specificMethod)
  
  # shuffle class labels within each study, recompute signatures
  perm_data_list <- lapply(tmp, function(x){
    ix <- sample(1:length(x$class), replace = F)
    x$class <- x$class[ix]
    names(x$class) <- rownames(x$pheno)
    return(x)
  })
  perm <- plotDatSignature(signature = current_sig, data_list = perm_data_list, performance_metric = metric, method = method, specific_method = specificMethod)
  perm <- perm %>%
    filter(Type == 'New') %>%
    mutate('Type' = 'Null')
  df <- df %>%
    bind_rows(perm)
  
  out[[signature_ID]] <- df
}
plot_dat <- bind_rows(out)
