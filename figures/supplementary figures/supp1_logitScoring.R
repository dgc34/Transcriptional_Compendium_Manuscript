library(openxlsx)
library(ggplot2)
library(dplyr)
library(MetaIntegrator)
library(stringr)
library(ggh4x)
library(tidyr)
library(patchwork)

source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))

# visualize signature performance in geomMean approach vs logit scoring
# 0 load data
# 1 score both ways
# 2 generate plot

# LOAD DATA & SIGNATURES
data_list <- readRDS('~/Documents/darpa-manuscript-data/data_list_final/data_list_v1.RDS')
data_list <- lapply(data_list, function(x){
  ix <- x$pheno$Class %in% c('Bact+Vir', 'Other.Infectious;Virus') # remove co-infections and immunomodulatory drugs
  x <- filterMIobj(x, !ix)
})

# remove time series studies from evaluation
time_series_df <- read.xlsx(paste0(Sys.getenv('DARPA'), 'Input/time_series.xlsx')) 
time_series_df <- time_series_df %>%
  mutate(Series = as.logical(Series)) %>%
  filter(Series)

time_ix <- sapply(names(data_list), parseAccession) %in% time_series_df$Accession
data_list <- data_list[!time_ix]

# read in signatures
signatures_file <- '~/Dropbox/compendium_manuscript/tables/table_1a.xlsx'
signatures <- read.xlsx(signatures_file, sheet = 2)

# settings for figure generation
generateTable <- function(data_list, signatures, metric, method, specificMethod){
  # generate plotDat table
  out <- list()
  for (signature_ID in 1:nrow(signatures)){
    
    if(!signatures[signature_ID,'Positive.Class'] %in% c("Bacteria", "Virus", "Non-Infectious", "Resp.Virus", "Resp.Bact")){
      print('skipping')
      next
    } else {
      signatures[signature_ID,'Positive.Class'] <- gsub(pattern = 'Resp.Virus', replacement = 'Virus', x = signatures[signature_ID,'Positive.Class'])
      signatures[signature_ID,'Positive.Class'] <- gsub(pattern = 'Resp.Bact', replacement = 'Bacteria', x = signatures[signature_ID,'Positive.Class'])
    }
    
    tmp <- lapply(X = data_list, function(x){
      x$pheno$Class[x$pheno$Class == 'Convalescent'] <- 'Healthy' # treat recovered subjects as healthy
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
    ix <- sapply(X = tmp, filterGenes, signature = sig)
    tmp <- tmp[ix]
    if(length(tmp) == 0){
      next
    }
    
    
    # evaluate non-time-series data
    sapply(tmp, function(x){
      print(x$formattedName);
      return(checkDataObject(x, objectType = 'Dataset'))
    })
    
    current_sig <- sig2Meta(signatures[signature_ID,])
    print(current_sig$filterDescription)
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
  
  if(nrow(plot_dat) == 0){
    return(NULL)
  } else {
    plot_dat <- plot_dat %>%
      mutate('Type' = factor(Type, levels = c('Discovery', 'Validation', 'New', 'Null'), ordered = T)) %>%
      mutate('Comparison' = sapply(X = Signature, FUN = function(x){
        strsplit(x, split = '\n')[[1]][2]
      }))  %>%
      arrange(Comparison) %>%
      mutate(Signature = factor(Signature, levels = unique(Signature), ordered = T))
    
    plot_dat <- plot_dat %>%
      filter(!grepl(pattern = 'NA|DARPA', x = Signature)) %>% # remove the single gene testing variants and DARPA signatures
      mutate('Exposure.Class' = sapply(Signature, str_extract, pattern = '\n.*vs')) %>%
      mutate(Exposure.Class = gsub(x = Exposure.Class, pattern = '\\n| |vs', replacement = '')) %>%
      mutate('Author' = str_extract(pattern = '^.* [:digit:]+', string = Signature)) %>%
      #mutate(Author = str_replace(pattern = '', replacement = '', string = Author)) %>%
      mutate(Author = trimws(Author)) %>%
      mutate('Method' = paste0(method, '_', specificMethod))
    return(plot_dat)
  }
  
}
date <- format(Sys.Date(), '%m%d%y') 
output_dir <- paste0('~/Dropbox/figS1/', date, '/')
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

metric <- 'AUROC'
method <- 'zScore'
specificMethod = 'geomMean'
#p1 <- generateTable(data_list, signatures, metric, method, specificMethod)
#saveRDS(object = p1, file = paste0(output_dir, 'p1.RDS'))
for (i in nrow(signatures):1){
  print(i)
  par_output_file <- paste0(output_dir, 'p2_', i, '.RDS')
  if(!file.exists(par_output_file)){
    p2_par <- generateTable(data_list, signatures[i,], metric = 'AUROC', method = 'logit', specificMethod = NULL)
    saveRDS(object = p2_par, file = par_output_file)
  } else {
    print(paste0(par_output_file, ' already exists!'))
  }
  
}
