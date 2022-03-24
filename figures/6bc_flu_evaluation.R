library(dplyr)
library(ggplot2)
library(openxlsx)

source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))

input_dir <- '~/Documents/darpa-manuscript-data/flu_data_v3/'
setwd(input_dir)

### Edit contrast as needed
# input contrast
contrast <- 'healthy'

# healthy and noninf specify controls, original uses the author-supplied contrast
if(!contrast %in% c('healthy', 'noninf', 'original')){
  stop('contrast must be either "healthy", "noninf", or "original"')
}

# current contrast is a dummy variable for modifying the input table before parsing
current_contrast <- 'original'
if(contrast == 'healthy'){
  current_contrast <- 'Healthy'
}
if(contrast == 'noninf'){
  current_contrast <- 'Non-Infectious'
}
signatures <- read.xlsx('~/Dropbox/compendium_manuscript/tables/supp_1.xlsx') %>%
  filter(Type == 'Influenza' | Signature.Label == 'V10')

if(current_contrast != 'original'){
  signatures <- signatures %>%
    mutate(Comparison = current_contrast)
}

flu_test <- readRDS('flu_test.RDS')
other_files <- list.files(pattern = 'other')
other_files <- other_files[!grepl(pattern = 'other_data', x = other_files)]
other_list <- lapply(other_files, readRDS) %>%
  unlist(recursive = F)

out <- list()
for (i in 1:nrow(signatures)){
  current_sig <- sig2Meta(signatures[i,])
  flu_df <- plotDatSignature(signature = current_sig, 
                             data_list = flu_test,
                             pval_calc = F)
  other_df <- plotDatSignature(current_sig, 
                               data_list = other_list, 
                               pval_calc = F)
  
  flu_df <- flu_df %>%
    mutate(Pathogen = 'Influenza')
  other_df <- other_df %>%
    mutate(Pathogen = 'Non-Influenza Virus')
  
  plot_dat <- bind_rows(flu_df, other_df)
  out[[i]] <- plot_dat
}

saveRDS(object = out, file = paste0('~/Documents/darpa-manuscript-data/flu_eval_heatmap_v3_withV10_', contrast, '.RDS'))
