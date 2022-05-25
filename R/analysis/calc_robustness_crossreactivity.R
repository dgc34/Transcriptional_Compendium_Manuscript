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

# this script calculates generalizability and specificity for bacterial/viral infection data, 
# as well as aging, obesity, and parasite datasets

# specify generalizability or specificity
TASK <- 'generalizability'
if(!TASK %in% c('generalizability', 'specificity')){
  stop('TASK must be either generalizability or specificity')
}
print(TASK)

# initialization
data_dir <-  '../../data_list_final/'
data_file <- 'data_list_v1.RDS'

source(paste0('../analysis/', 'helperFunctions.R'))

# read in data / phenoData / time series annotations
time_series_df <- read.xlsx(paste0(Sys.getenv('DARPA'), 'Input/time_series.xlsx')) 
signatures_file <- '~/Dropbox/compendium_manuscript/tables/supp_1.xlsx'
control_type <- str_extract('healthy|noninf', string = signatures_file)
signatures <- read.xlsx(signatures_file)

data_list <- readRDS(paste0(data_dir, data_file))
# output file label depends on input file string
if(grepl(pattern = 'hiv|bp|male',  x = data_file)){
  data_label <- str_extract(pattern = 'hiv|bp|female|male', string = data_file)
} else {
  data_label <- 'big'
}

# remove co-infections and immunomodulatory drugs
data_list <- lapply(data_list, function(x){
  ix <- x$pheno$Class %in% c('Bact+Vir', 'Other.Infectious;Virus', "Vir+Other.Infectious") 
  x <- filterMIobj(x, !ix)
})

# format signatures for evaluation
signatures <- signatures %>%
  mutate(Positive.Class = str_extract(pattern = '.*vs.', string = Comparison)) %>% 
  mutate(Negative.Class = str_extract(pattern = 'vs.*', string = Comparison)) %>%
  mutate(Positive.Class = gsub(pattern = ' vs.', replacement = '', x = Positive.Class)) %>% 
  mutate(Negative.Class = gsub(pattern = 'vs. ', replacement = '', x = Negative.Class)) %>%
  mutate(Positive.Class = ifelse(Positive.Class == 'Influenza', 'Virus', Positive.Class)) # this will let us compute bacterial specificity for influenza signatures

# signature positive/negative class labels are modified for cross-reactivity (aka specificity) calculations
if(TASK == 'specificity'){
  signatures <- specificityConversion(signatures)
}

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
  # cross-reactivity not calculated for V/B signatures
  if(TASK == 'specificity' & signatures[signature_ID, 'Type'] == 'VxB'){
    next
  }
  
  # cross-reactivity not calculated for flu vs. bact signatures
  if(TASK == 'specificity' & signatures[signature_ID, 'Comparison'] == 'Influenza vs. Bacteria'){
    next
  }
  
  # only consider v, b, v/b, i signatures
  if(!signatures[signature_ID,'Type'] %in% c("Bacteria", "Virus", "VxB", "Influenza")){
  	next
  } else {
    
  }

  # find datasets with classes matching signature
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
  
  # remove studies that don't have 50% of  signature genes
  sig <- sig2Meta(signatures[signature_ID,])
  ix <- sapply(X = tmp, filterGenes, signature = sig, threshold = 0.5)
  tmp <- tmp[ix]
  if(length(tmp) == 0){
    next
  }

  # identify time-series study designs
  time_ix <- sapply(names(tmp), parseAccession) %in% time_series_df$Accession
  time_tmp <- tmp[time_ix]
  tmp <- tmp[!time_ix]

  # evaluate non-time-series data
  sapply(tmp, function(x){
    print(x$formattedName);
    return(checkDataObject(x, objectType = 'Dataset'))
  })

  # parse signature into MetaIntegrator filterObject
  current_sig <- sig2Meta(signatures[signature_ID,])
  
  # calculate AUROCs for cross-sectional studies
  df <- plotDatSignature(signature = current_sig, 
  		data_list = tmp, 
  		performance_metric = metric, 
  		method = method, 
  		specific_method = specificMethod, pval_calc = F, parallel = T)

  df <- df %>%
    mutate('TS' = F)
  
  # evaluate time-series data
  if(length(time_tmp) > 0){
    ts <- plotDatSignatureTimeSeriesAltPval(signature = current_sig, 
                                               data_list = time_tmp, 
                                               #signatures = signatures, 
                                               #signature_ID = signature_ID, 
                                               pval_calc = F)
    ts <- ts %>%
      mutate('TS' = T)
    
    # merge
    df <- bind_rows(df, ts)
  }
  
  # save to output list
  out[[signature_ID]] <- df
}
plot_dat <- bind_rows(out)

date <- format(Sys.Date(), '%m%d%y') 
saveRDS(paste0('~/Dropbox/', data_label, '_eval_heatmap_list_', tolower(TASK), '_', date, '_', control_type, '_', 'sbj_ts.RDS'), object = out)

