library(patchwork)

source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))
source('~/Dropbox/time_series_fxn_manuscript.R')

time_series_df <- read.xlsx(paste0(Sys.getenv('DARPA'), 'Input/time_series.xlsx')) 
signatures_file <- '~/Dropbox/compendium_manuscript/tables/table_1a.xlsx'
signatures <- read.xlsx(signatures_file, sheet = 2)

# load imputed sex list object
if(sum(grepl(pattern = 'imputed_sex_list', x = ls())) != 1){
  imputed_sex_list <- readRDS('~/Documents/darpa-manuscript-data/data_list_final/imputed_sex_data_list.RDS')
}

filterSex <- function(MIobj, sex = 'male'){
  keep_ix <- MIobj$pheno$Imputed.Sex == sex
  MIobj %>%
    filterMIobj(index = keep_ix) %>%
    return()
}

# compute AUCs for male and female separately
male_list <- imputed_sex_list %>%
  lapply(filterSex, 'sex' = 'male')
#saveRDS(object = male_list, file = '~/Documents/darpa-manuscript-data/imputed_male_list.RDS')

female_list <- imputed_sex_list %>%
  lapply(filterSex, 'sex' = 'female')
#saveRDS(object = female_list, file = '~/Documents/darpa-manuscript-data/imputed_female_list.RDS')
#female_list <- readRDS('~/Documents/darpa-manuscript-data/imputed_female_list.RDS')

rm(imputed_sex_list)

# make sure to record N along the way
calculateSensitivity <- function(data_list, out_file, data_dir = '~/Documents/darpa-manuscript-data/'){
  # read in data / phenoData / time series annotations
  
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
    sig <- sig2Meta(signatures[signature_ID,], format_option = 'v1')
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
    
    current_sig <- sig#2Meta(signatures[signature_ID,])
    df <- plotDatSignature(signature = current_sig, 
                           data_list = tmp, 
                           performance_metric = metric, 
                           method = method, 
                           specific_method = specificMethod, pval_calc = T, parallel = T)
    
    # evaluate time-series data
    ts <- plotDatSignatureTimeSeries(signature = current_sig, 
                                     data_list = time_tmp, 
                                     signatures = signatures, 
                                     signature_ID = signature_ID,
                                     pval_calc = F)
    ts <- ts %>%
      mutate('TS' = T)
    df <- df %>%
      mutate('TS' = F)
    df <- bind_rows(df, ts)
    
    out[[signature_ID]] <- df
  }
  
  saveRDS(out_file, object = out)
  return(out)
}

male <- calculateSensitivity(data_list = male_list, out_file = '~/Dropbox/male_AUCs.RDS')
female <- calculateSensitivity(data_list = female_list, out_file = '~/Dropbox/female_AUCs.RDS')

male <- readRDS('~/Dropbox/male_AUCs.RDS')
female <- readRDS('~/Dropbox/female_AUCs.RDS')
# plot male vs female along with correlations?

male_df <- bind_rows(male) %>%
  #mutate('Imputed.Sex' = 'male') %>%
  dplyr::rename('Male.Scores' = Scores, 'Male.Pval' = P.Value, 'Male.N' = N)
female_df <- bind_rows(female) %>%
  #mutate('Imputed.Sex' = 'female') %>%
  dplyr::rename('Female.Scores' = Scores, 'Female.Pval' = P.Value, 'Female.N' = N)

plot_df <- male_df %>%
  left_join(female_df) %>%
  mutate('Comparison' = sapply(X = Signature, FUN = function(x){
    strsplit(x, split = '\n')[[1]][2]
  })) %>%
  mutate(Comparison = ifelse(grepl(pattern = 'Virus vs. Non-Infectious; Bacteria', x = Comparison), 'Virus vs. Bacteria', Comparison))  %>%
  mutate(Comparison = ifelse(grepl(pattern = 'Bacteria vs. H|Bacteria vs. N', x = Comparison), 'Bacteria', Comparison)) %>%
  mutate(Comparison = ifelse(grepl(pattern = 'Virus vs. H|Virus vs. N', x = Comparison), 'Virus', Comparison)) 


comparisons <- c("Virus vs. Bacteria", "Virus", "Bacteria")

head(plot_df)
sum(is.na(plot_df$Male.Scores))
sum(is.na(plot_df$Female.Scores))

label_dat <- plot_df %>%
  group_by(Signature, Comparison) %>%
  summarize('corr' = round(cor.test(x = Male.Scores, y = Female.Scores)$estimate,3),
            'pval' = round(cor.test(x = Male.Scores, y = Female.Scores)$p.value,3))

plot_list <- list()
for (comparison in comparisons){
  plot_list[[comparison]] <- 
    ggplot(plot_df %>% 
           filter(Comparison %in% comparison), 
         aes(x = Male.Scores, y = Female.Scores, color = Type)) + 
    geom_point() + 
    facet_grid(Comparison ~ Signature) + 
    xlim(c(0,1)) + 
    ylim(c(0,1)) + 
    geom_abline(slope = 1, intercept = 0, col = 'red') + 
    geom_text(data = label_dat %>% 
                filter(Comparison %in% comparison), 
              aes(x = 0.5, y = 0.2, color = NULL,
                  label = paste0('r = ', corr,'\np.value = ', pval))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}
plot_list[[1]]
plot_list[[2]]
plot_list[[3]]

layout <- "
AAAAA#####
BBBBBBBBBB
CCCCCCC###
"
p4 <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_layout(design = layout)
ggsave(p4, filename = '~/Dropbox/sex_supplement.png', width = 18, height = 9, units = 'in')
