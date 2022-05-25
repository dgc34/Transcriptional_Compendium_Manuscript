library(openxlsx)
library(dplyr)
library(GEOquery)
library(parallel)
library(tidyr)
library(ggplot2)
library(limma)
library(stringr)
library(affy)
library(oligo)
library(magrittr)
library(Biobase)

download_raw_flag <- F
annotation_overwrite_flag <- F # this must be true if download_raw_flag is true

# pipeline
setwd(paste0(Sys.getenv('DARPA'), '/Input/'))
source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))
source(paste0(Sys.getenv('DARPA'), 'Input/platform_specific/functions/preProcessingFunctions.R'))

# pull GEO raw files
source('download_RAW.R')

sheet_num = 4

if (sheet_num == 1){
  out_dir <- '~/Data/DARPA/' 
  phenoData <- read.xlsx(paste0(Sys.getenv('DARPA'), 'Input/classLabels/phenoData_convalescent_infectiousOnly.xlsx'))
} else if (sheet_num == 2){
  out_dir <- '~/Data/DARPA/Specificity/Obesity/'
  phenoData <- read.csv('~/darpa-echo/Input/classLabels/Obesity_phenoData.csv', stringsAsFactors = F)
} else if (sheet_num == 3){
  out_dir <- '~/Data/DARPA/Specificity/Aging/'
  phenoData <- read.csv('~/darpa-echo/Input/classLabels/Aging_phenoData.csv', stringsAsFactors = F)
} else if (sheet_num == 4){
  out_dir <- '~/Data/DARPA/Specificity/COVID19/'
  phenoData <- read.xlsx('~/Dropbox/phenoData_C19_comorbidity.xlsx') %>%
    mutate('Include' = as.logical(Include)) %>%
    filter(Include)
}

phenoData <- phenoData %>%
  mutate('Include' = as.logical(Include))

GEO_dl_dir <- paste0(out_dir, 'GEO_dl/')
master_list <- read.xlsx(paste0(Sys.getenv('DARPA'), 'DARPA.xlsx'), sheet = sheet_num)

if(sheet_num == 4){
  master_list <- phenoData %>%
    mutate(Suitable = ifelse(Include, 'Yes', 'No'))
}

# move these to helperFunctions once debugging is complete
source('download.R')
source('sortPlatformsFromRAW.R')
source('updateAnnotationsFxn.R')

# set up directories for intermediate and processed raw data
if(!dir.exists(paste0(out_dir, 'GEO_intermediate/'))){
  dir.create(paste0(out_dir, 'GEO_intermediate/'))
}

if(!dir.exists(paste0(out_dir, 'GEO_processed/'))){
  dir.create(paste0(out_dir, 'GEO_processed/'))
}

# get a list of all platforms 
platforms <- list.files(path = paste0(out_dir, 'GEO_raw'))
platforms <- c(platforms, 'missing_raw')
platform_types <- read.xlsx(paste0(Sys.getenv('DARPA'), 'Input/platform_specific/platform_types.xlsx'))

if(download_raw_flag){
  # pull GEO exprsets w phenodata
  download(output_dir = paste0(out_dir, 'GEO_dl'), master_list)
  
  # move raw files to corresponding directories
  download_RAW(output_dir = out_dir, master_list)
  sortPlatformsFromRaw(out_dir, master_list)
  
  platforms <- trimws(list.files(path = paste0(out_dir, 'GEO_raw')))
  platform_types <- read.xlsx(paste0(Sys.getenv('DARPA'), 'Input/platform_specific/platform_types.xlsx'))
  
  for (platform in platforms){
    print(paste0('Platform: ', platform))
    accessions <- list.files(path = paste0(out_dir, 'GEO_raw/', platform))
    
    platform_type <- platform_types %>%
      dplyr::filter(Platform == platform) %>%
      dplyr::select(Manufacturer)
    
    data_dir <- paste0(out_dir)
    intermediate_dir <- paste0(data_dir, 'GEO_intermediate/')
    output_dir <- paste0(data_dir, 'GEO_processed/', platform, '/')
    rds_dir <- paste0(data_dir, 'GEO_dl/')
    setwd(data_dir)
    
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    
    for (accession in accessions){
      if(length(list.files(path = output_dir, pattern = accession)) > 0){
        print(paste0(accession, ': already processed'))
        next
      }
      
      raw_dir <- paste0(data_dir, 'GEO_raw/', platform, '/', accession, '/')
      rds_file <- list.files(path = rds_dir, pattern = accession)
      result_file <- NULL
      print(accession)
      
      ## run raw files through normalization
      
      # throw in a check for ilmn vs affy
      if (('ilmn' %in% platform_type)){
        if(accession %in% c('GSE41233', 'GSE35846', 'GSE65219')){
        #if(accession %in% c('GSE35846')){
          # GSE35846 is missing a column from its raw data
          # GSE65219 required heavy munging corrections, but still fails neqc
          next
        }
        print('ilmn')
        # load raw file and apply preprocessing
        raw_file <- list.files(path = raw_dir, pattern = "txt$")
        print('pre-process')
        preProcessILMN(raw_file, intermediate_dir)
        # read intermediate preprocessed file, apply expr fixes and normalization
        intermediate_file <- list.files(path = intermediate_dir, pattern = accession, full.names = T)
        result_file <- read.ilmn(intermediate_file)
        options(warn = 1)
        print('fix expr')
        result_file <- fixExpr(result_file)
        print('neqc')
        result_file <- neqc(result_file)
        options(warn = 0)
      }
      
      if ('affy' %in% platform_type){
        print('affy')
        GEO_dl_files <- list.files(GEO_dl_dir)
        
        # if a single study has multiple affy platforms, this check prevents the data from being quantile normalized togeter
        cel_ix <- T
        if (sum(grepl(pattern = accession, x = GEO_dl_files) & grepl(pattern = platform, x = GEO_dl_files)) > 0){
          print(paste0('multiple platforms detected: ', accession))
          print('parsing GSMs from GEO download')
          file_ix <- which(grepl(pattern = accession, x = GEO_dl_files) & grepl(pattern = platform, x = GEO_dl_files))
          GSM_IDs <- colnames(readRDS(paste0(GEO_dl_dir, GEO_dl_files[file_ix])))
          
          cel_files <- list.files(path = raw_dir, pattern = 'CEL|cel')
          cel_ix <- sapply(cel_files, parseGSM) %in% GSM_IDs
          cel_files <- list.files(path = raw_dir, pattern = 'CEL|cel', full.names = T)[cel_ix]
          
        } 
        
        cel_files <- list.files(path = raw_dir, pattern = 'CEL|cel', full.names = T)[cel_ix]
        compress_flag <- sum(grepl(pattern = 'gz', x = cel_files)) > 0
        
        if(platform %in% c('GPL11532', 'GPL5175')){
          affybatch <- read.celfiles(cel_files)
          result_file <- oligo::rma(affybatch)
        } else {
          affybatch <- read.affybatch(filenames = cel_files, compress = compress_flag)
          result_file <- affy::rma(affybatch)
        }
        # run RMA
       
      }
      
      if(!is.null(result_file)){
        saveRDS(object = result_file, file = paste0(output_dir, accession, '.RDS'))
      } else {
        print(paste0(accession, ': not processed'))
      }
      
    }
    setwd(paste0(Sys.getenv('DARPA'), '/Input/'))
  }

  
  ## match phenodata with normalized data
  ## output exprset
  #source('updateAnnotations.R')
  if(!dir.exists(paste0(out_dir, 'GEO_processed_exprsets'))){
    dir.create(paste0(out_dir, 'GEO_processed_exprsets'))
  }

}


for (platform in platforms){
  print(platform)
  updateAnnotations(processed_dir = paste0(out_dir, 'GEO_processed/', platform, '/'), 
                    dl_dir = paste0(out_dir, 'GEO_dl/'), 
                    output_dir = paste0(out_dir, 'GEO_processed_exprsets/'), 
                    platform = platform,
                    phenoData_xlsx = phenoData,
                    overwrite = annotation_overwrite_flag)
}

# # munge exprset into metaintegrator object
source('metaIntegratorMunging.R')
source('makingDataList.R')
