library(illuminaHumanv4.db)
library(Biobase)
library(dplyr)
library(MetaIntegrator)
library(parallel)
library(aveytoolkit)
library(openxlsx)

library(illuminaHumanv4.db)
library(illuminaHumanv3.db)
library(illuminaHumanv2.db)
library(hgu133plus2.db)
library(hgu133a2.db)
library(hgu133a.db)
library(hgu133b.db)
library(hgu133plusPM.db)
library(hgu219.db)
library(hgfocus.db)
library(primeview.db)
library(huex10sttranscriptcluster.db)
library(hugene11sttranscriptcluster.db)
library(hugene10sttranscriptcluster.db)
library(hta20transcriptcluster.db)
library(clariomdhumantranscriptcluster.db)

# YOU NEED TO COME UP WITH A SYSTEM TO SPECIFY PLATFORM FOR EACH FILE
source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))

output_dir <- Sys.getenv('GEO_DATA')
output_dir <- '~/Data/DARPA/Specificity/Aging/'
output_dir <- '~/Documents/DARPA/Specificity/COVID19/'
output_dir <- '~/Data/DARPA/'
working_dir <- paste0(output_dir, 'GEO_processed_exprsets/')
files <- list.files(working_dir)
processed_files <- list.files(paste0(output_dir, 'metaIntegratorObjects'))

files_to_process <- files
#files_to_process <- files[!files %in% processed_files]

platforms <- sort(unique(sapply(files_to_process, function(x){strsplit(x, split = '_|.RDS')[[1]][2]})))
annotation_file <- read.xlsx(paste0(Sys.getenv('DARPA'), 'Input/platform_specific/key_columns.xlsx'))
rownames(annotation_file) <- annotation_file$PLATFORM

metaIntegratorMunge <- function(file, input_dir){
  print(paste0(input_dir, file))
  eset <- readRDS(paste0(input_dir, file))
  
  platform <- strsplit(file, split = '_|.RDS')[[1]][2]

  # check for missing expression values (should only be in data taken directly from GEO)
  eset <- fixEsetExpr(eset)
  
  annotation = strsplit(file, split = '_|.RDS')[[1]][2]
  #print(annotation)
  
  db = NULL
  switch(annotation, 
         'GPL10558' = {db = illuminaHumanv4.db},
         'GPL6947' = {db = illuminaHumanv3.db},
         'GPL571' = {db = hgu133a2.db},
         'GPL570' = {db = hgu133plus2.db},
         'GPL97' = {db = hgu133b.db},
         'GPL96' = {db = hgu133a.db},
         'GPL11532' = {db =  hugene11sttranscriptcluster.db},
         'GPL13158' = {db = hgu133plusPM.db},
         'GPL13667' = {db = hgu219.db},
         'GPL201' = {db = hgfocus.db},
         'GPL5175' = {db = huex10sttranscriptcluster.db},
         'GPL6102' = {db = illuminaHumanv2.db},
         'GPL6244' = {db = hugene10sttranscriptcluster.db},
         'GPL6883' = {db = illuminaHumanv3.db},
         'GPL6884' = {db = illuminaHumanv3.db},
         'GPL17586' = {db = hta20transcriptcluster.db},
         'GPL23126' = {db = clariomdhumantranscriptcluster.db},
         'GPL15207' = {db = primeview.db},
         NULL
         )
  
  if (is.null(db)){
      if (is.na(annotation_file[platform, 'PDAT_COLUMN'])){
        print('NO KEY MAPPING--SKIPPING STUDY')
        return(NULL)
      }
      col_ix <- grep(pattern = annotation_file[platform, 'PDAT_COLUMN'], x = names(fData(eset)), value = T)
      if(length(col_ix) != 1){
        stop('key column does not match')
      }
      keytype <- annotation_file[platform, 'KEYTYPE']
      
      map_keys <- fData(eset)[, col_ix]
      
      if (keytype != 'ENTREZID'){
        map_keys <- sapply(strsplit(fData(eset)[,col_ix], "\\."), "[", 1)
        comma_ix <- sapply(map_keys, function(x){sum(grepl(pattern = ',', x))})
        if(sum(comma_ix) > 0){
          map_keys <- sapply(map_keys, function(x){strsplit(x, split = ',')[[1]][1]})
        }
        
        # if(keytype == 'parseACC'){
        #   map_keys <- sapply(map_keys, function(x){
        #     tmp_output <- stringr::str_extract(string = x, pattern = 'HGNC Symbol;Acc:\\d+')
        #     tmp_output <- gsub(pattern = 'HGNC Symbol;Acc:', replacement = '', x = tmp_output)
        #     return(tmp_output)
        #   })
        #   keytype <- 'ACCNUM'
        # }
        
        map_keys <- mapIds(x = org.Hs.eg.db, 
                           keys = map_keys, 
                           column = 'ENTREZID', 
                           keytype = keytype)
#<<<<<<< HEAD
        debug_file <- '~/Dropbox/changed_studies.tsv'
        debug_output <- data.frame('Study' = file, stringsAsFactors = F)
        
        if(!file.exists(debug_file)){
          write.table(file = debug_file, x = debug_output, row.names = F, col.names = F, append = F)
        } else {
          write.table(file = debug_file, x = debug_output, row.names = F, col.names = F, append = T)
        }
        
        length_key_list <- length(map_keys)
        #map_keys <- sapply(map_keys, as.character)
#=======
        map_keys <- sapply(map_keys, as.character)
        map_keys <- as.character(map_keys)
        map_keys[map_keys == 'character(0)'] <- NA
        length_key_vec <- length(map_keys)
        if(length_key_vec != length_key_list){
          stop('key length mapping error')
        }
#>>>>>>> a889a49f532a52b3c6b7a496083d84a829a4858e
      }
      
      names(map_keys) <- rownames(exprs(eset))
  } else {
    map_keys = tryCatch({mapIds(x = db,
                    keys = rownames(exprs(eset)), 
                    keytype = 'PROBEID', 
                    column = 'ENTREZID')},
                    error = function(err){
                      print(err)
                      return(NULL)
                    })
  }

  
  # if(is.null(map_keys)){
  #   obj <- list()
  #   obj$expr <- out$exprsVals
  #   return(obj)
  # }
  
  #names(map_keys) <- rownames(exprs(eset))

  # summarize probes: max
  # build metaIntegrator obj
  obj <- list()
  obj$expr <- exprs(eset)
  
  # summarize genes in expr
  tmp <- data.frame('probe' = as.character(names(map_keys)), 'gene' = as.character(map_keys), exprs(eset), stringsAsFactors = F)
  tmp <- tmp %>%
    filter(!is.na(gene))
  tmp$avg_expr <- apply(tmp %>% select_if(is.numeric), 1, FUN = mean, na.rm = T)
  tmp <- tmp %>%
    group_by(gene) %>%
    filter(avg_expr == max(avg_expr)) %>%
    ungroup()
  
  map_keys <- tmp$gene
  names(map_keys) <- tmp$gene
  
  tmp <- tmp %>%
    dplyr::select(-avg_expr) %>%
    select_if(is.numeric)
  
  obj$expr <- fixExprLog(as.matrix(tmp))
  rownames(obj$expr) <- as.character(map_keys)
  
  rm(tmp)
  
    # if(max(table(map_keys)) == 1){# | is.null(db)){
    #   out <- list()
    #   out$exprsVals <- exprs(eset)
    #   out$probes <- names(map_keys)
    #   names(out$probes) <- map_keys
    # } else {
    #   out <- collapseDataset(exprsVals = exprs(eset), mapVector = map_keys, oper = 'max', returnProbes = T)
    #   # remove blank probes
    #   #rm_ix <- trimws(out$probes) == '' | is.na(out$probes)
    #   #out$probes <- out$probes[!rm_ix]
    #   #obj$expr <- out$exprsVals[!rm_ix, ]
    #   
    #   obj$expr <- out$exprsVals
    #   #rownames(obj$expr) <- names(out$probes)
    #   #ix <- !is.na(rownames(obj$expr))
    #   #obj$expr <- obj$expr[ix,]
    #   
    #   # keys
    #   map_keys <- names(out$probes)
    #   names(map_keys) <- names(out$probes)
    #   #obj$keys <- obj$keys[ix]
    #   rm(out)
    # }
  
  
  
  
  obj$keys <- map_keys
  
  
  # pheno
  obj$pheno <- pData(eset)
  
  # class
  obj$class <- rep(0, nrow(obj$pheno))
  names(obj$class) <- rownames(obj$pheno)
  colnames(obj$expr) <- rownames(obj$pheno)
  
  # formattedName
  obj$formattedName <- tools::file_path_sans_ext(file)
  
  if(MetaIntegrator::checkDataObject(obj, objectType = 'Dataset')){
    return(obj)
  } else {
    stop('error')
    return(NULL)
  }
}



for (platform in platforms){
  print(platform)
  
  platform_files_to_process <- grep(pattern = platform, x = files_to_process, value = T)
  
  #MI_objs <- mclapply(platform_files_to_process, metaIntegratorMunge, working_dir, mc.cores = 4)
  MI_objs <- lapply(platform_files_to_process, function(input_file, working_dir){
  #MI_objs <- lapply(platform_files_to_process, function(input_file, working_dir){
    MIobj <- metaIntegratorMunge(input_file, working_dir)
    flag <- tryCatch({checkDataObject(MIobj,  objectType = 'Dataset')}, error = function(err){
      print(err)
      return(F)
    })
    if(flag){
      return(MIobj)
    } else {
      return(NULL)
    }
  #}, working_dir = working_dir, mc.cores = 4)
  }, working_dir = working_dir)
  
  if(sum(sapply(MI_objs, is.null)) == 0){
    print(':)')
  } else {
    print(':(')
  }
  
  lapply(MI_objs, 
           function(x, output_dir){
             saveRDS(object = x, 
                     file = paste0(output_dir, 
                                   x$formattedName, 
                                   '.RDS')); 
             return(x$formattedName)}, 
           output_dir = paste0(output_dir, 'metaIntegratorObjects/'))
  
  
}

