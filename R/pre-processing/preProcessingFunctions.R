preProcessILMN <- function(raw_file, intermediate_dir){
  print(raw_file)
  # make sure that every row has the right number of elements... trim the top of files as needed
  length_check <- F
  print('parsing start point')
  i <- 0
  while(!length_check){
    tmp <- read.delim(file = paste0(raw_dir, raw_file), sep = '\t', header = F, nrows = 1, skip = i)
    
    if (length(tmp) > 10 & sum(grepl(pattern = '^ref_id$|^id_ref$|^idref$|^id-ref$|^probe|^ilmn_probe', x = sapply(sapply(tmp, tolower), trimws))) > 0){
      length_check <- T
    } else {
      i <- i + 1
    }
    
    if (i > 0 & i %% 10 == 0){
      print(i)
    }
  }
  
  print('reading file')
  tmp <- read.delim(file = paste0(raw_dir, raw_file), sep = '\t', header = T, skip = i, stringsAsFactors = F)
  
  # filter out blank columns
  na_rm <- colSums(is.na(tmp)) > 0
  tmp <- tmp[, !na_rm]
  
  # make sure all numeric columns were read in as numeric
  detectNum <- function(df_col){
    return((mean(is.na(suppressWarnings(as.numeric(df_col)))) < 0.1 & !is.numeric(df_col)))
  }
  num_conv_ix <- sapply(tmp, detectNum)
  if (sum(num_conv_ix) > 0){
    tmp[, num_conv_ix] <- sapply(FUN = as.numeric, X = tmp[, num_conv_ix])
  }
  
  
  # find and standardize the probe column
  print('standardizing probe column')
  probe_ix <- grep(pattern = 'id_ref|idref|ref_id|id-ref|id.ref', x = tolower(names(tmp)))
  if(length(probe_ix) == 0){
    probe_rm <- F
    probe_ix <- grep(pattern = '^Probe_Id$|^PROBE_ID$|^ProbeID$', x = names(tmp), ignore.case = T)
    probe_ix <- probe_ix[probe_ix < 10]
    if(length(probe_ix) == 0){
      probe_ix <- grep(pattern = '^Probe$', x = names(tmp))
      probe_ix <- probe_ix[probe_ix < 10]
      if(length(probe_ix) == 0){
        probe_ix <- grep(pattern = '^ILMN_Probe$', x = names(tmp))
        probe_ix <- probe_ix[probe_ix < 10]
        if (length(probe_ix) == 0){
          print(raw_file)
          error()
        }
      }
    }
  } else {
    # if we found an ID_REF, we can get rid of the other probe id columns
    probe_rm <- grepl(pattern = '^Probe_Id$|^PROBE_ID$|^ProbeID$', x = names(tmp)) 
  }
  
  print('removing extra columns')
  #print(names(tmp)[probe_ix])
  names(tmp)[probe_ix] <- 'Probe'
  tmp <- tmp[, !probe_rm]
  #rm_ix <- 1:ncol(tmp) < probe_ix
  #tmp <- tmp[, !rm_ix]
  
  if(raw_file == "GSE61821_non-normalized.txt.gz"){
    tmp[,2] <- as.numeric(tmp[,2])
  }
  
  # make sure you have an avg_signal and detection column for every sample
  probe_ix <- grep(pattern = '^Probe$', x = names(tmp))
  num_ix <- which(sapply(tmp, is.numeric))
  tmp <- tmp[, c(probe_ix, num_ix)]
  
  # remove BEAD_STDERR and Avg_NBEADS columns
  rm_ix <- grepl(pattern = 'BEAD_STDERR|Avg_NBEADS|^X$', x = names(tmp))
  if(sum(rm_ix) > 0){
    tmp <- tmp[, !rm_ix]
  }
  
  # remove misc offending columns
  rm_ix <- grepl(pattern = 'ProbeID|ENTREZ|PROBE|ARRAY|^GI$|CHROMOSOME|Entrez|Array|Probe_', x = names(tmp))
  if(sum(rm_ix) > 0){
    tmp <- tmp[, !rm_ix]
  }
  
  print('checking column names')
  # check that 'AVG_SIGNAL' and 'Detection' are written where needed
  if(sum(grepl(pattern = 'AVG_Signal', x = names(tmp))) > 0){
    # we need to check that samples and pvalues have matching names...
    sample_name <- gsub(pattern = 'AVG_Signal', replacement = '', x = names(tmp)[2])
    sample_name <- gsub(pattern = '.', replacement = '', x = sample_name, fixed = T)
    pval_name <- gsub(pattern = '.', replacement = '', x = names(tmp), fixed = T)[-2]
    if (sum(grepl(pattern = sample_name, x = pval_name)) > 0){
      print(T)
    } else {
      print(paste0('relabeling: ', raw_file))
      col_names <- vector(mode = 'character', length = ncol(tmp))
      col_names[1] <- 'Probe'
      for (i in 2:ncol(tmp)){
        if (i %% 2 == 0){
          sample_name <- gsub(pattern = '$X|$X_|AVG_Signal.|.AVG_Signal|AVG_SIGNAL', replacement = '', names(tmp)[i])
          col_names[i] <- paste0('AVG_Signal.', sample_name)
          col_names[i+1] <- paste0('Detection.Pval.', sample_name)
        }
      }
      names(tmp) <- col_names
    }
    
  } else {
    #print('renaming')
    col_names <- vector(mode = 'character', length = ncol(tmp))
    col_names[1] <- 'Probe'
    
    sample_name <- gsub(pattern = 'AVG_Signal', replacement = '', x = names(tmp)[2])
    sample_name <- gsub(pattern = '.', replacement = '', x = sample_name, fixed = T)
    pval_name <- gsub(pattern = '.', replacement = '', x = names(tmp), fixed = T)[-2]
    if (sum(grepl(pattern = sample_name, x = pval_name)) > 0){
      print(T)
      pval_ix <- grepl(pattern = 'detection', x = names(tmp), ignore.case = T)
      sample_ix <- !pval_ix
      sample_ix[1] <- F
      names(tmp)[sample_ix] <- paste0('AVG_Signal.', names(tmp)[sample_ix])
      
    } else {
      for (i in 2:ncol(tmp)){
        if (i %% 2 == 0){
          sample_name <- gsub(pattern = '$X|$X_|', replacement = '', names(tmp)[i])
          col_names[i] <- paste0('AVG_Signal.', sample_name)
          col_names[i+1] <- paste0('Detection.Pval.', sample_name)
        }
      }
      names(tmp) <- col_names
    }
    
  }
  
  if(raw_file == 'GSE86434_non-normalized_data.txt.gz'){
    rm_ix <- grepl(pattern = 'P008902.TTCgP.01', names(tmp))
    tmp <- tmp[, !rm_ix]
  }
  
  print('writing file')
  write.table(file = paste0(intermediate_dir, tools::file_path_sans_ext(raw_file)), x = tmp, sep = '\t', quote = F, row.names = F)
  return(0)
}

fixExpr <- function(lim_obj){
  ix <- which(rowSums(is.na(lim_obj$other$Detection)) > 0)
  if(length(ix) > 0){
    if (length(ix) < 100){
      lim_obj$E <- lim_obj$E[-ix, ]
      lim_obj$other$Detection <- lim_obj$other$Detection[-ix, ]
    } else {
      print('ERROR')
    }
  }
  
  ix <- rowSums(is.na(lim_obj$E)) > 0
  if(sum(ix) > 0){
    if (mean(ix) < 0.1){  # if more than 10% of probes have an issue, print ERROR
      lim_obj$E <- lim_obj$E[-ix, ]
      lim_obj$other$Detection <- lim_obj$other$Detection[-ix, ]
    } else {
      print('ERROR')
    }
  }
  return(lim_obj)
}


