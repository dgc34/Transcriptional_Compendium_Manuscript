library(dplyr)
library(MetaIntegrator)
library(ggplot2)
library(openxlsx)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)
library(RColorBrewer)
library(stringr)
library(caret)
library(doParallel)

createColorScale <- function(){
  myColors <- brewer.pal(4,"Set1")
  names(myColors) <- c('Discovery', 'Validation', 'New', 'Null')
  
  dat <- data.frame('x' = c(1,2,3,4), 'y' = c(1,2,3,4), color = myColors)
  ggplot(dat, aes(x = x, y = y, color = color)) + 
    geom_point()
  
  return(myColors)
}

parseSets <- function(str){
  str <- trimws(str)
  if(str == '-' | str == '—' | str  == '--' | is.na(str)){
    set <- ''
  } else {
    set <- trimws(unlist(strsplit(str, split = ';')))
    blank_ix <- set == ''
    set <- set[!blank_ix]
  }
  
  return(set)
}

mapGenes <- function(genes, key_type = 'SYMBOL', column_type = 'ENTREZID'){
  #print(genes)
  #genes <- alias2SymbolTable(genes) # this is wrong, do not treat all symbols as aliases
  
  na_ix <- is.na(genes)
  genes[na_ix] <- ''
  #print(genes)
  
  if(length(genes) == 1){
    if(genes == ''){
      entrez_ids <- ''
    } else {
      entrez_ids <- tryCatch({
        #genes <- HGNChelper::checkGeneSymbols(genes)$Suggested.Symbol
        #genes <- genes[!is.na(genes)]
        tmp_genes <- mapIds(x = org.Hs.eg.db, keys = genes, keytype = key_type, column = column_type)
        tmp_genes <- tmp_genes[!sapply(tmp_genes, is.null)]
        return(tmp_genes)
      }, error = function(err){
        print(err)
        return('')
      })
    }
  } else {
    entrez_ids <- tryCatch({
      #genes <- HGNChelper::checkGeneSymbols(genes)$Suggested.Symbol
      #genes <- genes[!is.na(genes)]
      tmp_genes <- mapIds(x = org.Hs.eg.db, keys = genes, keytype = key_type, column = column_type)
      tmp_genes <- tmp_genes[!sapply(tmp_genes, is.null)]
      return(tmp_genes)
    }, error = function(err){
      print(err)
      return('')
    })
  }
  return(entrez_ids)
}

parseAccession <- function(formatted_name){
  if(!is.character(formatted_name)){
  	stop("'formatted_name' should be a character")
  }
  accession <- str_extract(string = formatted_name, pattern = 'GSE\\d+')
  return(accession)
  #accession <- strsplit(formatted_name, split = '_')[[1]][1]
}

parsePlatform <- function(formatted_name){
  str <- strsplit(formatted_name, split = '_')[[1]]
  platform_ix <- grepl(pattern = 'GPL', str)
  if(sum(platform_ix) == 1){
    return(str[platform_ix])
  } else {
    return(NULL)
  }
}

parseClass <- function(MIobj, signature_row){
  #signature <- sig2Meta(signature_row)
  
  positive_class <- parseSets(signature_row$Positive.Class)
  negative_class <- parseSets(signature_row$Negative.Class)
  pos_str <- vector(mode = 'character', length = length(positive_class))
  for (i in 1:length(positive_class)){
    pos_str[i] <- switch(positive_class[i], 
                         'Bacterial' = 'Bacteria',
                         'Bacteria' = 'Bacteria',
                         'Viral' = 'Virus',
                         'Virus' = 'Virus',
                         'Influenza' = 'Virus',
                         'Flu' = 'Virus',
                         'Resp.Viral' = 'Virus',
                         'Non-Infectious' = 'Other.NonInfectious',
                         'Sepsis' = 'Other.Infectious', 
                         'HIV' = 'Virus',
                         '')
  }
  neg_str <- vector(mode = 'character', length = length(negative_class))
  for (i in 1:length(negative_class)){
    neg_str[i] <- switch(negative_class[i], 
                         'Bacterial' = 'Bacteria',
                         'Bacteria' = 'Bacteria',
                         'Non-Bacterial' = 'Non',
                         'Non-Bacteria' = 'Non',
                         'Viral' = 'Virus',
                         'Virus' = 'Virus',
                         'Non-Viral' = 'Non',
                         'Non-Virus' = 'Non',
                         'Non-Influenza' = 'Non',
                         'Resp.Viral' = 'Virus',
                         'Non-Infectious' = 'Other.NonInfectious', 
                         'Healthy' = 'Healthy',
                         '')
  }
  
  MIobj <- replaceClass(MIobj, positive_str = pos_str, negative_str = neg_str)
  return(MIobj)
}

filterClass <- function(MIobj){
  # return a logical vector
  # T if there are two classes in MIobj$class
  MIobj$class <- MIobj$class[!is.na(MIobj$class)]
  if(length(unique(MIobj$class)) == 2){
    return(T)
  } else {
    return(F)
  }
}

filterGenes <- function(MIobj, signature, threshold = 0){
  # need to calculate fraction of positive and negative signature genes found in data set
  
  # remove empty '' genes
  positive_genes <- setdiff(signature$posGeneNames, '')
  negative_genes <- setdiff(signature$negGeneNames, '')
  
  pos_flag <- T
  if (length(positive_genes) > 0){
    # if we have positive genes, check that we meet the minimum fraction
    pos_fraction = mean(positive_genes %in% MIobj$keys)
    if (pos_fraction <= threshold){
      pos_flag <- F
    }
  }
  
  neg_flag <- T
  if (length(negative_genes) > 0){
    # if we have negative genes, check that we meet the minimum fraction
    neg_fraction = mean(negative_genes %in% MIobj$keys)
    if (neg_fraction <= threshold){
      neg_flag <- F
    }
  }
  
  if (length(positive_genes) == 0 & length(negative_genes) == 0){
    return(F)
  }
  
  if (pos_flag & neg_flag){ # if we're above the threshold for both pos & neg genes, we can proceed
    return(T)
  } else {
    return(F)
  }
}


replaceClass <- function(MIobj, positive_str, negative_str){
  
  print(paste(MIobj$formattedName, positive_str, negative_str))
  new_class <- MIobj$pheno$Class
  
  positive_str <- paste(positive_str, collapse = '|')
  negative_str <- paste(negative_str, collapse = '|')
  
  pos_ix <- grepl(pattern = positive_str, x = new_class)
  if(negative_str == 'Non'){
    neg_ix <- !pos_ix
  } else {
    neg_ix <- grepl(pattern = negative_str, x = new_class)
  }
  
  new_class[pos_ix] <- 1
  new_class[neg_ix] <- 0
  na_ix <- !(pos_ix | neg_ix)
  new_class[na_ix] <- NA
  new_class <- as.numeric(new_class)
  
  MIobj$expr <- MIobj$expr[, !na_ix]
  MIobj$pheno <- MIobj$pheno[!na_ix,]
  new_class <- new_class[!na_ix]
  
  rownames(MIobj$pheno) <- colnames(MIobj$expr)
  MIobj$class <- new_class
  names(MIobj$class) <- rownames(MIobj$pheno)
  
  return(MIobj)
}

formatLabel <- function(signature_row, nchar = 7){
  sig <- signature_row
  sig <- sig %>%
    mutate_at(c('Author', 'Positive.Class', 'Negative.Class'), 
              substr, start = 0, stop = nchar)
  if (is.na(sig$Classifier.Label)){
    lab <- paste0(paste(sig$Author, sig$Year, sep = ' '), 
                  '\n', 
                  paste(sig$Positive.Class, 'vs.', sig$Negative.Class, sep = ' '))
  } else {
    lab <- paste0(paste(sig$Author, sig$Year, sig$Classifier.Label, sep = ' '), 
                  '\n', 
                  paste(sig$Positive.Class, 'vs.', sig$Negative.Class, sep = ' '))
  }
  return(lab)
}

formatLabel2 <- function(signature_row, nchar = 7){
  sig <- signature_row
  sig <- sig %>%
    mutate_at(c('Author', 'Positive.Class', 'Negative.Class'), 
              substr, start = 0, stop = nchar)
  if (is.na(sig$Print.Label)){
    lab <- paste0(paste(sig$Author, sig$Year, sep = ' '))
  } else {
    lab <- paste0(paste(sig$Author, sig$Year, sep = ' '), sig$Print.Label)
  }
  return(lab)
}

parseSignatureAccessions <- function(signature_row, type){
  if(! type %in% c('discovery', 'validation')){
    stop('type must be either "discovery" or "validation"')
  }
  if(type == 'discovery'){
    str <- parseSets(signature_row$Discovery.Accessions)
  }
  if(type == 'validation'){
    str <- parseSets(signature_row$Validation.Accessions)
  }
  return(str)
}

sig2Meta <- function(signature_row, format_option = 'v1', ...){
  signature <- list()
  signature$posGeneNames <- as.character(mapGenes(parseSets(signature_row$Positive.Genes)))
  signature$negGeneNames <- as.character(mapGenes(parseSets(signature_row$Negative.Genes)))
  signature$FDRThresh <- 1
  signature$effectSizeThresh <- 0
  signature$numberStudiesThresh <- 0
  signature$heterogeneityPvalThresh <- 1
  signature$isLeaveOneOut <- F
  #signature$filterDescription <- switch(format_option, 
  #                                      'v2' = formatLabel2(signature_row, nchar = 25, ...),
  #                                      'v1' = formatLabel(signature_row, nchar = 25, ...))
  signature$filterDescription <- signature_row$Signature.Label
  signature$timestamp <- Sys.time()
  
  signature$discovery <- parseSignatureAccessions(signature_row, type = 'discovery')
  signature$validation <- parseSignatureAccessions(signature_row, type = 'validation')
  flag <- checkDataObject(signature, objectType = 'MetaFilter')
  if(flag){
    return(signature)
  } else {
    stop('invalid signature')
  }
}

scoreStudy <- function(MIobj, signature, performance_metric, method = 'zScore', n_perm = 1000, specific_method){
  if (method == 'zScore'){
    score <- switch(performance_metric, 
                    'AUROC' = calculateAUROC(MIobj, signature, specific_method),
                    NULL)
  }
  
  if (method == 'singleSample'){
    print(MIobj$formattedName)
    score <- calculateSingleSampleAUROC(MIobj, signature)
  }
  
  if (method == 'logit'){
    print(paste0(MIobj$formattedName, '_logit'))
    score <- calculateLogit(MIobj, signature, performance_metric, n_perm)
  }
  
  return(score)
  
}

scoreStudySignificant <- function(MIobj, signature, performance_metric, method = 'zScore', n_perm = 1000, specific_method, suppressMessages = F){
  # NOTE: ONLY Z-SCORE IS TESTED
  
  # need to rewrite functions to calculate scores so we can permute them before computing auc / pvalue
  if (method == 'zScore'){
    labels <- MIobj$class
    scores <- calculateScoreRobust(filterObject = signature, datasetObject = MIobj, method = 'geomMean', zScore = T, suppressMessages = suppressMessages)
  }
  
  if (method == 'singleSample'){
    scores <- calculateSingleSampleScore(MIobj, signature)
    labels = MIobj$class
    if (length(scores) == 0){
      return(NA)
    }
  }
  
  if (method == 'logit'){
    print(paste0(MIobj$formattedName, '_logit'))
    tmp <- calculateLogitScores(MIobj, signature, performance_metric, n_perm)
    scores <- tmp$scores
    labels <- tmp$labels
  }
  
  auc <- auroc(bool = labels, score = scores)
  pval <- calculatePval(labels = labels, values = scores, true_value = auc, n_perm)
  
  return(data.frame('Scores' = auc, 'P.Value' = pval, stringsAsFactors = F))
}

calculateAUROC <- function(MIobj, signature, method, suppressMessages = F){
  labels <- MIobj$class
  scores <- calculateScoreRobust(filterObject = signature, datasetObject = MIobj, method = method, suppressMessages = suppressMessages)
  if(all(sapply(scores, is.na))){
    roc <- NA
  } else {
    roc <- calculateROC(labels = labels, predictions = scores, AUConly = T)
  }
  return(as.numeric(roc))
}

calculateSingleSampleScore <- function(MIobj, signature){
  expr_dat <- collapseDataset(exprsVals = MIobj$expr, mapVector = MIobj$keys, oper = max)
  rankData <- rankGenes(expr_dat)
  upSet <- GeneSet(unique(signature$posGeneNames))
  dnSet <- GeneSet(unique(signature$negGeneNames))
  scoredf <- simpleScore(rankData, upSet = upSet, downSet = dnSet)
  
  output <- scoredf$TotalScore
  names(output) <- rownames(scoredf)
  return(output)
  
}

calculateSingleSampleAUROC <- function(MIobj, signature){
  pred <- calculateSingleSampleScore(MIobj, signature)
  if (length(pred) == 0){
    return(NA)
  }
  ROCobj <- roc(response = MIobj$class, predictor = pred)
  return(as.numeric(ROCobj$auc))
}

calculateLogit <- function(MIobj, sig, performance_metric, N_PERM = 1000, parallel = F){
  if(parallel){
    
    if (nrow(MIobj$pheno) < 50){
      cl <- makePSOCKcluster(6)
      #print('6')
    } else if (nrow(MIobj$pheno) < 100){
      cl <- makePSOCKcluster(2)
      #print('2')
    } else if(nrow(MIobj$pheno) < 500){
      cl <- makePSOCKcluster(1)
      #print('1')
    } else {
      cl <- makePSOCKcluster(1)
      #print('1')
    } 
      
    registerDoParallel(cl)
  } else {
    #cl <- makePSOCKcluster(1)
    #registerDoParallel(cl)
  }
  
  out <- tryCatch(
    {
      #print(performance_metric)
      # performance_metric should be one of {accuracy, acc_pval, AUROC, AUROC_pval, debug}
      
      # add a check that the labels are approximately balanced, otherwise we have to add that
      # if the labels are unbalanced, you need nested cross validation because you want to estimate the effect of downsampling (SMOTE)
      # on performance
      
      # considerations
      ## can you enforce feature signs in caret? let's start with just concat feature lists and go from there
      
      # performance metrics
      ## can report logit LOOCV accuracy at a cutoff of 0.5
      ## can report p-value of LOOCV accuracy against permuted class labels 
      ## can report AUROC from hold out samples; if you're generating probabilities for each sample, these are directly comparable for generating a ROC
      
      # set a flag in case we can avoid permutation
      perm_flag <- T
      
      # filter MIobj expr based on signature genes
      sig_genes <- setdiff(c(sig$posGeneNames, sig$negGeneNames), '') # concat pos and neg genes, rm empty set
      sig_gene_ix <- rownames(MIobj$expr) %in% sig_genes
      MIobj$mat <- MIobj$expr[sig_gene_ix, ]
      MIobj$mat <- t(MIobj$mat)
      
      # check for NA expression values
      if(sum(is.na(MIobj$mat)) > 0){
        sample_rm_ix <- apply(X = MIobj$mat, MARGIN = 1, FUN = function(x){return(sum(is.na(x)))}) > 0
        
        if(mean(sample_rm_ix) == 1){
          print(paste0('all samples have NA values ', MIobj$formattedName))
          return(NA)
        }
        
        MIobj <- filterMIobj(MIobj, !sample_rm_ix)
        MIobj$mat <- MIobj$mat[!sample_rm_ix,]
      }
      
      # you cannot fit a linear model with more predictors than observations, and we do LOOCV so this is N + 1 > M
      if(ncol(MIobj$mat) > (nrow(MIobj$mat) - 1)){
        print(paste0('More params than observations for ', MIobj$formattedName))
        return(NA)
      }
      
      if(min(table(MIobj$class)) == 1){
        print(paste0('Not enough observations for ', MIobj$formattedName))
        return(NA)
      }
      
      # train LOOCV logit model
      original_model <- train(x = MIobj$mat, 
                              y = as.factor(paste0('class_', MIobj$class)), 
                              method = 'vglmAdjCat', 
                              trControl = trainControl(method = 'LOOCV', 
                                                       savePredictions = T, 
                                                       classProbs = T))
      
      orig_results_df <- original_model$pred %>%
        dplyr::filter(parallel == F)
      original_accuracy <- original_model$results$Accuracy[1]
      
      # pull out AUROC
      original_AUROC <- as.numeric(pROC::roc(response = MIobj$class, predictor = orig_results_df$class_1)$auc)
      rm(original_model) # free memory 
      
      acc_vect <- vector(mode = 'double', length = N_PERM)
      AUROC_vect <- vector(mode = 'double', length = N_PERM)
      
      if (performance_metric %in% c('accuracy', 'AUROC')){
        perm_flag <- F
      }
      
      if (perm_flag){
        already_analyzed <- list.files("~/Dropbox/DARPA_ECHO/logit_permutations/")
        current_study_ix <- grepl(pattern = MIobj$formattedName, x = already_analyzed)
        current_author <- strsplit(x = sig$filterDescription, split = '\n')[[1]][1]
        current_filter <- strsplit(x = sig$filterDescription, split = '\n')[[1]][2]
        current_author_ix <- grepl(pattern = current_author, x = already_analyzed)
        current_filter_ix <- grepl(pattern = current_author, x = already_analyzed)
        if(max(sum(current_study_ix & current_author_ix & current_filter_ix)) == 1){
          print(paste0('already analyzed!'))
          return(NA)
        }
        
        # permutation testing
        for (i in 1:N_PERM){
          if (i %% ifelse(N_PERM < 101, 10, 25) == 0){
            print(paste0('iteration', i, ' / ', N_PERM))
          }
          
          current_MIobj <- permuteLabels(MIobj)
          current_model <- train(x = current_MIobj$mat, 
                                 y = as.factor(paste0('class_', current_MIobj$class)), 
                                 method = 'vglmAdjCat', 
                                 trControl = trainControl(method = 'LOOCV', 
                                                          savePredictions = T, 
                                                          classProbs = T))
          results_df <- current_model$pred %>%
            dplyr::filter(parallel == F) # what does this do?
          
          acc_vect[i] <- current_model$results$Accuracy[1]
          AUROC_vect[i] <- as.numeric(pROC::roc(response = current_MIobj$class, predictor = results_df$class_1)$auc)
        }
        rm(current_model) # free memory
        acc_pval <- mean(acc_vect >= original_accuracy)
        auroc_pval <- mean(AUROC_vect >= original_AUROC)
        return_df <- data.frame('AUROC' = original_AUROC, 
                                'AUROC_pval' = auroc_pval, 
                                'Accuracy' = original_accuracy, 
                                'Accuracy_pval' = acc_pval)
        
        debug_df <- data.frame('Perm_accuracy' = acc_vect, 
                               'Perm_AUROC' = AUROC_vect, 
                               'AUROC' = rep(original_AUROC, N_PERM), 
                               'Accuracy' = rep(original_accuracy, N_PERM), 
                               'AUROC_pval' = rep(auroc_pval, N_PERM),
                               'Accuracy_pval' = rep(acc_pval, N_PERM))
        sig$filterDescription <- gsub(sig$filterDescription, pattern = ' ', replacement = '')
        sig$filterDescription <- gsub(sig$filterDescription, pattern = '.', replacement = '', fixed = T)
        filename_out <- paste0('~/Dropbox/DARPA_ECHO/logit_permutations/', MIobj$formattedName, '_', sig$filterDescription, '_logit.csv')
        print(filename_out)
        write.csv(file = filename_out, 
                  x = debug_df, 
                  row.names = F, 
                  quote = F)
      } else {
        acc_pval <- NULL
        auroc_pval <- NULL
        return_df <- NULL
        debug_df <- NULL
      }
      
      output <- switch(performance_metric, 
                       'accuracy' = original_accuracy,
                       'AUROC' = original_AUROC,
                       'acc_pval' = acc_pval, 
                       'AUROC' = auroc_pval, 
                       'all' = return_df,
                       'debug' = debug_df,
                       NULL
                       )
     return(output)
    }, error = function(err){
      print(paste0(sig$filterDescription, ' ', MIobj$formattedName, ' ', err))
      return(NA)
    })
}

permuteLabels <- function(MIobj){
  new_class <- sample(x = MIobj$class, size = length(MIobj$class), replace = F)
  names(new_class) <- names(MIobj$class)
  MIobj$class <- new_class
  if(checkDataObject(MIobj, 'Dataset')){
    return(MIobj)
  } else {
    stop('label permutation error')  
  }
  
}

getGeneFraction <- function(datasetObject, filterObject, positive = T){
  expr_genes <- datasetObject$keys
  filterObject$posGeneNames <- setdiff(filterObject$posGeneNames,'')
  filterObject$negGeneNames <- setdiff(filterObject$negGeneNames,'')
  pos_genes = intersect(filterObject$posGeneNames, expr_genes)
  neg_genes = intersect(filterObject$negGeneNames, expr_genes)
  
  if(positive){
    if(length(filterObject$posGeneNames) > 0){
      frac <- length(pos_genes) / length(filterObject$posGeneNames)
    } else {
      frac <- NA
    }
    return(frac)
  } else {
    if(length(filterObject$negGeneNames) > 0){
      frac <- length(neg_genes) / length(filterObject$negGeneNames)
    } else {
      frac <- NA
    }
    return(frac)
  }
}

plotDatSignature <- function(signature, data_list, performance_metric = 'AUROC', method = 'zScore', number_permutations = 1000, specific_method = 'geomMean', pval_calc = F, pval_perm = 1000, parallel = F){
  # return a data frame giving each study containing a class of interest a score for discovery, validation, and new
  studies <- sapply(data_list, function(x){return(x$formattedName)})
  
  
  
  plot_dat <- data.frame('Study' = studies, 
                         'Accession' = unlist(lapply(studies, parseAccession)), 
                         'Platform' = unlist(lapply(studies, parsePlatform)), 
                         'N' = unlist(lapply(data_list, function(x){min(table(x$class), na.rm = T)})),
                         'N_neg' = unlist(lapply(data_list, function(x){sum(x$class == 0, na.rm = T)})),
                         'N_pos' = unlist(lapply(data_list, function(x){sum(x$class == 1, na.rm = T)})),
                         'pos_genes' = paste(signature$posGeneNames, collapse = ' '),
                         'neg_genes' = paste(signature$negGeneNames, collapse = ' '),
                         'frac_pos' = sapply(data_list, getGeneFraction, filterObject = signature, positive = T),
                         'frac_neg' = sapply(data_list, getGeneFraction, filterObject = signature, positive = F),
                         'Type' = "New",
                         'Signature' = signature$filterDescription,
                         stringsAsFactors = F)
  if(!pval_calc){ # use the original method that returns the specified metric
    score_df <- data.frame('Study' = names(data_list), 
                           'Scores' = sapply(data_list, scoreStudy, 
                                             signature = signature, 
                                             performance_metric = performance_metric, 
                                             method = method, 
                                             n_perm = number_permutations,
                                             specific_method = specific_method), 
                           stringsAsFactors = F)
  } else { # use use an alternate function call that returns p-values as well
    if(parallel == T){
      score_list <- mclapply(data_list, 
                           scoreStudySignificant, 
                           signature = signature, 
                           performance_metric = performance_metric, 
                           method = method, 
                           n_perm = number_permutations,
                           specific_method = specific_method, 
                           mc.cores = detectCores()-1)
    } else {
      score_list <- lapply(data_list, 
                           scoreStudySignificant, 
                           signature = signature, 
                           performance_metric = performance_metric, 
                           method = method, 
                           n_perm = number_permutations,
                           specific_method = specific_method)
    }
    
    score_df <- bind_rows(score_list) %>%
      mutate('Study' = names(data_list))
  }
  
  # assign each study a label in {discovery, validation, new}
  validation_ix <- sapply(plot_dat$Accession, function(x){
    if(x %in% signature$validation){
      return(T)
    } else {
      return(F)
    }})
  discovery_ix <- sapply(plot_dat$Accession, function(x){
    if(x %in% signature$discovery){
      return(T)
    } else {
      return(F)
    }})
  plot_dat$Type[validation_ix] <- "Validation"
  plot_dat$Type[discovery_ix] <- "Discovery"
  out_df <- left_join(plot_dat, score_df, by = 'Study')
}


checkIndices <- function(ix = healthy_ix, df = p, column_ix = col_ix){
  p[ix, column_ix]
}


# checkIndices(healthy_ix)
# checkIndices(bacteria_ix)
# checkIndices(virus_ix)
# checkIndices(co_inf_ix)
# checkIndices(other_inf_ix)
# checkIndices(other_noninf_ix)
listFiles <- function(platform, directory){
  
  return(file_vector)
}

reorderE <- function(expr_mat, Key.Column, phenodat, accession = accession){
  current_order <- data.frame('label' = colnames(expr_mat), 'ix' = 1:ncol(expr_mat), 'expr_label' = colnames(expr_mat), stringsAsFactors = F)
  output_order <- data.frame('label' = phenodat[, Key.Column], 'pheno_label' = phenodat[, Key.Column], stringsAsFactors = F)
  
  if(Key.Column == 'geo_accession'){
    print('processing GSMs')
    current_order %<>% fixGSMs()
    output_order %<>% fixGSMs()
    out <- left_join(output_order, current_order, by = 'label') %>%
      filter(!is.na(ix))
  } else {
    if(sum(str_detect(pattern = '\\d{10}_[:alpha:]', string = current_order$label)) == nrow(current_order)){
      #print('ILMN chip lanes detected')
      current_order %<>% fixIlmnChannelIds()
      output_order %<>% fixIlmnChannelIds()
    } else {
      #print('attempting exact match')
      current_order %<>% fixCommonIssues()
      output_order %<>% fixCommonIssues()
    }
    
    out <- left_join(output_order, current_order, by = 'label') %>%
      filter(!is.na(ix))
    
    if(dim(out)[1] == 0){
      #print('attempting partial matches')
      out <- current_order
      #out$expr_label <- out$label
      for (i in 1:nrow(current_order)){
        #print(i)
        # appended $ or _ because we don't want e.g. 47 to match with 470
        new_ix <- grep(pattern = paste0(current_order$label[i], '$|^', current_order$label[i], '_|_', current_order$label[i], '_|\\[', current_order$label[i],'\\]'), x = output_order$label)
        if(length(new_ix) > 1){
          #print(paste(i, 'non-unique matches found'))
          out$ix[i] <- NA
          out$label[i] <- NA
        }
        if(length(new_ix) == 0){
          #print(paste(i, 'no match found'))
          out$ix[i] <- NA
          out$label[i] <- NA
        }
        if(length(new_ix) == 1){
          out$ix[i] <- new_ix
          out$label[i] <- output_order$label[new_ix]
        }
        
      }
    }
    
    out <- out %>%
      filter(!is.na(ix))
    
    if (dim(out)[1] == 0){
      #print('attempting field swap matching')
      out <- current_order
      #out$expr_label <- out$label
      out$matched <- NA
      fields <- strsplit(current_order$label, split = '_')
      
      # for (m in 1:length(fields)){
      #   current_field <- fields[[m]]
      #   keep_ix <- nchar(current_field) > 1
      #   current_field <- current_field[keep_ix]
      #   if(sum(!keep_ix) > 0){
      #     print('Fields are misbehaving')
      #   }
      #   fields[[m]] <- current_field
      # }
      
      for (i in 1:length(fields)){
        tmp_ix <- list(length = length(fields[[i]]))
        tmp_out <- rep(F, length(fields))
        current_field <- fields[[i]]
        for (j in 1:length(fields[[i]])){
          tmp_ix[[j]] <- grepl(pattern = paste0(current_field[j], '$|^', current_field[j], '_|_', current_field[j], '_'), x = output_order$label)
          
          if (j == 1){
            tmp_out <- tmp_ix[[j]]
          } else {
            tmp_out <- tmp_ix[[j]] & tmp_out
          }
          
          if(length(which(tmp_ix[[j]])) == 1){
            break
          }
        }
        new_ix <- which(tmp_out)
        if(length(new_ix) == 1){
          #print(paste(i, 'matched'))
          out$ix[i] <- new_ix
          out$matched[i] <- output_order$label[new_ix]
        } else {
          out$ix[i] <- NA
          #print('gdi')
        }
      }
      
      out <- out %>%
        filter(!is.na(ix))
      
      if(dim(out)[1] > 0){
        out$label <- out$matched
        out$matched <- NULL
      }
      
      if(max(table(out$ix)) > 1){
        print('non-unique matching error')
        error()
      }
      
    }
  }
  
  
  
  
  if (nrow(current_order) != nrow(output_order)){
    #print('warning: samples not aligned')
    out_mat <- expr_mat[, out$expr_label]
  } else {
    out_mat <- expr_mat[, out$ix]
  }
  
  colnames(out_mat) <- out$label
  
  if (ncol(out_mat) != ncol(expr_mat)){
    print(paste0("Couldn't match all columns: ", ncol(out_mat), " out of ", ncol(expr_mat), " returned."))
  } else {
    #print("Matched all samples")
  }
  
  if(ncol(out_mat) == 0){
    print('error: no samples matched')
    error()
  }
  
  if(length(out$pheno_label) == 0){
    out <- out %>% 
      left_join(output_order, by = 'label')
  }
  return(out)
}

fixCommonIssues <- function(set, accession = accession){
  str <- set$label
  str <- trimws(str)
  str <- str_replace(string = str, pattern = '^\\.X|^\\.|^\\.x|^X', replacement = '')
  #str <- str_replace(string = str, pattern = '\b', replacement = '_')
  str <- gsub(x = str, pattern = ' ', replacement = '_')
  str <- gsub(x = str, pattern = '-', replacement = '_')
  str <- gsub(pattern = '\\.', replacement = '_', x = str)
  str <- gsub(pattern = 'ctrl', replacement = 'control', x = str)
  #str <- str_replace(string = str, pattern = ' ', replacement = '_')
  str <- tolower(str)
  
  # this is specific to GSE 45919
  if(sum(grepl(pattern = 'pbmc_seronegative_subj', x = str)) > 0){
    print('GSE45919-specific fix')
    Sys.sleep(2)
    str <- gsub(pattern = 'pbmc_seronegative_subj', replacement = 'baseline', x = str)
  }
  
  # this is specific go GSE38900 GPL6884
  if(sum(grepl(pattern = '^#_|^__', x = str)) > 0){
    print('GSE38900-GPL6884-specific fix')
    Sys.sleep(2)
    str <- gsub(pattern = '^#', replacement = '', x = str)
    str <- gsub(pattern = '^_', replacement = '', x = str)
    str <- gsub(pattern = '^__', replacement = '', x = str)
    str <- gsub(pattern = '\\.', replacement = '', x = str)
    str <- gsub(pattern = '^X', replacement = '', x = str)
    str <- gsub(pattern = '^_', replacement = '', x = str)
  }
  
  if(sum(grepl(pattern = tolower('DefiniteB|ProbableB|DefiniteV|ProbableV'), x = str)) > 0){
    print('GSE72810-specific fix')
    Sys.sleep(2)
    str <- gsub(pattern = 'definiteb', replacement = 'definite_b', x = str)
    str <- gsub(pattern = 'definitev', replacement = 'definite_v', x = str)
    str <- gsub(pattern = 'probableb', replacement = 'probable_b', x = str)
    str <- gsub(pattern = 'probablev', replacement = 'probable_v', x = str)
  }
  
  set$label <- str
  return(set)
}

fixIlmnChannelIds <- function(set){
  set <- fixCommonIssues(set)
  set$label <- str_extract(string = set$label, pattern = '\\d{10}_[:alpha:]')
  return(set)
}

fixGSMs <- function(set){
  set$label <- sapply(set$label, parseGSM)
  return(set)
}

matchLabels <- function(matching_set = current_order, reference_set = output_order){
  matching_set$label <- str_extract(string = matching_set$label, pattern = '\\d{10}_[:alpha:]')
  reference_set$label <- str_extract(string = reference_set$label, pattern = '\\d{10}_[:alpha:]')
}

createAnnotatedDataFrame <- function(df){
  varmeta <- data.frame('var' = colnames(df), 'labelDescription' = '', stringsAsFactors = F)
  ann_df <- AnnotatedDataFrame(df, varmeta)
  return(ann_df)
}

parseGSM <- function(str){
  str <- gsub(pattern = '^\\.', replacement = '', x = str)
  str <- strsplit(x = str, split = '\\.|_')[[1]][1]
  return(str)
}

sortGEOdl <- function(platform, master_file){
  rds_file_dir <- '~/Data/DARPA/GEO_dl/'
  GEO_dl_files <- list.files('~/Data/DARPA/GEO_dl')
  platform_file <- master_file %>%
    filter(tolower(trimws(Suitable)) != 'no') %>%
    filter(grepl(pattern = platform, x = Platform))
  straightforward_studies <- platform_file %>%
    filter(tolower(trimws(Platform)) == tolower(platform))
  
  files_to_copy[[1]] <- sapply(straightforward_studies$Accession, function(accession){
    file_ix <- which(grepl(pattern = accession, x = GEO_dl_files))
    return(GEO_dl_files[file_ix])
  })
  
  non_straightforward_studies <- platform_file %>%
    filter(tolower(trimws(Platform)) != tolower(platform))
  
  files_to_copy[[2]] <- sapply(X = non_straightforward_studies$Accession, FUN = function(accession){
    file_ix <- which(grepl(pattern = accession, x = GEO_dl_files) & grepl(pattern = platform, x = GEO_dl_files))
    return(GEO_dl_files[file_ix])
  })
  
  length(unlist(files_to_copy)) == dim(platform_file)[1]
  for(rds_file in files_to_copy){
    file.rename(paste0(rds_file_dir, rds_file), to = paste0(rds_file_dir, platform, '/', rds_file))
  }
}

fixEsetExpr <- function(eset){
  exprMat <- exprs(eset)
  ix <- which(rowSums(is.na(exprMat)) > 0)
  if(length(ix) > 0){
    if (length(ix)/nrow(eset) < 0.5){
      print('fixing rows')
      exprMat <- exprMat[-ix, ]
      eset <- ExpressionSet(assayData = exprMat, 
                            phenoData = phenoData(eset), 
                            featureData = createAnnotatedDataFrame(fData(eset)[-ix,]), 
                            annotation = eset@annotation)
    } else {
      print('fixing columns rather than rows')
      ix <- which(colSums(is.na(exprMat)) > 0)
       if (length(ix) > 0 & length(ix) < (ncol(exprMat)/5)){
         print('removing columns')
         exprMat <- exprMat[, -ix]
         eset <- ExpressionSet(assayData = exprMat, 
                               phenoData = phenoData(eset)[-ix,], 
                               featureData = featureData(eset), 
                               annotation = eset@annotation)
       } else {
         print('ERROR, returning original eset')
         #ix <- which(rowSums(is.na(exprMat)) > 0)
         #exprMat <- exprMat[-ix, ]
         exprMat <- exprs(eset)
         eset <- ExpressionSet(assayData = exprMat, 
                               phenoData = phenoData(eset), 
                               featureData = createAnnotatedDataFrame(fData(eset)), 
                               annotation = eset@annotation)
       
       }
    }
  }
  
  return(eset)
}

fixExprLog <- function(expr_mat){
  if(max(max(expr_mat, na.rm = T), na.rm = T) > 100){
    print('non-log adjusted data detected. log transforming...')
    if(min(min(expr_mat, na.rm = T), na.rm = T) < 1){
      print('negative values detected. adjusting before log transform')
      expr_mat <- 1 + abs(min(min(expr_mat, na.rm = T), na.rm = T)) + expr_mat
    }
    expr_mat <- log2(expr_mat)
  }
  return(expr_mat)
}

MItoEset <- function(MIobj){
  eObj <- ExpressionSet(assayData = MIobj$expr,
                        phenoData = createAnnotatedDataFrame(MIobj$pheno), 
                        annotation = parsePlatform(MIobj$formattedName))
  return(eObj)
}

checkSignatureGenesInData <- function(dataObj, signature){
  pos_flag <- neg_flag <- T
  pos_genes <- signature$posGeneNames
  pos_genes <- pos_genes[pos_genes != '']
  neg_genes <- signature$negGeneNames
  neg_genes <- neg_genes[neg_genes != '']
  if (length(pos_genes) > 0){
    pos_flag <- mean(pos_genes %in% dataObj$keys) >= 0.5
  }
  if (length(neg_genes) > 0){
    neg_flag <- mean(neg_genes %in% dataObj$keys) >= 0.5
  }
  if (pos_flag & neg_flag){
    return (T)
  } else {
    return(F)
  }
}


filterMIobj <- function(MIobj, index){
  tmpMIobj <- MIobj
  tmpMIobj$expr <- tmpMIobj$expr[, index]
  tmpMIobj$pheno <- tmpMIobj$pheno[index, ]
  tmpMIobj$class <- tmpMIobj$class[index]
  
  if(checkDataObject(tmpMIobj, 'Dataset')){
    return(tmpMIobj)
  } else {
    stop()
  }
}

checkExpressionScaling <- function(MIobj){
  
}


calculateScoreRobust <- function(filterObject, datasetObject, suppressMessages = F, method = 'geomMean', zScore = T){
  if(method == 'original'){
    totalScore = calculateScore(filterObject = filterObject, datasetObject = datasetObject, suppressMessages = suppressMessages)
  } else {
    # method options: geomMean, geomMean_exp2, mean, mean_exp2
    # assumes data is being entered on a log2 scale
    
    datasetObjectmin <- min(datasetObject$expr, na.rm = TRUE)
    if (datasetObjectmin < 0) {
      datasetObject$expr <- datasetObject$expr + abs(datasetObjectmin) + 1
    }
    
    # pull list of genes present in data set
    expr_genes <- datasetObject$keys
    filterObject$posGeneNames <- setdiff(filterObject$posGeneNames,'')
    filterObject$negGeneNames <- setdiff(filterObject$negGeneNames,'')
    pos_genes = intersect(filterObject$posGeneNames, expr_genes)
    neg_genes = intersect(filterObject$negGeneNames, expr_genes)
    
    # for reporting fraction of genes used
    if (!suppressMessages){
      N_pos_str <- paste0(length(pos_genes), ' of ', length(filterObject$posGeneNames))
      N_neg_str <- paste0(length(neg_genes), ' of ', length(filterObject$negGeneNames))
      print(paste0('Used ', N_pos_str, ' pos genes and ', N_neg_str, ' neg genes'))
    }
    
    
    if(length(pos_genes) == 0 & length(neg_genes) == 0){
      if (!suppressMessages){
        print('no common genes between data and signature')
      }
      return(NA)
    }
    
    posScore = .calculate_scores(datasetObject, pos_genes, method = method)
    negScore = .calculate_scores(datasetObject, neg_genes, method = method)
    
    if(grepl(pattern = 'exp2', x = method)){
      totalScore = posScore/negScore
    } else {
      totalScore = posScore - negScore
    }
    
    
    if (sum(abs(totalScore), na.rm = T) != 0){
      if(zScore){
        totalScore = as.numeric(scale(totalScore))
      } else {
        totalScore = as.numeric(totalScore)
      }
      
    }
  }
  
  return(totalScore)
}


.calculate_scores = function(datasetObject, genes, method = 'geomMean'){
  
  # find the genes to average
  genes_idx = which(datasetObject$keys %in% genes)
  
  # if there aren't any genes, return 0
  if (length(genes_idx) == 0){
    if(method == 'mean'){
      output_val = 1
    } else {
      output_val = 0
    }
    
    return(rep(output_val, ncol(datasetObject$expr)))
  }
  
  # filter expr to sig genes
  gene_expr = datasetObject$expr[genes_idx,] 
  
  #in case only one gene is present, take the vector itself
  if (is.null(nrow(gene_expr))){
    sample_scores = gene_expr
  }
  
  if (!is.null(nrow(gene_expr))){
    sample_scores = apply(gene_expr, 2, function(method, x){switch(method, 
                                                                   geomMean = geom_mean(x),
                                                                   geomMean_exp2 = geom_mean_exp2(x),
                                                                   mean = mean(x, na.rm = T),
                                                                   mean_exp2 = mean_exp2(x))
    }, method = method)}

  return(sample_scores)
}

geom_mean = function(x){
  return(exp(mean(log(x), na.rm = T)))
}

geom_mean_exp2 = function(x){
  return(exp(mean(x, na.rm = T)))
}

mean_exp2 <- function(x){
  return(log(mean(exp(x), na.rm = T)))
}


auroc = function(bool, score) {
  if (length(unique(score)) == 1){
    print('only one value of score: setting AUC to NA')
    return(NA)
  } 
  n1 = sum(!bool)
  n2 = sum(bool)
  U  = sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

virtualFPR <- function(MIobj, sig, quantile_value = 0.95, plot = F, two_sided = F, fraction = T){
  df <- data.frame('Class' = as.factor(MIobj$class), 
                   'Score' = calculateScoreRobust(sig, MIobj, method = 'geomMean', suppressMessages = T), stringsAsFactors = F)
  background_score_df <- df %>% 
    filter(Class == 0)
  background_scores <- background_score_df$Score
  
  if(two_sided){
    # get two scores and put them in increasing order
    background_cutoff <- quantile(background_scores, probs = c(1-quantile_value, quantile_value), na.rm = T)
    background_cutoff <- sort(background_cutoff, decreasing = F)
  } else {
    background_cutoff <- quantile(background_scores, probs = quantile_value, na.rm = T)
  }
  
  
  exposure_score_df <- df %>%
    filter(Class == 1)
  exposure_scores <- exposure_score_df$Score
  if (fraction){
    if(two_sided){
      fpr <- mean(exposure_scores >= background_cutoff[2] | exposure_scores <= background_cutoff[1], na.rm = T)
    } else {
      fpr <- mean(exposure_scores >= background_cutoff, na.rm = T)
    }  
  } else {
    if(two_sided){
      fpr <- sum(exposure_scores >= background_cutoff[2] | exposure_scores <= background_cutoff[1], na.rm = T)
    } else {
      fpr <- sum(exposure_scores >= background_cutoff, na.rm = T)
    }
  }
  
  
  
  if(plot == T){
    if (two_sided){
      plot_dat <- df %>%
        mutate('Label' = '0') %>%
        mutate(Label = ifelse((Class == 1 & (Score >= background_cutoff[2] | Score <= background_cutoff[1])), '1', Label))
      plot(ggplot(plot_dat, aes(x = Class, y = Score)) + 
             geom_boxplot(outlier.shape = NULL) + 
             geom_point(aes(color = Label)) + 
             geom_hline(yintercept = background_cutoff) + 
             labs(title = MIobj$formattedName,
                  subtitle = paste0(gsub(pattern = '\n', 
                                         replacement = ' ', 
                                         x = sig$filterDescription),
                                    '\nPositive Class Rate: ', round(fpr,3)))
      )
    } else {
      plot_dat <- df %>%
        mutate('Label' = '0') %>%
        mutate(Label = ifelse((Class == 1 & Score >= background_cutoff), '1', Label))
      plot(ggplot(plot_dat, aes(x = Class, y = Score)) + 
             geom_boxplot(outlier.shape = NULL) + 
             geom_point(aes(color = Label)) + 
             geom_hline(yintercept = background_cutoff) + 
             labs(title = MIobj$formattedName,
                  subtitle = paste0(gsub(pattern = '\n', 
                                         replacement = ' ', 
                                         x = sig$filterDescription),
                                    '\nPositive Class Rate: ', round(fpr,3)))
      )
    }
    
  }
  
  return(fpr)
}

calculatePval <- function(labels, values, true_value, n_perm = 1000){
  aucs <- vector(mode = 'double', length = n_perm)
  
  if(sum(is.na(values)) > 0){
    print('missing scores')
    return(NA)
  }
  
  for (iteration in 1:n_perm){
    # if(iteration %% 1000 == 0){
    #   print(paste0(iteration, '/', n_perm))
    # }
    sample_shuffle_ix <- sample(x = 1:length(labels), size = length(labels), replace = F)
    #aucs[iteration] <- pROC::auc(response = labels[sample_shuffle_ix], predictor = values)
    aucs[iteration] <- calculateROC(labels = labels[sample_shuffle_ix], predictions = values, AUConly = T)
  }
  pval <- mean(aucs >= true_value)
  
  #hist(aucs, 25); abline(v = true_value, col = 'red')
  
  return(pval)
}


calculateLogitScores <- function(MIobj, sig, performance_metric, N_PERM = 1000, parallel = F){
  if(parallel){
    
    if (nrow(MIobj$pheno) < 50){
      cl <- makePSOCKcluster(6)
      #print('6')
    } else if (nrow(MIobj$pheno) < 100){
      cl <- makePSOCKcluster(2)
      #print('2')
    } else if(nrow(MIobj$pheno) < 500){
      cl <- makePSOCKcluster(1)
      #print('1')
    } else {
      cl <- makePSOCKcluster(1)
      #print('1')
    } 
    
    registerDoParallel(cl)
  } else {
    #cl <- makePSOCKcluster(1)
    #registerDoParallel(cl)
  }
  
  out <- tryCatch(
    {
      #print(performance_metric)
      # performance_metric should be one of {accuracy, acc_pval, AUROC, AUROC_pval, debug}
      
      # add a check that the labels are approximately balanced, otherwise we have to add that
      # if the labels are unbalanced, you need nested cross validation because you want to estimate the effect of downsampling (SMOTE)
      # on performance
      
      # considerations
      ## can you enforce feature signs in caret? let's start with just concat feature lists and go from there
      
      # performance metrics
      ## can report logit LOOCV accuracy at a cutoff of 0.5
      ## can report p-value of LOOCV accuracy against permuted class labels 
      ## can report AUROC from hold out samples; if you're generating probabilities for each sample, these are directly comparable for generating a ROC
      
      # set a flag in case we can avoid permutation
      perm_flag <- T
      
      # filter MIobj expr based on signature genes
      sig_genes <- setdiff(c(sig$posGeneNames, sig$negGeneNames), '') # concat pos and neg genes, rm empty set
      sig_gene_ix <- rownames(MIobj$expr) %in% sig_genes
      MIobj$mat <- MIobj$expr[sig_gene_ix, ]
      MIobj$mat <- t(MIobj$mat)
      
      # check for NA expression values
      if(sum(is.na(MIobj$mat)) > 0){
        sample_rm_ix <- apply(X = MIobj$mat, MARGIN = 1, FUN = function(x){return(sum(is.na(x)))}) > 0
        
        if(mean(sample_rm_ix) == 1){
          print(paste0('all samples have NA values ', MIobj$formattedName))
          return(NA)
        }
        
        MIobj <- filterMIobj(MIobj, !sample_rm_ix)
        MIobj$mat <- MIobj$mat[!sample_rm_ix,]
      }
      
      # you cannot fit a linear model with more predictors than observations, and we do LOOCV so this is N + 1 > M
      if(ncol(MIobj$mat) > (nrow(MIobj$mat) - 1)){
        print(paste0('More params than observations for ', MIobj$formattedName))
        return(NA)
      }
      
      if(min(table(MIobj$class)) == 1){
        print(paste0('Not enough observations for ', MIobj$formattedName))
        return(NA)
      }
      
      
      # train LOOCV logit model
      original_model <- train(x = MIobj$mat, 
                              y = as.factor(paste0('class_', MIobj$class)), 
                              method = 'vglmAdjCat', 
                              trControl = trainControl(method = 'LOOCV', 
                                                       savePredictions = T, 
                                                       classProbs = T))
      
      orig_results_df <- original_model$pred %>%
        dplyr::filter(parallel == F)
      
      output <- list('pred' = orig_results_df$class_1, 'labels' = MIobj$class)
      return(output)
    }, error = function(err){
      print(paste0(sig$filterDescription, ' ', MIobj$formattedName, ' ', err))
      return(NA)
    })
}

# this function takes the signatures data frame and swaps virus/bacteria labels 
# in both the positive.class and negative.class columns. 
specificityConversion <- function(signatures){
  #table(signatures$Positive.Class)
  #table(signatures$Negative.Class)
  
  tmp_pos_class <- sigClassSwap(signatures$Positive.Class)
  tmp_neg_class <- sigClassSwap(signatures$Negative.Class)
  
  #table(tmp_neg_class)
  #table(tmp_pos_class)
  
  # assign and return
  signatures$Positive.Class <- tmp_pos_class
  signatures$Negative.Class <- tmp_neg_class
  return(signatures)
}

sigClassSwap <- function(class_vector){ # take signatures $positive.class or $negative.class and swap bacteria and virus labels
  tmp_pos_class <- gsub(class_vector, pattern = 'Virus', replacement = 'placeholder')
  tmp_pos_class <- gsub(tmp_pos_class, pattern = 'Bacteria', replacement = 'Virus')
  tmp_pos_class <- gsub(tmp_pos_class, pattern = 'placeholder', replacement = 'Bacteria')
  return(tmp_pos_class)
}
