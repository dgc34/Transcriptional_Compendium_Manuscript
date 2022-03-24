library(dplyr)
library(tidyr)
library(ggplot2)
library(MetaIntegrator)
library(openxlsx)
library(scales)
library(ggh4x)

# load helper functions
source(paste0(Sys.getenv('DARPA'), 'helperFunctions.R'))

# load data
data_list <- readRDS('~/Documents/darpa-manuscript-data/noninfectious_data_list.RDS')
#data_list <- readRDS('~/Dropbox/aging_data_list.RDS')

exposure_df <- data.frame('Study' = names(data_list), 
                          'Exposure' = sapply(data_list, function(x){
                                                            return(unique(x$pheno$Standardized.Exposure))
                                                          }),
                          stringsAsFactors = F)

# filter down to just aging and obesity
exposure_ix <- exposure_df$Study[exposure_df$Exposure %in% c('Obesity', 'Aging')]
data_list <- data_list[exposure_ix]

signatures_file <- '~/Dropbox/compendium_manuscript/tables/supp_1.xlsx'
signatures <- read.xlsx(signatures_file)
signatures <- signatures %>%
  mutate(Positive.Class = str_extract(pattern = '.*vs.', string = Comparison)) %>% 
  mutate(Negative.Class = str_extract(pattern = 'vs.*', string = Comparison)) %>%
  mutate(Positive.Class = gsub(pattern = ' vs.', replacement = '', x = Positive.Class)) %>% 
  mutate(Negative.Class = gsub(pattern = 'vs. ', replacement = '', x = Negative.Class)) %>%
  mutate(Positive.Class = ifelse(Positive.Class == 'Influenza', 'Virus', Positive.Class)) # this will let us compute bacterial specificity for influenza signatures

#data_list <- readRDS('~/Dropbox/data_list.RDS')
metric <- 'AUROC'

data_list$GSE46449_GPL570 <- NULL

for (study in seq_along(data_list)){
	current_study <- data_list[[study]]
	print(current_study$formattedName)

	tmp_class <- rep(NA, length(current_study$pheno$Class))
	pos_ix <- current_study$pheno$Class == 'Other.NonInfectious'
	#current_study$pheno$Class[pos_ix] <- 'Viral'
	neg_ix <- grepl(pattern = 'Healthy', x = current_study$pheno$Class)
	#current_study$pheno$Class[neg_ix] <- 'Bacterial'

	tmp_class[pos_ix] <- 1
	tmp_class[neg_ix] <- 0
	names(tmp_class) <- rownames(current_study$pheno)

	ix <- !is.na(tmp_class)
	tmp_class <- tmp_class[ix]
	current_study$pheno <- current_study$pheno[ix,]
	current_study$expr <- current_study$expr[,ix]
	current_study$class <- tmp_class

	if(checkDataObject(current_study, objectType = 'Dataset')){
		data_list[[study]] <- current_study
	} else {
		print(paste("ERROR", current_study$formattedName))
	}
}

out <- list()
#for (signature_ID in c(1:2,3:10)){
for (signature_ID in 1:nrow(signatures)){
	print(paste0('signature', signature_ID))

  if(!signatures[signature_ID,'Positive.Class'] %in% c("Bacteria", "Virus", "Non-Infectious", "Resp.Virus", "Resp.Bact")){
    next
  } else {
    signatures[signature_ID,'Positive.Class'] <- gsub(pattern = 'Resp.Virus', replacement = 'Virus', x = signatures[signature_ID,'Positive.Class'])
    signatures[signature_ID,'Positive.Class'] <- gsub(pattern = 'Resp.Bact', replacement = 'Bacteria', x = signatures[signature_ID,'Positive.Class'])
  }
  #short_data_list <- data_list[1:25]
	tmp <- data_list
  #tmp <- lapply(X = data_list, function(x){
  #  parseClass(x, signature_row = signatures[signature_ID,])
  #})
  # remove studies that don't have multiple classes
  ix <- sapply(X = tmp, filterClass)
  if(mean(ix) != 1){
    stop('problem with class labeling')
  }
  tmp <- tmp[ix]
  if(length(tmp) == 0){
    next
  }
  # remove studies that don't have all signature genes
  sig <- sig2Meta(signatures[signature_ID,])
  ix <- sapply(X = tmp, filterGenes, signature = sig, threshold = 0.25)
  tmp <- tmp[ix]
  if(mean(ix) != 1){
    #stop('problem with genes')
  }
  if(length(tmp) == 0){
    next
  }
  
  sapply(tmp, function(x){
    print(x$formattedName);
    return(checkDataObject(x, objectType = 'Dataset'))
  })

  current_sig <- sig2Meta(signatures[signature_ID,])
  df <- plotDatSignature(signature = current_sig, data_list = tmp, performance_metric = metric, pval_calc = T)

  out[[signature_ID]] <- df
}
saveRDS('~/Documents/darpa-manuscript-data/figures/aging_obesity_heatmap_012522.RDS', object = out)
out <- readRDS('~/Documents/darpa-manuscript-data/figures/aging_obesity_heatmap_012522.RDS')

plot_dat <- bind_rows(out)
plot_dat[grepl(pattern = 'GSE65219', x = plot_dat$Study) & grepl(pattern = 'Sampson', x = plot_dat$Signature), 'Type'] <- 'Discovery'

plot_dat <- plot_dat %>%
  mutate('Type' = factor(Type, levels = c('Discovery', 'Validation', 'New', 'Null'), ordered = T)) %>%
  mutate('Comparison' = str_extract(pattern = 'VxB|V|B', string = Signature))  %>%
  arrange(Comparison) %>%
  mutate(Signature = factor(Signature, levels = unique(Signature), ordered = T))

table(plot_dat$Comparison)

plot_dat <- plot_dat %>%
  left_join(exposure_df)

N_threshold = 0
colors <- createColorScale()

# facet by pathogen and cluster rows/columns
plot_list <- list()
head(plot_dat)
table(plot_dat$Exposure.Class)
exposures <- unique(plot_dat$Comparison)

aodf <- getSignificantDf(plot_dat %>% dplyr::rename(Pathogen = Exposure)) %>%
  mutate(Comparison = str_extract(pattern = 'I|VxB|V|B', string = Signature)) %>%
  mutate(Comparison = factor(Comparison, levels = c('VxB', 'V', 'B'), ordered = T))

f5d <- ggplot(aodf %>%
         filter(Comparison %in% c('VxB', 'V', 'B')), aes(x = Signature, y = Pathogen, fill = Specific)) + 
  #geom_tile(alpha = 0.1, fill = 'white', color = 'black') +
  geom_tile(alpha = 0.1) +
  #geom_tile() +
  geom_text(aes(label = Slabel), size = 15, position = position_nudge(y = 0.15)) + 
  geom_vline(xintercept = (0:length(unique(aodf$Signature)))+0.5, color = 'gray95') + 
  geom_hline(yintercept = (0:length(unique(aodf$Pathogen)))+0.5) + 
  #scale_fill_manual(values = color_bar) + 
  #scale_fill_gradient2(low = scales::muted('blue'), high = scales::muted('red'), mid = 'white', midpoint = 0.5) + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) + 
  facet_grid(. ~ Comparison, scales = 'free', space = 'free') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

plot_list <- list(f5c, f5d)
for (i in 1:length(plot_list)){
  current_plot <- plot_list[[i]]
  current_plot <- egg::set_panel_size(current_plot, height=unit((length(unique(current_plot$data$Pathogen))/2), "cm"),
                                      width=unit((length(unique(current_plot$data$Signature))/1.5), "cm") )
  plot_list[[i]] <- current_plot
}
cowplot::plot_grid(f5c, f5d, nrow=2)


combined_plot_dat <- f5c$data %>%
  mutate(Signature = as.character(Signature)) %>%
  bind_rows(f5d$data %>% mutate(Signature = as.character(Signature))) %>%
  mutate('Infectious' = ifelse(Pathogen %in% c('Obesity', 'Aging'), F, T)) %>%
  mutate(Signature = factor(Signature, levels = rev(sig_order), ordered = T))

ggplot(combined_plot_dat %>%
         filter(Comparison %in% c('VxB', 'V', 'B')), aes(x = Signature, y = Pathogen, fill = Specific)) + 
  #geom_tile(alpha = 0.1, fill = 'white', color = 'black') +
  geom_tile(alpha = 0.1) +
  #geom_tile() +
  geom_text(aes(label = Slabel), size = 15, position = position_nudge(y = 0.75)) + 
  geom_vline(xintercept = (0:length(unique(sig_dat$Signature)))+0.5, color = 'gray95') + 
  geom_hline(yintercept = (0:length(unique(sig_dat$Pathogen)))+0.5) + 
  #scale_fill_manual(values = color_bar) + 
  #scale_fill_gradient2(low = scales::muted('blue'), high = scales::muted('red'), mid = 'white', midpoint = 0.5) + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) + 
  facet_grid(Infectious ~ Comparison, scales = 'free', space = 'free') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = '')
