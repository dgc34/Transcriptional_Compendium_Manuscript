
selectCols <- function(df){
  df %>%
    filter(N > 4) %>% 
    dplyr::select(Accession, Platform, Signature, Type, Scores, N_neg, N_pos) %>%
    dplyr::rename(AUROC = Scores, N_case = N_pos, N_control = N_neg) %>% 
    mutate(Type = ifelse(Type == 'New', 'Validation', Type)) %>% 
    return()
}


auroc_vb_healthy <- readRDS('~/Dropbox/big_eval_heatmap_list_generalizability_021322_healthy_sbj_ts.RDS') %>%
  bind_rows() %>%
  selectCols() %>% 
  filter(!grepl(pattern = '^I', x = Signature))

auroc_vb_noninf <- readRDS('~/Dropbox/big_eval_heatmap_list_generalizability_021322_noninf_sbj_ts.RDS') %>%
  bind_rows() %>%
  selectCols() %>% 
  filter(!grepl(pattern = '^I', x = Signature))

XR_vb_healthy <- readRDS('~/Dropbox/big_eval_heatmap_list_specificity_021322_healthy_sbj_ts.RDS') %>%
  bind_rows() %>%
  selectCols() %>% 
  filter(!grepl(pattern = '^I', x = Signature))

XR_vb_noninf <- readRDS('~/Dropbox/big_eval_heatmap_list_specificity_021322_noninf_sbj_ts.RDS') %>%
  bind_rows() %>%
  selectCols() %>% 
  filter(!grepl(pattern = '^I', x = Signature))

output_list <- list('robustness (vs. healthy)' = auroc_vb_healthy,
                    'robustness (vs. non-inf.)' = auroc_vb_noninf,
                    'cross-reactivity (vs. healthy)' = XR_vb_healthy,
                    'cross-reactivity (vs. non-inf.)' = XR_vb_noninf)

output_file <- '~/Dropbox/compendium_manuscript/tables/supp4.xlsx'
write.xlsx(x = output_list, file = output_file, overwrite = T)

auroc_flu <- readRDS('~/Documents/darpa-manuscript-data/flu_eval_heatmap_v3_withV10_healthy.RDS') %>%
  bind_rows() %>%
  filter(Pathogen == 'Influenza') %>%
  bind_rows() %>%
  selectCols()

auroc_nonflu <- readRDS('~/Documents/darpa-manuscript-data/flu_eval_heatmap_v3_withV10_healthy.RDS') %>%
  bind_rows() %>%
  filter(Pathogen == 'Non-Influenza Virus') %>%
  bind_rows() %>%
  selectCols()

output_file <- '~/Dropbox/compendium_manuscript/tables/supp5_flu.xlsx'
write.xlsx(x = list('influenza robustness' = auroc_flu, 'influenza cross-reactivity' = auroc_nonflu), file = output_file, overwrite = T)
