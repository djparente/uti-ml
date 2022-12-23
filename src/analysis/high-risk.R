getHighRiskStratification <- function(dataset, premut=NULL) {
  # Specify sort orders
  vkey_order <- c(
    'fever_subjective',
    'temp',
    'hr',
    'sbp',
    'ams',
    'flank_pain',
    'cva',
    'vomiting',
    'catheter',
    'immunocompromised',
    'pregnant'
  )
  vvar_order <- c('1', '0', '[18,25]', '(25,35]', '(35,45]', '(45,55]', '(55,65]', '(65,75]', '(75,999]', 'Male', 'Female', 'negative', 'positive', 'small', 'none', 'few', 'moderate', 'large', '4+', 'clear', 'not_clear', 'Yes', 'No', 'other', 'many', 'marked', 'not_reported')
  
  
  if( !is.null(premut) ) {
    dataset <- premut(dataset)
  }
  
  # Get the stratification table
  res <- dataset %>%
    dplyr::mutate(
      age = cut(dataset$age, breaks=c(18, 25, 35, 45, 55, 65, 75, 999), include.lowest=T)
      #race = races_final[match(dataset$race, races_initial)],
      #ethnicity = ethnicity_final[match(dataset$ethnicity, ethnicity_initial)]
    ) %>%
    pivot_longer(c(
      fever_subjective,
      temp,
      hr,
      sbp,
      ams,
      flank_pain,
      cva,
      vomiting,
      catheter,
      immunocompromised,
      pregnant
    ), names_to='vkey', values_to='vvar') %>%
    dplyr::group_by(split, UCX_abnormal, vkey, vvar) %>%
    dplyr::summarize(count=n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(split, UCX_abnormal, vkey) %>%
    dplyr::mutate(
      count = prettyNPct(count, count/sum(count)) 
    ) %>%
    ungroup() %>%
    pivot_wider(names_from=c('split', 'UCX_abnormal'), values_from=c('count'), values_fill='-') %>%
    arrange(
      match(vkey, vkey_order),
      match(vvar, vvar_order)
    ) %>%
    
    dplyr::ungroup()
  
  return(res)
}

simTrial <- function(pop, roc, prob_name) {
  threshValue <- 0.85
  thresh85Trial <- tail(dplyr::filter(coords(roc, "all", tr=F), sensitivity > threshValue), 1)$threshold
  
  above <- dplyr::filter(pop, !!sym(prob_name) >= thresh85Trial)
  below <- dplyr::filter(pop, !!sym(prob_name) < thresh85Trial)
  aboveRes <- above %>% group_by(rx_abx, UCX_abnormal) %>% summarize(count=n()) %>% arrange(rx_abx, desc(UCX_abnormal))
  belowRes <- below %>% group_by(rx_abx, UCX_abnormal) %>% summarize(count=n()) %>% arrange(rx_abx, desc(UCX_abnormal))
  
  return(list(above=aboveRes, below=belowRes, aboveSet=above, belowSet=below));
}

pc_hr <- getHighRiskStratification(primarycare_dat, premut=premutate.PC)
write.table(pc_hr, 'pc_highrisk_strat.txt', sep="\t", row.names=F, quote=F)


pc_nohrf <- dplyr::filter(primarycare_dat, anyhrf == F)
pc_nohrf_above <- dplyr::filter(pc_nohrf, prob_xgb >= thresh85)
pc_nohrf_below <- dplyr::filter(pc_nohrf, prob_xgb < thresh85)


pc_hrf <- dplyr::filter(primarycare_dat, anyhrf == T)
pc_hrf_above <- dplyr::filter(pc_hrf, prob_xgb >= thresh85)
pc_hrf_below <- dplyr::filter(pc_hrf, prob_xgb < thresh85)


trial_xgb <- simTrial(pc_nohrf, ext.rocs$pc$xgb, 'prob_xgb')
trial_rf  <- simTrial(pc_nohrf, ext.rocs$pc$rf,  'prob_rf')
trial_ann <- simTrial(pc_nohrf, ext.rocs$pc$ann, 'prob_ann')