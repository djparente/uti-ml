accumulateDeciles <- function(valid, truth, preds, cuts=seq(0, 1, 0.1)) {
  dat <- data.frame(truth=valid[[truth]], pred=valid[[preds]])
  
  # Quantize the predictions into the cut point bins
  dat$quantPred <- cut(dat$pred, seq(0, 1, 0.1), include.lowest=T, right=F)
  
  # Quantize the truth into 0 or 1. This has the useful property that mean(quantTruth) over some range
  # of interest (i.e. datapoints falling into a quantPred bins) is the percentage of deleterious variants
  # in the bin
  dat$quantTruth <- ifelse(dat$truth == 'yes', 1, 0)
  
  # Aggregate means and medians, stratified by quantPred
  ag <- aggregate(. ~ quantPred, dat, FUN=function(x) c(mean=mean(x), median=median(x)))
  
  
  simple <- data.frame(
    truth=ag$quantTruth[,1],
    pred=ag$pred[,1]
  )
  
  return(simple)
}

characterizeDeciles <- function(valid, preds) {
  linmod <- lm(accumulateDeciles(valid, 'UCX_abnormal', preds))
  slm <- summary(linmod)
  
  rsq <- slm$r.squared[1]
  b0 <- linmod$coefficients[1]
  b1 <- linmod$coefficients[2]
  
  return(c(b0=b0, b1=b1, rsq=rsq))
}


boot.deciles <- function(probe, data, indices, pb) {
    d <- data[indices,]
    pb$tick()
    return(characterizeDeciles(d, probe))
}

# Calculates the RMSE of a model including error from boostrap analysis (95% CI; lci = 2.5% bound, uci = 97.5% bound)
decileModelWithError <- function(dataset, probe, setname) {
  R=2000
  
  pb <- progress_bar$new(total=R+1)
  
  # Calculate bootstrap object
  bootobj <- boot(dataset, statistic=boot.deciles, R=R, strata=dataset$UCX_abnormal, probe=probe, parallel="multicore", pb=pb)
  
  resTibble <- tidy(bootobj, conf.int=T, conf.method="perc")
  
  res <- resTibble %>%
    mutate(estimate=fmt.ParameterErrorEstimates(statistic, conf.low, conf.high)) %>%
    dplyr::select(c('term', 'estimate'))
  
  
  metadata <- data.frame(rbind(
    c(term='dataset', estimate=setname),
    c(term='probe', estimate=probe)
  ))
  
  res <- rbind(metadata, res)
  
  
  res <- res %>%
    pivot_wider(names_from=term, values_from=estimate)
  
  return(res)
  
  
}


getDecileTable <- function() {
  deciletable <- decileModelWithError(urine_val, 'prob_xgb', 'ed')
  deciletable <- rbind(deciletable, decileModelWithError(urine_val, 'prob_ann', 'ed'))
  deciletable <- rbind(deciletable, decileModelWithError(urine_val, 'prob_rf', 'ed'))
  deciletable <- rbind(deciletable, decileModelWithError(urine_val, 'prob_micro', 'ed'))
  deciletable <- rbind(deciletable, decileModelWithError(urine_val, 'prob_rand', 'ed'))
  deciletable <- rbind(deciletable, decileModelWithError(primarycare_dat, 'prob_xgb', 'pc'))
  deciletable <- rbind(deciletable, decileModelWithError(primarycare_dat, 'prob_ann', 'pc'))
  deciletable <- rbind(deciletable, decileModelWithError(primarycare_dat, 'prob_rf', 'pc'))
  deciletable <- rbind(deciletable, decileModelWithError(primarycare_dat, 'prob_rand', 'pc'))
  return(deciletable)
}


