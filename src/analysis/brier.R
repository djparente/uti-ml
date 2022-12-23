brierpc_frame <- primarycare_dat %>% select(c('UCX_abnormal', 'prob_xgb', 'prob_rf', 'prob_ann', 'prob_rand'))
brierpc_frame$UCXabnl <- ifelse(brierpc_frame$UCX_abnormal == "yes", 1, 0)
briered_frame <- urine_val %>% select(c('UCX_abnormal', 'prob_xgb', 'prob_rf', 'prob_ann', 'prob_rand', 'prob_micro'))
briered_frame$UCXabnl <- ifelse(briered_frame$UCX_abnormal == "yes", 1, 0)

save(briered_frame, file='brier_pc.Rdata')
save(briered_frame, file='brier_ed.Rdata')

getBrier <- function(dataset, probe) {
  b <- mean((dataset$UCXabnl - dataset[[probe]])^2)
  
  mx <- mean(dataset$UCXabnl)
  
  bref <- mean((dataset$UCXabnl - mx)^2)
  
  res <- 1-(b/bref)
  return(res)
}

# Helper function for bootstrap analysis; calculates RMSE on indices from data passed by boot
brierModelBootstrap <- function(probe, data, indices, pb) {
  d <- data[indices,]
  pb$tick()
  return(getBrier(d, probe))
}

# Calculates the RMSE of a model including error from boostrap analysis (95% CI; lci = 2.5% bound, uci = 97.5% bound)
brierModelWithError <- function(dataset, probe, setname) {
  R=2000
  
  pb <- progress_bar$new(total=R+1)
  
  # Calculate bootstrap object
  bootobj <- boot(dataset, statistic=brierModelBootstrap, R=R, strata=dataset$UCX_abnormal, probe=probe, parallel="multicore", pb=pb)
  # Use the bootstrap object to get a 95% CI
  ci <- boot.ci(bootobj, type="perc")
  
  # Return a data frame showing the method, point estimate and 95% CI, for inclusion in a larger table
  return(data.frame(method=probe, setname=setname, brier=ci$t0, lci=ci$percent[1,4],  uci=ci$percent[1,5]))
}


getBrierTable <- function() {
  briertable <- brierModelWithError(briered_frame, 'prob_xgb', 'ed')
  briertable <- rbind(briertable, brierModelWithError(briered_frame, 'prob_ann', 'ed'))
  briertable <- rbind(briertable, brierModelWithError(briered_frame, 'prob_rf', 'ed'))
  briertable <- rbind(briertable, brierModelWithError(briered_frame, 'prob_micro', 'ed'))
  briertable <- rbind(briertable, brierModelWithError(briered_frame, 'prob_rand', 'ed'))
  briertable <- rbind(briertable, brierModelWithError(brierpc_frame, 'prob_xgb', 'pc'))
  briertable <- rbind(briertable, brierModelWithError(brierpc_frame, 'prob_ann', 'pc'))
  briertable <- rbind(briertable, brierModelWithError(brierpc_frame, 'prob_rf', 'pc'))
  briertable <- rbind(briertable, brierModelWithError(brierpc_frame, 'prob_rand', 'pc'))
  return(briertable)
}