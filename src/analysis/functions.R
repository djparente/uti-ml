library(broom)
library(PRROC)
library(progress)

dig3 <- function(x) {
	return(format(round(x, 3), nsmall=3))
}

dig2 <- function(x) {
	return(format(round(x, 2), nsmall=2))
}


aucstring <- function(x) {
	formatted <- sapply(ci(x), dig3)
	return(paste0(formatted[2], " (", formatted[1], "-", formatted[3], ')'))
}


prettyNPct <- function(n, freq) {
	
	prty <- c() # Empty vector to store pretty-formatted results
	for(i in 1:length(n)) {
		pct <- format(round(100*freq[i], digits=1), nsmall=1) # Calculate the percentage and round it
		prty = c(prty, paste0(n[i], " (", pct, ")")) # Add a formatted result like No (%) [e.g. "5 (3.2)"] to the result
	}
	
	return(prty) # Return the result
}


# Plot the primary care ROC curve
getTPRFPR <- function(roc_curve, method) {
	resframe <- data.frame("tpr" = roc_curve$sensitivities, "fpr" = 1-roc_curve$specificities)
	resframe <- resframe[order(resframe$fpr, resframe$tpr),]
	
	#resframe['method'] = paste0(method, ' AUC: ', aucstring(roc_curve))
	resframe['method'] <- method
	
	return(resframe)
}

getPR <- function(pr, method) {
	resframe <- data.frame("recall" = pr$curve[,1], "precision" = pr$curve[,2])
	resframe <- resframe[order(resframe$recall, resframe$precision),]
	
	resframe['method'] <- method
	
	#resframe['method'] = paste0(method, ' AUC: ', aucstring(roc_curve))
	
	return(resframe)
}

getPRBase <- function(dataset, prob_name, method) {
	fg <- dplyr::filter(dataset, UCX_abnormal == "yes")
	bg <- dplyr::filter(dataset, UCX_abnormal != "yes")
	
	pr_curve <- pr.curve(scores.class0 = fg[[prob_name]], scores.class1 = bg[[prob_name]], curve=T)
	
	return(getPR(pr_curve, method))
}

getPRCurve <- function(dataset, prob_name, method) {
	fg <- dplyr::filter(dataset, UCX_abnormal == "yes")
	bg <- dplyr::filter(dataset, UCX_abnormal != "yes")
	
	pr_curve <- pr.curve(scores.class0 = fg[[prob_name]], scores.class1 = bg[[prob_name]], curve=T)
	
	return(pr_curve)
}

createPRPlot <- function(combinedframe, graphtitle, baseline) {
	g <- ggplot(combinedframe, aes(x=recall, y=precision)) +
		geom_line(aes(linetype=method, color=method, shape=method), size=1) + 
		geom_point(aes(shape=method, color=method), size=2) + 
		scale_linetype_manual(values=c("solid", "dotted", "dotdash", "twodash", "solid")) + 
		scale_color_manual(values=c("#FF0000", "#0000FF", "#00DD00", "#DD00DD", "#00DDDD")) + 
		scale_shape_manual(values=c(NA, NA, NA, NA, 18)) + 
		ggtitle(graphtitle) + 
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
		#theme_bw() + 
		theme(plot.title = element_text(hjust = 0.5)) +
		theme(legend.position = c(0.35, 0.5)) +
		
		xlab('Recall (Sensitivity, True Positive Rate)') +
		ylab('Precision (Positive Predictive Value)') +
		geom_abline(intercept = baseline, slope = 0, size=1.2, lty=2) +
		labs(color="Model", linetype="Model", shape="Model")
		
	return(g)
}




# Pulled from the PolyBoost analysis and modified for use here
getCalibrationData <- function(dat, pred, truth, method=NULL) {
	dat <- data.frame(dat)
	
	# Quantize the predictions into the cut point bins
	dat$quantPred <- cut(dat[[pred]], seq(0, 1, 0.1), include.lowest=T, right=F)
	
	# Quantize the truth into 0 or 1. This has the useful property that mean(quantTruth) over some range
	# of interest (i.e. datapoints falling into a quantPred bins) is the percentage of deleterious variants
	# in the bin
	dat$quantTruth <- ifelse(dat[[truth]] == 'yes', 1, 0)
	
	# Aggregate means and medians, stratified by quantPred
	ag <- aggregate(. ~ quantPred, dat, FUN=function(x) c(mean=mean(x), median=median(x)))
	
	# quantTruth[,1] is the mean of quantTruth within each bin (i.e. the fraction deleterious) while pred[,1] is the mean and pred[,2] the median of the predictions within the bin
	return(data.frame(calib=ag$quantTruth[,1], means=ag[[pred]][,1], 1=ag[[pred]][,2], method=method))
}

getContingency <- function(x) {
res <- x %>%
	group_by(UCX_abnormal, rx_abx) %>%
	summarize(count=n()) %>% 
	pivot_wider(names_from=UCX_abnormal, values_from=count, names_prefix='UCX_abnl_')
	
	return(res)
}

# Calculates the RMSE of a data frame given cutpoints
rmsdDF <- function(dat, cuts) {
	dat <- data.frame(dat)
	
	# Quantize the predictions into the cut point bins
	dat$quantPred <- cut(dat$pred, seq(0, 1, 0.1), include.lowest=T, right=F)
	
	# Quantize the truth into 0 or 1. This has the useful property that mean(quantTruth) over some range
	# of interest (i.e. datapoints falling into a quantPred bins) is the percentage of deleterious variants
	# in the bin
	dat$quantTruth <- ifelse(dat$truth == 'yes', 1, 0)
	
	# Aggregate means and medians, stratified by quantPred
	ag <- aggregate(. ~ quantPred, dat, FUN=function(x) c(mean=mean(x), median=median(x)))
	
	# quantTruth[,1] is the mean of quantTruth within a bin (i.e. the percent deleterious) while pred[,2] is the 
	# median of the predictions (i.e. the characteristic prediction of the bin). The sqrt of the mean of the square
	# difference is the RMSD, the quantity we want
	rmsd <- sqrt(mean((ag$quantTruth[,1] - ag$pred[,2])^2))
	return(rmsd)
}

# Helper function for bootstrap analysis; calculates RMSE on indices from data passed by boot
rmsdModelBootstrap <- function(cuts, data, indices, pb) {
	d <- data[indices,]
	pb$tick()
	return(rmsdDF(d, cuts))
}

# Calculates the RMSE of a model including error from boostrap analysis (95% CI; lci = 2.5% bound, uci = 97.5% bound)
rmsdModelWithError <- function(valid, truth, preds, method, cuts=seq(0, 1, 0.1)) {
	# Set up a data frame relating ground truth to predictions
	predframe <- data.frame(truth=valid[[truth]], pred=valid[[preds]])
	
	R=2000
	
	pb <- progress_bar$new(total=R+1)
	
	# Calculate bootstrap object
	bootobj <- boot(predframe, statistic=rmsdModelBootstrap, R=R, cuts=cuts, parallel="multicore", pb=pb)
	# Use the bootstrap object to get a 95% CI
	ci <- boot.ci(bootobj, type="perc")
	
	# Return a data frame showing the method, point estimate and 95% CI, for inclusion in a larger table
	return(data.frame(method=method, rmsd=ci$t0, lci=ci$percent[1,4],  uci=ci$percent[1,5]))
}

createROCPlot <- function(combinedframe, graphtitle) {
  recomb <- combinedframe
  recomb$method <- factor(recomb$method, rev(sort(levels(recomb$method))))
  
  g <- ggplot(recomb, aes(x=fpr, y=tpr)) +
    geom_line(aes(linetype=method, color=method, shape=method), size=1) + 
    geom_point(aes(shape=method, color=method), size=2) + 
    scale_linetype_manual(values=c("solid", "solid", "dotdash", "twodash", "dotted")) + 
    #scale_color_manual(values=c("#FF0000", "#0000FF", "#00DD00", "#DD00DD", "#00DDDD")) + 
    scale_color_brewer(palette="Set1") + 
    scale_shape_manual(values=c(NA, NA, NA, NA, NA)) + 
    ggtitle(graphtitle) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    #theme_bw() + 
    #theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = c(0.65, 0.2)) +
    
    xlab('False Positive Rate') +
    ylab('True Positive Rate') +
    geom_abline(intercept = 0, slope = 1, size=1.2, lty=2) +
    labs(color="Model", linetype="Model", shape="Model")
  
  return(g)
}


calibFigure <- function(calibDF, legendtitle, graphtitle="Calibration") {
  recomb <- calibDF
  recomb$method <- factor(recomb$method, rev(sort(levels(recomb$method))))
  
	res <- ggplot(recomb, aes(x=means, y=calib, group=method)) +
		xlim(0, 1) +
		ylim(0, 1) +
		geom_line(aes(linetype=method, color=method), size=1) + 
		geom_point(aes(shape=method, color=method), size=3) + 
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
		ggtitle(graphtitle) + 
		xlab('Abnormal culture probability (predicted)') +
		ylab('Abnormal culture rate (observed)') +
		geom_abline(intercept = 0, slope = 1, size=0.5) +
		scale_linetype_manual(values=c("dashed", "dashed", "dashed", "dashed", "dashed")) + 
		#scale_color_manual(values=c("#000000", "#FF0000", "#00FF00", "#0000FF", "#00DDDD")) + 
	  #scale_color_manual(values=c("#FF0000", "#0000FF", "#00DD00", "#DD00DD", "#00DDDD")) +
	  scale_color_brewer(palette="Set1") + 
		scale_shape_manual(values=c(15, 16, 17, 18, 19)) + 
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
		theme(legend.position = c(0.3, .8)) +
		
		labs(color=legendtitle, linetype=legendtitle, shape=legendtitle)
		
	return(res)
}

getContTable <- function(dataset, prob_name, threshold) {
	return(table(dataset[[prob_name]] > threshold, dataset[['UCX_abnormal']]))
}

boot.getParameters <- function(data, indices, prob_name, threshold) {
	# Get contingency table
	ct <- getContTable(data[indices,], prob_name, threshold)
	
	# Calculate sensitivity, specificity, Youden index, positive/negative predictive values
	sen <- ct[2,2] / (ct[1, 2] + ct[2,2])
	spec <- ct[1, 1] / (ct[1,1] + ct[2,1])
	youden <- sen + spec - 1
	ppv <- ct[2,2] / (ct[2,1] + ct[2,2])
	npv <- ct[1,1] / (ct[1,1] + ct[1,2])
	
	# Calculate positive and negative likelihood ratios
	lrpos <- sen / (1-spec)
	lrneg <- (1-sen)/spec
	
	# Calculate diagnostic odds ratio
	dor <- lrpos / lrneg
	
	# Calculate F1 score
	f1 <- harmonic.mean(c(ppv, sen))
	
	# Calculate ROC-AUC
	auc <- roc(data[indices,][['UCX_abnormal']], data[indices,][[prob_name]], levels=c("no", "yes"), direction="<")$auc
	
	# Calculate PR-AUC
	fg <- dplyr::filter(data[indices,], UCX_abnormal == "yes")
	bg <- dplyr::filter(data[indices,], UCX_abnormal != "yes")
	prauc_curve <- pr.curve(scores.class0 = fg[[prob_name]], scores.class1 = bg[[prob_name]], curve=T)
	prauc <- prauc_curve$auc.integral
	
	return(c(auc=auc, sen=sen*100, spec=spec*100, ppv=ppv*100, npv=npv*100, youden=youden, lrpos=lrpos, lrneg=lrneg, dor=dor, f1=f1, prauc=prauc))
}

fmt.ParameterErrorEstimates <- function(statistic, conf.low, conf.high) {
	return(paste0(prettyNum(statistic, digits=3), " (", prettyNum(conf.low, digits=3), '-', prettyNum(conf.high, digits=3),')'))
}

getParametersWithError <- function(dataset, prob_name, threshold, roc, dataset_name='none', classifier_name='none', threshold_type='none') {
	boot.out <- boot(data=dataset, statistic = boot.getParameters, R=2000, strata=dataset$UCX_abnormal, threshold=threshold, prob_name=prob_name)
	resTibble <- tidy(boot.out, conf.int=T, conf.method="perc")
	
	res <- resTibble %>%
		mutate(estimate=fmt.ParameterErrorEstimates(statistic, conf.low, conf.high)) %>%
		dplyr::select(c('term', 'estimate'))
	
	metadata <- data.frame(rbind(
		c(term='dataset', estimate=dataset_name),
		c(term='classifier', estimate=classifier_name),
		c(term='threshold', estimate=threshold_type) ))
	
	res <- rbind(metadata, res)
	
	#res[[classifier_name]] <- res$estimate
	#res <- res %>%
	#	dplyr::select(-c('estimate')) %>%
	res <- res %>%
		pivot_wider(names_from=term, values_from=estimate)
	
	return(res)
}

getParameters <- function(dataset, prob_name, threshold, dataset_name='none', classifier_name='none', threshold_type='none') {
	ct <- getContTable(dataset, prob_name, threshold)
	sen <- ct[2,2] / (ct[1, 2] + ct[2,2])
	spec <- ct[1, 1] / (ct[1,1] + ct[2,1])
	youden <- sen + spec - 1
	ppv <- ct[2,2] / (ct[2,1] + ct[2,2])
	npv <- ct[1,1] / (ct[1,1] + ct[1,2])
	
	lrpos <- sen / (1-spec)
	lrneg <- (1-sen)/spec
	
	dor <- lrpos / lrneg
	f1 <- harmonic.mean(c(ppv, sen))
	
	res <- data.frame(dataset_name=dataset_name, classifier_name=classifier_name, threshold_type=threshold_type, thresh=threshold, sen=sen, spec=spec, ppv=ppv, npv=npv, youden=youden, lrpos=lrpos, lrneg=lrneg, dor=dor, f1=f1)
	return(res)
}

getParametersWithThreshold <- function(dataset, prob_name, roc, dataset_name=NULL, classifier_name=NULL, threshold_type="best") {
	thresh <- -1
	res <- NULL
	if( threshold_type == "best") {
		thresh <- coords(roc, "best", tr=F)$threshold
		threshold_type <- -1
	} else {
		thresh <- tail(dplyr::filter(coords(roc, "all", tr=F), sensitivity > threshold_type), 1)$threshold
	}
	
	res <- getParametersWithError(dataset, prob_name, thresh, roc, dataset_name, classifier_name, threshold_type=threshold_type)
	
	return(res)
		
}

mc.getParametersWithThreshold <- function(para) {
	return(do.call(getParametersWithThreshold, para))
}