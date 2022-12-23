# Daniel J. Parente, MD PhD
# 2021-03-26
# UTI ML Statistical Analysis
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(boot)
library(psych)

# Garbage collect
gc()

set.seed(1000)

# Load the training file
load('all-pc-models.RData')

models$pc$svm <- NULL
models$pc$glm <- NULL
attributes(models$pc)$fullnames <- c('XGBoost', 'Random forests', 'Artificial neural networks')

# And the file with the model trained using microscopic data
load('initialTraining-micro.Rdata')

# Read in the internal validation data
urine_full <- read.csv("urine_full3.csv")
urine_train <- filter(urine_full, split == "training")
urine_val <- filter(urine_full, split == "validation")

# Load all subroutines
source('functions.R')

# Load the primary data set, perform external validation
primarycare_dat <- read.csv('processed.txt', sep="\t", header=T)

# Pull in race and ethnicity demograhpic data
pc_race_eth <- read.csv('pc-race-eth.txt', sep="\t", head=T)
primarycare_dat <- left_join(primarycare_dat, pc_race_eth)

# Annotate high risk features
primarycare_dat$anyhrf <- with(primarycare_dat, temp==3 | hr==2 | sbp == 1 | ams == 1 | flank_pain == 1 | cva == 1 | vomiting == 1 | catheter == 1 | immunocompromised == 1 | pregnant == 1)

# For each primary care model
ext.rocs <- list(pc=list())
int.rocs <- list(pc=list())
ext.tprfpr <- list(pc=list())
int.tprfpr <- list(pc=list())
ext.pr <- list(pc=list())
int.pr <- list(pc=list())
for(i in 1:length(models$pc)) {
	# Get the model name
	model_name = attributes(models$pc)$names[i]
	prob_name = paste0('prob_', model_name) # call the row 'prob_<modelname>'
	model_pretty_name = attributes(models$pc)$fullnames[i]
	
	# External validation on primary care dataset
	primarycare_dat[[prob_name]] <- predict(models$pc[[model_name]], primarycare_dat, type="prob")[,2]
	ext.rocs$pc[[model_name]] <- roc(primarycare_dat$UCX_abnormal, primarycare_dat[[prob_name]])
	ext.tprfpr$pc[[model_name]] <- getTPRFPR(ext.rocs$pc[[model_name]], paste0('NoMicro: ', model_pretty_name))
	ext.pr$pc[[model_name]] <- getPRBase(primarycare_dat, prob_name, paste0('NoMicro: ', model_pretty_name))
	
	# Internal validation on Emergency Department dataset
	urine_val[[prob_name]] <- predict(models$pc[[model_name]], urine_val, type="prob")[,2]
	int.rocs$pc[[model_name]] <- roc(urine_val$UCX_abnormal, urine_val[[prob_name]])
	int.tprfpr$pc[[model_name]] <- getTPRFPR(int.rocs$pc[[model_name]], paste0('NoMicro: ', model_pretty_name))
	int.pr$pc[[model_name]] <- getPRBase(urine_val, prob_name, paste0('NoMicro: ', model_pretty_name))
}

int.allframe <- data.frame(int.tprfpr$pc[[1]])
ext.allframe <- data.frame(ext.tprfpr$pc[[1]])

int.prallframe <- data.frame(int.pr$pc[[1]])
ext.prallframe <- data.frame(ext.pr$pc[[1]])
for(i in 2:length(int.tprfpr$pc)) {
	int.allframe <- rbind(int.allframe, data.frame(int.tprfpr$pc[[i]]))
	ext.allframe <- rbind(ext.allframe, data.frame(ext.tprfpr$pc[[i]]))
	
	int.prallframe <- rbind(int.prallframe, data.frame(int.pr$pc[[i]]))
	ext.prallframe <- rbind(ext.prallframe, data.frame(ext.pr$pc[[i]]))
}

# Predict the probabilities and construct an ROC curve for micro validation
micro_probs <- predict(xgb_tune_micro, urine_val, type="prob")[,2]
micro_roc <- roc(urine_val$UCX_abnormal, micro_probs)
urine_val$prob_micro <- micro_probs

microframe <- getTPRFPR(micro_roc, 'NeedMicro-required: XGBoost')
micropr <- getPRBase(urine_val, 'prob_micro', 'NeedMicro-required: XGBoost')

# Construct probability for 'random' model
primarycare_dat$prob_rand <- runif(nrow(primarycare_dat))
random.roc.pc <- roc(primarycare_dat$UCX_abnormal, primarycare_dat$prob_rand, levels=c("no", "yes"), direction="<")
random.tprfpr.pc <- getTPRFPR(random.roc.pc, "Random prediction (unskilled) model")
random.pr.pc <- getPRBase(primarycare_dat, 'prob_rand', 'Random prediction (unskilled) model')

urine_val$prob_rand <- runif(nrow(urine_val))
random.roc.ed <- roc(urine_val$UCX_abnormal, urine_val$prob_rand, levels=c("no", "yes"), direction="<")
random.tprfpr.ed <- getTPRFPR(random.roc.ed, "Random prediction (unskilled) model")
random.pr.ed <- getPRBase(urine_val, 'prob_rand', 'Random prediction (unskilled) model')


# AUC
int.allframe <- rbind(int.allframe, microframe, random.tprfpr.ed)
int.allframe$method <- as.factor(int.allframe$method)
int.allframe$method <- ordered(int.allframe$method, levels=levels(int.allframe$method)[c(5,4,3,1,2)])

ext.allframe <- rbind(ext.allframe, random.tprfpr.pc)
ext.allframe$method <- as.factor(ext.allframe$method)
ext.allframe$method <- ordered(ext.allframe$method, levels=levels(ext.allframe$method)[c(4,3,2,1)])

# PR
int.prallframe <- rbind(int.prallframe, micropr, random.pr.ed)
int.prallframe$method <- as.factor(int.prallframe$method)
int.prallframe$method <- ordered(int.prallframe$method, levels=levels(int.prallframe$method)[c(5,4,3,1,2)])

ext.prallframe <- rbind(ext.prallframe, random.pr.pc)
ext.prallframe$method <- as.factor(ext.prallframe$method)
ext.prallframe$method <- ordered(ext.prallframe$method, levels=levels(ext.prallframe$method)[c(4,3,2,1)])

source('figures.R')

# Garbage collect
gc()

# Obtain statistics for each classifier
classifierStats <- getParametersWithThreshold(primarycare_dat, 'prob_xgb', ext.rocs$pc$xgb, 'PC', 'XGB', 'best')
classifierStats <- rbind(classifierStats, getParametersWithThreshold(primarycare_dat, 'prob_rf', ext.rocs$pc$rf, 'PC', 'RF', 'best'))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(primarycare_dat, 'prob_ann', ext.rocs$pc$ann, 'PC', 'ANN', 'best'))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(primarycare_dat, 'prob_xgb', ext.rocs$pc$xgb, 'PC', 'XGB', 0.85))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(primarycare_dat, 'prob_rf', ext.rocs$pc$rf, 'PC', 'RF', 0.85))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(primarycare_dat, 'prob_ann', ext.rocs$pc$ann, 'PC', 'ANN', 0.85))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(urine_val, 'prob_xgb', int.rocs$pc$xgb, 'ED', 'XGB', 'best'))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(urine_val, 'prob_rf', int.rocs$pc$rf, 'ED', 'RF', 'best'))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(urine_val, 'prob_ann', int.rocs$pc$ann, 'ED', 'ANN', 'best'))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(urine_val, 'prob_micro', micro_roc, 'ED', 'Micro/XGB', 'best'))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(urine_val, 'prob_xgb', int.rocs$pc$xgb, 'ED', 'XGB', 0.85))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(urine_val, 'prob_rf', int.rocs$pc$rf, 'ED', 'RF', 0.85))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(urine_val, 'prob_ann', int.rocs$pc$ann, 'ED', 'ANN', 0.85))
classifierStats <- rbind(classifierStats, getParametersWithThreshold(urine_val, 'prob_micro', micro_roc, 'ED', 'Micro/XGB', 0.85))
print(classifierStats)

# Determine a threshold value
threshValue <- 0.85
thresh85 <- tail(dplyr::filter(coords(ext.rocs$pc$rf, "all", tr=F), sensitivity > threshValue), 1)$threshold

# Get the data set that would have been below this threshold
primarycare_dat_belowThresh <- dplyr::filter(primarycare_dat, prob_xgb < thresh85)
primarycare_dat_aboveThresh <- dplyr::filter(primarycare_dat, prob_xgb >= thresh85)

primarycare_dat$xgb_above <- with(primarycare_dat,prob_xgb >= thresh85)

print("Data below threshold")
getContingency(primarycare_dat_belowThresh)	

print("Data above threshold")
getContingency(primarycare_dat_aboveThresh)

# Get the calibration figures
calibDF.ext <- getCalibrationData(primarycare_dat, 'prob_xgb', 'UCX_abnormal', 'NoMicro: XGBoost')
calibDF.ext <- rbind(calibDF.ext, getCalibrationData(primarycare_dat, 'prob_rf', 'UCX_abnormal', 'NoMicro: Random Forest'))
calibDF.ext <- rbind(calibDF.ext, getCalibrationData(primarycare_dat, 'prob_ann', 'UCX_abnormal', 'NoMicro: Artificial Neural Network'))
calibDF.ext <- rbind(calibDF.ext, getCalibrationData(primarycare_dat, 'prob_rand', 'UCX_abnormal', 'Random prediction (unskilled) model'))

calibDF.int <- getCalibrationData(urine_val, 'prob_xgb', 'UCX_abnormal', 'NoMicro: XGBoost')
calibDF.int <- rbind(calibDF.int, getCalibrationData(urine_val, 'prob_rf', 'UCX_abnormal', 'NoMicro: Random Forest'))
calibDF.int <- rbind(calibDF.int, getCalibrationData(urine_val, 'prob_ann', 'UCX_abnormal', 'NoMicro: Artificial Neural Network'))
calibDF.int <- rbind(calibDF.int, getCalibrationData(urine_val, 'prob_micro', 'UCX_abnormal', 'NeedMicro: XGBoost'))
calibDF.int <- rbind(calibDF.int, getCalibrationData(urine_val, 'prob_rand', 'UCX_abnormal', 'Random prediction (unskilled) model'))

calibFig.ext <- calibFigure(calibDF.ext, 'Model', 'External (Primary Care) Calibration')
ggsave('Calibration-Figure-Ext.png', dpi=600)

calibFig.int <- calibFigure(calibDF.int, 'Model', 'Internal (Emergency Department) Calibration')
ggsave('Calibration-Figure-Int.png', dpi=600)

table2 <- classifierStats %>% 
	dplyr::filter(threshold=="-1") %>% 
	select(c('dataset', 'classifier', 'auc', 'prauc'))
	
write.table(table2, file='table2.txt', sep="\t", quote=F, row.names=F)

write.table(classifierStats, file='class_stats.txt', sep="\t", quote=F, row.names=F)