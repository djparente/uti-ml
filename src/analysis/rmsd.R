
gc()

# Calculate RMSDs
rmsdTable <- rmsdModelWithError(urine_val, 'UCX_abnormal', 'prob_xgb', method="Int/NoMicro/XGB")
rmsdTable <- rbind(rmsdTable, rmsdModelWithError(urine_val, 'UCX_abnormal', 'prob_rf', method="Int/NoMicro/RF"))
rmsdTable <- rbind(rmsdTable, rmsdModelWithError(urine_val, 'UCX_abnormal', 'prob_ann', method="Int/NoMicro/ANN"))
rmsdTable <- rbind(rmsdTable, rmsdModelWithError(urine_val, 'UCX_abnormal', 'prob_micro', method="Int/NeedMicro/XGB"))
rmsdTable <- rbind(rmsdTable, rmsdModelWithError(urine_val, 'UCX_abnormal', 'prob_rand', method="Int/Random"))
rmsdTable <- rbind(rmsdTable, rmsdModelWithError(primarycare_dat, 'UCX_abnormal', 'prob_xgb', method="Ext/NoMicro/XGB"))
rmsdTable <- rbind(rmsdTable, rmsdModelWithError(primarycare_dat, 'UCX_abnormal', 'prob_rf', method="Ext/NoMicro/RF"))
rmsdTable <- rbind(rmsdTable, rmsdModelWithError(primarycare_dat, 'UCX_abnormal', 'prob_ann', method="Ext/NoMicro/ANN"))
rmsdTable <- rbind(rmsdTable, rmsdModelWithError(primarycare_dat, 'UCX_abnormal', 'prob_rand', method="Ext/Random"))

write.table(rmsdTable, 'calib-rmsds.txt', sep="\t", quote=F, row.names=F)