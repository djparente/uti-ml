source('functions.R')

# Create ROC Figures
fig.roc.int <- createROCPlot(int.allframe, 'Internal (Emergency Department) Validation ROC')
#ggsave('InternalValid.png', dpi=600)
fig.roc.ext <- createROCPlot(ext.allframe, 'External (Primary Care) Validation ROC')
#ggsave('ExternalValid.png', dpi=600)


calibFig.ext <- calibFigure(calibDF.ext, 'Model', 'External (Primary Care) Calibration')
ggsave('Calibration-Figure-Ext.png', dpi=600)

calibFig.int <- calibFigure(calibDF.int, 'Model', 'Internal (Emergency Department) Calibration')
ggsave('Calibration-Figure-Int.png', dpi=600)


plot_grid(
  fig.roc.int,
  calibFig.int, 
  fig.roc.ext,
  calibFig.ext,
  labels=c('A', 'B', 'C', 'D'), ncol=2)

ggsave('all.png', height=12, width=12, unit="in", dpi=600)

# Write out these data files
write.table(int.allframe, 'Internal ED ROC data.txt', sep="\t", quote=F, row.names=F)
write.table(ext.allframe, 'External PC ROC data.txt', sep="\t", quote=F, row.names=F)
write.table(calibDF.int, 'Internal ED Calibration data.txt', sep="\t", quote=F, row.names=F)
write.table(calibDF.ext, 'External PC Calibration data.txt', sep="\t", quote=F, row.names=F)

