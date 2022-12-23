# Create ROC Figures
fig.roc.int <- createROCPlot(int.allframe, 'Internal (Emergency Department) Validation ROC')
ggsave('InternalValid.png', dpi=600)
fig.roc.ext <- createROCPlot(ext.allframe, 'External (Primary Care) Validation ROC')
ggsave('ExternalValid.png', dpi=600)
