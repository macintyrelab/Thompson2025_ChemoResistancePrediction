
## Here we are going to test the prediction capacity of HRDetect for plat and taxane
# Clean environment
freshr::freshr()

# Libraries
library(survival)
library(survminer)
library(this.path)
library(data.table)

# Paths
pancancer_dir <- dirname(this.path())
base_dir <- dirname(dirname(pancancer_dir))
input_dir <- file.path(base_dir, 'Input_Data')
scripts_dir <- file.path(base_dir, 'Analysis/Helper_Scripts')
xml_dir <- file.path(input_dir, 'TCGA_XMLs')
figs_dir <- file.path(base_dir, 'Figures')
supp_figs_dir <- file.path(figs_dir, 'Supplementary')
supp_figs_flow_dir <- file.path(supp_figs_dir, 'Flowcharts')
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_figs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_figs_flow_dir, showWarnings = FALSE, recursive = TRUE)

# Files
source(file.path(scripts_dir, 'Helper_functions.R'))
source(file.path(scripts_dir, 'power_calc.R'))

## HRDetect - TCGA
# Link of ids
LINK=readRDS(file.path(input_dir, "TCGA_PCAWG_links_plusCancer.rds"))
# From Degasperi et al. 2020
HRDETECT=file.path(input_dir, "Degasperi2020_HRDetect_PCAWG.txt")
# From Knijnenburg et al. 2018, downloaded and preprocessed by Jordan Griffin (Gerke Lab, github repo: TCGAhrd)
HRDSCORES=file.path(input_dir, "Knijnenburg2018_HRD_Scores_TCGA_GerkeLabGithub.txt")
# From Marquard et al. 2015 (Has more OV data than Knijnenburg et al., 2018)
HRDSCORES2=file.path(input_dir, "Marquard2015_HRD_Scores_TCGA.txt")


##### HRDETECT PROCESS #####
# Transfer HRDetect status
dtHRDetect = fread(HRDETECT)
dtHRDetect$TCGA = LINK$TCGA[ match(dtHRDetect$V1, LINK$PCAWG) ]
dtHRDetect$HRDetectClass = dtHRDetect$`HRDetect single score` >= 0.7
dtHRDetect$HRDetectClass = factor(dtHRDetect$HRDetectClass, levels = c(TRUE,FALSE), 
                                labels = c("Positive","Negative"))

# Transfer HRD score
dtHRDScores = fread(HRDSCORES)
dtHRDetect$HRDScore = dtHRDScores$HRD_Score[ match(dtHRDetect$TCGA, substr(dtHRDScores$patient_id,1,12)) ]

# From Marquard et al. 2015
dtHRDScores2 = fread(HRDSCORES2)
dtHRDScores2$HRDScore = dtHRDScores2$NtAI + dtHRDScores2$LST + dtHRDScores2$`HRD-LOH`
dtHRDetect$HRDScore2 = dtHRDScores2$HRDScore[ match(dtHRDetect$TCGA, substr(dtHRDScores2$Tumor,1,12)) ]

# Combine both scores
dtHRDetect$HRDScore[ is.na(dtHRDetect$HRDScore)  ] = dtHRDetect$HRDScore2[ is.na(dtHRDetect$HRDScore)  ]
dtHRDetect$HRDScore2 = NULL

# Use Myriad myChoice classifier
dtHRDetect$MyriadmyChoice = ifelse(dtHRDetect$HRDScore >= 42, "Positive", "Negative")
dtHRDetect$MyriadmyChoice = factor(dtHRDetect$MyriadmyChoice, levels = c("Positive","Negative"))

dtHRDetect=dtHRDetect[,13:16]


##### PHASE II #####

### TCGA-BRCA Taxane
tcga_brca_tax_p2 <- phase_2_func('BRCA', 'taxane', 'tcga')
tcga_brca_tax_p2$HRDetectClass <- dtHRDetect$HRDetectClass[ match(tcga_brca_tax_p2$bcr_patient_barcode,dtHRDetect$TCGA) ]
tcga_brca_tax_p2$MyriadmyChoice <- dtHRDetect$MyriadmyChoice[ match(tcga_brca_tax_p2$bcr_patient_barcode,dtHRDetect$TCGA) ]

## HRDClass
# KM
cols=c("#E7B800","#2E9FDF")
km = survfit(Surv(TTF, censoring) ~ HRDetectClass, data = tcga_brca_tax_p2)
svg(file.path(supp_figs_dir, "KM_HRDetect_TCGA-BRCA_tax.svg"), width = 6, height = 4)
pKM = ggsurvplot(km, data = tcga_brca_tax_p2, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                 font.tickslab = c(5), fontsize = 3, size = 0.5, pval.size = 3,
                 ggtheme = theme(panel.background = element_rect(fill='transparent'),
                                 panel.grid.major = element_line(colour = "grey90"),
                                 panel.grid.minor = element_line(colour = "grey90"),
                                 axis.line = element_line(size = 0.5), 
                                 axis.ticks = element_line(size = 0.5),
                                 axis.ticks.length = unit(.1, "cm"),
                                 axis.text.x = element_text(size = 8),
                                 axis.text.y = element_text(size = 8)),
                 risk.table.height = 0.3,
                 palette = cols,
                 surv.median.line = "hv", risk.table.y.text.col = TRUE, 
                 risk.table.y.text = FALSE) + labs(x = "Time (days)") 
pKM
dev.off()

# Cox
summary(coxph(Surv(TTF, censoring) ~ HRDetectClass, data=tcga_brca_tax_p2)) #HR=1.843 [0.3067-11.08] #p=0.504 
table(tcga_brca_tax_p2$HRDetectClass)
#4P & 4N

## MyriadmyChoice
# KM
cols=c("#E7B800","#2E9FDF")
km = survfit(Surv(TTF, censoring) ~ MyriadmyChoice, data = tcga_brca_tax_p2)
svg(file.path(supp_figs_dir, "KM_Myriad_TCGA-BRCA_tax.svg"), width = 6, height = 4)
pKM = ggsurvplot(km, data = tcga_brca_tax_p2, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                 font.tickslab = c(5), fontsize = 3, size = 0.5, pval.size = 3,
                 ggtheme = theme(panel.background = element_rect(fill='transparent'),
                                 panel.grid.major = element_line(colour = "grey90"),
                                 panel.grid.minor = element_line(colour = "grey90"),
                                 axis.line = element_line(size = 0.5), 
                                 axis.ticks = element_line(size = 0.5),
                                 axis.ticks.length = unit(.1, "cm"),
                                 axis.text.x = element_text(size = 8),
                                 axis.text.y = element_text(size = 8)),
                 risk.table.height = 0.3,
                 palette = cols,
                 surv.median.line = "hv", risk.table.y.text.col = TRUE, 
                 risk.table.y.text = FALSE) + labs(x = "Time (days)") 
pKM
dev.off()
# Cox
summary(coxph(Surv(TTF, censoring) ~ MyriadmyChoice, data=tcga_brca_tax_p2)) #HR=1.843 [0.3067-11.08] #p=0.504 
table(tcga_brca_tax_p2$MyriadmyChoice)
#9P & 4N

### TCGA-OV Taxane
tcga_ov_tax_p2 <- phase_2_func('OV', 'taxane', 'tcga')
tcga_ov_tax_p2$HRDetectClass <- dtHRDetect$HRDetectClass[ match(tcga_ov_tax_p2$bcr_patient_barcode,dtHRDetect$TCGA) ]
tcga_ov_tax_p2$MyriadmyChoice <- dtHRDetect$MyriadmyChoice[ match(tcga_ov_tax_p2$bcr_patient_barcode,dtHRDetect$TCGA) ]

## HRDClass
# KM
cols=c("#E7B800","#2E9FDF")
km = survfit(Surv(TTF, censoring) ~ HRDetectClass, data = tcga_ov_tax_p2)
svg(file.path(supp_figs_dir, "KM_HRDetect_TCGA-OV_tax.svg"), width = 6, height = 4)
pKM = ggsurvplot(km, data = tcga_ov_tax_p2, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                 font.tickslab = c(5), fontsize = 3, size = 0.5, pval.size = 3,
                 ggtheme = theme(panel.background = element_rect(fill='transparent'),
                                 panel.grid.major = element_line(colour = "grey90"),
                                 panel.grid.minor = element_line(colour = "grey90"),
                                 axis.line = element_line(size = 0.5), 
                                 axis.ticks = element_line(size = 0.5),
                                 axis.ticks.length = unit(.1, "cm"),
                                 axis.text.x = element_text(size = 8),
                                 axis.text.y = element_text(size = 8)),
                 risk.table.height = 0.3,
                 palette = cols,
                 surv.median.line = "hv", risk.table.y.text.col = TRUE, 
                 risk.table.y.text = FALSE) + labs(x = "Time (days)") 
pKM
dev.off()

# Cox
summary(coxph(Surv(TTF, censoring) ~ HRDetectClass, data=tcga_ov_tax_p2)) #HR=0.6231 [0.1013-3.833] #p=0.61 
table(tcga_ov_tax_p2$HRDetectClass)
# 3N & 3P


## MyriadmyChoice
# KM
cols=c("#E7B800","#2E9FDF")
km = survfit(Surv(TTF, censoring) ~ MyriadmyChoice, data = tcga_ov_tax_p2)
svg(file.path(supp_figs_dir, "KM_Myriad_TCGA-OV_tax.svg"), width = 6, height = 4)
pKM = ggsurvplot(km, data = tcga_ov_tax_p2, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                 font.tickslab = c(5), fontsize = 3, size = 0.5, pval.size = 3,
                 ggtheme = theme(panel.background = element_rect(fill='transparent'),
                                 panel.grid.major = element_line(colour = "grey90"),
                                 panel.grid.minor = element_line(colour = "grey90"),
                                 axis.line = element_line(size = 0.5), 
                                 axis.ticks = element_line(size = 0.5),
                                 axis.ticks.length = unit(.1, "cm"),
                                 axis.text.x = element_text(size = 8),
                                 axis.text.y = element_text(size = 8)),
                 risk.table.height = 0.3,
                 palette = cols,
                 surv.median.line = "hv", risk.table.y.text.col = TRUE, 
                 risk.table.y.text = FALSE) + labs(x = "Time (days)") 
pKM
dev.off()
# Cox
summary(coxph(Surv(TTF, censoring) ~ MyriadmyChoice, data=tcga_ov_tax_p2)) #All positive!
table(tcga_ov_tax_p2$MyriadmyChoice)
#6P & 0N


### TCGA-OV Platinum
tcga_ov_plat_p2 <- phase_2_func('OV', 'platinum', 'tcga-ov')
tcga_ov_plat_p2$HRDetectClass <- dtHRDetect$HRDetectClass[ match(tcga_ov_plat_p2$bcr_patient_barcode,dtHRDetect$TCGA) ]
tcga_ov_plat_p2$MyriadmyChoice <- dtHRDetect$MyriadmyChoice[ match(tcga_ov_plat_p2$bcr_patient_barcode,dtHRDetect$TCGA) ]

## HRDClass
# KM
cols=c("#E7B800","#2E9FDF")
km = survfit(Surv(TTF, censoring) ~ HRDetectClass, data = tcga_ov_plat_p2)
svg(file.path(supp_figs_dir, "KM_HRDetect_TCGA-OV_plat.svg"), width = 6, height = 4)
pKM = ggsurvplot(km, data = tcga_ov_plat_p2, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                 font.tickslab = c(5), fontsize = 3, size = 0.5, pval.size = 3,
                 ggtheme = theme(panel.background = element_rect(fill='transparent'),
                                 panel.grid.major = element_line(colour = "grey90"),
                                 panel.grid.minor = element_line(colour = "grey90"),
                                 axis.line = element_line(size = 0.5), 
                                 axis.ticks = element_line(size = 0.5),
                                 axis.ticks.length = unit(.1, "cm"),
                                 axis.text.x = element_text(size = 8),
                                 axis.text.y = element_text(size = 8)),
                 risk.table.height = 0.3,
                 palette = cols,
                 surv.median.line = "hv", risk.table.y.text.col = TRUE, 
                 risk.table.y.text = FALSE) + labs(x = "Time (days)") 
pKM
dev.off()

# Cox
summary(coxph(Surv(TTF, censoring) ~ HRDetectClass, data=tcga_ov_plat_p2)) #HR=4.548 [1.896-10.91] #p=0.000691 
summary(coxph(Surv(TTF, censoring) ~ HRDetectClass + `Tumor Stage` + `Age at diagnosis`, data=tcga_ov_plat_p2)) #HR=3.8742 [1.3589-11.045] #p=0.0113  
table(tcga_ov_plat_p2$HRDetectClass)
# 18N & 17P


## MyriadmyChoice
# KM
cols=c("#E7B800","#2E9FDF")
km = survfit(Surv(TTF, censoring) ~ MyriadmyChoice, data = tcga_ov_plat_p2)
svg(file.path(supp_figs_dir, "KM_Myriad_TCGA-OV_plat.svg"), width = 6, height = 4)
pKM = ggsurvplot(km, data = tcga_ov_plat_p2, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                 font.tickslab = c(5), fontsize = 3, size = 0.5, pval.size = 3,
                 ggtheme = theme(panel.background = element_rect(fill='transparent'),
                                 panel.grid.major = element_line(colour = "grey90"),
                                 panel.grid.minor = element_line(colour = "grey90"),
                                 axis.line = element_line(size = 0.5), 
                                 axis.ticks = element_line(size = 0.5),
                                 axis.ticks.length = unit(.1, "cm"),
                                 axis.text.x = element_text(size = 8),
                                 axis.text.y = element_text(size = 8)),
                 risk.table.height = 0.3,
                 palette = cols,
                 surv.median.line = "hv", risk.table.y.text.col = TRUE, 
                 risk.table.y.text = FALSE) + labs(x = "Time (days)") 
pKM
dev.off()
# Cox
summary(coxph(Surv(TTF, censoring) ~ MyriadmyChoice , data=tcga_ov_plat_p2)) #HR=1.019 [0.3877-2.68] #p=0.969 
summary(coxph(Surv(TTF, censoring) ~ MyriadmyChoice + `Tumor Stage` + `Age at diagnosis`, data=tcga_ov_plat_p2)) #HR=0.7571 [0.2803-2.045] #p=0.5831   
table(tcga_ov_plat_p2$MyriadmyChoice)
#28P & 7N




