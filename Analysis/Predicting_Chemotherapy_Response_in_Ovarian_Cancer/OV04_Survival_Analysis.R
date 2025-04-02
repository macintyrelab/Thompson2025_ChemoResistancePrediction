# This script is for creating the platinum, paclitaxel, and doxorubicin hazard ratio plots for Figure 2

# Clean environment
freshr::freshr()

## Libraries
library(readr)
library(dplyr)
library(survival)
library(ggplot2)
library(this.path)
library(UpSetR)

## Paths
ovarian_dir <- dirname(this.path())
base_dir <- dirname(dirname(ovarian_dir))
input_dir <- file.path(base_dir, "Input_Data")
helper_dir <- file.path(dirname(ovarian_dir), 'Helper_Scripts')
figs_dir <- file.path(base_dir, "Figures")
supp_figs_dir <- file.path(figs_dir, 'Supplementary')
tabs_dir <- file.path(base_dir, "Tables")
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_figs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabs_dir, showWarnings = FALSE, recursive = TRUE)

## Files
# IDs of samples
id_mapping <- read_csv(file.path(input_dir, 'paclitaxel_IDmapping.csv'))
patient_info <- read_csv(file.path(input_dir, 'patient_info.csv'))
id_mapping <- rbind(patient_info %>% select(OV04_ID=`patient_id...1`, samplename=primary_sample), id_mapping)

# Copy number profiles => used to compute wgii & signatures
old_cohort_copynumbers <- readRDS(file.path(input_dir, 'downsampled_41_OV04_absolute_copyNumbers.rds'))
new_cohort_copynumbers <- readRDS(file.path(input_dir, 'IUK_cohort_downsampled_SLX-24424_absolute_copyNumbers_50kb.rds'))

# Signature activities => used to predict response
old_cohort_exposures <- readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_OV_tissue_new_50kb.rds'))
new_cohort_exposures <- readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_OV_tissue_old_50kb.rds'))

# Clinical data & input for cox models => these tables have all needed for analyses => only included patients that meet the inclusion criteria
platinum_data <- readRDS(file.path(input_dir, 'OV_tissue_platinum_predictions.rds'))
taxane_data <- readRDS(file.path(input_dir, 'OV_tissue_paclitaxel_predictions.rds'))
doxo_data <- readRDS(file.path(input_dir, 'OV_tissue_doxorubicin_predictions.rds'))

# TCGA activities
tcga_exposures <- readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_TCGA.rds'))


# Computing wgii => this is already included in the input data frame for cox analyses
old_segments=c()
samples<-colnames(old_cohort_copynumbers)
for(sample in samples){
  tab<-as.data.frame(CNpare::getSegTable(old_cohort_copynumbers[,sample]))
  old_segments=rbind(old_segments,cbind(tab,sample=rep(sample,nrow(tab))))
}
old_segments[,2:4]=sapply(2:4,function(x)as.numeric(old_segments[,x]))
tab=NULL
new_segments=c()
samples<-colnames(new_cohort_copynumbers)
for(sample in samples){
  tab<-as.data.frame(CNpare::getSegTable(new_cohort_copynumbers[,sample]))
  new_segments=rbind(new_segments,cbind(tab,sample=rep(sample,nrow(tab))))
}
new_segments[,2:4]=sapply(2:4,function(x)as.numeric(new_segments[,x]))
all_segments <- rbind(old_segments, new_segments)

wgii <- all_segments %>% mutate(length=end-start) %>% filter(segVal!=2) %>% group_by(sample) %>% summarise(wgii=sum(length/3e9))
colnames(wgii)[1] <- 'study_subject_id'
wgii$study_subject_id <- unlist(lapply(strsplit(wgii$study_subject_id, '-'), FUN= function(x) gsub('downsampled_', '', paste(x[1:2], collapse='-'))))
wgii$study_subject_id <- id_mapping$OV04_ID[match(wgii$study_subject_id, id_mapping$samplename)]
wgii$wgii <- wgii$wgii * 100


#### PLATINUM ANALYSES ####
## Preparation of input data => platinum_data
# Patients that do not meet the inclusion criteria are excluded => n=45 patients
# Clinical data has been curated
# We also annotated maintenance with a treatment
platinum_data$maintenance <- 'No'
platinum_data$maintenance[platinum_data$study_subject_id %in% c(645, 668, 713, 771, 828, 853, 875)] <- 'Yes'
# And we grouped samples by age_at_diagnosis in 4 clinically relevant categories
platinum_data$age_recoded <- as.character(car::recode(as.numeric(platinum_data$age_at_diagnosis), "lo:64=1; 65:hi=2"))


# # Predictions => CX3 and CX2 are scaled within the cohort & CX3>CX2 is predicted as sensitive
# # This is already included in the platinum_data input
# exposures <- rbind(old_cohort_exposures, new_cohort_exposures)
# rownames(exposures) <- unlist(lapply(strsplit(rownames(exposures), '-'), FUN= function(x) gsub('downsampled_', '', paste(x[1:2], collapse='-'))))
# exposures <- as.data.frame(exposures[rownames(exposures) %in% id_mapping$samplename[match(platinum_data$study_subject_id, id_mapping$OV04_ID)],])
# exposures$study_subject_id <- id_mapping$OV04_ID[match(rownames(exposures), id_mapping$samplename)]
# 
# mREF = as.matrix(exposures[ , c("CX3", "CX2")])
# lModel = list(mean = attributes(scale(mREF))$`scaled:center`,
#               scale = attributes(scale(mREF))$`scaled:scale`)
# exposures <- exposures %>%
#   mutate(sCX2=(CX2-lModel$mean[['CX2']])/lModel$scale[['CX2']],
#          sCX3=(CX3-lModel$mean[['CX3']])/lModel$scale[['CX3']])
# exposures <- exposures %>% mutate(prediction=ifelse(sCX3>sCX2,'Sensitive','Resistant'))
# exposures$prediction <- factor(exposures$prediction, levels=c('Sensitive', 'Resistant'))

# Predictions as.numeric
platinum_data$prediction <- as.numeric(platinum_data$prediction)

cox <- coxph(Surv(PFS_since_start, censoring) ~ prediction + tumour_stage + age_recoded + maintenance + wgii, data=platinum_data)
cox.zph(cox) #all ok!
plat_cox_df <- as.data.frame(summary(coxph(Surv(PFS_since_start, censoring) ~ prediction + tumour_stage + age_recoded + maintenance + wgii, data=platinum_data))$conf.int)
plat_cox_df$variate <- factor(rownames(plat_cox_df), levels=c('wgii', 'maintenanceYes', 'age_recoded2','tumour_stage4', 'tumour_stage3', 'prediction'))

plot_2a <- ggplot(plat_cox_df, aes(x = variate)) +
  geom_errorbar(aes(ymin = `lower .95`, ymax = `upper .95`), width = 0.15, linewidth = 1) +
  geom_point(aes(y = `exp(coef)`), size = 10) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1) +
  scale_y_log10(position = 'right') + #,
  #breaks = (c(5, 50, 100, 500, 5000, 50000))) +
  coord_flip() +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, vjust = 1))
ggsave(file.path(figs_dir, '2a_hr.svg'), plot_2a,
       device = 'svg', width = 2000, height = 2226, units = 'px')


#### TAXANE ANALYSES ####
## Preparation of input data => taxane_data
# Patients that do not meet the inclusion criteria are excluded => n=29 patients
# Clinical data has been curated
taxane_data$age_recoded <- as.character(car::recode(as.numeric(taxane_data$age_at_diagnosis), "lo:64=1; 65:hi=2"))
table(taxane_data$treatment_line)
taxane_data$treatment_line_recoded="2/3"
taxane_data$treatment_line_recoded[taxane_data$treatment_line>3]="4+"
taxane_data$treatment_line_recoded[taxane_data$treatment_line_recoded=="2" | taxane_data$treatment_line_recoded=="3"]="2/3"

# # Predictions => CX5 is scaled to the TCGA cohort & CX5>0 is predicted as sensitive
# # This is already included in the taxane_data input
# exposures <- rbind(old_cohort_exposures, new_cohort_exposures)
# rownames(exposures) <- unlist(lapply(strsplit(rownames(exposures), '-'), FUN= function(x) gsub('downsampled_', '', paste(x[1:2], collapse='-'))))
# exposures <- as.data.frame(exposures[rownames(exposures) %in% id_mapping$samplename[match(taxane_data$study_subject_id, id_mapping$OV04_ID)],])
# exposures$study_subject_id <- id_mapping$OV04_ID[match(rownames(exposures), id_mapping$samplename)]
# 
# mTCGA = as.matrix(tcga_exposures[ , c("CX3", "CX5")])
# lModel = list(mean = attributes(scale(mTCGA))$`scaled:center`,
#               scale = attributes(scale(mTCGA))$`scaled:scale`)
# exposures <- exposures %>% mutate(sCX5=(CX5-lModel$mean[['CX5']])/lModel$scale[['CX5']])
# exposures <- exposures %>% mutate(prediction=ifelse(sCX5>0,'Sensitive','Resistant'))
# exposures$prediction <- factor(exposures$prediction, levels=c('Sensitive', 'Resistant'))
# Predictions as.numeric
taxane_data$prediction <- as.numeric(taxane_data$prediction)

cox <- coxph(Surv(PFS, censoring) ~ prediction + tumour_stage + treatment_line_recoded + age_recoded + maintenance + wgii, data=taxane_data)
cox.zph(cox) #all ok!
tax_cox_df <- as.data.frame(summary(coxph(Surv(PFS, censoring) ~ prediction + tumour_stage + treatment_line_recoded + age_recoded +  maintenance + wgii, data=taxane_data))$conf.int)
tax_cox_df$variate <- factor(rownames(tax_cox_df), levels=c('wgii', 'maintenanceTRUE', 'age_recoded2', 'treatment_line_recoded4+', 'tumour_stage4', 'prediction'))

plot_2b <- ggplot(tax_cox_df, aes(x = variate)) +
  geom_errorbar(aes(ymin = `lower .95`, ymax = `upper .95`), width = 0.15, linewidth = 1) +
  geom_point(aes(y = `exp(coef)`), size = 10) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1) +
  scale_y_log10(position = 'right') + #,
  #breaks = (c(5, 50, 100, 500, 5000, 50000))) +
  coord_flip() +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, vjust = 1))
ggsave(file.path(figs_dir, '2b_hr.svg'), plot_2b,
       device = 'svg', width = 2000, height = 1934, units = 'px')


#### DOXORUBICIN ANALYSES ####
## Preparation of input data => doxo_data
# Patients that do not meet the inclusion criteria are excluded => n=36 patients
# Clinical data has been curated
doxo_data$age_recoded <- as.character(car::recode(as.numeric(doxo_data$age_at_diagnosis), "lo:64=1; 65:hi=2"))
table(doxo_data$treatment_line)
doxo_data$treatment_line_recoded="2"
doxo_data$treatment_line_recoded[doxo_data$treatment_line>2]="3+"
doxo_data=doxo_data[!is.na(doxo_data$primary_PFS),] #double-check

# # Predictions => CX8>0.01 or CX9>0.009 or CX13>0.009 is predicted as resistant
# # This is already included in the doxo_data input
# exposures <- rbind(old_cohort_exposures, new_cohort_exposures)
# rownames(exposures) <- unlist(lapply(strsplit(rownames(exposures), '-'), FUN= function(x) gsub('downsampled_', '', paste(x[1:2], collapse='-'))))
# exposures <- as.data.frame(exposures[rownames(exposures) %in% id_mapping$samplename[match(doxo_data$study_subject_id, id_mapping$OV04_ID)],])
# exposures$study_subject_id <- id_mapping$OV04_ID[match(rownames(exposures), id_mapping$samplename)]
# exposures <- exposures %>% mutate(prediction=ifelse(CX8 > 0.01 | CX9 > 0.009 | CX13 > 0.009, 'Resistant', 'Sensitive'))
# exposures$prediction <- factor(exposures$prediction, levels=c('Sensitive', 'Resistant'))
# Predictions as.numeric
doxo_data$prediction <- as.numeric(doxo_data$prediction)

# This interaction can be included in the model & the HR is superhigh
cox <- coxph(Surv(PFS, censoring) ~ prediction*primary_PFS + tumour_stage + age_recoded + maintenance + wgii, data=doxo_data)
cox.zph(cox) #all ok!

# We used intEST function of the interactionRCS package to obtain the HR of the prediction over the different primary_PFS at 6 months
# A relapse at 6 months is considered as resistant patients for 1st line platinum
# HR at 6 months = 20.020 (1.059-378.635)
# We plot this HR for the main figure
dox_cox_df <- as.data.frame(summary(coxph(Surv(PFS, censoring) ~ prediction*primary_PFS + tumour_stage + age_recoded + maintenance + wgii, data=doxo_data))$conf.int)
dox_cox_df$variate <- factor(rownames(dox_cox_df), levels=c('prediction:primary_PFS', 'wgii', 'maintenanceTRUE', 'age_recoded2', 'tumour_stage4', 'primary_PFS', 'prediction'))
dox_cox_df$`exp(coef)`[dox_cox_df$variate=="prediction"]=20.020
dox_cox_df$`lower .95`[dox_cox_df$variate=="prediction"]=1.059
dox_cox_df$`upper .95`[dox_cox_df$variate=="prediction"]=378.635

plot_2c <- ggplot(dox_cox_df, aes(x = variate)) +
  geom_errorbar(aes(ymin = `lower .95`, ymax = `upper .95`), width = 0.15, linewidth = 1) +
  geom_point(aes(y = `exp(coef)`), size = 10) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1) +
  scale_y_log10(position = 'right') + #,
  #breaks = (c(5, 50, 100, 500, 5000, 50000))) +
  coord_flip() +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, vjust = 1))
ggsave(file.path(figs_dir, '2c_hr.svg'), plot_2c,
       device = 'svg', width = 2000, height = 2810, units = 'px')


#### SUPPLEMENTARY ANALYSES #####
# Which samples overlap or do not overlap between the three cohorts?
# Generate UpSet plot for extended data figure
overlap <- data.frame(study_subject_id=sort(unique(c(platinum_data$study_subject_id, taxane_data$study_subject_id, doxo_data$study_subject_id))))
overlap$Platinum <- 0
overlap$Taxane <- 0
overlap$Anthracycline <- 0
overlap$Platinum[overlap$study_subject_id %in% platinum_data$study_subject_id] <- 1
overlap$Taxane[overlap$study_subject_id %in% taxane_data$study_subject_id] <- 1
overlap$Anthracycline[overlap$study_subject_id %in% doxo_data$study_subject_id] <- 1
png(filename=file.path(figs_dir, 'OV04_overlap_plat.png'), res = 300, width=2000, 1500)
upset(overlap, point.size = 4, text.scale=2, sets=c('Anthracycline', 'Taxane', 'Platinum'), keep.order=TRUE)
dev.off()


# Are feature components predictive?
# Generate cox models using components
old_components <- readRDS(file.path(input_dir, '3_SxC_OV_tissue_new_50kb.rds'))
new_components <- readRDS(file.path(input_dir, '3_SxC_OV_tissue_old_50kb.rds'))
old_components <- as.data.frame(old_components)
new_components <- as.data.frame(new_components)
new_components$study_subject_id <- id_mapping$OV04_ID[match(unlist(lapply(strsplit(rownames(new_components), '-'),
                                                                FUN= function(x) gsub('downsampled_', '', paste(x[1:2], collapse='-')))),
                                                  id_mapping$samplename)]
old_components$study_subject_id <- id_mapping$OV04_ID[match(unlist(lapply(strsplit(rownames(old_components), '-'),
                                                                FUN= function(x) gsub('downsampled_', '', paste(x[1:2], collapse='-')))),
                                                  id_mapping$samplename)]
new_components <- new_components %>% filter(!(study_subject_id %in% c(648, 788)))
components <- rbind(old_components, new_components)
rownames(components) <- unlist(lapply(strsplit(rownames(components), '-'), FUN= function(x) gsub('downsampled_', '', paste(x[1:2], collapse='-'))))
components$study_subject_id <- id_mapping$OV04_ID[match(rownames(components), id_mapping$samplename)]

platinum_components <- left_join(platinum_data, components, by="study_subject_id")
taxane_components <- left_join(taxane_data, components, by="study_subject_id")
doxo_components <- left_join(doxo_data, components,by="study_subject_id")
component_names <- colnames(components)[1:ncol(components)-1]
all_comp_df <- data.frame(chemo=c(), component=c(), HR=c(), p=c())
for(i in 1:length(component_names)){
  plat_cox <- summary(coxph(as.formula(paste0('Surv(PFS_since_start, censoring) ~ ', component_names[i])), data=platinum_components))
  tax_cox <- summary(coxph(as.formula(paste0('Surv(PFS, censoring) ~ ', component_names[i])), data=taxane_components))
  anth_cox <- summary(coxph(as.formula(paste0('Surv(PFS, censoring) ~ ', component_names[i])), data=doxo_components))
  tdf <- data.frame(chemo=c('Platinum', 'Taxane', 'Anthracycline'),
                    component=component_names[i],
                    HR=c(plat_cox$coefficients[1,2], tax_cox$coefficients[1,2], anth_cox$coefficients[1,2]),
                    p=c(plat_cox$coefficients[1,5], tax_cox$coefficients[1,5], anth_cox$coefficients[1,5]))
  all_comp_df <- rbind(all_comp_df, tdf)
}
all_comp_df <- all_comp_df %>% mutate(adj_p=p.adjust(p, method='BH'))


# Generate all signature corrected data frame
all_sig_df <- data.frame(chemo=c(), signature=c(), HR=c(), p=c())
for(i in 1:17){
  plat_cox <- summary(coxph(as.formula(paste0('Surv(PFS_since_start, censoring) ~ ', paste0('CX', i))), data=platinum_data))
  tax_cox <- summary(coxph(as.formula(paste0('Surv(PFS, censoring) ~ ', paste0('CX', i))), data=taxane_data))
  anth_cox <- summary(coxph(as.formula(paste0('Surv(PFS, censoring) ~ ', paste0('CX', i))), data=doxo_data))
  tdf <- data.frame(chemo=c('Platinum', 'Taxane', 'Anthracycline'),
                    signature=paste0('CX', i),
                    HR=c(plat_cox$coefficients[1,2], tax_cox$coefficients[1,2], anth_cox$coefficients[1,2]),
                    p=c(plat_cox$coefficients[1,5], tax_cox$coefficients[1,5], anth_cox$coefficients[1,5]))
  all_sig_df <- rbind(all_sig_df, tdf)
}
all_sig_df <- all_sig_df %>% mutate(adj_p=p.adjust(p, method='BH'))


# save results
write.table(x = all_comp_df,
            file = file.path(tabs_dir,
                             "Table_CoxModel_FeatureComponents_OV04.csv"),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

write.table(x = all_sig_df,
            file = file.path(tabs_dir,
                             "Table_CoxModel_CINSignatures_OV04.csv"),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
