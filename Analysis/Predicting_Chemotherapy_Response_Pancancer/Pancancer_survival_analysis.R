
# This script is for performing phaseII and phaseII clinical trials in the 
# different TCGA and HMF cohorts with statistical power.

# Clean environment
freshr::freshr()

# Libraries
library(survival)
library(survminer)
library(tidyr)
library(this.path)
library(rsvg)
library(DiagrammeR)
library(DiagrammeRsvg)
library(data.table)
library(doMC)
library(WeightIt)
library(cobalt)
library(ggplot2)

# Paths
pancancer_dir <- dirname(this.path())
base_dir <- dirname(dirname(pancancer_dir))
input_dir <- file.path(base_dir, 'Input_Data')
scripts_dir <- file.path(base_dir, 'Analysis/Helper_Scripts')
xml_dir <- file.path(input_dir, 'TCGA_XMLs')
figs_dir <- file.path(base_dir, 'Figures')
supp_figs_dir <- file.path(figs_dir, 'Supplementary')
supp_figs_flow_dir <- file.path(supp_figs_dir, 'Flowcharts')
ext_figs_dir <- file.path(figs_dir, 'Extended')
tabs_dir <- file.path(base_dir, 'Tables')
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_figs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_figs_flow_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabs_dir, showWarnings = FALSE, recursive = TRUE)

# Files
source(file.path(scripts_dir, 'Helper_functions.R'))
source(file.path(scripts_dir, 'power_calc.R'))


################################################################################

#### PROCESS INPUT DATA ####
# HMF
hmf_clin <- readRDS(paste0(input_dir,'/HMF_clinical_data.RDS'))
metadata <- read_delim(paste0(input_dir,'/metadata.tsv'),
                       delim = "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)
metadata$patientIdentifier <- gsub('TI*V*$', '', metadata$sampleId, perl=TRUE)
hmf_segments <- readRDS(file.path(input_dir, 'hmf_cnp_smooth.rds'))
tcga_segments_nocin <- readRDS(file.path(input_dir, '0_TCGA_Segments_noCIN.rds'))
tcga_segments <- readRDS(file.path(input_dir, '0_TCGA_Segments_dCIN.rds'))

WIGGLE=0.1
IGNOREDELS=FALSE

# Set everything very close to 2 to 2
hmf_segments$segVal[hmf_segments$segVal > (2-WIGGLE) & hmf_segments$segVal < (2+WIGGLE)] = 2
# Merge segments only when two normal follow each other -> SMOOTHINGFACTOR = 0
hmf_segments = idSmoothingTargets(hmf_segments, WIGGLE = 0, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = IGNOREDELS)
# Split by sample name
lRaw = split(hmf_segments, hmf_segments$sample)
# Smooth segments by taking the weighted average of the segVal and their lengths
lhmf_segments = smoothSegments(lRaw, CORES=1, SMOOTHINGFACTOR = 0, colNameMerge = "segVal", colNameChr = "chromosome",
                               colNameStart = "start", colNameEnd = "end", IGNOREDELS = IGNOREDELS, asDf = FALSE)

hmf_segments = rbindlist(lhmf_segments)
wgii_hmf <- hmf_segments %>% mutate(length=end-start) %>% filter(segVal!=2) %>% group_by(sample) %>% summarise(wgii=sum(length/3e9))
wgii_hmf$wgii[wgii_hmf$wgii > 1] <- 1
colnames(wgii_hmf)[1] <- 'sampleId'
wgii_hmf$wgii <- wgii_hmf$wgii * 100

tcga_segments <- rbind(tcga_segments, tcga_segments_nocin)
wgii_tcga <- tcga_segments %>% mutate(length=end-start) %>% filter(segVal!=2) %>% group_by(sample) %>% summarise(wgii=sum(length/3e9))
colnames(wgii_tcga)[1] <- 'bcr_patient_barcode'
wgii_tcga$wgii <- wgii_tcga$wgii * 100


#### POWER CALCULATIONS PHASE III ####
# Phase III is performed only in those cohorts in which we have enough power
Power_res_df <- data.frame(cohort=c(),
                       cancer=c(),
                       chemo=c(),
                       res_powered_hr=c(),
                       res_required_size=c(),
                       res_actual_size=c())
Power_sen_df <- data.frame(cohort=c(),
                           cancer=c(),
                           chemo=c(),
                           sen_powered_hr=c(),
                           sen_required_size=c(),
                           sen_actual_size=c())

chemos <- c('platinum', 'taxane', 'doxorubicin')

# HMF
hmf_cohorts <- unique(hmf_clin$cancer_type)
hmf_cohorts <- hmf_cohorts[-c(13,36)]

for(i in 1:length(hmf_cohorts)){
  for(j in 1:length(chemos)){
    if(hmf_cohorts[i]=="Breast"){
      a <- try(phase_3_func(cancer=hmf_cohorts[i], chemo=chemos[j], cohort='hmf', SA_exp=FALSE, SA_control=FALSE, control_pattern='citabine'))
    }else{
      a <- try(phase_3_func(cancer=hmf_cohorts[i], chemo=chemos[j], cohort='hmf', SA_exp=FALSE, SA_control=FALSE))
    }
    if(class(a)=='list'){
      # Resistant
      res <- a[[1]]
      P <- mean(res$censoring)
      pred_prop <- table(res$Exp_treatment)[1]/nrow(res)
      d <- nrow(res)
      if(j==1){
        hr <- 3 #OV04 plat
      }
      if(j==2){
        hr <- 5.5 #OV04 taxane
      }
      if(j==3){
        hr <- 5.5 #HR at 6 months in OV04 doxo is 22 => we limit to a putative expected HR
      }
      res_powered_hr <- power_calc_inv(P=P, Rsquared=0.05, pred_prop=pred_prop, d=d)[2]
      res_required_size <- ceiling(power_calc(P=P, Rsquared=0.05, pred_prop=pred_prop, hr=hr))
      res_actual_size=table(res$Exp_treatment)[1]
      
      # Output
      out <- data.frame(cohort='hmf',
                        cancer=hmf_cohorts[i],
                        chemo=chemos[j],
                        res_powered_hr=res_powered_hr,
                        res_required_size=res_required_size,
                        res_actual_size=res_actual_size)
      Power_res_df <- rbind(Power_res_df, out)
    }
  }
}

for(i in 1:length(hmf_cohorts)){
  for(j in 1:length(chemos)){
    if(hmf_cohorts[i]=="Breast"){
      a <- try(phase_3_func(cancer=hmf_cohorts[i], chemo=chemos[j], cohort='hmf', SA_exp=FALSE, SA_control=FALSE, control_pattern='citabine'))
    }else{
      a <- try(phase_3_func(cancer=hmf_cohorts[i], chemo=chemos[j], cohort='hmf', SA_exp=FALSE, SA_control=FALSE))
    }
    if(class(a)=='list'){
      # Sensitive
      sen <- a[[2]]
      P <- mean(sen$censoring)
      pred_prop <- table(sen$Exp_treatment)[1]/nrow(sen)
      d <- nrow(sen)
      if(j==1){
        hr <- 1/3 #OV04 plat
      }
      if(j==2){
        hr <- 1/5.5 #OV04 taxane
      }
      if(j==3){
        hr <- 1/5.5 #HR at 6 months in OV04 doxo is 22 => we limit to a putative expected HR
      }
      sen_powered_hr <- power_calc_inv(P=P, Rsquared=0.05, pred_prop=pred_prop, d=d)[1]
      sen_required_size <- ceiling(power_calc(P=P, Rsquared=0.05, pred_prop=pred_prop, hr=hr))
      sen_actual_size=table(sen$Exp_treatment)[1]
      
      # Output
      out <- data.frame(cohort='hmf',
                        cancer=hmf_cohorts[i],
                        chemo=chemos[j],
                        sen_powered_hr=sen_powered_hr,
                        sen_required_size=sen_required_size,
                        sen_actual_size=sen_actual_size)
      Power_sen_df <- rbind(Power_sen_df, out)
    }
  }
}


# TCGA
tcga_clin <- readRDS(paste0(input_dir,'/TCGA_clinical_data.RDS'))
tcga_cohorts <- unique(tcga_clin$cancer_type)

for(i in 1:length(tcga_cohorts)){
  for(j in 1:length(chemos)){
    a <- try(phase_3_func(cancer=tcga_cohorts[i], chemo=chemos[j], cohort='tcga', SA_exp=FALSE, SA_control=FALSE))
    if(class(a)=='list'){
      # Resistant
      res <- a[[1]]
      if(tcga_cohorts[i]=="BRCA"){res=res[res$type=="TNBC" | sen$type=="Other",]}
      P <- mean(res$censoring)
      pred_prop <- table(res$Exp_treatment)[1]/nrow(res)
      d <- nrow(res)
      if(j==1){
        hr <- 3 #OV04 plat
      }
      if(j==2){
        hr <- 5.5 #OV04 taxane
      }
      if(j==3){
        hr <- 5.5 #HR at 6 months in OV04 doxo is 22 => we limit to a putative expected HR
      }
      res_powered_hr <- power_calc_inv(P=P, Rsquared=0.05, pred_prop=pred_prop, d=d)[2]
      res_required_size <- ceiling(power_calc(P=P, Rsquared=0.05, pred_prop=pred_prop, hr=hr))
      res_actual_size=table(res$Exp_treatment)[1]
      
      # Output
      out <- data.frame(cohort='tcga',
                        cancer=tcga_cohorts[i],
                        chemo=chemos[j],
                        res_powered_hr=res_powered_hr,
                        res_required_size=res_required_size,
                        res_actual_size=res_actual_size)
      Power_res_df <- rbind(Power_res_df, out)
    }
  }
}

for(i in 1:length(tcga_cohorts)){
  for(j in 1:length(chemos)){
    a <- try(phase_3_func(cancer=tcga_cohorts[i], chemo=chemos[j], cohort='tcga', SA_exp=FALSE, SA_control=FALSE))
    if(class(a)=='list'){
      # Sensitive
      sen <- a[[2]]
      if(tcga_cohorts[i]=="BRCA"){sen=sen[sen$type=="TNBC" | sen$type=="Other",]}
      P <- mean(sen$censoring)
      pred_prop <- table(sen$Exp_treatment)[1]/nrow(sen)
      d <- nrow(sen)
      if(j==1){
        hr <- 1/3 #OV04 plat
      }
      if(j==2){
        hr <- 1/5.5 #OV04 taxane
      }
      if(j==3){
        hr <- 1/5.5 #HR at 6 months in OV04 doxo is 22 => we limit to a putative expected HR
      }
      sen_powered_hr <- power_calc_inv(P=P, Rsquared=0.05, pred_prop=pred_prop, d=d)[1]
      sen_required_size <- ceiling(power_calc(P=P, Rsquared=0.05, pred_prop=pred_prop, hr=hr))
      sen_actual_size=table(sen$Exp_treatment)[1]
      
      # Output
      out <- data.frame(cohort='tcga',
                        cancer=tcga_cohorts[i],
                        chemo=chemos[j],
                        sen_powered_hr=sen_powered_hr,
                        sen_required_size=sen_required_size,
                        sen_actual_size=sen_actual_size)
      Power_sen_df <- rbind(Power_sen_df, out)
    }
  }
}

# All cohorts we can perform phase III analyses
Power_Phase3_res_df <- Power_res_df[Power_res_df$res_actual_size>Power_res_df$res_required_size,] 
Power_Phase3_sen_df <- Power_sen_df[Power_sen_df$sen_actual_size>Power_sen_df$sen_required_size,] #8 cohorts
# save in tables!
Power = full_join(Power_res_df,Power_sen_df,by=c("cohort","cancer","chemo"))
# save results
write.table(x = Power,
            file = file.path(tabs_dir,
                             "Power_Analyses_Results.csv"),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


#### POWER CALCULATIONS PHASE II ####
Power_p2_df <- data.frame(cohort=c(),
                           cancer=c(),
                           chemo=c(),
                           powered_hr=c(),
                           required_size=c(),
                           actual_size=c())

chemos <- c('platinum', 'taxane', 'doxorubicin')

# HMF
hmf_cohorts <- unique(hmf_clin$cancer_type)
hmf_cohorts <- hmf_cohorts[-c(13,36)]

for(i in 1:length(hmf_cohorts)){
  for(j in 1:length(chemos)){
    a <- try(phase_2_func(cancer=hmf_cohorts[i], chemo=chemos[j], cohort='hmf'))
    
    P <- mean(a$censoring)
    pred_prop <- table(a$prediction)[1]/nrow(a)
    d <- nrow(a)
    if(j==1){
      hr <- 3 #OV04 plat
    }
    if(j==2){
      hr <- 5.5 #OV04 taxane
    }
    if(j==3){
      hr <- 5.5 #HR at 6 months in OV04 doxo is 22 => we limit to a putative expected HR
    }
    powered_hr <- power_calc_inv(P=P, Rsquared=0.05, pred_prop=pred_prop, d=d)[2]
    required_size <- ceiling(power_calc(P=P, Rsquared=0.05, pred_prop=pred_prop, hr=hr))
    actual_size=sum(table(a$prediction))
    
    # Output
    out <- data.frame(cohort='hmf',
                      cancer=hmf_cohorts[i],
                      chemo=chemos[j],
                      powered_hr=powered_hr,
                      required_size=required_size,
                      actual_size=actual_size)
    Power_p2_df <- rbind(Power_p2_df, out)
    
  }
}

# TCGA
tcga_clin <- readRDS(paste0(input_dir,'/TCGA_clinical_data.RDS'))
tcga_cohorts <- unique(tcga_clin$cancer_type)

for(i in 1:length(tcga_cohorts)){
  for(j in 1:length(chemos)){
    a <- try(phase_2_func(cancer=tcga_cohorts[i], chemo=chemos[j], cohort='tcga'))
    if(tcga_cohorts[i]=="BRCA"){a=a[a$type=="TNBC" | a$type=="Other",]}
    if(nrow(a)>0){
      P <- mean(a$censoring)
      pred_prop <- table(a$prediction)[1]/nrow(a)
      d <- nrow(a)
      if(j==1){
        hr <- 3 #OV04 plat
      }
      if(j==2){
        hr <- 5.5 #OV04 taxane
      }
      if(j==3){
        hr <- 5.5 #HR at 6 months in OV04 doxo is 22 => we limit to a putative expected HR
      }
      powered_hr <- power_calc_inv(P=P, Rsquared=0.05, pred_prop=pred_prop, d=d)[2]
      required_size <- ceiling(power_calc(P=P, Rsquared=0.05, pred_prop=pred_prop, hr=hr))
      actual_size=sum(table(a$prediction))
      
      # Output
      out <- data.frame(cohort='tcga',
                        cancer=tcga_cohorts[i],
                        chemo=chemos[j],
                        powered_hr=powered_hr,
                        required_size=required_size,
                        actual_size=actual_size)
      Power_p2_df <- rbind(Power_p2_df, out)
    }
  }
}


# All cohorts we can perform phase II analyses
Power_Phase2_df <- Power_p2_df[Power_p2_df$actual_size>Power_p2_df$required_size,] 
# save in tables!
write.table(x = Power_p2_df,
            file = file.path(tabs_dir,
                             "Power_Analyses_Results_PhaseII.csv"),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


#### PHASE II ####
## SURVIVAL ANALYSES
# We are doing phase II for TCGA-OV Platinum and TCGA-SARC Anthracyclines
# We have power to do phase III. However, we cannot apply this design. This is
# because platinum is the 1st treatment for all OV tumours, so there is not alternative SoC
# For Sarcoma we don't have the numbers to run a phase III analysis given the standard clinical treatment
# pathways.

# TCGA-OV Platinum => we are powered for 1st line
tcga_ov_plat_p2 <- phase_2_func('OV', 'platinum', 'tcga-ov')
tcga_ov_plat_p2 <- left_join(tcga_ov_plat_p2, wgii_tcga)
tcga_ov_plat_p2$`Tumor Stage` <-  gsub('[ABC]', '', tcga_ov_plat_p2$`Tumor Stage`)
colnames(tcga_ov_plat_p2)[colnames(tcga_ov_plat_p2)=='Tumor Stage'] <- 'Tumour stage'
colnames(tcga_ov_plat_p2)[colnames(tcga_ov_plat_p2)=='prediction'] <- 'Prediction'
colnames(tcga_ov_plat_p2)[colnames(tcga_ov_plat_p2)=='wgii'] <- 'wGII'
tcga_ov_plat_p2$age_recoded <- as.character(car::recode(as.numeric(tcga_ov_plat_p2$'Age at diagnosis'), "lo:59=1; 60:69=2; 70:hi=3"))
tcga_ov_plat_p2_cox <- coxph(Surv(TTF, censoring) ~ Prediction + `Tumour stage` + strata(age_recoded), data=tcga_ov_plat_p2) #all are 1st regimen, so treatment_line excluded as covariate
cox.zph(tcga_ov_plat_p2_cox) #all ok!
render_graph(phase_2_plot_func('OV', 'platinum', 'tcga')) %>%
  export_svg %>%
  charToRaw %>%
  rsvg_pdf(file.path(supp_figs_flow_dir, 'TCGA_OV_Plat_filtering_phase_2.pdf'))

# TCGA-SARC Dox => we are powered for phase II
tcga_sarc_dox_p2 <- phase_2_func('SARC', 'doxorubicin', 'tcga')
tcga_sarc_dox_p2 <- left_join(tcga_sarc_dox_p2 %>% select(-wgii), wgii_tcga)
tcga_sarc_dox_p2$wgii[is.na(tcga_sarc_dox_p2$wgii)]=0 #all segments are diploid in TCGA-MO-A47P
colnames(tcga_sarc_dox_p2)[colnames(tcga_sarc_dox_p2)=='prediction'] <- 'Prediction'
colnames(tcga_sarc_dox_p2)[colnames(tcga_sarc_dox_p2)=='wgii'] <- 'wGII'
colnames(tcga_sarc_dox_p2)[colnames(tcga_sarc_dox_p2)=='age'] <- 'Age at diagnosis'
tcga_sarc_dox_p2$age_recoded <- as.character(car::recode(as.numeric(tcga_sarc_dox_p2$'Age at diagnosis'), "lo:59=1; 60:hi=2"))
colnames(tcga_sarc_dox_p2)[colnames(tcga_sarc_dox_p2)=='treatment_line'] <- 'Treatment line'
tcga_sarc_dox_p2$`Treatment line`[tcga_sarc_dox_p2$`Treatment line`>1]="2+"
tcga_sarc_dox_p2$Ifosfamide <- grepl('ifo', tcga_sarc_dox_p2$drug_name)
tcga_sarc_dox_p2_cox <- coxph(Surv(TTF, censoring) ~ Prediction + Ifosfamide, data=tcga_sarc_dox_p2)
cox.zph(tcga_sarc_dox_p2_cox) #all ok
render_graph(phase_2_plot_func('SARC', 'doxorubicin', 'tcga')) %>%
  export_svg %>%
  charToRaw %>%
  rsvg_pdf(file.path(supp_figs_flow_dir, 'TCGA_SARC_Dox_filtering_phase_2.pdf'))


# TCGA-BRCA TNBC Taxane => we have not power for phase II
brca_comb_tax_p2 <- phase_2_func('BRCA', 'taxane', 'tcga', SA_exp = FALSE)
brca_comb_tax_p2 <- brca_comb_tax_p2[brca_comb_tax_p2$type=="TNBC"|brca_comb_tax_p2$type=="Other",]
powered_hr <- power_calc_inv(P=mean(brca_comb_tax_p2$censoring), Rsquared=0.05,
                             pred_prop=table(brca_comb_tax_p2$prediction)[2]/nrow(brca_comb_tax_p2),
                             d=nrow(brca_comb_tax_p2))[2]
required_size <- ceiling(power_calc(P=mean(brca_comb_tax_p2$censoring), Rsquared=0.05,
                                    pred_prop=table(brca_comb_tax_p2$prediction)[2]/nrow(brca_comb_tax_p2),
                                    hr=5.5))

# HMF Colorectum Platinum => we have not power for phase II (unbalanced data)
col_comb_plat_p2 <- phase_2_func('Colorectum', 'platinum', 'hmf')
powered_hr <- power_calc_inv(P=mean(col_comb_plat_p2$censoring), Rsquared=0.05,
                             pred_prop=table(col_comb_plat_p2$prediction)[2]/nrow(col_comb_plat_p2),
                             d=nrow(col_comb_plat_p2))[2]
required_size <- ceiling(power_calc(P=mean(col_comb_plat_p2$censoring), Rsquared=0.05,
                                    pred_prop=table(col_comb_plat_p2$prediction)[2]/nrow(col_comb_plat_p2),
                                    hr=3))


# HMF Lung taxane => we have not power for phase II 
lung_comb_tax_p2 <- phase_2_func('Lung', 'taxane', 'hmf')
powered_hr <- power_calc_inv(P=mean(lung_comb_tax_p2$censoring), Rsquared=0.05,
                             pred_prop=table(lung_comb_tax_p2$prediction)[2]/nrow(lung_comb_tax_p2),
                             d=nrow(lung_comb_tax_p2))[2]
required_size <- ceiling(power_calc(P=mean(col_comb_plat_p2$censoring), Rsquared=0.05,
                                    pred_prop=table(col_comb_plat_p2$prediction)[2]/nrow(col_comb_plat_p2),
                                    hr=5.5))


#### PHASE III #################################################################

##### PRIMARY TUMOUR COMBINATION #####
# TCGA-OV Taxane
ov_comb_tax_p3 <- phase_3_func('OV', 'taxane', 'tcga-ov', SA_exp=FALSE, SA_control=FALSE)
ov_comb_tax_p3_res <- ov_comb_tax_p3[[1]]
ov_comb_tax_p3_res$treatment_line <- as.numeric(substr(ov_comb_tax_p3_res$treatment_line, 1, 1))
ov_comb_tax_p3_res <- left_join(ov_comb_tax_p3_res, wgii_tcga)
ov_comb_tax_p3_res$`Tumor Stage` <- gsub('[ABC]', '', ov_comb_tax_p3_res$`Tumor Stage`)
colnames(ov_comb_tax_p3_res)[colnames(ov_comb_tax_p3_res)=='Exp_treatment'] <- 'Treatment arm'
ov_comb_tax_p3_res$`Treatment arm` <- ifelse(ov_comb_tax_p3_res$`Treatment arm`, 'Experimental', 'Control')
colnames(ov_comb_tax_p3_res)[colnames(ov_comb_tax_p3_res)=='TTF 1st line (months)'] <- '1st line TTF'
colnames(ov_comb_tax_p3_res)[colnames(ov_comb_tax_p3_res)=='Tumor Stage'] <- 'Tumour stage'
colnames(ov_comb_tax_p3_res)[colnames(ov_comb_tax_p3_res)=='treatment_line'] <- 'Treatment line'
colnames(ov_comb_tax_p3_res)[colnames(ov_comb_tax_p3_res)=='wgii'] <- 'wGII'
ov_comb_tax_p3_res$age_recoded <- as.character(car::recode(as.numeric(ov_comb_tax_p3_res$'Age at diagnosis'), "lo:64=1; 65:hi=2"))
ov_comb_tax_p3_res$year_of_treatment <- apply(ov_comb_tax_p3_res,
                                              MARGIN=1,
                                              FUN=function(x) get_year_of_treatment(x[which(colnames(ov_comb_tax_p3_res)=='bcr_patient_barcode')], x[which(colnames(ov_comb_tax_p3_res)=='Treatment line')]))
ov_comb_tax_p3_res_weights <- weightit(`Treatment arm` ~ year_of_treatment,
                                       data=ov_comb_tax_p3_res,
                                       method='glm',
                                       estimand = 'ATT')
bal.tab(ov_comb_tax_p3_res_weights,
        un=TRUE)
ov_comb_tax_p3_res_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm`*`1st line TTF` + strata(age_recoded), data=ov_comb_tax_p3_res, weights=ov_comb_tax_p3_res_weights$weights)
cox.zph(ov_comb_tax_p3_res_cox) #global ok!
# TCGA OV tax resistant point estimate (with intEST function): 7.435 (2.782-20.458)

ov_comb_tax_p3_sen <- ov_comb_tax_p3[[2]]
ov_comb_tax_p3_sen$treatment_line <- as.numeric(substr(ov_comb_tax_p3_sen$treatment_line, 1, 1))
ov_comb_tax_p3_sen <- left_join(ov_comb_tax_p3_sen, wgii_tcga)
ov_comb_tax_p3_sen$`Tumor Stage` <- gsub('[ABC]', '', ov_comb_tax_p3_sen$`Tumor Stage`)
colnames(ov_comb_tax_p3_sen)[colnames(ov_comb_tax_p3_sen)=='Exp_treatment'] <- 'Treatment arm'
ov_comb_tax_p3_sen$`Treatment arm` <- ifelse(ov_comb_tax_p3_sen$`Treatment arm`, 'Experimental', 'Control')
colnames(ov_comb_tax_p3_sen)[colnames(ov_comb_tax_p3_sen)=='TTF 1st line (months)'] <- '1st line TTF'
colnames(ov_comb_tax_p3_sen)[colnames(ov_comb_tax_p3_sen)=='Tumor Stage'] <- 'Tumour stage'
colnames(ov_comb_tax_p3_sen)[colnames(ov_comb_tax_p3_sen)=='treatment_line'] <- 'Treatment line'
colnames(ov_comb_tax_p3_sen)[colnames(ov_comb_tax_p3_sen)=='wgii'] <- 'wGII'
ov_comb_tax_p3_sen$age_recoded <- as.character(car::recode(as.numeric(ov_comb_tax_p3_sen$'Age at diagnosis'), "lo:64=1; 65:hi=2"))
ov_comb_tax_p3_sen$year_of_treatment <- apply(ov_comb_tax_p3_sen,
                                              MARGIN=1,
                                              FUN=function(x) get_year_of_treatment(x[which(colnames(ov_comb_tax_p3_sen)=='bcr_patient_barcode')], x[which(colnames(ov_comb_tax_p3_sen)=='Treatment line')]))
ov_comb_tax_p3_sen_weights <- weightit(`Treatment arm` ~ year_of_treatment,
                                       data=ov_comb_tax_p3_sen,
                                       method='glm',
                                       estimand = 'ATT')
bal.tab(ov_comb_tax_p3_sen_weights,
        un=TRUE)
ov_comb_tax_p3_sen_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm`*`1st line TTF` + strata(age_recoded), data=ov_comb_tax_p3_sen, weights=ov_comb_tax_p3_sen_weights$weights)
cox.zph(ov_comb_tax_p3_sen_cox) #all ok!
# TCGA OV tax sensitive point estimate (with intEST function): 1.335 (0.450-3.967)

render_graph(phase_3_plot_func('OV', 'taxane', 'tcga-ov', SA_exp=FALSE, SA_control=FALSE)) %>%
  export_svg %>%
  charToRaw %>%
  rsvg_pdf(file.path(supp_figs_flow_dir, 'TCGA_OV_Comb_Tax_filtering_phase_3.pdf'))


# TCGA-OV Doxorubicin
ov_comb_dox_p3 <- phase_3_func('OV', 'doxorubicin', 'tcga-ov', SA_exp=FALSE, SA_control=FALSE)
saveRDS(ov_comb_dox_p3,file.path(figs_dir, 'TCGA_OV_Doxo_phase_3.rds'))
ov_comb_dox_p3_res <- ov_comb_dox_p3[[1]]
ov_comb_dox_p3_res <- left_join(ov_comb_dox_p3_res, wgii_tcga)
ov_comb_dox_p3_res$`Tumor Stage` <- gsub('[ABC]', '', ov_comb_dox_p3_res$`Tumor Stage`)
colnames(ov_comb_dox_p3_res)[colnames(ov_comb_dox_p3_res)=='Exp_treatment'] <- 'Treatment arm'
ov_comb_dox_p3_res$`Treatment arm` <- ifelse(ov_comb_dox_p3_res$`Treatment arm`, 'Experimental', 'Control')
colnames(ov_comb_dox_p3_res)[colnames(ov_comb_dox_p3_res)=='TTF 1st line (months)'] <- '1st line TTF'
colnames(ov_comb_dox_p3_res)[colnames(ov_comb_dox_p3_res)=='Tumor Stage'] <- 'Tumour stage'
colnames(ov_comb_dox_p3_res)[colnames(ov_comb_dox_p3_res)=='treatment_line'] <- 'Treatment line'
ov_comb_dox_p3_res$`Treatment line` <- as.numeric(substr(ov_comb_dox_p3_res$`Treatment line`, 1, 1))
ov_comb_dox_p3_res$treatment_line_recoded <- ov_comb_dox_p3_res$`Treatment line`
ov_comb_dox_p3_res$treatment_line_recoded[ov_comb_dox_p3_res$`Treatment line`>3]="4+"
colnames(ov_comb_dox_p3_res)[colnames(ov_comb_dox_p3_res)=='wgii'] <- 'wGII'
ov_comb_dox_p3_res$TTF_1st_recoded <- as.character(car::recode(as.numeric(ov_comb_dox_p3_res$'1st line TTF'), "lo:6=1; 6.01:hi=2"))
ov_comb_dox_p3_res$age_recoded <- as.character(car::recode(as.numeric(ov_comb_dox_p3_res$'Age at diagnosis'), "lo:64=1; 65:hi=2"))
ov_comb_dox_p3_res$year_of_treatment <- apply(ov_comb_dox_p3_res,
                                              MARGIN=1,
                                              FUN=function(x) get_year_of_treatment(x[which(colnames(ov_comb_dox_p3_res)=='bcr_patient_barcode')], x[which(colnames(ov_comb_dox_p3_res)=='Treatment line')]))
ov_comb_dox_p3_res_weights <- weightit(`Treatment arm` ~ year_of_treatment,
                                       data=ov_comb_dox_p3_res,
                                       method='glm',
                                       estimand = 'ATT')
bal.tab(ov_comb_dox_p3_res_weights,
        un=TRUE)
ov_comb_dox_p3_res_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + TTF_1st_recoded + strata(age_recoded), data=ov_comb_dox_p3_res, weights=ov_comb_dox_p3_res_weights$weights)
cox.zph(ov_comb_dox_p3_res_cox) #all ok! => 1st line cannot be added as interaction term because ph assumption & needs to be recoded

ov_comb_dox_p3_sen <- ov_comb_dox_p3[[2]]
ov_comb_dox_p3_sen$treatment_line <- as.numeric(substr(ov_comb_dox_p3_sen$treatment_line, 1, 1))
ov_comb_dox_p3_sen <- left_join(ov_comb_dox_p3_sen, wgii_tcga)
ov_comb_dox_p3_sen$`Tumor Stage` <- gsub('[ABC]', '', ov_comb_dox_p3_sen$`Tumor Stage`)
colnames(ov_comb_dox_p3_sen)[colnames(ov_comb_dox_p3_sen)=='Exp_treatment'] <- 'Treatment arm'
ov_comb_dox_p3_sen$`Treatment arm` <- ifelse(ov_comb_dox_p3_sen$`Treatment arm`, 'Experimental', 'Control')
colnames(ov_comb_dox_p3_sen)[colnames(ov_comb_dox_p3_sen)=='TTF 1st line (months)'] <- '1st line TTF'
colnames(ov_comb_dox_p3_sen)[colnames(ov_comb_dox_p3_sen)=='Tumor Stage'] <- 'Tumour stage'
colnames(ov_comb_dox_p3_sen)[colnames(ov_comb_dox_p3_sen)=='treatment_line'] <- 'Treatment line'
ov_comb_dox_p3_sen$`Treatment line` <- as.numeric(substr(ov_comb_dox_p3_sen$`Treatment line`, 1, 1))
ov_comb_dox_p3_sen$treatment_line_recoded<-ov_comb_dox_p3_sen$`Treatment line`
ov_comb_dox_p3_sen$treatment_line_recoded[ov_comb_dox_p3_sen$`Treatment line`>3]="4+"
colnames(ov_comb_dox_p3_sen)[colnames(ov_comb_dox_p3_sen)=='wgii'] <- 'wGII'
ov_comb_dox_p3_sen$TTF_1st_recoded <- as.character(car::recode(as.numeric(ov_comb_dox_p3_sen$'1st line TTF'), "lo:6=1; 6.01:hi=2"))
ov_comb_dox_p3_sen$age_recoded <- as.character(car::recode(as.numeric(ov_comb_dox_p3_sen$'Age at diagnosis'), "lo:64=1; 65:hi=2"))
ov_comb_dox_p3_sen$year_of_treatment <- apply(ov_comb_dox_p3_sen,
                                              MARGIN=1,
                                              FUN=function(x) get_year_of_treatment(x[which(colnames(ov_comb_dox_p3_sen)=='bcr_patient_barcode')], x[which(colnames(ov_comb_dox_p3_sen)=='Treatment line')]))
ov_comb_dox_p3_sen_weights <- weightit(`Treatment arm` ~ year_of_treatment,
                                       data=ov_comb_dox_p3_sen,
                                       method='glm',
                                       estimand = 'ATT')
bal.tab(ov_comb_dox_p3_sen_weights,
        un=TRUE)
ov_comb_dox_p3_sen_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + TTF_1st_recoded + strata(age_recoded), data=ov_comb_dox_p3_sen, weights=ov_comb_dox_p3_sen_weights$weights)
cox.zph(ov_comb_dox_p3_sen_cox) #all ok!

render_graph(phase_3_plot_func('OV', 'doxorubicin', 'tcga-ov', SA_exp=FALSE, SA_control=FALSE)) %>%
  export_svg %>%
  charToRaw %>%
  rsvg_pdf(file.path(supp_figs_flow_dir, 'TCGA_OV_Comb_Anth_filtering_phase_3.pdf'))


##### PRIMARY TUMOUR SINGLE-AGENT #####
# TCGA-OV Taxane
ov_sa_tax_p3 <- phase_3_func('OV', 'taxane', 'tcga-ov', SA_exp=TRUE, SA_control=FALSE)
ov_sa_tax_p3_res <- ov_sa_tax_p3[[1]]
ov_sa_tax_p3_res$treatment_line <- as.numeric(substr(ov_sa_tax_p3_res$treatment_line, 1, 1))
ov_sa_tax_p3_res <- left_join(ov_sa_tax_p3_res, wgii_tcga)
ov_sa_tax_p3_res$`Tumor Stage` <- gsub('[ABC]', '', ov_sa_tax_p3_res$`Tumor Stage`)
colnames(ov_sa_tax_p3_res)[colnames(ov_sa_tax_p3_res)=='Exp_treatment'] <- 'Treatment arm'
ov_sa_tax_p3_res$`Treatment arm` <- ifelse(ov_sa_tax_p3_res$`Treatment arm`, 'Experimental', 'Control')
colnames(ov_sa_tax_p3_res)[colnames(ov_sa_tax_p3_res)=='TTF 1st line (months)'] <- '1st line TTF'
colnames(ov_sa_tax_p3_res)[colnames(ov_sa_tax_p3_res)=='Tumor Stage'] <- 'Tumour stage'
colnames(ov_sa_tax_p3_res)[colnames(ov_sa_tax_p3_res)=='treatment_line'] <- 'Treatment line'
colnames(ov_sa_tax_p3_res)[colnames(ov_sa_tax_p3_res)=='wgii'] <- 'wGII'
ov_sa_tax_p3_res$age_recoded <- as.character(car::recode(as.numeric(ov_sa_tax_p3_res$'Age at diagnosis'), "lo:64=1; 65:hi=2"))
ov_sa_tax_p3_res$year_of_treatment <- apply(ov_sa_tax_p3_res,
                                              MARGIN=1,
                                              FUN=function(x) get_year_of_treatment(x[which(colnames(ov_sa_tax_p3_res)=='bcr_patient_barcode')], x[which(colnames(ov_sa_tax_p3_res)=='Treatment line')]))
ov_sa_tax_p3_res_weights <- weightit(`Treatment arm` ~ year_of_treatment,
                                       data=ov_sa_tax_p3_res,
                                       method='glm',
                                       estimand = 'ATT')
bal.tab(ov_sa_tax_p3_res_weights,
        un=TRUE)
ov_sa_tax_p3_res_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm`*`1st line TTF` + strata(age_recoded), data=ov_sa_tax_p3_res, weights=ov_sa_tax_p3_res_weights$weights)
cox.zph(ov_sa_tax_p3_res_cox) #ok!

ov_sa_tax_p3_sen <- ov_sa_tax_p3[[2]]
ov_sa_tax_p3_sen$treatment_line <- as.numeric(substr(ov_sa_tax_p3_sen$treatment_line, 1, 1))
ov_sa_tax_p3_sen <- left_join(ov_sa_tax_p3_sen, wgii_tcga)
ov_sa_tax_p3_sen$`Tumor Stage` <- gsub('[ABC]', '', ov_sa_tax_p3_sen$`Tumor Stage`)
colnames(ov_sa_tax_p3_sen)[colnames(ov_sa_tax_p3_sen)=='Exp_treatment'] <- 'Treatment arm'
ov_sa_tax_p3_sen$`Treatment arm` <- ifelse(ov_sa_tax_p3_sen$`Treatment arm`, 'Experimental', 'Control')
colnames(ov_sa_tax_p3_sen)[colnames(ov_sa_tax_p3_sen)=='TTF 1st line (months)'] <- '1st line TTF'
colnames(ov_sa_tax_p3_sen)[colnames(ov_sa_tax_p3_sen)=='Tumor Stage'] <- 'Tumour stage'
colnames(ov_sa_tax_p3_sen)[colnames(ov_sa_tax_p3_sen)=='treatment_line'] <- 'Treatment line'
colnames(ov_sa_tax_p3_sen)[colnames(ov_sa_tax_p3_sen)=='wgii'] <- 'wGII'
ov_sa_tax_p3_sen$age_recoded <- as.character(car::recode(as.numeric(ov_sa_tax_p3_sen$'Age at diagnosis'), "lo:64=1; 65:hi=2"))
ov_sa_tax_p3_sen$year_of_treatment <- apply(ov_sa_tax_p3_sen,
                                            MARGIN=1,
                                            FUN=function(x) get_year_of_treatment(x[which(colnames(ov_sa_tax_p3_sen)=='bcr_patient_barcode')], x[which(colnames(ov_sa_tax_p3_sen)=='Treatment line')]))
ov_sa_tax_p3_sen_weights <- weightit(`Treatment arm` ~ year_of_treatment,
                                     data=ov_sa_tax_p3_sen,
                                     method='glm',
                                     estimand = 'ATT')
bal.tab(ov_sa_tax_p3_sen_weights,
        un=TRUE)
ov_sa_tax_p3_sen_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm`*`1st line TTF` + strata(age_recoded), data=ov_sa_tax_p3_sen, weights=ov_sa_tax_p3_sen_weights$weights)
cox.zph(ov_sa_tax_p3_sen_cox) #all ok!



# TCGA-OV Doxorubicin
ov_sa_dox_p3 <- phase_3_func('OV', 'doxorubicin', 'tcga-ov', SA_exp=TRUE, SA_control=FALSE)
ov_sa_dox_p3_res <- ov_sa_dox_p3[[1]]
ov_sa_dox_p3_res <- left_join(ov_sa_dox_p3_res, wgii_tcga)
ov_sa_dox_p3_res$`Tumor Stage` <- gsub('[ABC]', '', ov_sa_dox_p3_res$`Tumor Stage`)
colnames(ov_sa_dox_p3_res)[colnames(ov_sa_dox_p3_res)=='Exp_treatment'] <- 'Treatment arm'
ov_sa_dox_p3_res$`Treatment arm` <- ifelse(ov_sa_dox_p3_res$`Treatment arm`, 'Experimental', 'Control')
colnames(ov_sa_dox_p3_res)[colnames(ov_sa_dox_p3_res)=='TTF 1st line (months)'] <- '1st line TTF'
colnames(ov_sa_dox_p3_res)[colnames(ov_sa_dox_p3_res)=='Tumor Stage'] <- 'Tumour stage'
colnames(ov_sa_dox_p3_res)[colnames(ov_sa_dox_p3_res)=='treatment_line'] <- 'Treatment line'
ov_sa_dox_p3_res$`Treatment line` <- as.numeric(substr(ov_sa_dox_p3_res$`Treatment line`, 1, 1))
ov_sa_dox_p3_res$treatment_line_recoded<-ov_sa_dox_p3_res$`Treatment line`
ov_sa_dox_p3_res$treatment_line_recoded[ov_sa_dox_p3_res$`Treatment line`>3]="4+"
colnames(ov_sa_dox_p3_res)[colnames(ov_sa_dox_p3_res)=='wgii'] <- 'wGII'
ov_sa_dox_p3_res$TTF_1st_recoded <- as.character(car::recode(as.numeric(ov_sa_dox_p3_res$'1st line TTF'), "lo:6=1; 6.01:hi=2"))
ov_sa_dox_p3_res$age_recoded <- as.character(car::recode(as.numeric(ov_sa_dox_p3_res$'Age at diagnosis'), "lo:64=1; 65:hi=2"))
ov_sa_dox_p3_res$year_of_treatment <- apply(ov_sa_dox_p3_res,
                                              MARGIN=1,
                                              FUN=function(x) get_year_of_treatment(x[which(colnames(ov_sa_dox_p3_res)=='bcr_patient_barcode')], x[which(colnames(ov_sa_dox_p3_res)=='Treatment line')]))
ov_sa_dox_p3_res_weights <- weightit(`Treatment arm` ~ year_of_treatment,
                                       data=ov_sa_dox_p3_res,
                                       method='glm',
                                       estimand = 'ATT')
bal.tab(ov_sa_dox_p3_res_weights,
        un=TRUE)
ov_sa_dox_p3_res_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + TTF_1st_recoded + strata(age_recoded), data=ov_sa_dox_p3_res, weights=ov_sa_dox_p3_res_weights$weights)
cox.zph(ov_sa_dox_p3_res_cox) #global model!

ov_sa_dox_p3_sen <- ov_sa_dox_p3[[2]]
ov_sa_dox_p3_sen$treatment_line <- as.numeric(substr(ov_sa_dox_p3_sen$treatment_line, 1, 1))
ov_sa_dox_p3_sen <- left_join(ov_sa_dox_p3_sen, wgii_tcga)
ov_sa_dox_p3_sen$`Tumor Stage` <- gsub('[ABC]', '', ov_sa_dox_p3_sen$`Tumor Stage`)
colnames(ov_sa_dox_p3_sen)[colnames(ov_sa_dox_p3_sen)=='Exp_treatment'] <- 'Treatment arm'
ov_sa_dox_p3_sen$`Treatment arm` <- ifelse(ov_sa_dox_p3_sen$`Treatment arm`, 'Experimental', 'Control')
colnames(ov_sa_dox_p3_sen)[colnames(ov_sa_dox_p3_sen)=='TTF 1st line (months)'] <- '1st line TTF'
colnames(ov_sa_dox_p3_sen)[colnames(ov_sa_dox_p3_sen)=='Tumor Stage'] <- 'Tumour stage'
colnames(ov_sa_dox_p3_sen)[colnames(ov_sa_dox_p3_sen)=='treatment_line'] <- 'Treatment line'
ov_sa_dox_p3_sen$`Treatment line` <- as.numeric(substr(ov_sa_dox_p3_sen$`Treatment line`, 1, 1))
ov_sa_dox_p3_sen$treatment_line_recoded<-ov_sa_dox_p3_sen$`Treatment line`
ov_sa_dox_p3_sen$treatment_line_recoded[ov_sa_dox_p3_sen$`Treatment line`>3]="4+"
colnames(ov_sa_dox_p3_sen)[colnames(ov_sa_dox_p3_sen)=='wgii'] <- 'wGII'
ov_sa_dox_p3_sen$TTF_1st_recoded <- as.character(car::recode(as.numeric(ov_sa_dox_p3_sen$'1st line TTF'), "lo:6=1; 6.01:hi=2"))
ov_sa_dox_p3_sen$age_recoded <- as.character(car::recode(as.numeric(ov_sa_dox_p3_sen$'Age at diagnosis'), "lo:64=1; 65:hi=2"))
ov_sa_dox_p3_sen$year_of_treatment <- apply(ov_sa_dox_p3_sen,
                                            MARGIN=1,
                                            FUN=function(x) get_year_of_treatment(x[which(colnames(ov_sa_dox_p3_sen)=='bcr_patient_barcode')], x[which(colnames(ov_sa_dox_p3_sen)=='Treatment line')]))
ov_sa_dox_p3_sen_weights <- weightit(`Treatment arm` ~ year_of_treatment,
                                     data=ov_sa_dox_p3_sen,
                                     method='glm',
                                     estimand = 'ATT')
bal.tab(ov_sa_dox_p3_sen_weights,
        un=TRUE)
ov_sa_dox_p3_sen_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + TTF_1st_recoded + strata(age_recoded), data=ov_sa_dox_p3_sen, weights=ov_sa_dox_p3_sen_weights$weights)
cox.zph(ov_sa_dox_p3_sen_cox) #all ok!


##### METASTATIC TUMOUR COMBINATION THERAPY #####

# HMF prostate Taxane
pro_comb_tax_p3 <- phase_3_func('Prostate', 'taxane', 'hmf', SA_exp=FALSE, SA_control=FALSE)
pro_comb_tax_p3_res <- pro_comb_tax_p3[[1]]
pro_comb_tax_p3_res <- left_join(pro_comb_tax_p3_res, wgii_hmf)
colnames(pro_comb_tax_p3_res)[colnames(pro_comb_tax_p3_res)=='Exp_treatment'] <- 'Treatment arm'
pro_comb_tax_p3_res$`Treatment arm` <- ifelse(pro_comb_tax_p3_res$`Treatment arm`, 'Experimental', 'Control')
colnames(pro_comb_tax_p3_res)[colnames(pro_comb_tax_p3_res)=='age'] <- 'Age at diagnosis'
pro_comb_tax_p3_res$age_recoded <- as.character(car::recode(as.numeric(pro_comb_tax_p3_res$'Age at diagnosis'), "lo:69=1; 70:hi=2"))
colnames(pro_comb_tax_p3_res)[colnames(pro_comb_tax_p3_res)=='treatment_line'] <- 'Treatment line'
colnames(pro_comb_tax_p3_res)[colnames(pro_comb_tax_p3_res)=='wgii'] <- 'wGII'
pro_comb_tax_p3_res$biopsyYear <- as.numeric(substr(pro_comb_tax_p3_res$biopsyDate, 1, 4))
pro_comb_tax_p3_res_weights <- weightit(`Treatment arm` ~ biopsyYear,
                                        data=pro_comb_tax_p3_res,
                                        method='glm',
                                        estimand = 'ATT')
bal.tab(pro_comb_tax_p3_res_weights,
        un=TRUE)
pro_comb_tax_p3_res_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + `Age at diagnosis`, data=pro_comb_tax_p3_res, weights=pro_comb_tax_p3_res_weights$weights)
cox.zph(pro_comb_tax_p3_res_cox) #all ok!

pro_comb_tax_p3_sen <- pro_comb_tax_p3[[2]]
pro_comb_tax_p3_sen <- left_join(pro_comb_tax_p3_sen, wgii_hmf)
colnames(pro_comb_tax_p3_sen)[colnames(pro_comb_tax_p3_sen)=='Exp_treatment'] <- 'Treatment arm'
pro_comb_tax_p3_sen$`Treatment arm` <- ifelse(pro_comb_tax_p3_sen$`Treatment arm`, 'Experimental', 'Control')
colnames(pro_comb_tax_p3_sen)[colnames(pro_comb_tax_p3_sen)=='age'] <- 'Age at diagnosis'
pro_comb_tax_p3_sen$age_recoded <- as.character(car::recode(as.numeric(pro_comb_tax_p3_sen$'Age at diagnosis'), "lo:54=1; 55:hi=2"))
colnames(pro_comb_tax_p3_sen)[colnames(pro_comb_tax_p3_sen)=='treatment_line'] <- 'Treatment line'
colnames(pro_comb_tax_p3_sen)[colnames(pro_comb_tax_p3_sen)=='wgii'] <- 'wGII'
pro_comb_tax_p3_sen$biopsyYear <- as.numeric(substr(pro_comb_tax_p3_sen$biopsyDate, 1, 4))
pro_comb_tax_p3_sen_weights <- weightit(`Treatment arm` ~ biopsyYear,
                                        data=pro_comb_tax_p3_sen,
                                        method='glm',
                                        estimand = 'ATT')
bal.tab(pro_comb_tax_p3_sen_weights,
        un=TRUE)
pro_comb_tax_p3_sen_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + `Age at diagnosis`, data=pro_comb_tax_p3_sen, weights=pro_comb_tax_p3_sen_weights$weights)
cox.zph(pro_comb_tax_p3_sen_cox) #global!

render_graph(phase_3_plot_func('Prostate', 'taxane', 'hmf', SA_exp=FALSE, SA_control=FALSE)) %>%
  export_svg %>%
  charToRaw %>%
  rsvg_pdf(file.path(supp_figs_flow_dir, 'HMF_Prostate_Comb_Tax_filtering_phase_3.pdf'))


# HMF breast Taxane
bre_comb_tax_p3 <- phase_3_func('Breast', 'taxane', 'hmf', SA_exp=FALSE, SA_control=FALSE, control_pattern='citabine')
bre_comb_tax_p3_res <- bre_comb_tax_p3[[1]]
bre_comb_tax_p3_res <- bre_comb_tax_p3_res %>% filter(name!='Eftilagimod alpha or placebo/Paclitaxel')
bre_comb_tax_p3_res <- left_join(bre_comb_tax_p3_res, wgii_hmf)
colnames(bre_comb_tax_p3_res)[colnames(bre_comb_tax_p3_res)=='Exp_treatment'] <- 'Treatment arm'
bre_comb_tax_p3_res$`Treatment arm` <- ifelse(bre_comb_tax_p3_res$`Treatment arm`, 'Experimental', 'Control')
colnames(bre_comb_tax_p3_res)[colnames(bre_comb_tax_p3_res)=='age'] <- 'Age at diagnosis'
bre_comb_tax_p3_res$age_recoded <- as.character(car::recode(as.numeric(bre_comb_tax_p3_res$'Age at diagnosis'), "lo:59=1; 60:hi=2"))
colnames(bre_comb_tax_p3_res)[colnames(bre_comb_tax_p3_res)=='primaryTumorSubType'] <- 'Tumour subtype'
bre_comb_tax_p3_res$`Tumour subtype` <- gsub('-negative', '-', gsub('-positive', '+', bre_comb_tax_p3_res$`Tumour subtype`))
bre_comb_tax_p3_res$`Tumour subtype`[bre_comb_tax_p3_res$`Tumour subtype`=='Triple negative'] <- 'TNBC'
colnames(bre_comb_tax_p3_res)[colnames(bre_comb_tax_p3_res)=='treatment_line'] <- 'Treatment line'
colnames(bre_comb_tax_p3_res)[colnames(bre_comb_tax_p3_res)=='wgii'] <- 'wGII'
bre_comb_tax_p3_res$biopsyYear <- as.numeric(substr(bre_comb_tax_p3_res$biopsyDate, 1, 4))
bre_comb_tax_p3_res_weights <- weightit(`Treatment arm` ~ biopsyYear,
                                        data=bre_comb_tax_p3_res,
                                        method='glm',
                                        estimand = 'ATT')
bal.tab(bre_comb_tax_p3_res_weights,
        un=TRUE)
bre_comb_tax_p3_res_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + `Age at diagnosis`, data=bre_comb_tax_p3_res, weights=bre_comb_tax_p3_res_weights$weights)
cox.zph(bre_comb_tax_p3_res_cox) #all ok!

bre_comb_tax_p3_sen <- bre_comb_tax_p3[[2]]
bre_comb_tax_p3_sen <- bre_comb_tax_p3_sen %>% filter(name!='Mucopolysaccharidosis I (MPS I) inhibitor/Paclitaxel' & name!='Carboplatin/Gemcitabine/Paclitaxel')
bre_comb_tax_p3_sen <- left_join(bre_comb_tax_p3_sen, wgii_hmf)
colnames(bre_comb_tax_p3_sen)[colnames(bre_comb_tax_p3_sen)=='Exp_treatment'] <- 'Treatment arm'
bre_comb_tax_p3_sen$`Treatment arm` <- ifelse(bre_comb_tax_p3_sen$`Treatment arm`, 'Experimental', 'Control')
colnames(bre_comb_tax_p3_sen)[colnames(bre_comb_tax_p3_sen)=='age'] <- 'Age at diagnosis'
bre_comb_tax_p3_sen$age_recoded <- as.character(car::recode(as.numeric(bre_comb_tax_p3_sen$'Age at diagnosis'), "lo:59=1; 60:hi=2"))
colnames(bre_comb_tax_p3_sen)[colnames(bre_comb_tax_p3_sen)=='primaryTumorSubType'] <- 'Tumour subtype'
bre_comb_tax_p3_sen$`Tumour subtype` <- gsub('-negative', '-', gsub('-positive', '+', bre_comb_tax_p3_sen$`Tumour subtype`))
bre_comb_tax_p3_sen$`Tumour subtype`[bre_comb_tax_p3_sen$`Tumour subtype`=='Triple negative'] <- 'TNBC'
colnames(bre_comb_tax_p3_sen)[colnames(bre_comb_tax_p3_sen)=='treatment_line'] <- 'Treatment line'
colnames(bre_comb_tax_p3_sen)[colnames(bre_comb_tax_p3_sen)=='wgii'] <- 'wGII'
bre_comb_tax_p3_sen$biopsyYear <- as.numeric(substr(bre_comb_tax_p3_sen$biopsyDate, 1, 4))
bre_comb_tax_p3_sen_weights <- weightit(`Treatment arm` ~ biopsyYear,
                                        data=bre_comb_tax_p3_sen,
                                        method='glm',
                                        estimand = 'ATT')
bal.tab(bre_comb_tax_p3_sen_weights,
        un=TRUE)
bre_comb_tax_p3_sen_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + `Age at diagnosis`, data=bre_comb_tax_p3_sen, weights=bre_comb_tax_p3_sen_weights$weights)
cox.zph(bre_comb_tax_p3_sen_cox) #all ok!

render_graph(phase_3_plot_func('Breast', 'taxane', 'hmf', SA_exp=FALSE, SA_control=FALSE, control_pattern='citabine')) %>%
  export_svg %>%
  charToRaw %>%
  rsvg_pdf(file.path(supp_figs_flow_dir, 'HMF_Breast_Comb_Tax_filtering_phase_3.pdf'))


# HMF breast Doxo
bre_comb_dox_p3 <- phase_3_func('Breast', 'doxorubicin', 'hmf', SA_exp=FALSE, SA_control=FALSE, control_pattern='citabine')
bre_comb_dox_p3_res <- bre_comb_dox_p3[[1]]
bre_comb_dox_p3_res <- bre_comb_dox_p3_res %>% filter(name!='Eftilagimod alpha or placebo/Paclitaxel')
bre_comb_dox_p3_res <- left_join(bre_comb_dox_p3_res, wgii_hmf)
colnames(bre_comb_dox_p3_res)[colnames(bre_comb_dox_p3_res)=='Exp_treatment'] <- 'Treatment arm'
bre_comb_dox_p3_res$`Treatment arm` <- ifelse(bre_comb_dox_p3_res$`Treatment arm`, 'Experimental', 'Control')
colnames(bre_comb_dox_p3_res)[colnames(bre_comb_dox_p3_res)=='age'] <- 'Age at diagnosis'
bre_comb_dox_p3_res$age_recoded <- as.character(car::recode(as.numeric(bre_comb_dox_p3_res$'Age at diagnosis'), "lo:59=1; 60:hi=2"))
colnames(bre_comb_dox_p3_res)[colnames(bre_comb_dox_p3_res)=='primaryTumorSubType'] <- 'Tumour subtype'
bre_comb_dox_p3_res$`Tumour subtype` <- gsub('-negative', '-', gsub('-positive', '+', bre_comb_dox_p3_res$`Tumour subtype`))
bre_comb_dox_p3_res$`Tumour subtype`[bre_comb_dox_p3_res$`Tumour subtype`=='Triple negative'] <- 'TNBC'
colnames(bre_comb_dox_p3_res)[colnames(bre_comb_dox_p3_res)=='treatment_line'] <- 'Treatment line'
colnames(bre_comb_dox_p3_res)[colnames(bre_comb_dox_p3_res)=='wgii'] <- 'wGII'
bre_comb_dox_p3_res$biopsyYear <- as.numeric(substr(bre_comb_dox_p3_res$biopsyDate, 1, 4))
bre_comb_dox_p3_res_weights <- weightit(`Treatment arm` ~ biopsyYear,
                                        data=bre_comb_dox_p3_res,
                                        method='glm',
                                        estimand = 'ATT')
bal.tab(bre_comb_dox_p3_res_weights,
        un=TRUE)
bre_comb_dox_p3_res_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + `Age at diagnosis`, data=bre_comb_dox_p3_res, weights=bre_comb_dox_p3_res_weights$weights)
cox.zph(bre_comb_dox_p3_res_cox) #all ok!

bre_comb_dox_p3_sen <- bre_comb_dox_p3[[2]]
bre_comb_dox_p3_sen <- bre_comb_dox_p3_sen %>% filter(name!='Eftilagimod alpha or placebo/Paclitaxel')
bre_comb_dox_p3_sen <- left_join(bre_comb_dox_p3_sen, wgii_hmf)
colnames(bre_comb_dox_p3_sen)[colnames(bre_comb_dox_p3_sen)=='Exp_treatment'] <- 'Treatment arm'
bre_comb_dox_p3_sen$`Treatment arm` <- ifelse(bre_comb_dox_p3_sen$`Treatment arm`, 'Experimental', 'Control')
colnames(bre_comb_dox_p3_sen)[colnames(bre_comb_dox_p3_sen)=='age'] <- 'Age at diagnosis'
bre_comb_dox_p3_sen$age_recoded <- as.character(car::recode(as.numeric(bre_comb_dox_p3_sen$'Age at diagnosis'), "lo:59=1; 60:hi=2"))
colnames(bre_comb_dox_p3_sen)[colnames(bre_comb_dox_p3_sen)=='primaryTumorSubType'] <- 'Tumour subtype'
bre_comb_dox_p3_sen$`Tumour subtype` <- gsub('-negative', '-', gsub('-positive', '+', bre_comb_dox_p3_sen$`Tumour subtype`))
bre_comb_dox_p3_sen$`Tumour subtype`[bre_comb_dox_p3_sen$`Tumour subtype`=='Triple negative'] <- 'TNBC'
colnames(bre_comb_dox_p3_sen)[colnames(bre_comb_dox_p3_sen)=='treatment_line'] <- 'Treatment line'
colnames(bre_comb_dox_p3_sen)[colnames(bre_comb_dox_p3_sen)=='wgii'] <- 'wGII'
bre_comb_dox_p3_sen_weights <- weightit(`Treatment arm` ~ `Age at diagnosis` + `Treatment line` + `Tumour subtype`,
                                        data=bre_comb_dox_p3_sen,
                                        method='glm',
                                        estimand = 'ATT')
bal.tab(bre_comb_dox_p3_sen_weights,
        un=TRUE)
bre_comb_dox_p3_sen_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + `Age at diagnosis` + wGII, data=bre_comb_dox_p3_sen, weights=bre_comb_dox_p3_sen_weights$weights)
cox.zph(bre_comb_dox_p3_sen_cox) #all ok!

render_graph(phase_3_plot_func('Breast', 'doxorubicin', 'hmf', SA_exp=FALSE, SA_control=FALSE, control_pattern='citabine')) %>%
  export_svg %>%
  charToRaw %>%
  rsvg_pdf(file.path(supp_figs_flow_dir, 'HMF_Breast_Comb_Anth_filtering_phase_3.pdf'))


# HMF ovary Doxo => we are not powered when platinum is excluded in the experimental and control arm (following the same than in TCGA)
ovary_comb_dox_p3 <- phase_3_func('Ovary', 'doxorubicin', 'hmf', SA_exp=FALSE, SA_control=FALSE)
ovary_comb_dox_p3_res <- ovary_comb_dox_p3[[1]]
ovary_comb_dox_p3_res <- ovary_comb_dox_p3_res %>% filter(!grepl('platin', name))

# TCGA UCEC Doxo => there are only 8 patients treated with dox and the model is not proportional => excluded
ucec_comb_dox_p3 <- phase_3_func('UCEC', 'doxorubicin', 'tcga', SA_exp=FALSE, SA_control=FALSE)
ucec_comb_dox_p3_res <- ucec_comb_dox_p3[[1]]
ucec_comb_dox_p3_res <- left_join(ucec_comb_dox_p3_res, wgii_hmf)
colnames(ucec_comb_dox_p3_res)[colnames(ucec_comb_dox_p3_res)=='Exp_treatment'] <- 'Treatment arm'
ucec_comb_dox_p3_res$`Treatment arm` <- ifelse(ucec_comb_dox_p3_res$`Treatment arm`, 'Experimental', 'Control')
table(ucec_comb_dox_p3_res$`Treatment arm`)
ucec_comb_dox_p3_res_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm`, data=ucec_comb_dox_p3_res)
cox.zph(ucec_comb_dox_p3_res_cox) #not proportional! not good representation of both arms


##### METASTATIC SINGLE-AGENT #####
# HMF prostate Taxane
pro_sa_tax_p3 <- phase_3_func('Prostate', 'taxane', 'hmf', SA_exp=TRUE, SA_control=FALSE)
pro_sa_tax_p3_res <- pro_sa_tax_p3[[1]]
pro_sa_tax_p3_res <- left_join(pro_sa_tax_p3_res, wgii_hmf)
colnames(pro_sa_tax_p3_res)[colnames(pro_sa_tax_p3_res)=='Exp_treatment'] <- 'Treatment arm'
pro_sa_tax_p3_res$`Treatment arm` <- ifelse(pro_sa_tax_p3_res$`Treatment arm`, 'Experimental', 'Control')
colnames(pro_sa_tax_p3_res)[colnames(pro_sa_tax_p3_res)=='age'] <- 'Age at diagnosis'
pro_sa_tax_p3_res$age_recoded <- as.character(car::recode(as.numeric(pro_sa_tax_p3_res$'Age at diagnosis'), "lo:69=1; 70:hi=2"))
colnames(pro_sa_tax_p3_res)[colnames(pro_sa_tax_p3_res)=='treatment_line'] <- 'Treatment line'
colnames(pro_sa_tax_p3_res)[colnames(pro_sa_tax_p3_res)=='wgii'] <- 'wGII'
pro_sa_tax_p3_res_weights <- weightit(`Treatment arm` ~ `Age at diagnosis` + `Treatment line`,
                                        data=pro_sa_tax_p3_res,
                                        method='glm',
                                        estimand = 'ATT')
bal.tab(pro_sa_tax_p3_res_weights,
        un=TRUE)
pro_sa_tax_p3_res_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm` + `Age at diagnosis`, data=pro_sa_tax_p3_res, weights=pro_sa_tax_p3_res_weights$weights)
cox.zph(pro_sa_tax_p3_res_cox) #all ok!

pro_sa_tax_p3_sen <- pro_sa_tax_p3[[2]]
pro_sa_tax_p3_sen <- left_join(pro_sa_tax_p3_sen, wgii_hmf)
colnames(pro_sa_tax_p3_sen)[colnames(pro_sa_tax_p3_sen)=='Exp_treatment'] <- 'Treatment arm'
pro_sa_tax_p3_sen$`Treatment arm` <- ifelse(pro_sa_tax_p3_sen$`Treatment arm`, 'Experimental', 'Control')
colnames(pro_sa_tax_p3_sen)[colnames(pro_sa_tax_p3_sen)=='age'] <- 'Age at diagnosis'
pro_sa_tax_p3_sen$age_recoded <- as.character(car::recode(as.numeric(pro_sa_tax_p3_sen$'Age at diagnosis'), "lo:54=1; 55:hi=2"))
colnames(pro_sa_tax_p3_sen)[colnames(pro_sa_tax_p3_sen)=='treatment_line'] <- 'Treatment line'
colnames(pro_sa_tax_p3_sen)[colnames(pro_sa_tax_p3_sen)=='wgii'] <- 'wGII'
pro_sa_tax_p3_sen_weights <- weightit(`Treatment arm` ~ `Age at diagnosis` + `Treatment line`,
                                      data=pro_sa_tax_p3_sen,
                                      method='glm',
                                      estimand = 'ATT')
bal.tab(pro_sa_tax_p3_sen_weights,
        un=TRUE)
pro_sa_tax_p3_sen_cox <- coxph(Surv(TTF, censoring) ~ `Treatment arm`+ `Age at diagnosis` , data=pro_sa_tax_p3_sen, weights=pro_sa_tax_p3_sen_weights$weights)
cox.zph(pro_sa_tax_p3_sen_cox) #all ok!


##### PLOTTING #################################################################

# Total number of samples
tcga_ids=unique(c(tcga_ov_plat_p2$bcr_patient_barcode,tcga_sarc_dox_p2$bcr_patient_barcode,
      ov_comb_dox_p3_res$bcr_patient_barcode,ov_comb_dox_p3_sen$bcr_patient_barcode,
      ov_comb_tax_p3_res$bcr_patient_barcode,ov_comb_tax_p3_sen$bcr_patient_barcode))

hmf_ids=unique(c(pro_comb_tax_p3_res$patientIdentifier,pro_comb_tax_p3_res$patientIdentifier,
             bre_comb_dox_p3_res$patientIdentifier,bre_comb_dox_p3_sen$patientIdentifier,
             bre_comb_tax_p3_res$patientIdentifier,bre_comb_tax_p3_sen$patientIdentifier))
length(c(tcga_ids,hmf_ids))

#### FIGURE 4 PHASE II RESULTS ####

ov_plat_upper_hr <- summary(tcga_ov_plat_p2_cox)$conf.int[1, 4]
ov_plat_lower_hr <- summary(tcga_ov_plat_p2_cox)$conf.int[1, 3]
ov_plat_hr <- summary(tcga_ov_plat_p2_cox)$coefficients[1, 2]
ov_plat_pval <- summary(tcga_ov_plat_p2_cox)$coefficients[1, 5]
ov_plat_n <- summary(tcga_ov_plat_p2_cox)$n
ov_plat_df <- data.frame(test = "PhaseII",
                             lower_hr = ov_plat_lower_hr,
                             upper_hr = ov_plat_upper_hr,
                             hr = ov_plat_hr)
ov_plat_2 <- ggplot(ov_plat_df, aes(x = test)) +
  geom_errorbar(aes(ymin = lower_hr, ymax = upper_hr),
                width = 0.15, linewidth = 1) +
  geom_point(aes(y = hr), size = 10) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1) +
  scale_y_log10(position = 'right') +
  coord_flip() +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, vjust = 1))
ggsave(file.path(figs_dir, 'ov_plat_p2_hr.svg'),
       ov_plat_2, device = 'svg', width = 225, height = 150, units = 'mm')


sarc_dox_upper_hr <- summary(tcga_sarc_dox_p2_cox)$conf.int[1, 4]
sarc_dox_lower_hr <- summary(tcga_sarc_dox_p2_cox)$conf.int[1, 3]
sarc_dox_hr <- summary(tcga_sarc_dox_p2_cox)$coefficients[1, 2]
sarc_dox_pval <- summary(tcga_sarc_dox_p2_cox)$coefficients[1, 5]
sarc_dox_n <- summary(tcga_sarc_dox_p2_cox)$n
sarc_dox_df <- data.frame(test = "PhaseII",
                         lower_hr = sarc_dox_upper_hr,
                         upper_hr = sarc_dox_lower_hr,
                         hr = sarc_dox_hr)
sarc_dox_2 <- ggplot(sarc_dox_df, aes(x = test)) +
  geom_errorbar(aes(ymin = lower_hr, ymax = upper_hr),
                width = 0.15, linewidth = 1) +
  geom_point(aes(y = hr), size = 10) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1) +
  scale_y_log10(position = 'right') +
  coord_flip() +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, vjust = 1))
ggsave(file.path(figs_dir, 'sarc_dox_p2_hr.svg'),
       sarc_dox_2, device = 'svg', width = 225, height = 150, units = 'mm')


#### FIGURE 3 PHASE III RESULTS ####

ov_comb_dox_upper_hr <- c(summary(ov_comb_dox_p3_res_cox)$conf.int[1, 4],
                  summary(ov_comb_dox_p3_sen_cox)$conf.int[1, 4])
ov_comb_dox_lower_hr <- c(summary(ov_comb_dox_p3_res_cox)$conf.int[1, 3],
                  summary(ov_comb_dox_p3_sen_cox)$conf.int[1, 3])
ov_comb_dox_hr <- c(summary(ov_comb_dox_p3_res_cox)$coefficients[1, 2],
            summary(ov_comb_dox_p3_sen_cox)$coefficients[1, 2])
ov_comb_dox_pval <- c(summary(ov_comb_dox_p3_res_cox)$coefficients[1, 5],
              summary(ov_comb_dox_p3_sen_cox)$coefficients[1, 5])
ov_comb_dox_n <- c(summary(ov_comb_dox_p3_res_cox)$n,
           summary(ov_comb_dox_p3_sen_cox)$n)
ov_comb_dox_df <- data.frame(test = factor(c('Resistant',
                                     'Sensitive'),
                                   levels = c('Sensitive',
                                              'Resistant')),
                     lower_hr = ov_comb_dox_lower_hr,
                     upper_hr = ov_comb_dox_upper_hr,
                     hr = ov_comb_dox_hr)
ov_comb_dox_3 <- ggplot(ov_comb_dox_df, aes(x = test)) +
  geom_errorbar(aes(ymin = lower_hr, ymax = upper_hr),
                width = 0.15, linewidth = 1) +
  geom_point(aes(y = hr), size = 10) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1) +
  scale_y_log10(position = 'right') +
  coord_flip() +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, vjust = 1))
ggsave(file.path(figs_dir, 'ov_comb_dox_p3_hr.svg'),
       ov_comb_dox_3, device = 'svg', width = 225, height = 150, units = 'mm')


ov_comb_tax_upper_hr <- c(20.458,3.967) #95CI with intEST function at 6months
ov_comb_tax_lower_hr <- c(2.782,0.450) #95CI with intEST function at 6months
ov_comb_tax_hr <- c(7.435,1.335) #HR with intEST function at 6months
ov_comb_tax_pval <- c(summary(ov_comb_tax_p3_res_cox)$coefficients[1, 5],
              summary(ov_comb_tax_p3_sen_cox)$coefficients[1, 5])
ov_comb_tax_n <- c(summary(ov_comb_tax_p3_res_cox)$n,
           summary(ov_comb_tax_p3_sen_cox)$n)
ov_comb_tax_df <- data.frame(test = factor(c('Resistant',
                                     'Sensitive'),
                                   levels = c('Sensitive',
                                              'Resistant')),
                     lower_hr = ov_comb_tax_lower_hr,
                     upper_hr = ov_comb_tax_upper_hr,
                     hr = ov_comb_tax_hr)
ov_comb_tax_3 <- ggplot(ov_comb_tax_df, aes(x = test)) +
  geom_errorbar(aes(ymin = lower_hr, ymax = upper_hr),
                width = 0.15, linewidth = 1) +
  geom_point(aes(y = hr), size = 10) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1) +
  scale_y_log10(position = 'right') +
  coord_flip() +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, vjust = 1))
ggsave(file.path(figs_dir, 'ov_comb_tax_p3_hr.svg'),
       ov_comb_tax_3, device = 'svg', width = 225, height = 150, units = 'mm')


pro_comb_tax_upper_hr <- c(summary(pro_comb_tax_p3_res_cox)$conf.int[1, 4],
                     summary(pro_comb_tax_p3_sen_cox)$conf.int[1, 4])
pro_comb_tax_lower_hr <- c(summary(pro_comb_tax_p3_res_cox)$conf.int[1, 3],
                     summary(pro_comb_tax_p3_sen_cox)$conf.int[1, 3])
pro_comb_tax_hr <- c(summary(pro_comb_tax_p3_res_cox)$coefficients[1, 2],
               summary(pro_comb_tax_p3_sen_cox)$coefficients[1, 2])
pro_comb_tax_pval <- c(summary(pro_comb_tax_p3_res_cox)$coefficients[1, 5],
                 summary(pro_comb_tax_p3_sen_cox)$coefficients[1, 5])
pro_comb_tax_n <- c(summary(pro_comb_tax_p3_res_cox)$n,
              summary(pro_comb_tax_p3_sen_cox)$n)
pro_comb_tax_df <- data.frame(test = factor(c('Resistant',
                                        'Sensitive'),
                                      levels = c('Sensitive',
                                                 'Resistant')),
                        lower_hr = pro_comb_tax_lower_hr,
                        upper_hr = pro_comb_tax_upper_hr,
                        hr = pro_comb_tax_hr)
pro_comb_tax_3 <- ggplot(pro_comb_tax_df, aes(x = test)) +
  geom_errorbar(aes(ymin = lower_hr, ymax = upper_hr),
                width = 0.15, linewidth = 1) +
  geom_point(aes(y = hr), size = 10) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1) +
  scale_y_log10(position = 'right') +
  coord_flip() +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, vjust = 1))
ggsave(file.path(figs_dir, 'pro_comb_tax_p3_hr.svg'),
       pro_comb_tax_3, device = 'svg', width = 225, height = 150, units = 'mm')


bre_comb_tax_upper_hr <- c(summary(bre_comb_tax_p3_res_cox)$conf.int[1, 4],
                           summary(bre_comb_tax_p3_sen_cox)$conf.int[1, 4])
bre_comb_tax_lower_hr <- c(summary(bre_comb_tax_p3_res_cox)$conf.int[1, 3],
                           summary(bre_comb_tax_p3_sen_cox)$conf.int[1, 3])
bre_comb_tax_hr <- c(summary(bre_comb_tax_p3_res_cox)$coefficients[1, 2],
                     summary(bre_comb_tax_p3_sen_cox)$coefficients[1, 2])
bre_comb_tax_pval <- c(summary(bre_comb_tax_p3_res_cox)$coefficients[1, 5],
                       summary(bre_comb_tax_p3_sen_cox)$coefficients[1, 5])
bre_comb_tax_n <- c(summary(bre_comb_tax_p3_res_cox)$n,
                    summary(bre_comb_tax_p3_sen_cox)$n)
bre_comb_tax_df <- data.frame(test = factor(c('Resistant',
                                              'Sensitive'),
                                            levels = c('Sensitive',
                                                       'Resistant')),
                              lower_hr = bre_comb_tax_lower_hr,
                              upper_hr = bre_comb_tax_upper_hr,
                              hr = bre_comb_tax_hr)
bre_comb_tax_3 <- ggplot(bre_comb_tax_df, aes(x = test)) +
  geom_errorbar(aes(ymin = lower_hr, ymax = upper_hr),
                width = 0.15, linewidth = 1) +
  geom_point(aes(y = hr), size = 10) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1) +
  scale_y_log10(position = 'right') +
  coord_flip() +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, vjust = 1))
ggsave(file.path(figs_dir, 'bre_comb_tax_p3_hr.svg'),
       bre_comb_tax_3, device = 'svg', width = 225, height = 150, units = 'mm')


bre_comb_dox_upper_hr <- c(summary(bre_comb_dox_p3_res_cox)$conf.int[1, 4],
                           summary(bre_comb_dox_p3_sen_cox)$conf.int[1, 4])
bre_comb_dox_lower_hr <- c(summary(bre_comb_dox_p3_res_cox)$conf.int[1, 3],
                           summary(bre_comb_dox_p3_sen_cox)$conf.int[1, 3])
bre_comb_dox_hr <- c(summary(bre_comb_dox_p3_res_cox)$coefficients[1, 2],
                     summary(bre_comb_dox_p3_sen_cox)$coefficients[1, 2])
bre_comb_dox_pval <- c(summary(bre_comb_dox_p3_res_cox)$coefficients[1, 5],
                       summary(bre_comb_dox_p3_sen_cox)$coefficients[1, 5])
bre_comb_dox_n <- c(summary(bre_comb_dox_p3_res_cox)$n,
                    summary(bre_comb_dox_p3_sen_cox)$n)
bre_comb_dox_df <- data.frame(test = factor(c('Resistant',
                                              'Sensitive'),
                                            levels = c('Sensitive',
                                                       'Resistant')),
                              lower_hr = bre_comb_dox_lower_hr,
                              upper_hr = bre_comb_dox_upper_hr,
                              hr = bre_comb_dox_hr)
bre_comb_dox_3 <- ggplot(bre_comb_dox_df, aes(x = test)) +
  geom_errorbar(aes(ymin = lower_hr, ymax = upper_hr),
                width = 0.15, linewidth = 1) +
  geom_point(aes(y = hr), size = 10) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1) +
  scale_y_log10(position = 'right') +
  coord_flip() +
  theme_light() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, vjust = 1))
ggsave(file.path(figs_dir, 'bre_comb_dox_p3_hr.svg'),
       bre_comb_dox_3, device = 'svg', width = 225, height = 150, units = 'mm')

#### GGFOREST SUPPL PLOTS ####
# Make supplementary HR plots

# Phase II
# OV Plat
png(file.path(supp_figs_dir, 'OV_Comb_Plat_PhaseII_ggforest.png'), width=8, height=4, units = "in", res = 300)
ggforest_altered(tcga_ov_plat_p2_cox, data = NULL, main = "Ovarian Platinum Phase II", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

# SARC Dox
png(file.path(supp_figs_dir, 'SARC_Comb_Anth_PhaseII_ggforest.png'), width=8, height=4, units = "in", res = 300)
ggforest_altered(tcga_sarc_dox_p2_cox, data = NULL, main = "Sarcoma Anthracycline Phase II", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

# Phase III Combination
# TCGA-OV Taxane
svg(file.path(supp_figs_dir, 'OV_Comb_Tax_PhaseIII_resistant_ggforest.svg'), width=8, height=4)
ggforest_altered(ov_comb_tax_p3_res_cox, data = NULL, main = "Relapsed Ovarian Taxane Phase III - Resistant", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

svg(file.path(supp_figs_dir, 'OV_Comb_Tax_PhaseIII_sensitive_ggforest.svg'), width=8, height=4)
ggforest_altered(ov_comb_tax_p3_sen_cox, data = NULL, main = "Relapsed Ovarian Taxane Phase III - Sensitive", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()


# TCGA-OV Doxorubicin
svg(file.path(supp_figs_dir, 'OV_Comb_Anth_PhaseIII_resistant_ggforest.svg'), width=8, height=4)
ggforest_altered(ov_comb_dox_p3_res_cox, data = NULL, main = "Relapsed Ovarian Anthracycline Phase III - Resistant", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

svg(file.path(supp_figs_dir, 'OV_Comb_Anth_PhaseIII_sensitive_ggforest.svg'), width=8, height=4)
ggforest_altered(ov_comb_dox_p3_sen_cox, data = NULL, main = "Relapsed Ovarian Anthracycline Phase III - Sensitive", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()


# HMF prostate Taxane
png(file.path(supp_figs_dir, 'Pro_Comb_Tax_PhaseIII_resistant_ggforest.png'), width=8, height=4, units = "in", res = 300)
ggforest_altered(pro_comb_tax_p3_res_cox, data = NULL, main = "Metastatic Prostate Taxane Phase III - Resistant", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

png(file.path(supp_figs_dir, 'Pro_Comb_Tax_PhaseIII_sensitive_ggforest.png'), width=8, height=4, units = "in", res = 300)
ggforest_altered(pro_comb_tax_p3_sen_cox, data = NULL, main = "Metastatic Prostate Taxane Phase III - Sensitive", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()


# HMF breast Taxane
png(file.path(supp_figs_dir, 'Bre_Comb_Tax_PhaseIII_resistant_ggforest.png'), width=8, height=4, units = "in", res = 300)
ggforest_altered(bre_comb_tax_p3_res_cox, data = NULL, main = "Metastatic Breast Taxane Phase III - Resistant", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

png(file.path(supp_figs_dir, 'Bre_Comb_Tax_PhaseIII_sensitive_ggforest.png'), width=8, height=4, units = "in", res = 300)
ggforest_altered(bre_comb_tax_p3_sen_cox, data = NULL, main = "Metastatic Breast Taxane Phase III - Sensitive", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()


# HMF breast Doxo
png(file.path(supp_figs_dir, 'Bre_Comb_Dox_PhaseIII_resistant_ggforest.png'), width=8, height=4, units = "in", res = 300)
ggforest_altered(bre_comb_dox_p3_res_cox, data = NULL, main = "Metastatic Breast Anthracycline Phase III - Resistant", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

png(file.path(supp_figs_dir, 'Bre_Comb_Dox_PhaseIII_sensitive_ggforest.png'), width=8, height=4, units = "in", res = 300)
ggforest_altered(bre_comb_dox_p3_sen_cox, data = NULL, main = "Metastatic Breast Anthracycline Phase III - Sensitive", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()


# Phase III SA
# TCGA-OV Taxane
svg(file.path(supp_figs_dir, 'OV_SA_Tax_PhaseIII_resistant_ggforest.svg'), width=8, height=4)
ggforest_altered(ov_sa_tax_p3_res_cox, data = NULL, main = "Relapsed Ovarian SA Taxane Phase III - Resistant", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

svg(file.path(supp_figs_dir, 'OV_SA_Tax_PhaseIII_sensitive_ggforest.svg'), width=8, height=4)
ggforest_altered(ov_sa_tax_p3_sen_cox, data = NULL, main = "Relapsed Ovarian SA Taxane Phase III - Sensitive", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()


# TCGA-OV Doxorubicin
svg(file.path(supp_figs_dir, 'OV_SA_Anth_PhaseIII_resistant_ggforest.svg'), width=8, height=4)
ggforest_altered(ov_sa_dox_p3_res_cox, data = NULL, main = "Relapsed Ovarian SA Anthracycline Phase III - Resistant", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

svg(file.path(supp_figs_dir, 'OV_SA_Anth_PhaseIII_sensitive_ggforest.svg'), width=8, height=4)
ggforest_altered(ov_sa_dox_p3_sen_cox, data = NULL, main = "Relapsed Ovarian SA Anthracycline Phase III - Sensitive", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()


# HMF prostate Taxane
png(file.path(supp_figs_dir, 'Pro_SA_Tax_PhaseIII_resistant_ggforest.png'), width=8, height=4, units = "in", res = 300)
ggforest_altered(pro_sa_tax_p3_res_cox, data = NULL, main = "Metastatic Prostate SA Taxane Phase III - Resistant", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

png(file.path(supp_figs_dir, 'Pro_SA_Tax_PhaseIII_sensitive_ggforest.png'), width=8, height=4, units = "in", res = 300)
ggforest_altered(pro_sa_tax_p3_sen_cox, data = NULL, main = "Metastatic Prostate SA Taxane Phase III - Sensitive", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()


#### BARPLOTS CO-THERAPIES ####

# TCGA-OV Platinum
svg(file.path(supp_figs_dir, 'Cotherapies_TCGA-OV_Plat_PhaseII.svg'), width=5, height=6)
tcga_ov_plat_p2 %>%
  group_by(Prediction, drug_name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=drug_name, y=count, fill=Prediction)) +
  geom_col() + labs(subtitle="Ovarian 1st Line Platinum Phase II") +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), size = 3) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()


# TCGA-OV Taxane
svg(file.path(supp_figs_dir, 'Cotherapies_TCGA-OV_Tax_PhaseIII_resistant.svg'), width=5, height=6)
ov_comb_tax_p3_res %>%
  group_by(`Treatment arm`, drug_name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=drug_name, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Relapsed Ovarian Taxane Phase III - Resistant") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()

svg(file.path(supp_figs_dir, 'Cotherapies_TCGA-OV_Tax_PhaseIII_sensitive.svg'), width=5, height=6)
ov_comb_tax_p3_sen %>%
  group_by(`Treatment arm`, drug_name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=drug_name, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Relapsed Ovarian Taxane Phase III - Sensitive") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()

# TCGA-OV Anthracyclines
svg(file.path(supp_figs_dir, 'Cotherapies_TCGA-OV_Anth_PhaseIII_resistant.svg'), width=5, height=6)
ov_comb_dox_p3_res %>%
  group_by(`Treatment arm`, drug_name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=drug_name, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Relapsed Ovarian Anthracycline Phase III - Resistant") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()

svg(file.path(supp_figs_dir, 'Cotherapies_TCGA-OV_Anth_PhaseIII_sensitive.svg'), width=5, height=6)
ov_comb_dox_p3_sen %>%
  group_by(`Treatment arm`, drug_name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=drug_name, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Relapsed Ovarian Anthracycline Phase III - Sensitive") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()

# TCGA-SARC Anthracycline
svg(file.path(supp_figs_dir, 'Cotherapies_TCGA-SARC_Anth_PhaseII.svg'), width=5, height=6)
tcga_sarc_dox_p2 %>%
  group_by(Prediction, drug_name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=drug_name, y=count, fill=Prediction)) +
  geom_col() + labs(subtitle="Sarcoma Anthracycline Phase II") +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), size = 3) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()

# HMF-Prostate Taxane
svg(file.path(supp_figs_dir, 'Cotherapies_HMF-Pro_Tax_PhaseIII_resistant.svg'), width=5, height=6)
pro_comb_tax_p3_res %>%
  group_by(`Treatment arm`, name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=name, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Metastatic Prostate Taxane Phase III - Resistant") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()

svg(file.path(supp_figs_dir, 'Cotherapies_HMF-Pro_Tax_PhaseIII_sensitive.svg'), width=5, height=6)
pro_comb_tax_p3_sen %>%
  group_by(`Treatment arm`, name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=name, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Metastatic Prostate Taxane Phase III - Sensitive") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()


# HMF-Breast Taxane
svg(file.path(supp_figs_dir, 'Cotherapies_HMF-Breast_Tax_PhaseIII_resistant.svg'), width=5, height=6)
bre_comb_tax_p3_res %>%
  group_by(`Treatment arm`, name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=name, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Metastatic Breast Taxane Phase III - Resistant") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()

svg(file.path(supp_figs_dir, 'Cotherapies_HMF-Breast_Tax_PhaseIII_sensitive.svg'), width=5, height=6)
bre_comb_tax_p3_sen %>%
  group_by(`Treatment arm`, name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=name, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Metastatic Breast Taxane Phase III - Sensitive") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()

# HMF-Breast Dox
svg(file.path(supp_figs_dir, 'Cotherapies_HMF-Breast_Anth_PhaseIII_resistant.svg'), width=5, height=6)
bre_comb_dox_p3_res %>%
  group_by(`Treatment arm`, name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=name, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Metastatic Breast Anthracycline Phase III - Resistant") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()

svg(file.path(supp_figs_dir, 'Cotherapies_HMF-Breast_Anth_PhaseIII_sensitive.svg'), width=5, height=6)
bre_comb_dox_p3_sen %>%
  group_by(`Treatment arm`, name) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=name, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Metastatic Breast Anthracycline Phase III - Sensitive") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()


#### BARPLOTS TREATMENT LINES FOR OVARY ####
# TCGA-OV Taxane
dt=rbind(ov_comb_tax_p3_res,ov_comb_tax_p3_sen)
svg(file.path(supp_figs_dir, 'TreatLines_TCGA-OV_Tax_PhaseIII.svg'), width=5, height=5)
dt %>%
  group_by(`Treatment arm`, `Treatment line`) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=`Treatment line`, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Relapsed Ovarian Taxane Phase III") +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), size = 3) +
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5))
dev.off()

# TCGA-OV Anthracyclines
dt=rbind(ov_comb_dox_p3_res,ov_comb_dox_p3_sen)
svg(file.path(supp_figs_dir, 'TreatLines_TCGA-OV_Anth_PhaseIII.svg'), width=5, height=5)
dt %>%
  group_by(`Treatment arm`, `Treatment line`) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=`Treatment line`, y=count, fill=`Treatment arm`)) +
  geom_col() + labs(subtitle="Relapsed Ovarian Anthracycline Phase III") +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), size = 3) +
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5))
dev.off()

#### KAPLAN-MEIER PLOTS ####
# Phase 3 OV Comb Tax
colnames(ov_comb_tax_p3_res)[colnames(ov_comb_tax_p3_res)=='Treatment arm'] <- 'Treatment_arm'
ov_comb_tax_p3_res$TTF_days <- round(ov_comb_tax_p3_res$TTF * (365.25/12))
ov_comb_tax_p3_sen$TTF_days <- round(ov_comb_tax_p3_sen$TTF * (365.25/12))
ov_tax_max_x <- ceiling(max(c(ov_comb_tax_p3_res$TTF_days, ov_comb_tax_p3_sen$TTF_days)))
svg(file.path(ext_figs_dir, 'OV_Tax_Res_KM.svg'), width=5, height=5)
ggsurvplot(survfit(Surv(TTF_days, censoring) ~ Treatment_arm, data=ov_comb_tax_p3_res),
           data=ov_comb_tax_p3_res,
           pval=TRUE,
           pval.coord=c(ov_tax_max_x*0.8, 0.75),
           break.time.by=250,
           pval.size=4,
           xlim=c(0, ov_tax_max_x),
           conf.int=FALSE,
           risk.table=TRUE,
           risk.table.col='strata',
           linetype='strata',
           surv.median.line='hv',
           ggtheme=theme_bw(),
           legend.labs=c('Control', 'Experimental'),
           legend.title='',
           xlab='Time (days)',
           size=0.5,
           risk.table.y.text = FALSE,
           title='Relapsed Ovarian - Taxane Resistant')
dev.off()

colnames(ov_comb_tax_p3_sen)[colnames(ov_comb_tax_p3_sen)=='Treatment arm'] <- 'Treatment_arm'
svg(file.path(ext_figs_dir, 'OV_Tax_Sen_KM.svg'), width=5, height=5)
ggsurvplot(survfit(Surv(TTF_days, censoring) ~ Treatment_arm, data=ov_comb_tax_p3_sen),
           data=ov_comb_tax_p3_sen,
           pval=TRUE,
           pval.coord=c(ov_tax_max_x*0.8, 0.75),
           pval.size=4,
           break.time.by=250,
           xlim=c(0, ov_tax_max_x),
           conf.int=FALSE,
           risk.table=TRUE,
           risk.table.col='strata',
           linetype='strata',
           surv.median.line='hv',
           ggtheme=theme_bw(),
           legend.labs=c('Control', 'Experimental'),
           legend.title='',
           xlab='Time (days)',
           size=0.5,
           risk.table.y.text = FALSE,
           title='Relapsed Ovarian - Taxane Sensitive')
dev.off()


# Phase 3 OV Comb Dox
colnames(ov_comb_dox_p3_res)[colnames(ov_comb_dox_p3_res)=='Treatment arm'] <- 'Treatment_arm'
ov_comb_dox_p3_res$TTF_days <- round(ov_comb_dox_p3_res$TTF * (365.25/12))
ov_comb_dox_p3_sen$TTF_days <- round(ov_comb_dox_p3_sen$TTF * (365.25/12))
ov_dox_max_x <- ceiling(max(c(ov_comb_dox_p3_res$TTF_days, ov_comb_dox_p3_sen$TTF_days)))
svg(file.path(ext_figs_dir, 'OV_Dox_Res_KM.svg'), width=5, height=5)
ggsurvplot(survfit(Surv(TTF_days, censoring) ~ Treatment_arm, data=ov_comb_dox_p3_res),
           data=ov_comb_dox_p3_res,
           pval=TRUE,
           pval.coord=c(ov_dox_max_x*0.8, 0.75),
           pval.size=4,
           break.time.by=250,
           xlim=c(0, ov_dox_max_x),
           conf.int=FALSE,
           risk.table=TRUE,
           risk.table.col='strata',
           linetype='strata',
           surv.median.line='hv',
           ggtheme=theme_bw(),
           legend.labs=c('Control', 'Experimental'),
           legend.title='',
           xlab='Time (days)',
           size=0.5,
           risk.table.y.text = FALSE,
           title='Relapsed Ovarian - Anthracycline Resistant')
dev.off()

colnames(ov_comb_dox_p3_sen)[colnames(ov_comb_dox_p3_sen)=='Treatment arm'] <- 'Treatment_arm'
svg(file.path(ext_figs_dir, 'OV_Dox_Sen_KM.svg'), width=5, height=5)
ggsurvplot(survfit(Surv(TTF_days, censoring) ~ Treatment_arm, data=ov_comb_dox_p3_sen),
           data=ov_comb_dox_p3_sen,
           pval=TRUE,
           pval.coord=c(ov_dox_max_x*0.8, 0.75),
           pval.size=4,
           break.time.by=250,
           xlim=c(0, ov_dox_max_x),
           conf.int=FALSE,
           risk.table=TRUE,
           risk.table.col='strata',
           linetype='strata',
           surv.median.line='hv',
           ggtheme=theme_bw(),
           legend.labs=c('Control', 'Experimental'),
           legend.title='',
           xlab='Time (days)',
           size=0.5,
           risk.table.y.text = FALSE,
           title='Relapsed Ovarian - Anthracycline Sensitive')
dev.off()

# Phase 3 Prostate Comb Tax
colnames(pro_comb_tax_p3_res)[colnames(pro_comb_tax_p3_res)=='Treatment arm'] <- 'Treatment_arm'
pro_tax_max_x <- ceiling(max(c(pro_comb_tax_p3_res$TTF, pro_comb_tax_p3_sen$TTF)))
svg(file.path(ext_figs_dir, 'Pro_Tax_Res_KM.svg'), width=5, height=5)
ggsurvplot(survfit(Surv(TTF, censoring) ~ Treatment_arm, data=pro_comb_tax_p3_res),
           data=pro_comb_tax_p3_res,
           pval=TRUE,
           pval.coord=c(pro_tax_max_x*0.8, 0.75),
           pval.size=4,
           break.time.by=250,
           xlim=c(0, pro_tax_max_x),
           conf.int=FALSE,
           risk.table=TRUE,
           risk.table.col='strata',
           linetype='strata',
           surv.median.line='hv',
           ggtheme=theme_bw(),
           legend.labs=c('Control', 'Experimental'),
           legend.title='',
           xlab='Time (days)',
           size=0.5,
           risk.table.y.text = FALSE,
           title='Metastatic Prostate - Taxane Resistant')
dev.off()

colnames(pro_comb_tax_p3_sen)[colnames(pro_comb_tax_p3_sen)=='Treatment arm'] <- 'Treatment_arm'
svg(file.path(ext_figs_dir, 'Pro_Tax_Sen_KM.svg'), width=5, height=5)
ggsurvplot(survfit(Surv(TTF, censoring) ~ Treatment_arm, data=pro_comb_tax_p3_sen),
           data=pro_comb_tax_p3_sen,
           pval=TRUE,
           pval.coord=c(pro_tax_max_x*0.8, 0.75),
           pval.size=4,
           break.time.by=250,
           xlim=c(0, pro_tax_max_x),
           conf.int=FALSE,
           risk.table=TRUE,
           risk.table.col='strata',
           linetype='strata',
           surv.median.line='hv',
           ggtheme=theme_bw(),
           legend.labs=c('Control', 'Experimental'),
           legend.title='',
           xlab='Time (days)',
           size=0.5,
           risk.table.y.text = FALSE,
           title='Metastatic Prostate - Taxane Sensitive')
dev.off()

# Phase 3 Breast Comb Tax
colnames(bre_comb_tax_p3_res)[colnames(bre_comb_tax_p3_res)=='Treatment arm'] <- 'Treatment_arm'
bre_tax_max_x <- ceiling(max(c(bre_comb_tax_p3_res$TTF, bre_comb_tax_p3_sen$TTF)))
svg(file.path(ext_figs_dir, 'Bre_Tax_Res_KM.svg'), width=5, height=5)
ggsurvplot(survfit(Surv(TTF, censoring) ~ Treatment_arm, data=bre_comb_tax_p3_res),
           data=bre_comb_tax_p3_res,
           pval=TRUE,
           pval.coord=c(bre_tax_max_x*0.8, 0.75),
           pval.size=4,
           break.time.by=250,
           xlim=c(0, bre_tax_max_x),
           conf.int=FALSE,
           risk.table=TRUE,
           risk.table.col='strata',
           linetype='strata',
           surv.median.line='hv',
           ggtheme=theme_bw(),
           legend.labs=c('Control', 'Experimental'),
           legend.title='',
           xlab='Time (days)',
           size=0.5,
           risk.table.y.text = FALSE,
           title='Metastatic Breast - Taxane Resistant')
dev.off()

colnames(bre_comb_tax_p3_sen)[colnames(bre_comb_tax_p3_sen)=='Treatment arm'] <- 'Treatment_arm'
svg(file.path(ext_figs_dir, 'Bre_Tax_Sen_KM.svg'), width=5, height=5)
ggsurvplot(survfit(Surv(TTF, censoring) ~ Treatment_arm, data=bre_comb_tax_p3_sen),
           data=bre_comb_tax_p3_sen,
           pval=TRUE,
           pval.coord=c(bre_tax_max_x*0.8, 0.75),
           pval.size=4,
           break.time.by=250,
           xlim=c(0, bre_tax_max_x),
           conf.int=FALSE,
           risk.table=TRUE,
           risk.table.col='strata',
           linetype='strata',
           surv.median.line='hv',
           ggtheme=theme_bw(),
           legend.labs=c('Control', 'Experimental'),
           legend.title='',
           xlab='Time (days)',
           size=0.5,
           risk.table.y.text = FALSE,
           title='Metastatic Breast - Taxane Sensitive')
dev.off()

# Phase 3 Breast Comb Dox
colnames(bre_comb_dox_p3_res)[colnames(bre_comb_dox_p3_res)=='Treatment arm'] <- 'Treatment_arm'
bre_dox_max_x <- ceiling(max(bre_comb_dox_p3_res$TTF))
svg(file.path(ext_figs_dir, 'Bre_Dox_Res_KM.svg'), width=5, height=5)
ggsurvplot(survfit(Surv(TTF, censoring) ~ Treatment_arm, data=bre_comb_dox_p3_res),
           data=bre_comb_dox_p3_res,
           pval=TRUE,
           pval.coord=c(bre_dox_max_x*0.8, 0.75),
           pval.size=4,
           break.time.by=250,
           xlim=c(0, bre_dox_max_x),
           conf.int=FALSE,
           risk.table=TRUE,
           risk.table.col='strata',
           linetype='strata',
           surv.median.line='hv',
           ggtheme=theme_bw(),
           legend.labs=c('Control', 'Experimental'),
           legend.title='',
           xlab='Time (days)',
           size=0.5,
           risk.table.y.text = FALSE,
           title='Metastatic Breast - Anthracycline Resistant')
dev.off()
