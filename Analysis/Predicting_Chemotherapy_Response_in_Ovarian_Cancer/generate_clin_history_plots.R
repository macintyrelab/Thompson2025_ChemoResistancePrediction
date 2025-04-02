# This script generates the clinical history plots for the OV04 patients, Supplementary Figs 1-41

# Clean environment
freshr::freshr()

# Libraries
library(readr)
library(dplyr)
library(ggplot2)
library(this.path)

# Paths
ovarian_dir <- dirname(this.path())
base_dir <- dirname(dirname(ovarian_dir))
input_dir <- file.path(base_dir,"Input_Data")
figs_dir <- file.path(base_dir, 'Figures')
supp_figs_dir <- file.path(figs_dir, 'Supplementary')

# Files
clinical_history <- read_csv(file.path(input_dir, 'OV04_clinical_history_anonymised.csv'))
platinum_data_from_survival <- readRDS(file.path(input_dir, 'OV_tissue_platinum_predictions.rds'))
taxane_data_from_survival <- readRDS(file.path(input_dir, 'OV_tissue_paclitaxel_predictions.rds'))
doxorubicin_data_from_survival <- readRDS(file.path(input_dir, 'OV_tissue_doxorubicin_predictions.rds'))
platinum_data <- readRDS(file.path(input_dir, 'OV04_plat_PFS_anonymised.rds'))
taxane_data <- readRDS(file.path(input_dir, 'OV04_tax_PFS_anonymised.rds'))
doxorubicin_data <- readRDS(file.path(input_dir, 'OV04_doxo_PFS_anonymised.rds'))


# Split clinical history into CA125 readings and treatment administrations
clinical_history_readings <- clinical_history %>%
  mutate(is_ca125 = !is.na(as.numeric(event))) %>%
  filter(is_ca125) %>% select(!is_ca125)
clinical_history_treatments <- clinical_history %>%
  mutate(is_ca125 = !is.na(as.numeric(event))) %>%
  filter(!is_ca125) %>% select(!is_ca125)
clinical_history_treatments$event <- toupper(clinical_history_treatments$event)

# Treatments
treatments <- c("CARBOPLATIN","PACLITAXEL", "CISPLATIN", "DOXORUBICIN","GEMCITABINE" , "MEDROXYPROGESTERONE", 
                "RUCAPARIB",  "EPIRUBICIN", "CAPECITABINE","MYOCET", "BEVACIZUMAB", "ONX-0801", "OLAPARIB", 
                "DENOSUMAB", "FULVESTRANT", "NIRAPARIB",  "PEMBROLIZUMAB", "PLATINUM", "CLINICAL TRIAL",
                "LETROZOLE", "TRABECTEDIN","TAMOXIFEN",  "TOPOTECAN",  "ANASTROZOLE","PAZOPANIB")
treatment_colours <- c("#9d823e", "#b70cb2", "#acb839", "#573687", "#6eab42",
                       "#be72c9", "#5fc875", "#cb4883", "#47bb8a", "#c43d56",
                       "#33d4d1", "#ce5040", "#628bd5", "#d1972c", "#84295f",
                       "#9dbb69", "#d573b3", "#427228", "#ce5040", "#be5967",
                       "#9b9937", "#842f1a", "#d5a550", "#bd6825", "#cf7c53")
treatment_colours <- setNames(treatment_colours, treatments)

treatment_linetypes <- c('solid', 'longdash', 'dotted')
treatment_linetypes <- setNames(treatment_linetypes, c('Platinum', 'Taxane', 'Anthracycline'))

treatment_shapes <- c('square', 'triangle', 'circle')
treatment_shapes <- setNames(treatment_shapes, c('Platinum', 'Taxane', 'Anthracycline'))

for(patient in unique(clinical_history$study_subject_id)){
  clinical_history_readings_patient <- clinical_history_readings %>%
    filter(study_subject_id == patient)
  clinical_history_readings_patient$event <- as.numeric(clinical_history_readings_patient$event)
  
  clinical_history_treatments_patient <- clinical_history_treatments %>%
    filter(study_subject_id == patient) %>%
    mutate(multiple_treatment = duplicated(days_since_first_treatment)) %>%
    mutate(days_since_first_treatment = ifelse(multiple_treatment,
                                               days_since_first_treatment+5, # Prevent overplotting of same-day treatments
                                               days_since_first_treatment)) 
  # treat lines for background
  treatment_lines <- clinical_history %>%
    filter(study_subject_id==patient) %>%
    group_by(treat_line) %>%
    summarise(start=min(days_since_first_treatment), end=max(treat_line_duration)) %>%
    arrange(start) %>%
    filter(treat_line!='Followup')
  
  # PFS - take the ca125 reading that is closer to the progression time point
  plat_pfs <- platinum_data %>% filter(study_subject_id==patient)
  if(nrow(plat_pfs)==1){
    plat_pfs$event <- clinical_history_readings_patient$event[which.min(abs(clinical_history_readings_patient$days_since_first_treatment - plat_pfs$PFS_from_initial_treatment))]
  }
  tax_pfs <- taxane_data %>% filter(study_subject_id==patient) 
  tax_pfs <- tax_pfs %>% ungroup()
  if(nrow(tax_pfs)>1){
    tax_pfs <- tax_pfs[1,]
  }
  if(nrow(tax_pfs)==1){
    tax_pfs$event <- clinical_history_readings_patient$event[which.min(abs(clinical_history_readings_patient$days_since_first_treatment - tax_pfs$PFS_from_initial_treatment))]
  }
  dox_pfs <- doxorubicin_data %>% filter(study_subject_id==patient) 
  if(nrow(dox_pfs)==1){
    dox_pfs$event <- clinical_history_readings_patient$event[which.min(abs(clinical_history_readings_patient$days_since_first_treatment - dox_pfs$PFS_from_initial_treatment))]
  }
  
  pfsdat <- rbind(plat_pfs, tax_pfs, dox_pfs)
  pfsdat$therapy <- factor(pfsdat$therapy, levels=c('Platinum', 'Taxane', 'Anthracycline'))
  colnames(pfsdat)[colnames(pfsdat)=='therapy'] <-'Therapy'
  treatment_lines$treat_line <- factor(treatment_lines$treat_line, levels=unique(clinical_history$treat_line))  
  
  ymax <- max(clinical_history_readings_patient$event)
  y_break_seq_max <- ceiling(log2(ymax))
  
  clinical_history_plot_strict <-  ggplot() +
    geom_rect(data = treatment_lines,
              aes(xmin=start, xmax=end, ymin=0, ymax=ymax, fill=as.character(treat_line)), alpha=0.2) +
    geom_vline(data = clinical_history_treatments_patient, 
               aes(xintercept = days_since_first_treatment, color = event),
               linewidth = 0.8) +
    geom_line(data = clinical_history_readings_patient,
              aes(x = days_since_first_treatment, y = event)) +
    geom_point(data = clinical_history_readings_patient,
               aes(x = days_since_first_treatment, y = event)) +
    # progression point for primary treatment, doxo and tax
    geom_point(data = pfsdat, 
               aes(x = PFS_from_initial_treatment, y = event,shape=Therapy),size=4, color='red') +
    geom_hline(yintercept = 35,size=0.3, color='blue') +
    geom_hline(data=pfsdat, aes(yintercept=threshold.to.use, linetype=Therapy)) +
    scale_y_continuous(trans = 'log2', breaks=c(0, 2^(1:y_break_seq_max))) +
    scale_color_manual(values = treatment_colours) + 
    scale_linetype_manual(values=treatment_linetypes) +
    scale_shape_manual(values=treatment_shapes) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          title = element_text(size = 18),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12)) +
    labs(x = 'Days since start of clinical history',
         y = 'CA125 units/ml',
         title = paste('Patient', patient), color = 'Treatment',fill = 'Treatment Line') +
    guides(color=guide_legend(order=2),
           shape=guide_legend(order=1),
           fill=guide_legend(order=3),
           linetype=guide_legend(order=1))
  ggsave(paste0(supp_figs_dir, '/strict_OV04_clin/', patient, '_clin_history_strict.png'),
         clinical_history_plot_strict, device = 'png', width = 40, height = 20, units = 'cm')
}




#### SWIMMER PLOTS ####
library(swimplot)

# Functions
source(file.path(base_dir, 'Analysis/Helper_Scripts/Helper_functions.R'))

# Prepare files
clinical_history_readings$event = "CA125 reading"
clinical_history = rbind(clinical_history_readings,clinical_history_treatments)

# OV04 platinum data
platinum_data <- platinum_data[platinum_data$study_subject_id%in%platinum_data_from_survival$study_subject_id,]
platinum_data$prediction = platinum_data_from_survival$prediction[match(platinum_data_from_survival$study_subject_id,platinum_data_from_survival$study_subject_id)]
platinum_data$censoring = platinum_data_from_survival$censoring[match(platinum_data_from_survival$study_subject_id,platinum_data_from_survival$study_subject_id)]
plat_data_swimmer <- prepare_swimmer(chemo_history=platinum_data,clinical_history=clinical_history,chemo="Platinum")

treatment=plat_data_swimmer$treatment
treatment=as.data.frame(reshape2::dcast(treatment,study_subject_id~event, value.var = "days_since_first_treatment"))
treatment$event="Treatment_length"
plat_swimmer_plot <- swimmer_plot(df=as.data.frame(plat_data_swimmer$main),id='study_subject_id',end='PFS_from_initial_treatment',name_fill='prediction',width=.85,id_order='prediction') +
  swimmer_lines(df_lines=treatment,id='study_subject_id',start="Treatment_start",end="Treatment_end",name_col='event',size=1) + 
  swimmer_points(df_points=as.data.frame(plat_data_swimmer$next_treatment),id='study_subject_id',time='days_since_first_treatment',name_shape='event',size=1.5,fill='white',col=alpha('black',0.5)) + # add points 
  theme_bw() + labs(subtitle="OV04 platinum",
                    x = "Patient", y = "PFS (days)") + 
  scale_fill_manual(name="Prediction",values=c("Resistant" = "skyblue1", "Sensitive" = "lightpink1")) +
  scale_shape_manual(name="Event",values=c(17,8)) +
  scale_color_manual(name="Treatment length",values="red")+
  theme(axis.text.y = element_text(size=6))
ggsave(file.path(supp_figs_dir, 'OV04_platinum_swimmer.svg'), plot = plat_swimmer_plot, width=25, height=20, units='cm', dpi=300, device='svg')

# Differences in treatment time between sensitive and resistant patients
treatment$length=treatment$Treatment_end-treatment$Treatment_start
treat.res=treatment$length[treatment$study_subject_id%in%platinum_data$study_subject_id[platinum_data$prediction=="Resistant" & platinum_data$censoring!=0]]
treat.sen=treatment$length[!treatment$study_subject_id%in%platinum_data$study_subject_id[platinum_data$prediction=="Resistant" & platinum_data$censoring!=0]]
wilcox.test(treat.res,treat.sen) #p-val=0.4388


# OV04 paclitaxel data 
taxane_data <- taxane_data[taxane_data$study_subject_id%in%taxane_data_from_survival$study_subject_id,]
taxane_data$prediction = taxane_data_from_survival$prediction[match(taxane_data$study_subject_id,taxane_data_from_survival$study_subject_id)]
taxane_data$treatment_line = taxane_data_from_survival$treatment_line[match(taxane_data$study_subject_id,taxane_data_from_survival$study_subject_id)]
taxane_data$censoring = taxane_data_from_survival$censoring[match(taxane_data$study_subject_id,taxane_data_from_survival$study_subject_id)]

tax_data_swimmer <- prepare_swimmer(chemo_history=taxane_data,clinical_history=clinical_history,chemo="Taxane")

treatment=tax_data_swimmer$treatment
treatment=as.data.frame(reshape2::dcast(treatment,study_subject_id~event, value.var = "days_since_first_treatment"))
treatment$event="Treatment_length"
tax_swimmer_plot <- swimmer_plot(df=as.data.frame(tax_data_swimmer$main),id='study_subject_id',end='PFS',name_fill='prediction',width=.85,id_order='prediction') +
  swimmer_lines(df_lines=treatment,id='study_subject_id',start="Treatment_start",end="Treatment_end",name_col='event',size=1) + 
  swimmer_points(df_points=as.data.frame(tax_data_swimmer$next_treatment),id='study_subject_id',time='days_since_first_treatment',name_shape='event',size=1.5,fill='white',col=alpha('black',0.5)) + # add points 
  theme_bw() + labs(subtitle="OV04 paclitaxel",
                    x = "Patient", y = "PFS (days)") + 
  scale_fill_manual(name="Prediction",values=c("Resistant" = "skyblue1", "Sensitive" = "lightpink1")) +
  scale_shape_manual(name="Event",values=c(17,8)) +
  scale_color_manual(name="Treatment length",values="red")+
  theme(axis.text.y = element_text(size=6))
ggsave(file.path(supp_figs_dir, 'OV04_paclitaxel_swimmer.svg'), plot = tax_swimmer_plot, width=25, height=15, units='cm', dpi=300, device='svg')

# Differences in treatment time between sensitive and resistant patients
treatment$length=treatment$Treatment_end-treatment$Treatment_start
treat.res=treatment$length[treatment$study_subject_id%in%taxane_data$study_subject_id[taxane_data$prediction=="Resistant" & taxane_data$censoring!=0]]
treat.sen=treatment$length[!treatment$study_subject_id%in%taxane_data$study_subject_id[taxane_data$prediction=="Resistant"& taxane_data$censoring!=0]]
wilcox.test(treat.res,treat.sen) #p-val=0.5736


# OV04 doxorubicin data
doxorubicin_data <- doxorubicin_data[doxorubicin_data$study_subject_id%in%doxorubicin_data_from_survival$study_subject_id,]
doxorubicin_data$prediction = doxorubicin_data_from_survival$prediction[match(doxorubicin_data$study_subject_id,doxorubicin_data_from_survival$study_subject_id)]
doxorubicin_data$treatment_line = doxorubicin_data_from_survival$treatment_line[match(doxorubicin_data$study_subject_id,doxorubicin_data_from_survival$study_subject_id)]
doxorubicin_data$censoring = doxorubicin_data_from_survival$censoring[match(doxorubicin_data$study_subject_id,doxorubicin_data_from_survival$study_subject_id)]

dox_data_swimmer <- prepare_swimmer(chemo_history=doxorubicin_data,clinical_history=clinical_history,chemo="Doxorubicin")

treatment=dox_data_swimmer$treatment
treatment=as.data.frame(reshape2::dcast(treatment,study_subject_id~event, value.var = "days_since_first_treatment"))
treatment$event="Treatment_length"
dox_swimmer_plot <- swimmer_plot(df=as.data.frame(dox_data_swimmer$main),id='study_subject_id',end='PFS',name_fill='prediction',width=.85,id_order='prediction') +
  swimmer_lines(df_lines=treatment,id='study_subject_id',start="Treatment_start",end="Treatment_end",name_col='event',size=1) + 
  swimmer_points(df_points=as.data.frame(dox_data_swimmer$next_treatment),id='study_subject_id',time='days_since_first_treatment',name_shape='event',size=1.5,fill='white',col=alpha('black',0.5)) + # add points 
  theme_bw() + labs(subtitle="OV04 doxorubicin",
                    x = "Patient", y = "PFS (days)") + 
  scale_fill_manual(name="Prediction",values=c("Resistant" = "skyblue1", "Sensitive" = "lightpink1")) +
  scale_shape_manual(name="Event",values=c(17,8)) +
  scale_color_manual(name="Treatment length",values="red")+
  theme(axis.text.y = element_text(size=6))
ggsave(file.path(supp_figs_dir, 'OV04_doxorubicin_swimmer.svg'), plot = dox_swimmer_plot, width=25, height=15, units='cm', dpi=300, device='svg')

# Differences in treatment time between sensitive and resistant patients
treatment$length=treatment$Treatment_end-treatment$Treatment_start
treat.res=treatment$length[treatment$study_subject_id%in%doxorubicin_data$study_subject_id[doxorubicin_data$prediction=="Resistant"]]
treat.sen=treatment$length[!treatment$study_subject_id%in%doxorubicin_data$study_subject_id[doxorubicin_data$prediction=="Resistant"]]
wilcox.test(treat.res,treat.sen) #p-value=0.2078



