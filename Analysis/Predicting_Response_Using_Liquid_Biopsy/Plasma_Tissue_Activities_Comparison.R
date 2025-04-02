################################################################################
## This script is used for confirming that our predictor can be applied in plasma
## samples from 41 HGSOC samples
##
## We quantify CIN signatures in tissue and liquid biopsies from 41 ovarian tumours 
## included in the OV04 clinical study. Copy number profiles from both tissue 
## and ctDNA were derived from sWGS data.
################################################################################

# Clean environment
freshr::freshr()

# Libraries
library(this.path)
library(data.table)
library(ggplot2)
library(dplyr)

# Paths
liquid_dir <- dirname(this.path())
base_dir <- dirname(dirname(liquid_dir))
helper_dir <- file.path(dirname(liquid_dir), 'Helper_Scripts')
input_dir <- file.path(base_dir, "Input_Data")
figs_dir <- file.path(base_dir, "Figures")
extended_figs_dir <- file.path(figs_dir, 'Extended')
dir.create(extended_figs_dir, showWarnings = FALSE, recursive = TRUE)

# Files
source(file.path(helper_dir, 'Helper_functions.R'))
metadata <- readRDS(paste0(input_dir, "/patient_data.rds"))
plasma_sigs <- readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_OV_plasma.rds'))
tissue_sigs <- readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_OV_tissue_old_50kb.rds'))
tcga_exposures <- readRDS(file.path(input_dir,"4_Activities_CompendiumCINSigs_THRESH95_TCGA.rds"))
tcga_meta <- readRDS(file.path(input_dir,"Metadata_TCGA_ASCAT_penalty70.rds"))
thresh <- readRDS(file.path(input_dir, 'Jul2024_TCGA_Signature_Thresholds.rds'))
brca_status <- readRDS(file.path(input_dir,"TCGA_BRCA_MutStatus.rds"))
platinum_data <- readRDS(file.path(input_dir, 'OV_tissue_platinum_predictions.rds'))
paclitaxel_data <- readRDS(file.path(input_dir, 'OV_tissue_paclitaxel_predictions.rds'))
doxorubicin_data <- readRDS(file.path(input_dir, 'OV_tissue_doxorubicin_predictions.rds'))
patient_info <- read_csv(file.path(input_dir, 'patient_info.csv'))
ids <- fread(paste0(input_dir, "/OV04_JBLAB_id_mapping.csv"))
qc <- readRDS(paste0(input_dir,"/final_patients.RDS"))
ids <- ids[ids$OV04_ID%in%qc$study_subject_id,]
ids <- ids[!ids$samplename%in%c("JBLAB-17049","JBLAB-17058"),] #we use the sample from the previous cohort


## Apply signature-threshold in exposures to trust in low activities
tissue_sigs = do.call(cbind,lapply(1:ncol(tissue_sigs), function(thisSig){
  sig=colnames(tissue_sigs)[thisSig]
  acts=tissue_sigs[,thisSig]
  thrs=thresh[names(thresh)==sig]
  acts[acts<thrs]=0
  return(acts)
}))
colnames(tissue_sigs)=names(thresh)
tissue_sigs=as.data.frame(t(apply(tissue_sigs, 1, function(x) x/sum(x))))
rownames(tissue_sigs) <- unlist(lapply(strsplit(gsub('downsampled_', '', rownames(tissue_sigs), fixed = T), '-Tseq'), FUN = function(x) x[1]))

plasma_sigs = do.call(cbind,lapply(1:ncol(plasma_sigs), function(thisSig){
  sig=colnames(plasma_sigs)[thisSig]
  acts=plasma_sigs[,thisSig]
  thrs=thresh[names(thresh)==sig]
  acts[acts<thrs]=0
  return(acts)
}))
colnames(plasma_sigs)=names(thresh)
plasma_sigs=as.data.frame(t(apply(plasma_sigs, 1, function(x) x/sum(x))))
rownames(plasma_sigs) <- unlist(lapply(strsplit(gsub('downsampled_', '', rownames(plasma_sigs), fixed = T), '-Tseq'), FUN = function(x) x[1]))


################################################################################
########### COMPARE TISSUE & PLASMA EXPOSURES ##################################
signatures=c(paste0("CX",1:17))
metadata <- metadata[!is.na(metadata$plasma_id), ]
metadata <- metadata[metadata$study_subject_id != 78, ]

plasma_exps=cbind(sample=row.names(plasma_sigs),plasma_sigs)
plasma_exps=melt(plasma_exps)
plasma_exps$study_subject_id=metadata$study_subject_id[match(plasma_exps$sample, metadata$plasma_id)]
plasma_exps$condition="plasma"
ids=metadata$tissue_primary_sample[metadata$plasma_id%in%plasma_exps$sample]
tissue_exps=cbind(sample=row.names(tissue_sigs),tissue_sigs)
tissue_exps=melt(tissue_exps)
tissue_exps=tissue_exps[tissue_exps$sample%in%ids,]
tissue_exps$study_subject_id=metadata$study_subject_id[match(tissue_exps$sample, metadata$tissue_primary_sample)]
tissue_exps$condition="tissue"

pair_exps=rbind(plasma_exps,tissue_exps)
pair_exps=pair_exps[!is.na(pair_exps$study_subject_id),]
pair_exps$variable <- factor(pair_exps$variable, levels = paste0('CX', 1:17))
pair_exps$study_subject_id <- factor(pair_exps$study_subject_id,
                                     levels = as.character(sort(as.numeric(unique(pair_exps$study_subject_id)))))

custom_palette <- data.frame(colour = c("#4cbcd1", "#db9540",
                                         "#b16f45", "#4d75b2",
                                         "#73792f", "#d35138",
                                         "#ad5fd0", "#6cb543",
                                         "#ca48a0", "#3a814f",
                                         "#6861c4", "#b9b04b",
                                         "#8b99e4", "#ce82bc",
                                         "#d95070", "#a34e6a", "#5bbf8a"),
                             signature = paste0("CX", rep(1:17)))


# Plot pancancer signatures in platinum samples
sig_comparison <- ggplot(pair_exps, aes(fill = variable, y = value, x = condition)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = custom_palette$colour, name = "Signatures") +
  facet_wrap(~study_subject_id, nrow = 1) +
  ylab(c("Pancancer signature activity")) +
  xlab("") +
  ggtitle("")  +
  theme_bw() +
  theme(legend.key.height =  unit(0.3, 'cm'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())
ggsave(file.path(figs_dir, 'OV_plasma_tissue_signature_concordance.svg'),
       sig_comparison, width = 150/25.4, height = 80/25.4)


## Cosine similarities between tissue and plasma exposures
ids=unique(pair_exps$study_subject_id)
cosine <- c()
for(thisSamp in ids){
  tiss <- tissue_sigs[row.names(tissue_sigs)%in%pair_exps$sample[pair_exps$study_subject_id==thisSamp], ]
  plas <- plasma_sigs[row.names(plasma_sigs)%in%pair_exps$sample[pair_exps$study_subject_id==thisSamp], ]
  cos <- (sum(tiss * plas, na.rm = TRUE)) / 
    (sqrt(sum(tiss^2, na.rm = TRUE)) * sqrt(sum(plas^2, na.rm = TRUE)))
  cosine <- rbind(cosine, c(thisSamp, cos))
}
cosine <- as.data.frame(cosine)
colnames(cosine) <- c("id", "cosine")
cosine$cosine <- as.numeric(cosine$cosine)
mean(cosine$cosine) #0.9000742
sd(cosine$cosine) #0.07408859


################################################################################
########### MAP TISSUE & PLASMA PREDICTIONS OF PATIENTS ########################

### Platinum predictions
## Scaling model
lModel = readRDS(paste0(input_dir,"/OV04_plat_scaling.RDS"))

## Plasma
sc_platinum_plasma_exposures = sweep(sweep(plasma_sigs[,c("CX3", "CX2")], 2, lModel$mean, FUN = '-'), 2, lModel$scale, FUN = "/")
colnames(sc_platinum_plasma_exposures)=c("sCX3","sCX2")
sc_platinum_plasma_exposures=as.data.frame(sc_platinum_plasma_exposures)
sc_platinum_plasma_exposures$plasma_id <- row.names(sc_platinum_plasma_exposures)
sc_platinum_plasma_exposures$study_subject_id=metadata$study_subject_id[match(sc_platinum_plasma_exposures$plasma_id, metadata$plasma_id)]
sc_platinum_plasma_exposures$platinum_prediction_plasma <- ifelse(sc_platinum_plasma_exposures$sCX3 > sc_platinum_plasma_exposures$sCX2,
                                                                  'Plat sensitive',
                                                                  'Plat resistant')

## Compare tissue and plasma predictions
platinum_data=platinum_data[,c(1,26)]
platinum_data=rbind(platinum_data,c(71,"Sensitive")) #removed from cox because stage I
platinum_data=rbind(platinum_data,c(648,"Resistant")) #removed from cox because enrollment in clinical trial
platinum_data$study_subject_id<-as.numeric(platinum_data$study_subject_id)
platinum_data=left_join(platinum_data,sc_platinum_plasma_exposures,by="study_subject_id")
platinum_data_filt=platinum_data[!is.na(platinum_data$platinum_prediction_plasma),]
table(platinum_data_filt$prediction,platinum_data_filt$platinum_prediction_plasma)


### Paclitaxel predictions
## Scaling model
mTCGA = as.matrix(tcga_exposures[ , c("CX3", "CX5")])
lModel = list(mean = attributes(scale(mTCGA))$`scaled:center`, 
              scale = attributes(scale(mTCGA))$`scaled:scale`)

sc_paclitaxel_plasma_exposures = sweep(sweep(plasma_sigs[,c("CX3", "CX5")], 2, lModel$mean, FUN = '-'), 2, lModel$scale, FUN = "/")
colnames(sc_paclitaxel_plasma_exposures)=c("sCX3","sCX5")
sc_paclitaxel_plasma_exposures=as.data.frame(sc_paclitaxel_plasma_exposures)
sc_paclitaxel_plasma_exposures$plasma_id <- row.names(sc_paclitaxel_plasma_exposures)
sc_paclitaxel_plasma_exposures$study_subject_id=metadata$study_subject_id[match(sc_paclitaxel_plasma_exposures$plasma_id, metadata$plasma_id)]
sc_paclitaxel_plasma_exposures$taxane_prediction_plasma <- ifelse(sc_paclitaxel_plasma_exposures$sCX5<0,
                                                                  'Tax resistant',
                                                                  'Tax sensitive')
## Compare tissue and plasma predictions
paclitaxel_data=left_join(paclitaxel_data,sc_paclitaxel_plasma_exposures,by="study_subject_id")
paclitaxel_data_filt=paclitaxel_data[!is.na(paclitaxel_data$taxane_prediction_plasma),]
table(paclitaxel_data_filt$prediction, paclitaxel_data_filt$taxane_prediction_plasma)


### Doxorubicin predictions
## Plasma predictions
dox_plasma_exposures <- as.data.frame(plasma_sigs[,c("CX8","CX9","CX13")])
dox_plasma_exposures$plasma_id <- row.names(dox_plasma_exposures)
dox_plasma_exposures$study_subject_id=metadata$study_subject_id[match(dox_plasma_exposures$plasma_id, metadata$plasma_id)]
dox_plasma_exposures$doxorubicin_prediction_plasma <- ifelse(dox_plasma_exposures$CX8 > 0.01 | dox_plasma_exposures$CX9 > 0.009 | dox_plasma_exposures$CX13 > 0.009, 
                                                             'Dox resistant', 'Dox sensitive')
## Compare tissue and plasma predictions
doxorubicin_data=left_join(doxorubicin_data,dox_plasma_exposures,by="study_subject_id")
doxorubicin_data_filt=doxorubicin_data[!is.na(doxorubicin_data$doxorubicin_prediction_plasma),]
table(doxorubicin_data_filt$prediction, doxorubicin_data_filt$doxorubicin_prediction_plasma)


## TOTAL NUMBERS => 83.33%
# Sensitive tissue & sensitive plasma = 8
# Sensitive tissue & resistant plasma = 1
# Resistant tissue & sensitive plasma = 2
# Resistant tissue & resistant plasma = 7
