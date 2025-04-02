################################################################################
## This script is used for confirming that our predictor can be applied in profiles
## derived from TSO500
##
## We quantify CIN signatures in profiles coming from sWGS and TSO500 techs, and then
## compare classifier predictions
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
target_sigs <- readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_OV_tissue_TSO500.rds'))
tso500_power <- fread(paste0(input_dir,"/TSO500_power.tsv"))
tso500_power$sample=paste0("JBLAB-",gsub("JBLAB", "", tso500_power$sample))
swgs_sigs <- rbind(readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_OV_tissue_new_50kb.rds')),
                   readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_OV_tissue_old_50kb.rds')))
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

## Apply signature-threshold in exposures to trust in low activities => this is already applied but just to double-check
swgs_sigs = do.call(cbind,lapply(1:ncol(swgs_sigs), function(thisSig){
  sig=colnames(swgs_sigs)[thisSig]
  acts=swgs_sigs[,thisSig]
  thrs=thresh[names(thresh)==sig]
  acts[acts<thrs]=0
  return(acts)
}))
colnames(swgs_sigs)=names(thresh)
swgs_sigs=as.data.frame(t(apply(swgs_sigs, 1, function(x) x/sum(x))))
rownames(swgs_sigs) <- unlist(lapply(strsplit(gsub('downsampled_', '', rownames(swgs_sigs), fixed = T), '-Tseq'), FUN = function(x) x[1]))
rownames(swgs_sigs)=gsub('-sWGS-EXP09', '', rownames(swgs_sigs))
rownames(swgs_sigs)=gsub('-sWGS-EXP10', '', rownames(swgs_sigs))
swgs_sigs=swgs_sigs[row.names(swgs_sigs)%in%ids$samplename,] # include only those samples that have good profiles

target_sigs = do.call(cbind,lapply(1:ncol(target_sigs), function(thisSig){
  sig=colnames(target_sigs)[thisSig]
  acts=target_sigs[,thisSig]
  thrs=thresh[names(thresh)==sig]
  acts[acts<thrs]=0
  return(acts)
}))
colnames(target_sigs)=names(thresh)
row.names(target_sigs)=paste0("JBLAB-",gsub("JBLAB", "", row.names(target_sigs)))
target_sigs=as.data.frame(t(apply(target_sigs, 1, function(x) x/sum(x))))


################################################################################
########### COMPARE TISSUE & PLASMA EXPOSURES ##################################
signatures=c(paste0("CX",1:17))

target_exps=cbind(sample=row.names(target_sigs),target_sigs)
target_exps=melt(target_exps)
target_exps$OV04_ID=ids$OV04_ID[match(target_exps$sample, ids$samplename)]
target_exps$condition="tso500"
target_exps=target_exps[!is.na(target_exps$OV04_ID),] #remove samples that are not in sWGS set

swgs_exps=cbind(sample=row.names(swgs_sigs),swgs_sigs)
swgs_exps=melt(swgs_exps)
swgs_exps=swgs_exps[swgs_exps$sample%in%target_exps$sample,]
swgs_exps$OV04_ID=ids$OV04_ID[match(swgs_exps$sample, ids$samplename)]
swgs_exps$condition="sWGS"

pair_exps=rbind(target_exps,swgs_exps)
pair_exps$variable <- factor(pair_exps$variable, levels = paste0('CX', 1:17))
pair_exps$OV04_ID <- factor(pair_exps$OV04_ID,
                            levels = as.character(sort(as.numeric(unique(pair_exps$OV04_ID)))))

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
  facet_wrap(~OV04_ID, nrow = 1) +
  ylab(c("Pancancer signature activity")) +
  xlab("") +
  ggtitle("")  +
  theme_bw() +
  theme(legend.key.height =  unit(0.3, 'cm'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())
ggsave(file.path(figs_dir, 'OV_tso500_tissue_signature_concordance.svg'),
       sig_comparison, width = 150/25.4, height = 80/25.4)


## Cosine similarities between tissue and plasma exposures
id=unique(pair_exps$OV04_ID)
cosine <- c()
for(thisSamp in id){
  swgs <- swgs_sigs[row.names(swgs_sigs)%in%pair_exps$sample[pair_exps$OV04_ID==thisSamp], ]
  tso <- target_sigs[row.names(target_sigs)%in%pair_exps$sample[pair_exps$OV04_ID==thisSamp], ]
  cos <- (sum(swgs * tso, na.rm = TRUE)) / 
    (sqrt(sum(swgs^2, na.rm = TRUE)) * sqrt(sum(tso^2, na.rm = TRUE)))
  cosine <- rbind(cosine, c(thisSamp, cos))
}
cosine <- as.data.frame(cosine)
colnames(cosine) <- c("id", "cosine")
cosine$cosine <- as.numeric(cosine$cosine)
mean(cosine$cosine) #0.9165082
sd(cosine$cosine) #0.1026145


################################################################################
############# MAP sWGS & TSO500 PREDICTIONS OF PATIENTS ########################


### Platinum predictions
## Scaling model
mTSO500=target_sigs[,c("CX3", "CX2")]
lModel = list(mean = attributes(scale(mTSO500))$`scaled:center`,
              scale = attributes(scale(mTSO500))$`scaled:scale`)

## TSO500
sc_platinum_target_exposures = sweep(sweep(target_sigs[,c("CX3", "CX2")], 2, lModel$mean, FUN = '-'), 2, lModel$scale, FUN = "/")
colnames(sc_platinum_target_exposures)=c("sCX3","sCX2")
sc_platinum_target_exposures=as.data.frame(sc_platinum_target_exposures)
sc_platinum_target_exposures$tso_id <- row.names(sc_platinum_target_exposures)
sc_platinum_target_exposures$study_subject_id=ids$OV04_ID[match(sc_platinum_target_exposures$tso_id, ids$samplename)]
sc_platinum_target_exposures$platinum_prediction_tso <- ifelse(sc_platinum_target_exposures$sCX3 > sc_platinum_target_exposures$sCX2,
                                                               'Plat sensitive',
                                                               'Plat resistant')

## Compare sWGS and tso500 predictions
platinum_data=left_join(platinum_data,sc_platinum_target_exposures,by="study_subject_id")
platinum_data_filt=platinum_data[!is.na(platinum_data$platinum_prediction_tso),]
table(platinum_data_filt$prediction,platinum_data_filt$platinum_prediction_tso)
#+535 patient => sensitive in both => removed because of clinical data in cox analyses


### Paclitaxel predictions
## Scaling model
mTCGA = as.matrix(tcga_exposures[ , c("CX3", "CX5")])
lModel = list(mean = attributes(scale(mTCGA))$`scaled:center`, 
              scale = attributes(scale(mTCGA))$`scaled:scale`)

# TSO500
sc_paclitaxel_target_exposures = sweep(sweep(swgs_sigs[,c("CX3", "CX5")], 2, lModel$mean, FUN = '-'), 2, lModel$scale, FUN = "/")
colnames(sc_paclitaxel_target_exposures)=c("sCX3","sCX5")
sc_paclitaxel_target_exposures=as.data.frame(sc_paclitaxel_target_exposures)
sc_paclitaxel_target_exposures$tso_id <- row.names(sc_paclitaxel_target_exposures)
sc_paclitaxel_target_exposures$study_subject_id=ids$OV04_ID[match(sc_paclitaxel_target_exposures$tso_id, ids$samplename)]
sc_paclitaxel_target_exposures$taxane_prediction_tso <- ifelse(sc_paclitaxel_target_exposures$sCX5<0,
                                                               'Tax resistant',
                                                               'Tax sensitive')

## Compare sWGS and tso500 predictions
paclitaxel_data=left_join(paclitaxel_data,sc_paclitaxel_target_exposures,by="study_subject_id")
paclitaxel_data_filt=paclitaxel_data[!is.na(paclitaxel_data$taxane_prediction_tso),]
table(paclitaxel_data_filt$prediction, paclitaxel_data_filt$taxane_prediction_tso)
#+535 patient => sensitive in tso500 & resistant in sWGS => removed because of clinical data in cox analyses


### Doxorubicin predictions
## TSO500 predictions
dox_target_exposures <- as.data.frame(target_sigs[,c("CX8","CX9","CX13")])
dox_target_exposures$tso_id <- row.names(dox_target_exposures)
dox_target_exposures$study_subject_id=ids$OV04_ID[match(dox_target_exposures$tso_id, ids$samplename)]
dox_target_exposures$doxorubicin_prediction_tso <- ifelse(dox_target_exposures$CX8 > 0.01 | dox_target_exposures$CX9 > 0.009 | dox_target_exposures$CX13 > 0.009, 
                                                          'Dox resistant', 'Dox sensitive')
## Compare sWGS and tso500 predictions
doxorubicin_data=left_join(doxorubicin_data,dox_target_exposures,by="study_subject_id")
doxorubicin_data_filt=doxorubicin_data[!is.na(doxorubicin_data$doxorubicin_prediction_tso),]
table(doxorubicin_data_filt$prediction, doxorubicin_data_filt$doxorubicin_prediction_tso)
#+535 patient => resistant in both => removed because of clinical data in cox analyses


## TOTAL NUMBERS => 100% concordance
# Sensitive sWGS & sensitive tso500 = 6
# Sensitive sWGS & resistant tso500 = 0
# Resistant sWGS & sensitive tso500 = 1
# Resistant sWGS & resistant tso500 = 7
