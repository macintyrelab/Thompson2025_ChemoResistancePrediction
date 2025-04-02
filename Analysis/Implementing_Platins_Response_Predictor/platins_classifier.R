################################################################################
## This script is used to generalize the classifier for platins prediction to
## all different cohorts. This is crucial for moving the biomarker into the clinic
## Script adapted from Drews et al. Nature 2022
################################################################################

rm(list=ls(all=TRUE))

### LIBRARIES
suppressPackageStartupMessages({
    # Standard packages used throughout repo
    library(this.path)
    library(data.table)
    library(ggplot2)
    library(ggthemes)
    # # For survival analysis
    library(survival)
    library(survminer)
})

### PLOTTING THEMES
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(size = 0.5),
             axis.ticks = element_line(size = 0.5),
             axis.ticks.length = unit(.1, "cm"))

### PATHS
PLATDIR=dirname(this.path())
BASE=dirname(dirname(PLATDIR))
INDIR=paste0(BASE,"/Input_Data")
FIGDIR=paste0(BASE,"/Figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

### LOAD DATA
## Data - TCGA
ACTIVITIES=file.path(INDIR, "Signature_Compendium_v5_Cosine-0.74_Activities_THRESH95_NAMESAPRIL21.rds")
META=file.path(INDIR, "Metadata_TCGA_ASCAT_penalty70.rds")
BRCA.STATE=file.path(INDIR, "TCGA_BRCA_MutStatus.rds")
## From previous survival analysis (Compile_survival_data_onlyclin.R)
SURVIVALOV=file.path(INDIR, "Survival_and_BRCA_Status_TCGA-OV.rds")

LINK=file.path(INDIR,"TCGA_PCAWG_links_plusCancer.rds")

## Data - PCAWG
PCAWG=file.path(INDIR, "PCAWG_signature_activities_THRESH095_NAMESAPRIL21.rds")
METAPCAWG=file.path(INDIR, "PCAWG_1900_samples_CINSig_activities_metadata_plus_HRDetect.rds")
SURVPCAWG=file.path(INDIR, "pcawg_donor_clinical_August2016_v9.fixed.tsv")
TSPCAWG=file.path(INDIR, "pcawg_donor.tsv")

## Data - ICGC (560 Breast cancers)
ICGC560=file.path(INDIR, "ICGC_560BCs_signature_activities_THRESH095_NAMESAPRIL21.rds")
METAICGC560=file.path(INDIR, "Meta_560BreastCancers.tsv")

## Other essential variables
AGEBRACKETS=c(0, 53, 62, 76, 100)

## Load essential files
links = readRDS(LINK)


### FUNCTIONS
## Little function for applying the classifier to a data frame with columns "CX3", "CX2" and "Name"
applyClinClass = function(dtSurvOV, lModel) {
    
    ## Scale the activities
    mTCGAOV = as.matrix(dtSurvOV[ , c("CX3", "CX2")])
    rownames(mTCGAOV) = dtSurvOV$Name
    smTCGAOV = sweep(sweep(mTCGAOV, 2, lModel$mean, FUN = '-'), 2, lModel$scale, FUN = "/")
    
    dtSurvOV$sCX3 = smTCGAOV[,"CX3"]
    dtSurvOV$sCX2 = smTCGAOV[,"CX2"]
    
    ## Apply classifier
    if(identical(dtSurvOV$Name, rownames(smTCGAOV))) {
        dtSurvOV$Classifier = ifelse(dtSurvOV$sCX3 >= dtSurvOV$sCX2, "Predicted sensitive", "Predicted resistant")  
        dtSurvOV$Classifier = factor(dtSurvOV$Classifier, levels = c("Predicted sensitive","Predicted resistant"))
    } else {
        stop("Something went wrong.")
    }
    
    return(dtSurvOV)
    
}

## Bread and butter function of this markdown
plotKMandCoxOS = function(dtSurvOV, COHORTNAME = "TCGA-OV", CLASSIFIERNAME = "CX3/CX2 classifier", COXMODEL = "Full") {
    
    ## Argument: CLASSIFIER = "Classifier"
    
    ### Cannot get the automatic handling of variable CLASSIFIER from within the function. It is not recognised within the function despite being defined.
    # kmOV = survfit(Surv(OS.time, OS) ~ Classifier, data = dtSurvOV)
    # CLASSIFIER = noquote(COLNAME)
    # print(CLASSIFIER)
    # print(exists(CLASSIFIER))
    # kmOV = survfit(Surv(OS.time, OS) ~ get(CLASSIFIER), data = dtSurvOV)
    kmOV = survfit(reformulate(termlabels = CLASSIFIER, response = 'Surv(OS.time, OS)'), data = dtSurvOV)
    # kmOV = survfit(as.formula(paste("Surv(OS.time, OS)~", paste(TEST))), data = dtSurvOV)
    pKMOV = ggsurvplot(kmOV, data = dtSurvOV, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                       xlim = c(0,4000), break.time.by = 365, 
                       ggtheme = theme_tufte(base_size = 6, base_family = "ArialMT"),
                       surv.median.line = "hv", risk.table.y.text.col = TRUE, risk.table.y.text = FALSE) + 
        ggtitle(paste("OS in", COHORTNAME, "patients by", CLASSIFIERNAME))
    
    # Cox model
    # resCoxOV = coxph(Surv(OS.time, OS) ~ Classifier + AgeCategory + TS, data =  dtSurvOV)
    if(COXMODEL == "Full") {
        resCoxOV = coxph(reformulate(termlabels = paste(CLASSIFIER, "AgeCategory", "TS", sep = "+"), response = 'Surv(OS.time, OS)'), data = dtSurvOV)
    } else {
        resCoxOV = coxph(reformulate(termlabels = paste(CLASSIFIER, "AgeCategory", sep = "+"), response = 'Surv(OS.time, OS)'), data = dtSurvOV)
    }
    # Don't attempt to plot if Cox model failed
    nFailedCox = sum(is.na(resCoxOV$coefficients))
    if(nFailedCox != length(resCoxOV$coefficients)) {
        pCoxOV = ggforest(resCoxOV, data = dtSurvOV, main = paste("Cox HR model for", COHORTNAME))  
    } else {
        # Just an empty plot to not break the plotting later on.
        pCoxOV = ggplot() + theme_void() + xlab(NULL) +
            geom_text(aes(0,0,label = paste("Cox HR model for", COHORTNAME, "\nEmpty because Cox model failed to converge.")))
    }
    
    lOut = list(km = kmOV, kmPlot = pKMOV, cox = resCoxOV, coxPlot = pCoxOV)
    return(lOut)
}

### MODEL FOR SCALING ACTIVITIES BASED ON PAN-CANCER BRCAmut
## Get all BRCA1/2mut TCGA samples
brcaMuts = readRDS(BRCA.STATE)
brca1 = unique(as.character(brcaMuts$Sample[ grepl("[germline|somatic] BRCA1", 
                                                   brcaMuts$Status) ]))
brca2 = unique(as.character(brcaMuts$Sample[ grepl("[germline|somatic] BRCA2", 
                                                   brcaMuts$Status) ]))
# rad51c = unique(as.character(brcaMuts$Sample[ grepl("RAD51C", brcaMuts$Status) ]))
brcaMuts = unique(c(brca1,brca2))

# Model for scaling
mBRCA = readRDS(ACTIVITIES)
mBRCA = mBRCA[row.names(mBRCA)%in%brcaMuts,]
mBRCA = as.matrix(mBRCA[ , c("CX3", "CX2")])

## Model for scaling
lModel = list(mean = attributes(scale(mBRCA))$`scaled:center`, 
              scale = attributes(scale(mBRCA))$`scaled:scale`)


################################################################################

#### TCGA-OV ####
## Apply classifier to all TCGA-OV samples
dtSurvOV = readRDS(SURVIVALOV)

## Apply classifier to all TCGA-OV samples
dtSurvOV = applyClinClass(dtSurvOV, lModel)

## Age categorisation chosen from age histogram so that it works across all cohorts equally
dtSurvOV$AgeCategory = cut(dtSurvOV$AgeY, breaks = AGEBRACKETS, labels = 1:(length(AGEBRACKETS)-1))

## Survival analysis
CLASSIFIER = "Classifier"
lTCGAOVCX = plotKMandCoxOS(dtSurvOV, COHORTNAME = "TCGA-OV", CLASSIFIERNAME = "CX3/CX2 classifier")
summary(lTCGAOVCX$cox)
# HR=0.7248 & p-val=0.0067


#### PCAWG-OV ####
## Prepare signature activities
mPCAWG = readRDS(PCAWG)

## Add survival
dtSurv = fread(SURVPCAWG)
metaPCAWG = readRDS(METAPCAWG)
dtPcawg = metaPCAWG[, c("projectcode", "cancer_type", "samplename", "icgc_donor_id", "icgc_sample_id",
                        "donor_age_at_diagnosis", "tumour_type", "histology_tier4", "HRDetect")]

## Check that age at diagnosis is identical in both - just sanity check
dtPcawg$comp = dtSurv$donor_age_at_diagnosis[ match(dtPcawg$icgc_donor_id, dtSurv$icgc_donor_id) ]
summary(dtPcawg$comp - dtPcawg$donor_age_at_diagnosis)
dtPcawg$comp = NULL

dtPcawg$vital_status = dtSurv$donor_vital_status[ match(dtPcawg$icgc_donor_id, 
                                                        dtSurv$icgc_donor_id) ]
dtPcawg$survival_time = dtSurv$donor_survival_time[ match(dtPcawg$icgc_donor_id, 
                                                          dtSurv$icgc_donor_id) ]
dtPcawg$last_follow_up = dtSurv$donor_interval_of_last_followup[ match(dtPcawg$icgc_donor_id, 
                                                                       dtSurv$icgc_donor_id) ]

## Combine survival_time and last_follow_up into one OS.time column
dtPcawg$OS.time = dtPcawg$survival_time
dtPcawg$OS.time[ is.na(dtPcawg$OS.time) ] = dtPcawg$last_follow_up[ is.na(dtPcawg$OS.time) ]


## Create data table with all information
dtPcawg$CX3 = mPCAWG[ match(dtPcawg$samplename, rownames(mPCAWG)), "CX3"]
dtPcawg$CX2 = mPCAWG[ match(dtPcawg$samplename, rownames(mPCAWG)), "CX2"]

## Convert survival
dtPcawg$OS = ifelse(dtPcawg$vital_status == "alive", 0, 1)
dtPcawg$OS[ dtPcawg$vital_status == "" ] = NA

## Add staging
## PCAWG
dtTS = fread(TSPCAWG)
dtPcawg$TS = dtTS$donor_tumour_stage_at_diagnosis[ match(dtPcawg$icgc_donor_id, dtTS$icgc_donor_id) ]

## TCGA
dtPcawg$TCGA = links$TCGA[ match(dtPcawg$icgc_donor_id, links$ICGC) ]

## Stratify by age
## Look at histogram of age and you see three dips, neatly dividing the cohort into four age groups.
dtPcawg$AgeCategory = cut(dtPcawg$donor_age_at_diagnosis, breaks = AGEBRACKETS, labels = 1:(length(AGEBRACKETS)-1))

#### PCAWG-OV-AU
## Use only PCAWG-AU as PCAWG-US is the TCGA-OV cohort

## CX classifier => scaling within the cohort
mPCAWG.OV = as.matrix(dtPcawg[ dtPcawg$projectcode == "OV-AU", c("CX3", "CX2")])

## Scaling
dtOV = applyClinClass(dtPcawg[ dtPcawg$projectcode == "OV-AU", ], lModel)
dtOV$TS[ dtOV$TS == "" ] = NA

## Survival analysis
CLASSIFIER = "Classifier"
lPCAWGAUCX = plotKMandCoxOS(dtOV, COHORTNAME = "PCAWG-AU", CLASSIFIERNAME = "CX3/CX2 classifier")
summary(lPCAWGAUCX$cox)
# HR=0.6437 & p-value=0.0932


#### PCAWG-ESAD ####
mPCAWG.ESAD = as.matrix(dtPcawg[ dtPcawg$cancer_type == "ESAD", c("CX3", "CX2")])

## Scaling
dtESAD = applyClinClass(dtPcawg[ dtPcawg$cancer_type == "ESAD", ], lModel)

## Simpify staging
dtESAD$NewTS = NA
dtESAD$NewTS[ grepl("T1", dtESAD$TS) ] = "I"
dtESAD$NewTS[ grepl("T2", dtESAD$TS) ] = "II"
dtESAD$NewTS[ grepl("T3", dtESAD$TS) ] = "III"
dtESAD$NewTS[ grepl("T4", dtESAD$TS) ] = "IV"
dtESAD$NewTS = factor(dtESAD$NewTS, levels = c("I", "II", "III", "IV")) 
dtESAD$TS = dtESAD$NewTS
dtESAD$NewTS = NULL

## Survival analysis
CLASSIFIER = "Classifier"
lPCAWGESADCX = plotKMandCoxOS(dtESAD, COHORTNAME = "PCAWG-ESAD", CLASSIFIERNAME = "CX3/CX2 classifier")
summary(lPCAWGESADCX$cox)
# HR=0.3489 & p-value=0.0141


#### ICGC ####
## Prepare survival
metaIcgc = fread(METAICGC560)
metaIcgc$OS = ifelse(metaIcgc$donor_vital_status == "alive", 0, 1)
metaIcgc$OS.time = NA

## Combine column of OS.time
## First use interval last check up as base (data points for survived patients)
metaIcgc$OS.time = metaIcgc$`donor_interval_of_last_follow-up_in_DAYS`

## Then add data from deceased people
metaIcgc$OS.time[ metaIcgc$donor_vital_status == "deceased" ] = 
    metaIcgc$donor_survival_time_in_DAYS[ metaIcgc$donor_vital_status == "deceased" ]

## Convert
metaIcgc$OS.time[ metaIcgc$OS.time == "no_data_supplied" ] = NA
metaIcgc$OS.time = as.numeric(metaIcgc$OS.time)

## CX classifier
m560BC = readRDS(ICGC560)
rownames(m560BC) = substr(rownames(m560BC), 1, nchar(rownames(m560BC))-1)
metaIcgc$CX3 = m560BC[ match(metaIcgc$sample_name, rownames(m560BC)), "CX3" ]
metaIcgc$CX2 = m560BC[ match(metaIcgc$sample_name, rownames(m560BC)), "CX2" ]

## Scaling
metaIcgc = applyClinClass(metaIcgc, lModel)

## Convert age and put into categories
metaIcgc$Age = metaIcgc$donor_age_at_diagnosis
metaIcgc$Age[ metaIcgc$Age == "over_80" ] = 80
metaIcgc$Age = as.numeric(metaIcgc$Age)
metaIcgc$AgeCategory = cut(metaIcgc$Age, breaks = AGEBRACKETS, labels = 1:(length(AGEBRACKETS)-1))

## Has tumour stage info already but in a different column
metaIcgc$TS = metaIcgc$T_stage

## Survival analysis
CLASSIFIER = "Classifier"
l560BCCX = plotKMandCoxOS(metaIcgc, COHORTNAME = "ICGC 560 Breast cancers", CLASSIFIERNAME = "CX3/CX2 classifier")
summary(l560BCCX$cox)
# HR=2.6707 & p-value=0.00367


#### PLOTTING: COMPARE RESULTS WITH ORIGINAL PAPER #############################
## Load original results
dfSurv.org=readRDS(file.path(INDIR,"Classifier_SurvResults_Drews2022.rds"))
dfSurv.org$method="original"
dfSurv.org=dfSurv.org[dfSurv.org$test!="TCGA-OV BRCAmut",]

## Table with new results
dfSurv <- data.frame(test = factor(c('TCGA-OV',
                                     'PCAWG-OV',
                                     'PCAWG-ESAD',
                                     'ICGC-BRCA'),
                                   levels = c('TCGA-OV',
                                              'PCAWG-OV',
                                              'PCAWG-ESAD',
                                              'ICGC-BRCA')),
                     lower_hr = c(summary(lTCGAOVCX$cox)$conf.int[1, 3],
                                  summary(lPCAWGAUCX$cox)$conf.int[1, 3],
                                  summary(lPCAWGESADCX$cox)$conf.int[1, 3],
                                  summary(l560BCCX$cox)$conf.int[1, 3]),
                     upper_hr = c(summary(lTCGAOVCX$cox)$conf.int[1, 4],
                                  summary(lPCAWGAUCX$cox)$conf.int[1, 4],
                                  summary(lPCAWGESADCX$cox)$conf.int[1, 4],
                                  summary(l560BCCX$cox)$conf.int[1, 4]),
                     hr = c(summary(lTCGAOVCX$cox)$conf.int[1, 1],
                            summary(lPCAWGAUCX$cox)$conf.int[1, 1],
                            summary(lPCAWGESADCX$cox)$conf.int[1, 1],
                            summary(l560BCCX$cox)$conf.int[1, 1]))

dfSurv$method="new"

## Join both for plotting
dfSurv=rbind(dfSurv,dfSurv.org)
dfSurv$group=paste0(dfSurv$test,"+",dfSurv$method)
dfSurv$group=factor(dfSurv$group,levels = c("TCGA-OV+original","TCGA-OV+new",
                                            "PCAWG-OV+original","PCAWG-OV+new",
                                            "PCAWG-ESAD+original","PCAWG-ESAD+new",
                                            "ICGC-BRCA+original","ICGC-BRCA+new"))

p <- ggplot(dfSurv[dfSurv$test=="PCAWG-ESAD",], aes(x = group, color = method, fill = method)) +
    geom_errorbar(aes(ymin = lower_hr, ymax = upper_hr), width = 0.1, linewidth = 0.8) +
    geom_point(aes(y = hr), size = 3) +
    geom_hline(yintercept = 1, linetype = 3, linewidth = 0.5) +
    scale_y_log10(position = 'right',
                  breaks = (c(0.2, 0.5, 1, 2, 3)),
                  labels = c('0.2', '0.5', '1', '2', '3')) +
    scale_color_manual(values = c("original" = "sienna2", "new" = "skyblue2")) +
    coord_flip() +
    theme_light() +
    theme(#panel.grid = element_blank(),
          #axis.ticks = element_blank(),
          axis.title = element_blank(),
          #panel.border = element_blank(),
          axis.text.y = element_text(size = 6, vjust = 0.5),
          axis.text.x = element_text(size = 6, vjust = 1),
          legend.position = "none")
ggsave(file.path(FIGDIR, 'Platins_Classifier_ESAD_hr.svg'), p, device = 'svg', width = 50, height = 20, units = 'mm')
