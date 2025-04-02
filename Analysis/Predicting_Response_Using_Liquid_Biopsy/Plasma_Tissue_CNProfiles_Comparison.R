################################################################################
## This script is used for comparing copy number profiles derived from plasma
## and tissue samples from 41 HGSOC samples
##
## We use CNpare tool for doing this comparison. Copy number profiles from both  
## tissue and ctDNA were derived from sWGS data.
################################################################################

# Clean environment
freshr::freshr()

# Libraries
library(data.table)
library(this.path)

# Paths
liquid_dir <- dirname(this.path())
base_dir <- dirname(dirname(liquid_dir))
input_dir <- file.path(base_dir, "Input_Data")
figs_dir <- file.path(base_dir, "Figures")
supp_figs_dir <- file.path(figs_dir, 'Supplementary')

# Files
tissue_copynumbers <- readRDS(paste0(input_dir, "/downsampled_41_OV04_absolute_copyNumbers.rds"))
liquid_copynumbers <- readRDS(paste0(input_dir, "/plasma_absolute_copy_number_profiles_50kb_correct.rds"))
metadata <- readRDS(paste0(input_dir, "/patient_data.rds"))

ids <- fread(paste0(input_dir, "/OV04_JBLAB_id_mapping.csv"))
qc <- readRDS(paste0(input_dir,"/final_patients.RDS"))
ids <- ids[ids$OV04_ID%in%qc$study_subject_id,]
ids <- ids[!ids$samplename%in%c("JBLAB-17049","JBLAB-17058"),] #we use the sample from the original cohort


################################################################################
########### COMPARE CN PROFILES FROM LIQUID & TISSUE SAMPLES ###################

### 1: Get segTable
## Tissue profiles
tissue_cn_table <- do.call(rbind, lapply(colnames(tissue_copynumbers), function(thisSamp){
    tab <- as.data.frame(CNpare::getSegTable(tissue_copynumbers[, thisSamp]))
    tab <- cbind(tab, sample = rep(thisSamp, nrow(tab)))
    return(tab)
}))
tissue_cn_table[, 2:4] <- sapply(2:4, function(x) as.numeric(tissue_cn_table[, x]))


## Plasma profiles
liquid_cn_table <- do.call(rbind, lapply(colnames(liquid_copynumbers), function(thisSamp){
    tab <- as.data.frame(CNpare::getSegTable(liquid_copynumbers[, thisSamp]))
    tab <- cbind(tab, sample = rep(thisSamp, nrow(tab)))
    return(tab)
}))
liquid_cn_table[, 2:4] <- sapply(2:4, function(x) as.numeric(liquid_cn_table[, x]))



### 2: Plot differences
metadata <- metadata[!is.na(metadata$plasma_id), ]
#Remove patient 78 because plasma sample was wrongly labeled in the lab & to be on the safe side
metadata <- metadata[metadata$study_subject_id != 78, ]

percentage_diff <- c()
for(thisSamp in unique(liquid_cn_table$sample)){
    print(thisSamp)
    thisSamp <- gsub("downsampled_", "", thisSamp)
    thisTiss <- metadata$tissue_primary_samplename[metadata$plasma_samplename == thisSamp]
    
    if(length(thisTiss)==1){ #not all 12 plasma samples have tissue samples
        tiss <- unique(tissue_cn_table$sample[grepl(thisTiss, tissue_cn_table$sample)])
        plas <- unique(liquid_cn_table$sample[grepl(thisSamp, liquid_cn_table$sample)])
        diff<-CNpare::getDifference(events = tissue_cn_table[tissue_cn_table$sample == tiss, ], 
                                    events_2 = liquid_cn_table[liquid_cn_table$sample == plas, ],
                                    method = "non-normalized")
        percentage_diff <- rbind(percentage_diff,
                                 c(metadata$study_subject_id[metadata$tissue_primary_samplename == thisTiss],
                                   thisSamp, thisTiss, diff))
        
        pdf(file.path(supp_figs_dir, paste0(thisSamp, "_CNpare_profiles_plasma.pdf")),
            width = 7, height = 4)
        p <- CNpare::CNPlot_events(tissue_cn_table[tissue_cn_table$sample == tiss, ],
                                   liquid_cn_table[liquid_cn_table$sample == plas, ],
                                   method_diff = "non-normalized", plot_diff = FALSE)
        print(p)
        dev.off()
    }
    
}
percentage_diff <- as.data.frame(percentage_diff)
colnames(percentage_diff) <- c("subject_id", "plasma_id", "tissue_id", "diff")
percentage_diff$diff <- as.numeric(percentage_diff$diff)
mean(percentage_diff$diff) #20.83051
sd(percentage_diff$diff) #14.4484

