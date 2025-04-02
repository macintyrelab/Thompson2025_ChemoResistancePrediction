################################################################################
## This script is used for comparing copy number profiles derived from TSO500 and
## sWGS of tissue OV04 samples
##
## We use CNpare tool for doing this comparison. Copy number profiles from both
## technologies are derived from the same tissue sample
################################################################################

# Clean environment
freshr::freshr()

# Libraries
library(data.table)
library(this.path)

# Paths
target_dir <- dirname(this.path())
base_dir <- dirname(dirname(target_dir))
input_dir <- file.path(base_dir, "Input_Data")
figs_dir <- file.path(base_dir, "Figures")
supp_figs_dir <- file.path(figs_dir, 'Supplementary')

# Files
swgs_copynumbers_old <- readRDS(paste0(input_dir, "/downsampled_41_OV04_absolute_copyNumbers.rds"))
swgs_copynumbers_new <- readRDS(paste0(input_dir, "/IUK_cohort_downsampled_SLX-24424_absolute_copyNumbers_50kb.rds"))

target_copynumbers <- readRDS(paste0(input_dir,"/TSO500_OV04_absoluteCN.rds"))
target_copynumbers$sample=paste0("JBLAB-",gsub("JBLAB", "", target_copynumbers$sample))

ids <- fread(paste0(input_dir, "/OV04_JBLAB_id_mapping.csv"))
qc <- readRDS(paste0(input_dir,"/final_patients.RDS"))
ids <- ids[ids$OV04_ID%in%qc$study_subject_id,]
ids <- ids[!ids$samplename%in%c("JBLAB-17049","JBLAB-17058"),] #we use the sample from the previous cohort


################################################################################
########### COMPARE CN PROFILES FROM LIQUID & TISSUE SAMPLES ###################

### 1: Get segTable
## sWGS profiles
# old
swgs_cn_table <- do.call(rbind, lapply(colnames(swgs_copynumbers_old), function(thisSamp){
    tab <- as.data.frame(CNpare::getSegTable(swgs_copynumbers_old[, thisSamp]))
    tab <- cbind(tab, sample = rep(thisSamp, nrow(tab)))
    return(tab)
}))

# add new
for(thisSamp in colnames(swgs_copynumbers_new)){
    tab <- as.data.frame(CNpare::getSegTable(swgs_copynumbers_new[, thisSamp]))
    tab <- cbind(tab, sample = rep(thisSamp, nrow(tab)))
    swgs_cn_table <- rbind(swgs_cn_table, tab)
}
swgs_cn_table[, 2:4] <- sapply(2:4, function(x) as.numeric(swgs_cn_table[, x]))
swgs_cn_table$sample <- gsub("downsampled_", "", swgs_cn_table$sample)


### 2: Plot differences
samples=unique(target_copynumbers$sample)
samples=samples[samples%in%ids$samplename]

percentage_diff <- c()
for(thisSamp in unique(target_copynumbers$sample)){
    print(thisSamp)
    thisTiss <- unique(swgs_cn_table$sample[grepl(thisSamp,swgs_cn_table$sample)])
    thisID <- ids$OV04_ID[ids$samplename == thisSamp]
    if(length(thisID)>0){
        
        diff<-CNpare::getDifference(events = swgs_cn_table[swgs_cn_table$sample == thisTiss, ], 
                                    events_2 = target_copynumbers[target_copynumbers$sample == thisSamp, ],
                                    method = "non-normalized")
        percentage_diff <- rbind(percentage_diff,
                                 c(ids$OV04_ID[ids$samplename == thisSamp],
                                   thisSamp, thisTiss, diff))
        
        pdf(file.path(supp_figs_dir, paste0(thisSamp, "_CNpare_profiles_tso500.pdf")),
            width = 7, height = 4)
        p <- CNpare::CNPlot_events(swgs_cn_table[swgs_cn_table$sample == thisTiss, ],
                                   target_copynumbers[target_copynumbers$sample == thisSamp, ],
                                   method_diff = "non-normalized", plot_diff = FALSE)
        print(p)
        dev.off()
    }
}
percentage_diff <- as.data.frame(percentage_diff)
colnames(percentage_diff) <- c("subject_id", "tso500_id", "sWGS_id", "diff")
percentage_diff$diff <- as.numeric(percentage_diff$diff)
mean(percentage_diff$diff) #19.91916
sd(percentage_diff$diff) #11.29797

