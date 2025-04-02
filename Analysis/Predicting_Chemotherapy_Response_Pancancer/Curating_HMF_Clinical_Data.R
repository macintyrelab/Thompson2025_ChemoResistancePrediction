# This script is for generating the curated HMF clinical history file 'HMF_clinical_history.RDS'

# Clean environment
freshr::freshr()

# Libraries
library(this.path)
library(readr)
library(dplyr)

# Files
pancancer_dir <- dirname(this.path())
base_dir <- dirname(dirname(pancancer_dir))
input_dir <- file.path(base_dir, 'Input_Data')
helper_dir <- file.path(base_dir, 'Analysis/Helper_Scripts')
source(file.path(helper_dir, 'Helper_functions.R'))
treatment_responses <- read_delim(file.path(input_dir, "treatment_responses.tsv"),
                                  delim = "\t",
                                  escape_double = FALSE,
                                  trim_ws = TRUE)
pre_biopsy_drugs <- read_delim(file.path(input_dir, "pre_biopsy_drugs.tsv"),
                               delim = "\t",
                               escape_double = FALSE,
                               trim_ws = TRUE)
post_biopsy_drugs <- read_delim(file.path(input_dir, "post_biopsy_drugs.tsv"),
                                delim = "\t",
                                escape_double = FALSE,
                                trim_ws = TRUE)
metadata <- read_delim(file.path(input_dir, "metadata.tsv"),
                       delim = "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)
metadata$patientIdentifier <- gsub('TI*V*$', '', metadata$sampleId, perl=TRUE)
hmf_segments <- readRDS(file.path(input_dir, 'hmf_cnp_smooth.rds'))
  

# Combine pre and post-biopsy drugs for patients that received post-biopsy treatment
combined_drugs <- data.frame(patientIdentifier=c(),
                             sampleId=c(),
                             startDate=c(),
                             endDate=c(),
                             name=c(),
                             type=c(),
                             mechanism=c(),
                             biopsy_state=c())
for(i in unique(post_biopsy_drugs$patientIdentifier)){
  temp_pre_df <- NULL
  sampleIds <- unique(post_biopsy_drugs$sampleId[post_biopsy_drugs$patientIdentifier==i])
  temp_pre_df <- pre_biopsy_drugs %>% filter(patientIdentifier==i)
  if(length(sampleIds)>0){
    temp_post_df <- NULL
    temp_pre_df_sample <- NULL
    temp_post_df <- post_biopsy_drugs %>% filter(sampleId %in% sampleIds)
    temp_df <- NULL
    if(nrow(temp_pre_df)>0){
      temp_pre_df_sample <- temp_pre_df %>%
        mutate(sampleId=sampleIds[1]) %>%
        select('patientIdentifier',
               'sampleId',
               'startDate',
               'endDate',
               'name',
               'type',
               'mechanism')
      temp_df <- rbind(temp_pre_df_sample, temp_post_df)
      combined_drugs <- rbind(combined_drugs, temp_df)
    }else{
      combined_drugs <- rbind(combined_drugs, temp_post_df)
    }
  }else{
    temp_pre_df <- temp_pre_df %>%
      mutate(sampleId=j) %>%
      select('patientIdentifier',
             'sampleId',
             'startDate',
             'endDate',
             'name',
             'type',
             'mechanism')
    combined_drugs <- rbind(combined_drugs, temp_pre_df)
  }
}

combined_drugs <- left_join(combined_drugs,
                            metadata %>%
                              mutate(dup=duplicated(patientIdentifier)) %>%
                              filter(!dup) %>%
                              select(patientIdentifier, cancer_type=primaryTumorLocation))

combined_drugs_filtered <- combined_drugs %>%
  filter(endDate!='null') %>%
  mutate(endDate=as.Date(endDate)) %>%
  filter(endDate>1950 & startDate>1950 & (endDate>=startDate)) %>%
  group_by(sampleId) %>%
  arrange(startDate)


# Combine overlapping treatment lines
for(i in 1:length(unique(combined_drugs_filtered$patientIdentifier))){
  print(i)
  if(i==1){
    temp_df <- assign_treatment_lines(combined_drugs_filtered %>%
                                        filter(patientIdentifier==unique(combined_drugs_filtered$patientIdentifier)[i]))
  }else{
    temp_df <- rbind(temp_df,
                     assign_treatment_lines(combined_drugs_filtered %>%
                                              filter(patientIdentifier==unique(combined_drugs_filtered$patientIdentifier)[i])))
  }
}

combined_drugs_collapsed <- temp_df %>%
  group_by(patientIdentifier, treatment_line, cancer_type) %>%
  summarise(name=paste(sort(unique(name)), collapse='/'),
            startDate=min(startDate),
            endDate=max(endDate)) %>%
  mutate(treatment_length=endDate-startDate)


# Calculate TTF
for(i in 1:length(unique(combined_drugs_collapsed$patientIdentifier))){
  id <- unique(combined_drugs_collapsed$patientIdentifier)[i]
  temp_drugs <- combined_drugs_collapsed %>% filter(patientIdentifier==id)
  temp_responses <- treatment_responses %>% filter(patientIdentifier==id)
  if(i==1){
    combined_drugs_TTF <- calc_TTF_HMF(temp_drugs, temp_responses)
  }else{
    combined_drugs_TTF <- rbind(combined_drugs_TTF, calc_TTF_HMF(temp_drugs, temp_responses))
  }
}


# Annotate with biopsy status
combined_drugs_TTF$biopsyState <- 'Unknown'
for(i in unique(combined_drugs_TTF$patientIdentifier)){
  bdate <- unique(metadata$biopsyDate[metadata$patientIdentifier==i])[1]
  if(bdate!='null'){
    bdate <- as.Date(bdate)
    combined_drugs_indexes <- which(combined_drugs_TTF$patientIdentifier==i)
    pre_biopsy_dates <- which(combined_drugs_TTF$startDate[combined_drugs_indexes] < bdate & combined_drugs_TTF$endDate[combined_drugs_indexes] <= bdate)
    post_biopsy_dates <- which(combined_drugs_TTF$startDate[combined_drugs_indexes] >= bdate & combined_drugs_TTF$endDate[combined_drugs_indexes] > bdate)
    biopsy_dates <- which(combined_drugs_TTF$startDate[combined_drugs_indexes] <= bdate & combined_drugs_TTF$endDate[combined_drugs_indexes] >= bdate)
    combined_drugs_TTF$biopsyState[combined_drugs_indexes[pre_biopsy_dates]] <- 'Pre'
    combined_drugs_TTF$biopsyState[combined_drugs_indexes[post_biopsy_dates]] <- 'Post'
    combined_drugs_TTF$biopsyState[combined_drugs_indexes[biopsy_dates]] <- 'Biopsy'
  }
}

saveRDS(combined_drugs_TTF, file.path(input_dir, 'HMF_clinical_data.RDS'))
