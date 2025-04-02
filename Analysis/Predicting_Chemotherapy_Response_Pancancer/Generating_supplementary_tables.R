# Clean environment
freshr::freshr()

# Libraries
library(readr)
library(readxl)
library(survival)
library(tidyr)
library(this.path)
library(rsvg)
library(DiagrammeR)
library(DiagrammeRsvg)
library(data.table)
library(foreach)
library(dplyr)
library(openxlsx)

# Paths
pancancer_dir <- dirname(this.path())
base_dir <- dirname(dirname(pancancer_dir))
input_dir <- file.path(base_dir, 'Input_Data')
scripts_dir <- file.path(base_dir, 'Analysis/Helper_Scripts')
xml_dir <- file.path(input_dir, 'TCGA_XMLs')
tabs_dir <- file.path(base_dir, 'Tables')

# Files
source(file.path(scripts_dir, 'Helper_functions.R'))
source(file.path(scripts_dir, 'power_calc.R'))


################################################################################

# #### PROCESS INPUT DATA ####
# # HMF
hmf_clin <- readRDS(paste0(input_dir,'/HMF_clinical_data.RDS'))
metadata <- read_delim(paste0(input_dir,'/metadata.tsv'),
                       delim = "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)
metadata$patientIdentifier <- gsub('TI*V*$', '', metadata$sampleId, perl=TRUE)


# SUMMARIZE CLINICAL DATA FOR EACH COHORT

summarize_generic <- function (df, values, colname, unknown_str=NA) {
  freqs <- c()
  for( value in values) {
    filter <- grepl(value, df[[colname]])
    freqs <- append(freqs,nrow(df[filter, ]))
  }
  filterNA <- is.na(df[[colname]])
  if (!is.na(unknown_str)){
    filterNA <- filterNA | grepl(unknown_str, df[[colname]])
  }
  freqs <- append(freqs,nrow(df[filterNA, ]))

  return(data.frame(Count = freqs, row.names = append(values,"Unknown")))

}

summarize_generic_stage <- function(cohort, colname){
  stages = c("I","II","III","IV")
  freqs <- c()
  for( stage in stages) {
    filter <- stringr::str_count(cohort[[colname]], "I") == stringr::str_count(stage, "I") &
              stringr::str_count(cohort[[colname]], "V") == stringr::str_count(stage, "V") &
              !is.na(cohort[[colname]])
    freqs <- append(freqs,nrow(cohort[filter, ]))
  }
  filterNA <- is.na(cohort[[colname]])
  filterNA <- filterNA | grepl("X", cohort[[colname]])
  freqs <- append(freqs,nrow(cohort[filterNA, ]))

  return(data.frame(Count = freqs, row.names = append(stages,"Unknown")))
}

summarize_age_p2 <- function(c, resistant){
  cohort <- if (resistant) c[c$prediction == "Resistant",] else c[c$prediction == "Sensitive",]
  colname <- if ("age" %in% colnames(cohort)) "age" else "Age at diagnosis"

  df <- sapply(cohort[colname], function(x) {
    table(cut(x,
              breaks = c(0, 49, 59, 69, 79, Inf),
              labels = c("1. <50", "2. 50-59", "3. 60-69", "4. 70-79", "5. >80")), useNA = "ifany")
  })

  colnames(df)[colnames(df)==colname] <- "Count"
  return(df)
}

summarize_tumour_stage_p2 <- function(c, resistant){
  cohort <- if (resistant) c[c$prediction == "Resistant",] else c[c$prediction == "Sensitive",]
  colname <- "Tumor Stage"
  df <- summarize_generic_stage(cohort, colname)

  return(df)
}

summarize_tumour_grade_p2 <- function(c, resistant){
  cohort <- if (resistant) c[c$prediction == "Resistant",] else c[c$prediction == "Sensitive",]
  stages <- c("1","2","3","4")
  colname <- if(is.null(cohort$grade)) "Tumor Grade" else "grade"
  df <- summarize_generic(cohort, stages, colname, unknown_str="X")

  return(df)
}

summarize_treatment_lines_p2 <- function(c, resistant){
  cohort <- if (resistant) c[c$prediction == "Resistant",] else c[c$prediction == "Sensitive",]
  stages <- unique(cohort$treatment_line)
  colname <- "treatment_line"
  df <- summarize_generic(cohort, stages, colname)

  return(df)
}

summarize_censoring_p2 <- function(c, resistant){
  cohort <- if (resistant) c[c$prediction == "Resistant",] else c[c$prediction == "Sensitive",]
  stages = c(1,0)
  colname <- "censoring"
  df <- summarize_generic(cohort, stages, colname)

  return(df)
}

summarize_age_p3 <- function(c, resistant, arm){
  cohort <- if (resistant) c[[1]] else c[[2]]
  colname <- if ("age" %in% colnames(cohort)) "age" else "Age at diagnosis"
  row_names = c("1. <50", "2. 50-59", "3. 60-69", "4. 70-79", "5. >80")

  df <- sapply(cohort[cohort$SA_treatment == arm,][colname], function(x) {
    table(cut(x,
              breaks = c(0, 49, 59, 69, 79, Inf),
              labels = row_names), useNA = "always")
  })
  rownames(df) <- append(row_names, "Unknown")
  colnames(df)[colnames(df)==colname] <- "Count"
  return(df)
}

# Only for OV
summarize_tumour_stage_p3 <- function(c, resistant, arm){
  cohort <- if (resistant) c[[1]] else c[[2]]
  colname <- "Tumor Stage"
  df <- summarize_generic_stage(cohort[cohort$SA_treatment == arm,], colname)

  return(df)
}

summarize_tumour_grade_p3 <- function(c, resistant, arm){
  cohort <- if (resistant) c[[1]] else c[[2]]
  stages <- c("1","2","3","4")
  colname <- if(is.null(cohort$grade)) "Tumor Grade" else "grade"
  df <- summarize_generic(cohort[cohort$SA_treatment == arm,], stages, colname, unknown_str="X")

  return(df)
}

summarize_clinical_stage_p3 <- function(c, resistant, arm){
  cohort <- if (resistant) c[[1]] else c[[2]]
  colname <- "clin_stage"
  df <- summarize_generic_stage(cohort[cohort$SA_treatment == arm,], colname)

  return(df)
}

summarize_pathological_stage_p3 <- function(c, resistant, arm){
  cohort <- if (resistant) c[[1]] else c[[2]]
  colname <- "path_stage"
  df <- summarize_generic_stage(cohort[cohort$SA_treatment == arm,], colname)

  return(df)
}

summarize_brca_mut_status_p3 <- function(c, resistant, arm){
  cohort <- if (resistant) c[[1]] else c[[2]]
  stages = c("WT BRCA1/2 TCGA","WT BRCA1/2","WT BRCA1/2+LOH","somatic BRCA1","germline BRCA1",
             "somatic BRCA1+LOH","germline BRCA1+LOH","somatic BRCA2","germline BRCA2","somatic BRCA2+LOH",
             "germline BRCA2+LOH","BRCA1 Hypermethyl.","RAD51C Hypermethyl.")
  colname <- "brca_mut_status"
  df <- summarize_generic(cohort[cohort$SA_treatment == arm,], stages, colname)

  return(df)
}

summarize_brca_subtype_p3 <- function(c, resistant, arm){
  cohort <- if (resistant) c[[1]] else c[[2]]
  stages = c("ER","ER/HER2","ER/PR","ER/PR/HER2","HER2","PR","PR/HER2","TNBC")
  colname <- "type"
  df <- summarize_generic(cohort[cohort$SA_treatment == arm,], stages, colname, unknown_str="Other")

  return(df)
}

summarize_sex_p3 <- function(c, resistant, arm){
  cohort <- if (resistant) c[[1]] else c[[2]]
  stages = c("male","female")
  colname <- "gender"
  df <- summarize_generic(cohort[cohort$SA_treatment == arm,], stages, colname)

  return(df)
}

summarize_tumour_subtypes_p3 <- function(c, resistant, arm){
  cohort <- if (resistant) c[[1]] else c[[2]]
  stages <- unique(cohort$primaryTumorSubType)
  stages <- stages[!is.na(stages)]
  colname <- "primaryTumorSubType"
  df <- summarize_generic(cohort[cohort$SA_treatment == arm,], stages, colname)

  return(df)
}

summarize_biopsy_site_p3 <- function(c, resistant, arm){
  cohort <- if (resistant) c[[1]] else c[[2]]
  stages <- unique(cohort$biopsySite)
  stages <- stages[!is.na(stages)]
  colname <- "biopsySite"
  df <- summarize_generic(cohort[cohort$SA_treatment == arm,], stages, colname)

  return(df)
}


generate_joined_table_p2 <- function (func, cohort) {
  df <- merge(
      as.data.frame(func(cohort,TRUE)),
      as.data.frame(func(cohort,FALSE)),
      by=0,
      all.x=TRUE,
      all.y=TRUE
    )

  colnames(df)[colnames(df)=='Count.x'] <- "Resistant"
  colnames(df)[colnames(df)=='Count.y'] <- "Sensitive"
  return(df)
}

generate_joined_table_p3 <- function (func, cohort) {
  # print(func)
  df <- merge(
    merge(
      as.data.frame(func(cohort,TRUE,TRUE)),
      as.data.frame(func(cohort,TRUE,FALSE)),
      by=0,
      all.x=TRUE,
      all.y=TRUE
    ),
    merge(
      as.data.frame(func(cohort,FALSE,TRUE)),
      as.data.frame(func(cohort,FALSE,FALSE)),
      by=0,
      all.x=TRUE,
      all.y=TRUE
    ),
    by=0,
    all.x=TRUE
  )
  df$Row.names <- NULL
  df$Row.names.y <- NULL

  colnames(df)[colnames(df)=='Count.x.x'] <- "Resistant_Exp"
  colnames(df)[colnames(df)=='Count.y.x'] <- "Resistant_Control"
  colnames(df)[colnames(df)=='Count.x.y'] <- "Sensitive_Exp"
  colnames(df)[colnames(df)=='Count.y.y'] <- "Sensitive_Control"

  return(df)
}

make_worksheet <- function(cancer, chemo, cohort, phase, cohort_summary, wb){
    title_style <- createStyle(halign='center', textDecoration = 'bold', fontSize = 15)
    subtitle_style <- createStyle(halign='center', textDecoration = 'bold', fontSize = 12)
    characteristics_title_style <- createStyle(textDecoration = 'bold', fontSize = 11)
    nonzero_style <- createStyle(textDecoration = 'bold')
    sheet_name <- paste0(toupper(cohort), '_', toupper(cancer), '_', toupper(chemo), '_PHASE', phase)
    if(phase=='II'){
        columns <- 3:4
    }else{
        columns <- 3:6
    }
    addWorksheet(wb, sheet_name)
    setColWidths(wb, sheet_name, cols=columns, 22)
    writeData(wb, sheet_name, paste0(toupper(cohort), '-', toupper(cancer), ' ', toupper(chemo)), startCol=3, startRow=1)
    mergeCells(wb, sheet_name, cols=columns, rows=1)
    addStyle(wb, sheet_name, title_style, rows=1, cols=columns)
    setRowHeights(wb, sheet_name, 1, 20)
    writeData(wb, sheet_name, t(as.data.frame(gsub('_', ' ', colnames(cohort_summary[[1]])[-1]))), startCol=3, startRow=2, colNames=FALSE)
    addStyle(wb, sheet_name, subtitle_style, rows=2, cols=columns)
    row_index <- 3
    for(i in 1:length(cohort_summary)){
        writeData(wb, sheet_name, names(cohort_summary)[i], startCol=1, startRow=row_index)
        writeData(wb, sheet_name, cohort_summary[[i]], startCol=2, startRow=row_index+1, colNames=FALSE)
        row_index <- row_index + nrow(cohort_summary[[i]]) + 1
    }
    setColWidths(wb, sheet_name, cols=1, 19)
    setColWidths(wb, sheet_name, cols=2, 15)
    addStyle(wb, sheet_name, characteristics_title_style, rows=3:row_index, cols=1)
    conditionalFormatting(wb, sheet_name, cols=columns, rows=4:row_index, rule='!=0', style=nonzero_style)
    saveWorkbook(wb, file.path(tabs_dir, 'SupplementaryTable5_Cohort_clinical_characteristics.xlsx'), overwrite = TRUE)
}

##### GENERATE CLINICAL DATA

tcga_ov_plat_p2 <- phase_2_func('OV', 'platinum', 'tcga-ov')
tcga_sarc_dox_p2 <- phase_2_func('SARC', 'doxorubicin', 'tcga')

tcga_brca_tax_p3 <- phase_3_func('BRCA', 'taxane', 'tcga', SA_exp=FALSE, SA_control=FALSE)
tcga_ov_tax_p3 <- phase_3_func('OV', 'taxane', 'tcga-ov', SA_exp=FALSE, SA_control=FALSE)
tcga_ov_dox_p3 <- phase_3_func('OV', 'doxorubicin', 'tcga-ov', SA_exp=FALSE, SA_control=FALSE)
tcga_hnsc_plat_p3 <- phase_3_func('HNSC', 'platinum', 'tcga', SA_exp=FALSE, SA_control=FALSE)
tcga_luad_tax_p3 <- phase_3_func('LUAD', 'taxane', 'tcga', SA_exp=FALSE, SA_control=FALSE)
tcga_lusc_plat_p3 <- phase_3_func('LUSC', 'platinum', 'tcga', SA_exp=FALSE, SA_control=FALSE)
tcga_lusc_tax_p3 <- phase_3_func('LUSC', 'taxane', 'tcga', SA_exp=FALSE, SA_control=FALSE)
tcga_sarc_dox_p3 <- phase_3_func('SARC', 'doxorubicin', 'tcga', SA_exp=FALSE, SA_control=FALSE)
tcga_stad_plat_p3 <- phase_3_func('STAD', 'platinum', 'tcga', SA_exp=FALSE, SA_control=FALSE)
tcga_ucec_plat_p3 <- phase_3_func('UCEC', 'platinum', 'tcga', SA_exp=FALSE, SA_control=FALSE)
tcga_ucec_tax_p3 <- phase_3_func('UCEC', 'taxane', 'tcga', SA_exp=FALSE, SA_control=FALSE)

hmf_bone_dox_p3 <- phase_3_func('Bone/Soft tissue', 'doxorubicin', 'hmf', SA_exp=FALSE, SA_control=FALSE)
hmf_bre_tax_p3 <- phase_3_func('Breast', 'taxane', 'hmf', SA_exp=FALSE, SA_control=FALSE)
hmf_bre_dox_p3 <- phase_3_func('Breast', 'doxorubicin', 'hmf', SA_exp=FALSE, SA_control=FALSE)
hmf_ov_plat_p3 <- phase_3_func('Ovary', 'platinum', 'hmf', SA_exp=FALSE, SA_control=FALSE)
hmf_ov_tax_p3 <- phase_3_func('Ovary', 'taxane', 'hmf', SA_exp=FALSE, SA_control=FALSE)
hmf_ov_dox_p3 <- phase_3_func('Ovary', 'doxorubicin', 'hmf', SA_exp=FALSE, SA_control=FALSE)
hmf_ute_dox_p3 <- phase_3_func('Uterus', 'doxorubicin', 'hmf', SA_exp=FALSE, SA_control=FALSE)
hmf_pro_tax_p3 <- phase_3_func('Prostate', 'taxane', 'hmf', SA_exp=FALSE, SA_control=FALSE)




##### SUMMARIZE

# TCGA phase II

tcga_ov_plat_p2_summary <- list( Age = generate_joined_table_p2(summarize_age_p2,tcga_ov_plat_p2),
                                 Tumour_Stage = generate_joined_table_p2(summarize_tumour_stage_p2,tcga_ov_plat_p2),
                                 Tumour_grade = generate_joined_table_p2(summarize_tumour_grade_p2,tcga_ov_plat_p2),
                                 Treatment_lines = generate_joined_table_p2(summarize_treatment_lines_p2,tcga_ov_plat_p2),
                                 Censoring =  generate_joined_table_p2(summarize_censoring_p2,tcga_ov_plat_p2))
tcga_sarc_dox_p2_summary <- list( Age = generate_joined_table_p2(summarize_age_p2,tcga_sarc_dox_p2),
                                  Treatment_lines = generate_joined_table_p2(summarize_treatment_lines_p2,tcga_sarc_dox_p2),
                                  Censoring =  generate_joined_table_p2(summarize_censoring_p2,tcga_sarc_dox_p2))

# TCGA

tcga_brca_tax_p3_summary <- list(Age = generate_joined_table_p3(summarize_age_p3,tcga_brca_tax_p3),
                                  Clinical_Stage = generate_joined_table_p3(summarize_clinical_stage_p3,tcga_brca_tax_p3),
                                  Pathological_Stage = generate_joined_table_p3(summarize_pathological_stage_p3,tcga_brca_tax_p3),
                                  Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_brca_tax_p3),
                                  BRCA_Mutation_Status = generate_joined_table_p3(summarize_brca_mut_status_p3,tcga_brca_tax_p3),
                                  BRCA_Subtype = generate_joined_table_p3(summarize_brca_subtype_p3,tcga_brca_tax_p3))
tcga_ov_tax_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,tcga_ov_tax_p3),
                                Tumour_Stage = generate_joined_table_p3(summarize_tumour_stage_p3,tcga_ov_tax_p3),
                                Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_ov_tax_p3))
tcga_ov_dox_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,tcga_ov_dox_p3),
                                Tumour_Stage = generate_joined_table_p3(summarize_tumour_stage_p3,tcga_ov_dox_p3),
                                Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_ov_dox_p3))
tcga_hnsc_plat_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,tcga_hnsc_plat_p3),
                                   Clinical_Stage = generate_joined_table_p3(summarize_clinical_stage_p3,tcga_hnsc_plat_p3),
                                   Pathological_Stage = generate_joined_table_p3(summarize_pathological_stage_p3,tcga_hnsc_plat_p3),
                                   Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_hnsc_plat_p3),
                                   BRCA_Mutation_Status = generate_joined_table_p3(summarize_brca_mut_status_p3,tcga_hnsc_plat_p3))
tcga_luad_tax_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,tcga_luad_tax_p3),
                                  Clinical_Stage = generate_joined_table_p3(summarize_clinical_stage_p3,tcga_luad_tax_p3),
                                  Pathological_Stage = generate_joined_table_p3(summarize_pathological_stage_p3,tcga_luad_tax_p3),
                                  Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_luad_tax_p3),
                                  BRCA_Mutation_Status = generate_joined_table_p3(summarize_brca_mut_status_p3,tcga_luad_tax_p3))
tcga_lusc_plat_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,tcga_lusc_plat_p3),
                                   Clinical_Stage = generate_joined_table_p3(summarize_clinical_stage_p3,tcga_lusc_plat_p3),
                                   Pathological_Stage = generate_joined_table_p3(summarize_pathological_stage_p3,tcga_lusc_plat_p3),
                                   Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_lusc_plat_p3),
                                   BRCA_Mutation_Status = generate_joined_table_p3(summarize_brca_mut_status_p3,tcga_lusc_plat_p3))
tcga_lusc_tax_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,tcga_lusc_tax_p3),
                                  Clinical_Stage = generate_joined_table_p3(summarize_clinical_stage_p3,tcga_lusc_tax_p3),
                                  Pathological_Stage = generate_joined_table_p3(summarize_pathological_stage_p3,tcga_lusc_tax_p3),
                                  Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_lusc_tax_p3),
                                  BRCA_Mutation_Status = generate_joined_table_p3(summarize_brca_mut_status_p3,tcga_lusc_tax_p3))
tcga_sarc_dox_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,tcga_sarc_dox_p3),
                                  Clinical_Stage = generate_joined_table_p3(summarize_clinical_stage_p3,tcga_sarc_dox_p3),
                                  Pathological_Stage = generate_joined_table_p3(summarize_pathological_stage_p3,tcga_sarc_dox_p3),
                                  Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_sarc_dox_p3),
                                  BRCA_Mutation_Status = generate_joined_table_p3(summarize_brca_mut_status_p3,tcga_sarc_dox_p3))
tcga_stad_plat_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,tcga_stad_plat_p3),
                                   Clinical_Stage = generate_joined_table_p3(summarize_clinical_stage_p3,tcga_stad_plat_p3),
                                   Pathological_Stage = generate_joined_table_p3(summarize_pathological_stage_p3,tcga_stad_plat_p3),
                                   Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_stad_plat_p3),
                                   BRCA_Mutation_Status = generate_joined_table_p3(summarize_brca_mut_status_p3,tcga_stad_plat_p3))
tcga_ucec_plat_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,tcga_ucec_plat_p3),
                                   Clinical_Stage = generate_joined_table_p3(summarize_clinical_stage_p3,tcga_ucec_plat_p3),
                                   Pathological_Stage = generate_joined_table_p3(summarize_pathological_stage_p3,tcga_ucec_plat_p3),
                                   Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_ucec_plat_p3),
                                   BRCA_Mutation_Status = generate_joined_table_p3(summarize_brca_mut_status_p3,tcga_ucec_plat_p3))
tcga_ucec_tax_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,tcga_ucec_tax_p3),
                                  Clinical_Stage = generate_joined_table_p3(summarize_clinical_stage_p3,tcga_ucec_tax_p3),
                                  Pathological_Stage = generate_joined_table_p3(summarize_pathological_stage_p3,tcga_ucec_tax_p3),
                                  Tumour_grade = generate_joined_table_p3(summarize_tumour_grade_p3,tcga_ucec_tax_p3),
                                  BRCA_Mutation_Status = generate_joined_table_p3(summarize_brca_mut_status_p3,tcga_ucec_tax_p3))

# HMF

hmf_bone_dox_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,hmf_bone_dox_p3),
                                 Sex = generate_joined_table_p3(summarize_sex_p3,hmf_bone_dox_p3),
                                 Tumor_subtype = generate_joined_table_p3(summarize_tumour_subtypes_p3,hmf_bone_dox_p3),
                                 Biopsy_site = generate_joined_table_p3(summarize_biopsy_site_p3,hmf_bone_dox_p3))
hmf_bre_tax_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,hmf_bre_tax_p3),
                                Tumour_subtype = generate_joined_table_p3(summarize_tumour_subtypes_p3,hmf_bre_tax_p3),
                                Biopsy_site = generate_joined_table_p3(summarize_biopsy_site_p3,hmf_bre_tax_p3))
hmf_bre_dox_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,hmf_bre_dox_p3),
                                Tumour_subtype = generate_joined_table_p3(summarize_tumour_subtypes_p3,hmf_bre_dox_p3),
                                Biopsy_site = generate_joined_table_p3(summarize_biopsy_site_p3,hmf_bre_dox_p3))
hmf_ov_plat_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,hmf_ov_plat_p3),
                                Sex = generate_joined_table_p3(summarize_sex_p3,hmf_ov_plat_p3),
                                Tumor_subtype = generate_joined_table_p3(summarize_tumour_subtypes_p3,hmf_ov_plat_p3),
                                Biopsy_site = generate_joined_table_p3(summarize_biopsy_site_p3,hmf_ov_plat_p3))
hmf_ov_tax_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,hmf_ov_tax_p3),
                               Sex = generate_joined_table_p3(summarize_sex_p3,hmf_ov_tax_p3),
                               Tumor_subtype = generate_joined_table_p3(summarize_tumour_subtypes_p3,hmf_ov_tax_p3),
                               Biopsy_site = generate_joined_table_p3(summarize_biopsy_site_p3,hmf_ov_tax_p3))
hmf_ov_dox_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,hmf_ov_dox_p3),
                               Sex = generate_joined_table_p3(summarize_sex_p3,hmf_ov_dox_p3),
                               Tumor_subtype = generate_joined_table_p3(summarize_tumour_subtypes_p3,hmf_ov_dox_p3),
                               Biopsy_site = generate_joined_table_p3(summarize_biopsy_site_p3,hmf_ov_dox_p3))
hmf_ute_dox_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,hmf_ute_dox_p3),
                                Sex = generate_joined_table_p3(summarize_sex_p3,hmf_ute_dox_p3),
                                Tumor_subtype = generate_joined_table_p3(summarize_tumour_subtypes_p3,hmf_ute_dox_p3),
                                Biopsy_site = generate_joined_table_p3(summarize_biopsy_site_p3,hmf_ute_dox_p3))
hmf_pro_tax_p3_summary <- list( Age = generate_joined_table_p3(summarize_age_p3,hmf_pro_tax_p3),
                                Tumour_subtype = generate_joined_table_p3(summarize_tumour_subtypes_p3,hmf_pro_tax_p3),
                                Biopsy_site = generate_joined_table_p3(summarize_biopsy_site_p3,hmf_pro_tax_p3))

# Generate excel sheet
wb <- createWorkbook()
make_worksheet(cancer = 'OV', cohort = 'tcga', chemo = 'Platinum', phase = 'II', cohort_summary = tcga_ov_plat_p2_summary, wb = wb)
make_worksheet(cancer = 'SARC', cohort = 'tcga', chemo = 'Anthracycline', phase = 'II', cohort_summary = tcga_sarc_dox_p2_summary, wb = wb)

make_worksheet(cancer='BRCA', cohort='tcga', chemo='Taxane', phase='III', cohort_summary = tcga_brca_tax_p3_summary, wb=wb)
make_worksheet(cancer='OV', cohort='tcga', chemo='Taxane', phase='III', cohort_summary = tcga_ov_tax_p3_summary, wb=wb)
make_worksheet(cancer='OV', cohort='tcga', chemo='Anthracycline', phase='III', cohort_summary = tcga_ov_dox_p3_summary, wb=wb)
make_worksheet(cancer='HNSC', cohort='tcga', chemo='Platinum', phase='III', cohort_summary = tcga_hnsc_plat_p3_summary, wb=wb)
make_worksheet(cancer='LUAD', cohort='tcga', chemo='Taxane', phase='III', cohort_summary = tcga_luad_tax_p3_summary, wb=wb)
make_worksheet(cancer='LUSC', cohort='tcga', chemo='Platinum', phase='III', cohort_summary = tcga_lusc_plat_p3_summary, wb=wb)
make_worksheet(cancer='LUSC', cohort='tcga', chemo='Taxane', phase='III', cohort_summary = tcga_lusc_tax_p3_summary, wb=wb)
make_worksheet(cancer='SARC', cohort='tcga', chemo='Anthracyclin', phase='III', cohort_summary = tcga_sarc_dox_p3_summary, wb=wb)
make_worksheet(cancer='STAD', cohort='tcga', chemo='Platinum', phase='III', cohort_summary = tcga_stad_plat_p3_summary, wb=wb)
make_worksheet(cancer='UCEC', cohort='tcga', chemo='Platinum', phase='III', cohort_summary = tcga_ucec_plat_p3_summary, wb=wb)
make_worksheet(cancer='UCEC', cohort='tcga', chemo='Taxane', phase='III', cohort_summary = tcga_ucec_tax_p3_summary, wb=wb)

make_worksheet(cancer='Bone', cohort='hmf', chemo='Anthracycline', phase='III', cohort_summary = hmf_bone_dox_p3_summary, wb=wb)
make_worksheet(cancer='Breast', cohort='hmf', chemo='Taxane', phase='III', cohort_summary = hmf_bre_tax_p3_summary, wb=wb)
make_worksheet(cancer='Breast', cohort='hmf', chemo='Anthracycli', phase='III', cohort_summary = hmf_bre_dox_p3_summary, wb=wb)
make_worksheet(cancer='Ovarian', cohort='hmf', chemo='Platinum', phase='III', cohort_summary = hmf_ov_plat_p3_summary, wb=wb)
make_worksheet(cancer='Ovarian', cohort='hmf', chemo='Taxane', phase='III', cohort_summary = hmf_ov_tax_p3_summary, wb=wb)
make_worksheet(cancer='Ovarian', cohort='hmf', chemo='Anthracycl', phase='III', cohort_summary = hmf_ov_dox_p3_summary, wb=wb)
make_worksheet(cancer='Uterine', cohort='hmf', chemo='Anthracycl', phase='III', cohort_summary = hmf_ute_dox_p3_summary, wb=wb)
make_worksheet(cancer='Prostate', cohort='hmf', chemo='Taxane', phase='III', cohort_summary = hmf_pro_tax_p3_summary, wb=wb)
