library(xml2)
library(readxl)
library(readr)
library(dplyr)
library(broom)
library(grid)
library(this.path)
library(rlang)

# Paths
helper_dir <- dirname(this.path())
base_dir <- dirname(dirname(helper_dir))
input_dir <- file.path(base_dir, 'Input_Data')
xml_dirs <- file.path(input_dir, 'TCGA_XMLs')

##### TCGA Parsing Functions #####

download_gdc <- function(tcga_type){
  tcga_str <- paste0('TCGA-', tcga_type)
  gdc_cases <- GenomicDataCommons::cases() %>%
    GenomicDataCommons::filter(project.project_id==tcga_str) %>%
    GenomicDataCommons::results_all() %>%
    GenomicDataCommons::ids() #GenomicDataCommons loaded separately as it masks dplyr::filter
  clinical_xmls <- GenomicDataCommons::files() %>%
    GenomicDataCommons::filter(cases.case_id==gdc_cases) %>%
    GenomicDataCommons::manifest()
  clinical_xmls <- clinical_xmls[grepl('clinical', clinical_xmls$filename),]
  clinical_xmls <- clinical_xmls[grepl('xml', clinical_xmls$filename),]
  GenomicDataCommons::gdcdata(clinical_xmls$id) # Download GDC files
  gdc_xmls <- list.files('~/.cache/GenomicDataCommons/', recursive=TRUE)
  xml_dir <- file.path(xml_dirs, tcga_type)
  system(paste('mkdir', xml_dir))
  for(i in gdc_xmls){
    system(paste0('mv ~/.cache/GenomicDataCommons/', i, ' ', xml_dir))
  }
  system('rm -rf ~/.cache/GenomicDataCommons/*')
}

get_brca_types <- function(){
  lehmann_tnbc <- read_xlsx(file.path(input_dir, 'pone.0157368.s008.xlsx'),
                            sheet = 'TCGA')
  thennavan_tnbc <- read_xlsx(file.path(input_dir, 'mmc2.xlsx'), 'TCGA Histologic type annotation')
  kalecky_tnbc <- read_xlsx(file.path(input_dir, '12885_2020_6600_MOESM1_ESM.xlsx'),
                            sheet = 'Sheet1')
  tblab_tnbc <- read_xlsx(file.path(input_dir, 'Table_S2_TNBCsubtype_clinical_information_and_signatures.xlsx'),
                          sheet = 'A-TCGA_TNBC_subtype')
  lehmann_tnbc$bcr_patient_barcode <- unlist(lapply(strsplit(lehmann_tnbc$TCGA_SAMPLE, '-'),
                                                    FUN = function(x) paste(x[1:3], collapse = '-')))
  thennavan_tnbc$bcr_patient_barcode <- unlist(lapply(strsplit(thennavan_tnbc$CLID, '-'),
                                                      FUN = function(x) paste(x[1:3], collapse = '-')))
  er_pr_her2 <- data.frame(bcr_patient_barcode = unique(c(lehmann_tnbc$bcr_patient_barcode,
                                                          thennavan_tnbc$bcr_patient_barcode)),
                           ER = 0,
                           PR = 0,
                           HER2 = 0)
  er_pr_her2$ER[er_pr_her2$bcr_patient_barcode %in% intersect(intersect(lehmann_tnbc$bcr_patient_barcode[as.numeric(lehmann_tnbc$ER_CALL) == 1],
                                                                        lehmann_tnbc$bcr_patient_barcode[as.numeric(lehmann_tnbc$ER_PATH) == 1]),
                                                              thennavan_tnbc$bcr_patient_barcode[thennavan_tnbc$er_status_by_ihc == 'Positive'])] <- 1
  er_pr_her2$PR[er_pr_her2$bcr_patient_barcode %in% intersect(lehmann_tnbc$bcr_patient_barcode[as.numeric(lehmann_tnbc$PR_PATH)>0],
                                                              thennavan_tnbc$bcr_patient_barcode[thennavan_tnbc$pr_status_by_ihc == 'Positive'])] <- 1
  er_pr_her2$HER2[er_pr_her2$bcr_patient_barcode %in% intersect(lehmann_tnbc$bcr_patient_barcode[lehmann_tnbc$HER2_PATH_CALL...42 == '1'],
                                                                thennavan_tnbc$bcr_patient_barcode[thennavan_tnbc$HER2.newly.derived == 'Positive'])] <- 1
  
  kalecky_xmls <- GenomicDataCommons::files() %>%
    GenomicDataCommons::filter(cases.case_id == kalecky_tnbc$`GDC Case UUID`) %>%
    GenomicDataCommons::manifest()
  kalecky_xmls <- kalecky_xmls[grepl('clinical', kalecky_xmls$filename), ]
  kalecky_xmls <- kalecky_xmls[grepl('xml', kalecky_xmls$filename), ]
  kalecky_bcr <- character(length = nrow(kalecky_xmls))
  for(i in 1:length(kalecky_bcr)){
    brca_xml <- file.path(xml_dirs,
                          'BRCA/',
                          kalecky_xmls$filename[i])
    kalecky_bcr[i] <- xml2::as_list(read_xml(brca_xml))$tcga_bcr$patient$bcr_patient_barcode[[1]]
  }
  kalecky_tnbc$bcr_patient_barcode <- kalecky_bcr
  lehmann_tnbc <- lehmann_tnbc[lehmann_tnbc$TNBC == 'YES', ]
  thennavan_tnbc <- thennavan_tnbc[thennavan_tnbc$`Triple Negative Status` == 'Yes', ]
  tblab_tnbc$bcr_patient_barcode <- unlist(lapply(strsplit(tblab_tnbc$Sample, '-'),
                                                  FUN = function(x) paste(x[1:3], collapse = '-')))
  
  tnbc <- data.frame(bcr_patient_barcode = unique(c(kalecky_tnbc$bcr_patient_barcode,
                                                    thennavan_tnbc$bcr_patient_barcode,
                                                    tblab_tnbc$bcr_patient_barcode)),
                     Kalecky = 0,
                     Thennavan = 0,
                     TBLab = 0)
  tnbc$Kalecky <- as.numeric(tnbc$bcr_patient_barcode %in% kalecky_tnbc$bcr_patient_barcode)
  tnbc$Thennavan <- as.numeric(tnbc$bcr_patient_barcode %in% thennavan_tnbc$bcr_patient_barcode)
  tnbc$TBLab <- as.numeric(tnbc$bcr_patient_barcode %in% tblab_tnbc$bcr_patient_barcode)
  tnbc$count <- rowSums(tnbc[, 2:4])
  
  brca_type <- merge(tnbc, er_pr_her2,
                     by = 'bcr_patient_barcode',
                     all.x = TRUE,
                     all.y = TRUE)
  brca_type$type <- 'Other'
  for(i in 1:nrow(brca_type)){
    if(ifelse(!is.na(brca_type$count[i]), brca_type$count[i] == 3, FALSE)){
      brca_type$type[i] <- 'TNBC'
    }else if(brca_type$ER[i] == 1 | brca_type$PR[i] == 1 | brca_type$HER2[i] == 1){
      brca_type$type[i] <- paste(colnames(brca_type[6:8])[brca_type[i, 6:8] == 1], collapse = '/')
    }
  }
  return(brca_type)
}

get_wgii_tcga <- function(){
  TCGA_segments <- readRDS(file.path(input_dir, '0_TCGA_Segments_dCIN.rds'))
  wgii <- TCGA_segments %>%
    mutate(length=end-start) %>%
    filter(segVal!=2) %>%
    group_by(sample) %>%
    summarise(wgii=sum(length/3e9) * 100)
  colnames(wgii)[1] <- 'bcr_patient_barcode'
  return(wgii)
}

parse_xmls <- function(xml_paths){
  results_df <- data.frame(bcr_patient_barcode=character(),
                           days_to_start=numeric(),
                           days_to_end=numeric(),
                           drug_name=character(),
                           cycles=numeric(),
                           treatment_line=numeric(),
                           overall_survival=numeric(),
                           overall_survival_reason=character(),
                           days_to_last_followup=numeric(),
                           clin_stage=character(),
                           path_stage=character(),
                           grade=character(),
                           age=numeric(),
                           neoadjuvant=character(),
                           tissue_source_site=character(),
                           days_to_next_treatment=numeric(),
                           TTF=numeric(),
                           progression_reason=character(),
                           TTF_months=numeric(),
                           censoring=numeric())
  for(i in xml_paths){
    print(strsplit(basename(i), '.', fixed=T)[[1]][3])
    xml_results <- parse_xml(i)
    results_df <- rbind(results_df, xml_results)
  }
  results_df
}

parse_xml <- function(xml_path){
  # Load xml
  clin_xml <- xml2::as_list(read_xml(xml_path))
  # Get overall_survival data
  os <- get_os(clin_xml)
  days_to_last_followup <- get_days_to_last_followup(clin_xml)
  clin_stage <- get_clin_stage(clin_xml)
  path_stage <- get_path_stage(clin_xml)
  grade <- get_grade(clin_xml)
  age <- get_age(clin_xml)
  neoadjuvant <- get_neoadjuvant(clin_xml)
  tissue_source_site <- get_tissue_source_site(clin_xml)
  # Make results df
  results_df <- data.frame(bcr_patient_barcode=clin_xml$tcga_bcr$patient$bcr_patient_barcode[[1]],
                           overall_survival=unname(os),
                           overall_survival_reason=names(os),
                           days_to_last_followup=days_to_last_followup,
                           clin_stage=clin_stage,
                           path_stage=path_stage,
                           grade=grade,
                           age=age,
                           neoadjuvant=neoadjuvant,
                           tissue_source_site=tissue_source_site)
  # Get drug dataframe
  drug_df <- parse_drugs(clin_xml)
  results_df <- merge(drug_df, results_df, by='bcr_patient_barcode')
  # Collapse adjuvant treatment
  results_df <- results_df %>% mutate(days_to_next_treatment=lead(days_to_start)-days_to_end)
  results_df <- join_overlapping(results_df)
  # Get NTE dataframe
  nte_df <- pull_days_to_new_tumour_event(clin_xml)
  # Calculate TTF
  results_df <- calc_TTF_TCGA(results_df, nte_df)
  results_df
}

get_clin_stage <- function(xml_obj){
  ifelse(length(xml_obj$tcga_bcr$patient$stage_event$clinical_stage)==1,
         xml_obj$tcga_bcr$patient$stage_event$clinical_stage[[1]],
         NA)
}

get_path_stage <- function(xml_obj){
  ifelse(length(xml_obj$tcga_bcr$patient$stage_event$pathologic_stage)==1,
         xml_obj$tcga_bcr$patient$stage_event$pathologic_stage[[1]],
         NA)
}

get_grade <- function(xml_obj){
  ifelse(length(xml_obj$tcga_bcr$patient$neoplasm_histologic_grade)==1,
         xml_obj$tcga_bcr$patient$neoplasm_histologic_grade[[1]],
         NA)
}

get_age <- function(xml_obj){
  age_days <- ifelse(length(xml_obj$tcga_bcr$patient$days_to_birth)==1,
                     as.numeric(xml_obj$tcga_bcr$patient$days_to_birth)*-1,
                     NA)
  round(age_days/365)
}

get_neoadjuvant <- function(xml_obj){
  ifelse(length(xml_obj$tcga_bcr$patient$history_of_neoadjuvant_treatment)==1,
         xml_obj$tcga_bcr$patient$history_of_neoadjuvant_treatment[[1]],
         NA)
}

get_tissue_source_site <- function(xml_obj){
  ifelse(length(xml_obj$tcga_bcr$patient$tissue_source_site)==1,
         xml_obj$tcga_bcr$patient$tissue_source_site[[1]],
         NA)
}

get_days_to_last_followup <- function(xml_obj){
  ifelse(length(xml_obj$tcga_bcr$patient$days_to_last_followup)==1,
         as.numeric(xml_obj$tcga_bcr$patient$days_to_last_followup[[1]]),
         NA)
}

get_os <- function(xml_obj){
  days_to_last_followup <- ifelse(length(xml_obj$tcga_bcr$patient$days_to_last_followup)==1,
                                  as.numeric(xml_obj$tcga_bcr$patient$days_to_last_followup[[1]]),
                                  NA)
  days_to_death <- ifelse(length(xml_obj$tcga_bcr$patient$days_to_death),
                          as.numeric(xml_obj$tcga_bcr$patient$days_to_death[[1]]),
                          NA)
  drug_nodes <- xml_obj$tcga_bcr$patient$drugs
  if(length(drug_nodes)>0){
    names(drug_nodes) <- paste0(names(drug_nodes), 1:length(names(drug_nodes)))
    last_drug_treatment_day <- unlist(lapply(drug_nodes, function(x) as.numeric(unlist(x$days_to_drug_therapy_end))))
    if(length(last_drug_treatment_day)>0){
      last_drug_treatment_day <- max(last_drug_treatment_day)
    }else{
      last_drug_treatment_day <- NA
    }
  }else{
    last_drug_treatment_day <- NA
  }
  last_date_vec <- c(days_to_last_followup, days_to_death, last_drug_treatment_day)
  names(last_date_vec) <- c('Last followup', 'Death', 'End of last treatment')
  if(all(is.na(last_date_vec))){
    os <- NA
    names(os) <- 'No event'
  }else{
    os <- max(last_date_vec, na.rm=T)
    names(os) <- sort(paste(names(which(last_date_vec==os)), collapse='/'))
    if(os==0){
      os <- NA
      names(os) <- 'No event'
    }
  }
  os
}

parse_drugs <- function(xml_obj){
  drugs_nodeset <- xml_obj$tcga_bcr$patient$drugs
  barcode <- xml_obj$tcga_bcr$patient$bcr_patient_barcode[[1]]
  drug_df <- data.frame(bcr_patient_barcode=character(),
                        drug_name=character(),
                        days_to_start=numeric(),
                        days_to_end=numeric(),
                        cycles=numeric())
  drug_names <- read_csv(file.path(input_dir, 'tcga_standardised_drug_names.csv'))
  if(length(drugs_nodeset)){
    names(drugs_nodeset) <- paste0(names(drugs_nodeset), 1:length(names(drugs_nodeset)))
    for(i in 1:length(drugs_nodeset)){
      drug_name <- tolower(ifelse(length(drugs_nodeset[[i]]$drug_name),
                                  drugs_nodeset[[i]]$drug_name[[1]],
                                  ""))
      drug_name <- gsub(' ', '', drug_name)
      drug_name <- unlist(strsplit(drug_name,'[/&,+]'))
      if(length(drug_name)){
        for(j in 1:length(drug_name)){
          if(drug_name[j] %in% drug_names$input_name){
            drug_name[j] <- drug_names$standardised[drug_names$input_name==drug_name[j]]
          }
        }
        days_to_start <- ifelse(length(drugs_nodeset[[i]]$days_to_drug_therapy_start)==1,
                                as.numeric(drugs_nodeset[[i]]$days_to_drug_therapy_start[[1]]),
                                NA)
        days_to_end <- ifelse(length(drugs_nodeset[[i]]$days_to_drug_therapy_end)==1,
                              as.numeric(drugs_nodeset[[i]]$days_to_drug_therapy_end[[1]]),
                              NA)
        cycles <- ifelse(length(drugs_nodeset[[i]]$number_cycles)==1,
                         as.numeric(drugs_nodeset[[i]]$number_cycles[[1]]),
                         NA)
        temp_df <- data.frame(bcr_patient_barcode=barcode,
                              drug_name=drug_name,
                              days_to_start=days_to_start,
                              days_to_end=days_to_end,
                              cycles=cycles)
        drug_df <- rbind(drug_df, temp_df)
      }
    }
    drug_df <- drug_df %>% arrange(days_to_start)
    as.data.frame(collapse_treatments(drug_df))
  }else{
    drug_df
  }
}

collapse_treatments <- function(drug_df){
  drug_df <- drug_df %>%
    group_by(bcr_patient_barcode, days_to_start, days_to_end) %>%
    summarise(drug_name=paste(sort(unique(drug_name)),collapse='/'),
              cycles=min(cycles))
  drug_df$treatment_line <- rank(drug_df$days_to_start)
  drug_df
}

pull_days_to_new_tumour_event <- function(xml_obj){
  tumour_event_df <- data.frame(bcr_patient_barcode=character(),
                                days_to_new_tumour_event=numeric())
  barcode <- xml_obj$tcga_bcr$patient$bcr_patient_barcode[[1]]
  new_tumour_events <- unlist(xml_obj$tcga_bcr$patient$follow_ups$follow_up$new_tumor_events)
  new_tumour_events_days_index <- sort(which(names(new_tumour_events)=='new_tumor_event.days_to_new_tumor_event_after_initial_treatment'))
  new_tumour_events2 <- unlist(xml_obj$tcga_bcr$patient$new_tumor_events)
  new_tumour_events_days_index2 <- sort(which(names(new_tumour_events2)=='new_tumor_event.days_to_new_tumor_event_after_initial_treatment'))
  new_tumour_events_days <- as.numeric(c(new_tumour_events[new_tumour_events_days_index], new_tumour_events2[new_tumour_events_days_index2]))
  if(length(new_tumour_events_days)>0){
    tumour_event_df <- rbind(tumour_event_df,
                             data.frame(bcr_patient_barcode=barcode,
                                        days_to_new_tumour_event=new_tumour_events_days))
    tumour_event_df <- tumour_event_df[!duplicated(tumour_event_df),] %>%
      arrange(days_to_new_tumour_event)
  }
  return(tumour_event_df)
}

calc_TTF_TCGA <- function(results_df, nte_df){
  results_df <- results_df %>%
    mutate(TTF=lead(days_to_start) - days_to_start) %>%
    mutate(progression_reason=ifelse(is.na(TTF), NA, 'Next treatment line'))
  nte_df <- nte_df %>% dplyr::filter(days_to_new_tumour_event>results_df$days_to_start[1])
  if(nrow(results_df)>0){
    if(nrow(nte_df)>0){
      for(i in 1:nrow(results_df)){
        days_to_start <- min(c(results_df$days_to_start[i],
                               results_df$days_to_end[i]), na.rm=TRUE)
        nte_after_start <- days_to_start < nte_df$days_to_new_tumour_event
        nte_after_start[is.na(nte_after_start)] <- FALSE
        nte_before_next_start <- results_df$days_to_start[i+1] > nte_df$days_to_new_tumour_event
        nte_comp <- apply(matrix(c(nte_after_start, nte_before_next_start), ncol=2),
                          MARGIN = 1,
                          FUN = xfun)
        if(any(nte_comp, na.rm=T)){
          nte_df_ind <- min(which(nte_comp))
          results_df$TTF[i] <- nte_df$days_to_new_tumour_event[nte_df_ind] - results_df$days_to_start[i]
          results_df$progression_reason[i] <- 'New tumour event'
        }
      }
    }
    if(any(is.na(results_df$TTF))){
      na_ind <- min(which(is.na(results_df$TTF)))
      last_day_vec <- c(results_df$overall_survival[na_ind], results_df$days_to_last_followup[na_ind])
      names(last_day_vec) <- c(results_df$overall_survival_reason[na_ind], 'Last followup')
      last_day_vec <- last_day_vec[last_day_vec>results_df$days_to_start[na_ind]]
      last_day_vec <- last_day_vec[!is.na(last_day_vec)]
      if(length(last_day_vec)>0){
        results_df$TTF[na_ind] <- min(last_day_vec) - results_df$days_to_start[na_ind]
        results_df$progression_reason[na_ind] <- paste(sort(unique(names(last_day_vec[which(last_day_vec==min(last_day_vec))]))), collapse='/')
      }else{
        results_df$TTF[na_ind] <- NA
      }
    } 
  }
  results_df <- results_df %>%
    mutate(TTF_months=TTF / (365/12),
           censoring=ifelse(grepl('Last followup', progression_reason) | grepl('End of last treatment', progression_reason) | is.na(progression_reason),
                            0,
                            1))
  results_df
}

xfun <- function(x){
  if(sum(is.na(x))==2){
    FALSE
  }else if(sum(is.na(x))==1){
    x[which(!is.na(x))]
  }else{
    x[1] & x[2]
  }
}

join_overlapping <- function(results_df){
  processed <- FALSE
  i <- 1
  if(nrow(results_df)>0){
    while(!processed){
      if(i >= nrow(results_df)){
        processed <- TRUE
      }else{
        # Check if the treatment times are within adjuvant time distance
        if(!is.na(results_df$days_to_next_treatment[i]) & results_df$days_to_next_treatment[i]<0){
          results_df <- join_rows(results_df, i)
        }else{
          i <- i+1
        }
      }
    }
    results_df$treatment_line <- 1:nrow(results_df)
    rownames(results_df) <- 1:nrow(results_df)
  }
  results_df
}

join_rows <- function(results_df, i){
  replace_df <- results_df[i,]
  replace_df$days_to_end <- results_df$days_to_end[i+1]
  replace_df$drug_name <- paste(sort(unique(unlist(strsplit(results_df$drug_name[i:(i+1)], '/')))), collapse='/')
  replace_df$days_to_next_treatment <- results_df$days_to_next_treatment[(i+1)]
  if(i>1){
    replace_df <- rbind(results_df[1:(i-1),],
                        replace_df)
  }
  if(nrow(results_df)>(i+2)){
    replace_df <- rbind(replace_df,
                        results_df[(i+2):nrow(results_df),])
  }
  replace_df
}

##### Survival Analysis Functions #####

.get_data <- function(fit, data = NULL, complain = FALSE){
  if(is.null(data)){
    if (complain)
      warning ("The `data` argument is not provided. Data will be extracted from model fit.")
    data <- eval(fit$call$data)
    if (is.null(data))
      stop("The `data` argument should be provided either to ggsurvfit or survfit.")
  }
  data
}

ggforest_altered <- function(model, data = NULL, main = "Hazard ratio", cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2){
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- .get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(tidy(model, conf.int = TRUE))
  gmodel <- glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    var <- gsub('`','', var)
    if(startsWith(var, 'strata')){
      cbind(var = "strata", Var1 = "", Freq = "", pos = 1)
    }
    else if (terms[i] %in% c("factor", "character") & !startsWith(var, 'strata')) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric" & var!="1st line TTF") {
      data.frame(var = var, Var1 = "", Freq = nrow(data), 
                 pos = 1)
    }
    else if (var == "1st line TTF"){
      data.frame(cbind(var = c("1st line TTF","Treatment armExperimental:1st line TTF"),
                       Var1 = c("",""),
                       Freq = c(table(data$`Treatment arm`)),
                       pos = c(1:2)))
    }
    else {
      vars = grep(paste0("^", var, "*."), coef$term, value = TRUE)
      data.frame(var = vars, Var1 = "", Freq = nrow(data), 
                 pos = seq_along(vars))
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  allTermsDF <- allTermsDF[allTermsDF$var!="strata",]
  allTermsDF <- allTermsDF[allTermsDF$var!="(weights)",]
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  allTermsDF$term <- inds
  coef$term <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- left_join(allTermsDF,coef,by="term")
  toShow <- toShow[, c("var", "level","N", "p.value","estimate", "conf.low","conf.high", "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 4],
                                                       noDigits + 1),
                               toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, noDigits + 1),
                                 " ",
                                 ifelse(toShowExpClean$p.value < 0.05, "*", ""),
                                 ifelse(toShowExpClean$p.value < 0.01, "*", ""),
                                 ifelse(toShowExpClean$p.value < 0.001, "*", ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], 
                              " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, ]
  if(any(toShowExpClean$conf.high.1==' Inf')){
    temp_index <- which(toShowExpClean$conf.high.1==' Inf')
    toShowExpClean$conf.high[temp_index] <- NA
    toShowExpClean$conf.low[temp_index] <- NA
  }
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, 
                  na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(convertX(unit(theme_get()$text$size, 
                                                       "pt"), "mm"))
  annot_size_mm2 <- (fontsize+2) * as.numeric(convertX(unit(theme_get()$text$size, 
                                                            "pt"), "mm"))
  annot_size_mm3 <- (fontsize+1) * as.numeric(convertX(unit(theme_get()$text$size, 
                                                            "pt"), "mm"))
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) + 
    geom_rect(aes(xmin = seq_along(var) - 0.5,
                  xmax = seq_along(var) + 0.5,
                  ymin = exp(rangeplot[1]),
                  ymax = exp(rangeplot[2]),
                  fill = ordered(seq_along(var) %% 2 + 1))) +
    scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
    geom_point(pch = 15, size = 3) +
    geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), width = 0.15) +
    geom_hline(yintercept = 1, linetype = 3) + 
    coord_flip(ylim = exp(rangeplot)) +
    ggtitle(main) +
    scale_y_log10(name = "", labels = sprintf("%g", breaks), expand = c(0.02, 0.02), breaks = breaks) +
    theme_light() +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none",
          panel.border = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          title=element_text(size=annot_size_mm2),
          axis.text.x=element_text(size=annot_size_mm3)) + 
    xlab("") +
    annotate(geom = "text", x = x_annotate, y = exp(y_variable), label = toShowExpClean$var, fontface = "bold", hjust = 0, size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = exp(y_nlevel), hjust = 0, label = toShowExpClean$level, vjust = -0.1, size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = exp(y_nlevel), label = toShowExpClean$N, fontface = "italic", hjust = 0, vjust = ifelse(toShowExpClean$level == "", 0.5, 1.1), size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring), label = toShowExpClean$estimate.1, size = annot_size_mm, vjust = ifelse(toShowExpClean$estimate.1 == "reference", 0.5, -0.1)) +
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring), label = toShowExpClean$ci, size = annot_size_mm, vjust = 1.1, fontface = "italic") + 
    annotate(geom = "text", x = x_annotate, y = exp(y_stars),  label = toShowExpClean$stars, size = annot_size_mm, hjust = -0.2, fontface = "italic") +
    annotate(geom = "text", x = 0.5, y = exp(y_variable), label = paste0("# Events: ", gmodel$nevent, "; Global p-value (Log-Rank): ", format.pval(gmodel$p.value.log, eps = ".001"), " \nAIC: ", round(gmodel$AIC, 2), "; Concordance Index: ", round(gmodel$concordance, 2)), size = annot_size_mm, hjust = 0, vjust = 1.2, fontface = "italic")
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}

##### Cell Line Functions #####

make_srb_long <- function(raw_srb_data, time){
  control_cols <- 2:4
  dose0.1_cols <- 5:7
  dose0.5_cols <- 8:10
  dose1_cols <- 11:13
  CIOV2_rows <- 1:2
  CIOV4_rows <- 3:4
  CIOV6_rows <- 5:6
  PEO23_rows <- 7:8
  CIOV1_rows <- 9:10
  OVCAR3_rows <- 11:12
  srb_data_long <- rbind(
    data.frame(cell_type = 'OVCAR3', time, rep = 1:6, dose = 0, value = as.numeric(as.matrix(raw_srb_data[OVCAR3_rows, control_cols]))),
    data.frame(cell_type = 'CIOV1', time, rep = 1:6, dose = 0, value = as.numeric(as.matrix(raw_srb_data[CIOV1_rows, control_cols]))),
    data.frame(cell_type = 'PEO23', time, rep = 1:6, dose = 0, value = as.numeric(as.matrix(raw_srb_data[PEO23_rows, control_cols]))),
    data.frame(cell_type = 'CIOV4', time, rep = 1:6, dose = 0, value = as.numeric(as.matrix(raw_srb_data[CIOV4_rows, control_cols]))),
    data.frame(cell_type = 'CIOV6', time, rep = 1:6, dose = 0, value = as.numeric(as.matrix(raw_srb_data[CIOV6_rows, control_cols]))),
    data.frame(cell_type = 'CIOV2', time, rep = 1:6, dose = 0, value = as.numeric(as.matrix(raw_srb_data[CIOV2_rows, control_cols]))),
    data.frame(cell_type = 'OVCAR3', time, rep = 1:6, dose = 0.1, value = as.numeric(as.matrix(raw_srb_data[OVCAR3_rows, dose0.1_cols]))),
    data.frame(cell_type = 'CIOV1', time, rep = 1:6, dose = 0.1, value = as.numeric(as.matrix(raw_srb_data[CIOV1_rows, dose0.1_cols]))),
    data.frame(cell_type = 'PEO23', time, rep = 1:6, dose = 0.1, value = as.numeric(as.matrix(raw_srb_data[PEO23_rows, dose0.1_cols]))),
    data.frame(cell_type = 'CIOV4', time, rep = 1:6, dose = 0.1, value = as.numeric(as.matrix(raw_srb_data[CIOV4_rows, dose0.1_cols]))),
    data.frame(cell_type = 'CIOV6', time, rep = 1:6, dose = 0.1, value = as.numeric(as.matrix(raw_srb_data[CIOV6_rows, dose0.1_cols]))),
    data.frame(cell_type = 'CIOV2', time, rep = 1:6, dose = 0.1, value = as.numeric(as.matrix(raw_srb_data[CIOV2_rows, dose0.1_cols]))),
    data.frame(cell_type = 'OVCAR3', time, rep = 1:6, dose = 0.5, value = as.numeric(as.matrix(raw_srb_data[OVCAR3_rows, dose0.5_cols]))),
    data.frame(cell_type = 'CIOV1', time, rep = 1:6, dose = 0.5, value = as.numeric(as.matrix(raw_srb_data[CIOV1_rows, dose0.5_cols]))),
    data.frame(cell_type = 'PEO23', time, rep = 1:6, dose = 0.5, value = as.numeric(as.matrix(raw_srb_data[PEO23_rows, dose0.5_cols]))),
    data.frame(cell_type = 'CIOV4', time, rep = 1:6, dose = 0.5, value = as.numeric(as.matrix(raw_srb_data[CIOV4_rows, dose0.5_cols]))),
    data.frame(cell_type = 'CIOV6', time, rep = 1:6, dose = 0.5, value = as.numeric(as.matrix(raw_srb_data[CIOV6_rows, dose0.5_cols]))),
    data.frame(cell_type = 'CIOV2', time, rep = 1:6, dose = 0.5, value = as.numeric(as.matrix(raw_srb_data[CIOV2_rows, dose0.5_cols]))),
    data.frame(cell_type = 'OVCAR3', time, rep = 1:6, dose = 1, value = as.numeric(as.matrix(raw_srb_data[OVCAR3_rows, dose1_cols]))),
    data.frame(cell_type = 'CIOV1', time, rep = 1:6, dose = 1, value = as.numeric(as.matrix(raw_srb_data[CIOV1_rows, dose1_cols]))),
    data.frame(cell_type = 'PEO23', time, rep = 1:6, dose = 1, value = as.numeric(as.matrix(raw_srb_data[PEO23_rows, dose1_cols]))),
    data.frame(cell_type = 'CIOV4', time, rep = 1:6, dose = 1, value = as.numeric(as.matrix(raw_srb_data[CIOV4_rows, dose1_cols]))),
    data.frame(cell_type = 'CIOV6', time, rep = 1:6, dose = 1, value = as.numeric(as.matrix(raw_srb_data[CIOV6_rows, dose1_cols]))),
    data.frame(cell_type = 'CIOV2', time, rep = 1:6, dose = 1, value = as.numeric(as.matrix(raw_srb_data[CIOV2_rows, dose1_cols])))
  )
  srb_data_long
}

make_incucyte_long <- function(raw_incucyte_data, time, time_row){
  control_cols <- 1:3
  dose0.025_cols <- 4:6
  dose0.05_cols <- 7:9
  dose1_cols <- 10:12
  CIOV2_lab <- 'B'
  CIOV4_lab <- 'D'
  CIOV6_lab <- 'E'
  PEO23_lab <- 'G'
  CIOV1_lab <- 'A'
  OVCAR3_lab <- 'F'
  
  d1 <- do.call('rbind', mapply(function(x, y, z){
    data.frame(cell_type = x, time, rep = 1:3, dose = 0, 
               value = as.numeric(as.matrix(raw_incucyte_data[time_row, paste0(y, control_cols)])))},
    c('OVCAR3', 'CIOV2', 'PEO23', 'CIOV4', 'CIOV6', 'CIOV1'),
    c(OVCAR3_lab, CIOV2_lab, PEO23_lab, CIOV4_lab, CIOV6_lab, CIOV1_lab),
    SIMPLIFY = F))
  
  d2 <- do.call('rbind', mapply(function(x, y, z){
    data.frame(cell_type = x, time, rep = 1:3, dose = 0.025,
               value = as.numeric(as.matrix(raw_incucyte_data[time_row, paste0(y, dose0.025_cols)])))},
    c('OVCAR3', 'CIOV2', 'PEO23', 'CIOV4', 'CIOV6', 'CIOV1'),
    c(OVCAR3_lab, CIOV2_lab, PEO23_lab, CIOV4_lab, CIOV6_lab, CIOV1_lab),
    SIMPLIFY = F))
  
  d3 <- do.call('rbind', mapply(function(x, y, z){
    data.frame(cell_type = x, time, rep = 1:3, dose = 0.05,
               value = as.numeric(as.matrix(raw_incucyte_data[time_row, paste0(y, dose0.05_cols)])))},
    c('OVCAR3', 'CIOV2', 'PEO23', 'CIOV4', 'CIOV6', 'CIOV1'),
    c(OVCAR3_lab, CIOV2_lab, PEO23_lab, CIOV4_lab, CIOV6_lab, CIOV1_lab),
    SIMPLIFY = F))
  
  d4 <- do.call('rbind', mapply(function(x, y, z){
    data.frame(cell_type = x, time, rep = 1:3, dose = 0.1,
               value = as.numeric(as.matrix(raw_incucyte_data[time_row, paste0(y, dose1_cols)])))},
    c('OVCAR3', 'CIOV2', 'PEO23', 'CIOV4', 'CIOV6', 'CIOV1'),
    c(OVCAR3_lab, CIOV2_lab, PEO23_lab, CIOV4_lab, CIOV6_lab, CIOV1_lab),
    SIMPLIFY = F))
  
  rbind(d1, d2, d3, d4)
}

compute_exp_mn_first_division <- function(f, MN0, d){
  MNt <- (1 - f) * MN0 + d * f * (MN0 * 0.5 + MN0 + (1 - MN0) / 2)
  MNt
}

adjust_exp_mn_subsequent_divisions <- function(f, MNt, MN0){
  MNt <- (1 - f) * MNt + 0.75 * f * MNt / 2 #it is missing background prop of MN formation
  MNt
}

cal_sens <- function(variable, response){
  threshold <- max(variable[response])
  sens <- sum(variable[!response] > threshold) / sum(!response)
  sens
}

get_thresh <- function(variable, response){
  threshold <- max(variable[response])
  threshold
}

get_sens_spec <- function(variable, thresh){
  response <- sph_resp$response == 'sensitive'
  sens <- sum(variable[!response] > thresh) / sum(!response)
  spec <- sum(variable[response] <= thresh) / sum(response)
  c(spec, sens)
}

##### HMF Parsing Functions #####

calc_TTF_HMF <- function(temp_drugs, temp_responses){
  # all treatment lines except the last have ttf as difference between start date and next treatment's start date
  temp_drugs <- temp_drugs %>%
    ungroup() %>%
    arrange(startDate) %>%
    mutate(TTF=lead(startDate) - startDate)
  temp_responses <- temp_responses %>%
    filter(responseDate>temp_drugs$startDate[1])
  if(nrow(temp_drugs)>0){
    temp_drugs$censoring <- FALSE
    temp_drugs$censoring[nrow(temp_drugs)] <- TRUE
    temp_drugs$reason <- 'Next treatment line'
    deathDate <- NA
    temp_responses <- temp_responses %>%
      filter(responseDate>temp_drugs$startDate[nrow(temp_drugs)])
    if(!is_empty(metadata$deathDate[metadata$patientIdentifier==temp_drugs$patientIdentifier[1]]) & sum(temp_responses$response=='PD')==0){
      if(unique(metadata$deathDate[metadata$patientIdentifier==temp_drugs$patientIdentifier[1]])!='null'){
        deathDate <- unique(as.Date(metadata$deathDate[metadata$patientIdentifier==temp_drugs$patientIdentifier[1]]))
        temp_drugs$TTF[nrow(temp_drugs)] <- deathDate - temp_drugs$startDate[nrow(temp_drugs)]
        temp_drugs$reason[nrow(temp_drugs)] <- 'Death'
      }
    }
    if(is.na(deathDate) & sum(temp_responses$response=='PD')==0 & sum(temp_responses$response!='PD')>0){
      temp_drugs$TTF[nrow(temp_drugs)] <- max(temp_responses$responseDate, na.rm=TRUE) - temp_drugs$startDate[nrow(temp_drugs)]
      temp_drugs$reason[nrow(temp_drugs)] <- 'RECIST'
    }
  }
  if(is.na(temp_drugs$TTF[nrow(temp_drugs)])){
    temp_drugs$TTF[nrow(temp_drugs)] <- temp_drugs$endDate[nrow(temp_drugs)] - temp_drugs$startDate[nrow(temp_drugs)]
    temp_drugs$reason[nrow(temp_drugs)] <- 'Length of treatment'
  }
  temp_drugs
}

assign_treatment_lines <- function(in_df){
  in_df <- in_df %>% ungroup()
  in_df$endDate <- as.Date(in_df$endDate)
  treatment_line <- numeric(length=nrow(in_df))
  treatment_line[1] <- 1
  if(nrow(in_df)>1){
    for(i in 2:nrow(in_df)){
      # First check for the case where the previous line of treatment only lasted 1 day
      if(in_df$startDate[i-1]==in_df$endDate[i-1]){
        if(in_df$startDate[i]==in_df$startDate[i-1]){
          #Because the lines are sorted by startdate, if the startdates are equal the lines can be joined
          in_df$endDate[c(i-1,i)] <- max(in_df$endDate[c(i-1,i)]) 
        }
      }
      # Now check for the cases where the previous line of treatment lasted multiple days
      # In these cases, the newer line startDate has to be >= previous startdate, but < previous enddate
      # If the new treatment starts on the day the old treatment ends they won't be joined
      else{
        if(in_df$startDate[i]>=in_df$startDate[i-1] & in_df$startDate[i]<in_df$endDate[i-1]){
          in_df$startDate[c(i-1,i)] <- min(in_df$startDate[c(i-1,i)])
          in_df$endDate[c(i-1,i)] <- max(in_df$endDate[c(i-1,i)])
        }
      }
    }
    in_df <- in_df %>% mutate(dateDiff=startDate-lag(endDate))
    # Again there is a special case when the previous treatment line lasted 1 day
    # In these cases we accept a dateDiff==0 for joining.
    # Otherwise we accept a dateDiff<0, not including 0
    for(i in 2:length(treatment_line)){
      if(in_df$startDate[i-1]==in_df$endDate[i-1]){
        if(in_df$dateDiff[i]==0){
          treatment_line[i] <- treatment_line[i-1]
        }else{
          treatment_line[i] <- treatment_line[i-1] + 1
        }
      }else{
        if(in_df$dateDiff[i]<0){
          treatment_line[i] <- treatment_line[i-1]
        }else{
          treatment_line[i] <- treatment_line[i-1] + 1
        }
      }
    }
  }else{
    in_df$dateDiff <- NA
  }
  in_df$treatment_line <- treatment_line
  return(in_df)
}

plot_lines <- function(patientID){
  pre_biopsy_temp <- pre_biopsy_drugs %>% filter(patientIdentifier==patientID)
  post_biopsy_temp <- post_biopsy_drugs %>% filter(patientIdentifier==patientID)
  biopsy_dates <- metadata$biopsyDate[metadata$patientIdentifier==patientID]
  biopsy_dates <- as.Date(biopsy_dates, format='%Y-%m-%d')
  death_date <- as.Date(unique(metadata$deathDate[metadata$patientIdentifier==patientID]), format='%Y-%m-%d')
  response_dates <- treatment_responses[treatment_responses$patientIdentifier==patientID,]
  pre_post_temp <- rbind(pre_biopsy_temp[,2:4], post_biopsy_temp[,3:5])
  pre_post_temp <- pre_post_temp %>%
    filter(endDate!='null') %>%
    mutate(endDate=as.Date(endDate)) %>%
    filter(endDate>1950 & startDate>1950 & (endDate>=startDate)) %>%
    arrange(startDate) %>%
    mutate(treatment_length=endDate-startDate)
  pre_post_temp <- assign_treatment_lines(pre_post_temp)
  x_range <- c(min(c(pre_post_temp$startDate, biopsy_dates, response_dates$responseDate, death_date), na.rm=TRUE),
               max(c(pre_post_temp$startDate, biopsy_dates, response_dates$responseDate, death_date), na.rm=TRUE))
  p <- pre_post_temp %>%
    mutate(name=factor(paste0(1:length(name), ' - ', name),
                       levels=rev(paste0(1:length(name), ' - ', name)))) %>%
    ggplot(aes(color=as.factor(treatment_line))) +
    xlim(x_range[1], x_range[2]) +
    geom_segment(aes(y=name, yend=name, x=startDate, xend=endDate), linewidth=2) +
    geom_point(aes(y=name, x=startDate)) +
    geom_vline(xintercept=biopsy_dates, linetype=5) +
    geom_vline(data=response_dates, aes(xintercept=responseDate, color=response)) +
    geom_vline(xintercept = death_date, linetype=5, color='red') +
    geom_text(aes(x=startDate, y=name, label=treatment_length), nudge_y=0.25) + 
    labs(color='Treatment line', title=patientID) +
    theme(axis.title.y=element_blank())
  return(p)
}

##### Phase 2 and 3 Analysis Functions #####

restrict_cancer <- function(in_df, cancer){
  return(in_df %>% filter(cancer_type==cancer))
}

immediate_post_biopsy_samples_only <- function(in_df, metadata){
  biopsy_dates <- metadata %>%
    select(patientIdentifier, biopsyDate, sampleId) %>%
    mutate(biopsyDate=as.Date(biopsyDate, format='%Y-%m-%d')) %>%
    filter(!is.na(biopsyDate) & patientIdentifier %in% in_df$patientIdentifier)
  for(i in 1:length(unique(in_df$patientIdentifier))){
    id <- unique(in_df$patientIdentifier)[i]
    temp_df <- in_df[in_df$patientIdentifier==id,]
    bdates <- biopsy_dates[biopsy_dates$patientIdentifier==id,]
    bdates <- bdates %>% arrange(biopsyDate)
    id_rows <- c()
    for(j in 1:nrow(bdates)){
      id_rows <- c(id_rows, which(temp_df$startDate >= bdates$biopsyDate[j])[1])
    }
    temp_df <- cbind(temp_df[id_rows,], bdates[,2:3])
    if(i==1){
      out_df <- temp_df
    }else{
      out_df <- rbind(out_df, temp_df)
    }
  }
  out_df <- out_df[!is.na(out_df$patientIdentifier),] # Why does this happen?
  return(out_df)
}

remove_na_ttf <- function(in_df){
  return(in_df[!is.na(in_df$TTF),])
}

prepare_tcga_ov <- function(){
  tcga_clin <- read_xlsx(file.path(input_dir, 'ds_CCI.17.00096-2_curated.xlsx'), sheet = 'Months')
  tcga_clin <- tcga_clin %>%
    pivot_longer(cols = c(`1st regimen`,
                          `2nd regimen`,
                          `3rd regimen`,
                          `4th regimen`,
                          `5th regimen`,
                          `6th regimen`),
                 names_to = 'treatment_line', values_to = 'drug_name')
  colnames(tcga_clin)[colnames(tcga_clin)=='TTF3rd line (months)'] <- 'TTF 3rd line (months)'
  colnames(tcga_clin)[colnames(tcga_clin)=='TTF  6th line (months)'] <- 'TTF 6th line (months)'
  TTF <- numeric(length = nrow(tcga_clin))
  censoring <- numeric(length=nrow(tcga_clin))
  for(i in 1:nrow(tcga_clin)){
    t_string <- tcga_clin$treatment_line[i]
    t_col <- which(grepl(paste0('TTF ', strsplit(t_string, ' ')[[1]][1], ' line'),
                         colnames(tcga_clin), fixed = T))[1]
    TTF[i] <- as.numeric(tcga_clin[i, t_col])
    censoring[i] <- tcga_clin[[paste0(gsub(' ', '_chemo_', t_string), '_days_outcome')]][i]
    
  }
  tcga_clin$TTF <- TTF
  tcga_clin$censoring <- censoring
  tcga_clin$`TTF 1st line (months)` <- as.numeric(tcga_clin$`TTF 1st line (months)`)
  return(tcga_clin)
}

remove_herceptin <- function(in_df, cohort){
  if(cohort=='hmf'){
    anti_her2 <- c('Trastuzumab', 'Pertuzumab',
                   'Ado-trastuzumab emtansine (T-DM1)', 'Lapatinib')
    herceptin_patients <- unique(in_df %>%
                                   filter(grepl(paste(anti_her2, collapse='|'),
                                                name)) %>%
                                   pull(patientIdentifier))
    in_df <- in_df %>% filter(!(patientIdentifier %in% herceptin_patients)) 
  }else if(cohort=='tcga'){
    herceptin_patients <- unique(in_df %>%
                                   filter(grepl('her2_inhib', drug_name)) %>%
                                   pull(bcr_patient_barcode))
    in_df <- in_df %>% filter(!(bcr_patient_barcode %in% herceptin_patients)) 
  }
  return(in_df)
}

remove_stage1 <- function(in_df){
  in_df <- in_df %>% filter(grepl('I{2,}|IV|X', path_stage, perl=TRUE) | is.na(path_stage))
  in_df <- in_df %>% filter(grepl('I{2,}|IV|X', clin_stage, perl=TRUE) | is.na(clin_stage))
  return(in_df)
}

predict_chemo_resistance <- function(in_df, chemo, cohort){
  tcga_sigs <- as.data.frame(readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_TCGA.rds')))
  tcga_sigs$bcr_patient_barcode <- rownames(tcga_sigs)
  if(cohort=='hmf'){
    hmf_sigs <- as.data.frame(readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_HMF.rds'))) #need to be computed from copy number profiles => use the script in /Analysis/Helper_Scripts/signature_quantification.R
    hmf_sigs$sample <- rownames(hmf_sigs)
    hmf_sigs$patientIdentifier <- gsub('TI*V*$', '', hmf_sigs$sample, perl=TRUE)
    sigs_to_predict <- hmf_sigs %>%
      filter(sample %in% in_df$sampleId)
  }else if(grepl('tcga', cohort)){
    sigs_to_predict <- tcga_sigs %>%
      filter(bcr_patient_barcode %in% in_df$bcr_patient_barcode)
  }
  if(chemo=='platinum'){
    brca_status <- readRDS(file.path(input_dir, 'TCGA_BRCA_MutStatus.rds'))
    brca1 = unique(as.character(brca_status$Sample[ grepl("[germline|somatic] BRCA1", 
                                                          brca_status$Status) ]))
    brca2 = unique(as.character(brca_status$Sample[ grepl("[germline|somatic] BRCA2", 
                                                          brca_status$Status) ]))
    brcaMuts = unique(c(brca1,brca2))
    mBRCA = tcga_sigs[row.names(tcga_sigs) %in% brcaMuts, ]
    mBRCA = as.matrix(mBRCA[ , c("CX3", "CX2")])
    plat_scaling_model = list(mean = attributes(scale(mBRCA))$`scaled:center`, 
                              scale = attributes(scale(mBRCA))$`scaled:scale`)
    sigs_to_predict <- sigs_to_predict %>%
      mutate(sCX2 = (CX2 - plat_scaling_model$mean[['CX2']]) / plat_scaling_model$scale[['CX2']],
             sCX3 = (CX3 - plat_scaling_model$mean[['CX3']]) / plat_scaling_model$scale[['CX3']]) %>%
      mutate(prediction = ifelse(sCX3>sCX2, 'Sensitive', 'Resistant'))
    if(cohort=='hmf'){
      sigs_to_predict <- sigs_to_predict %>%
        select(sampleId=sample, sCX2, sCX3, prediction)
      nocin_samps <- readRDS(file.path(input_dir, 'hmf_all_segs_nocin_samps.rds')) #download copy number profiles from HMF & apply dCIN threshold
      nocin_samps <- unique(nocin_samps$sample[nocin_samps$sample %in% in_df$sampleId])
      if(length(nocin_samps)>0){
        nocin_sigs <- as.data.frame(cbind(nocin_samps, NA, NA, 'Resistant'))
        colnames(nocin_sigs) <- c('sampleId', 'sCX2', 'sCX3', 'prediction')
        sigs_to_predict <- rbind(sigs_to_predict, nocin_sigs)
      }
    }else if(grepl('tcga', cohort)){
      sigs_to_predict <- sigs_to_predict %>%
        select(bcr_patient_barcode, sCX2, sCX3, prediction)
      metadata_cin_tcga <- readRDS(file.path(input_dir, 'Metadata_ALL_TCGA.rds'))
      nocin_samps <- metadata_cin_tcga$sample[metadata_cin_tcga$QC=="YES" & metadata_cin_tcga$CIN=="NO"]
      nocin_samps <- nocin_samps[nocin_samps %in% in_df$bcr_patient_barcode]
      nocin_sigs <- as.data.frame(cbind(nocin_samps, NA, NA, 'Resistant'))
      colnames(nocin_sigs) <- c('bcr_patient_barcode', 'sCX2', 'sCX3', 'prediction')
      sigs_to_predict <- rbind(sigs_to_predict, nocin_sigs)
    }
    in_df <- left_join(in_df, sigs_to_predict) %>% filter(!is.na(prediction))
  }
  if(chemo=='taxane'){
    tax_scaling_sigs <- as.matrix(tcga_sigs[, c('CX3', 'CX5')])
    tax_scaling_model <- list(mean = attributes(scale(tax_scaling_sigs))$`scaled:center`, 
                              scale = attributes(scale(tax_scaling_sigs))$`scaled:scale`)
    sigs_to_predict <- sigs_to_predict %>%
      mutate(sCX5 = (CX5 - tax_scaling_model$mean[['CX5']]) / tax_scaling_model$scale[['CX5']],
             sCX3 = (CX3 - tax_scaling_model$mean[['CX3']]) / tax_scaling_model$scale[['CX3']]) %>%
      mutate(prediction = ifelse(sCX5<0, 'Resistant', 'Sensitive'))
    if(cohort=='hmf'){
      sigs_to_predict <- sigs_to_predict %>%
        select(sampleId=sample, CX3, CX5, sCX3, sCX5, prediction)
      nocin_samps <- readRDS(file.path(input_dir, 'hmf_all_segs_nocin_samps.rds')) #download copy number profiles from HMF & apply dCIN threshold
      nocin_samps <- unique(nocin_samps$sample[nocin_samps$sample %in% in_df$sampleId])
      if(length(nocin_samps)>0){
        nocin_sigs <- as.data.frame(cbind(nocin_samps, NA, NA, NA, NA, 'Resistant'))
        colnames(nocin_sigs) <- c('sampleId', 'CX3', 'CX5', 'sCX3', 'sCX5', 'prediction')
        sigs_to_predict <- rbind(sigs_to_predict, nocin_sigs)
      }
    }else if(grepl('tcga', cohort=='tcga')){
      sigs_to_predict <- sigs_to_predict %>%
        select(bcr_patient_barcode, CX3, CX5, sCX3, sCX5, prediction)
      metadata_cin_tcga <- readRDS(file.path(input_dir, 'Metadata_ALL_TCGA.rds'))
      nocin_samps <- metadata_cin_tcga$sample[metadata_cin_tcga$QC=="YES" & metadata_cin_tcga$CIN=="NO"]
      nocin_samps <- nocin_samps[nocin_samps %in% in_df$bcr_patient_barcode]
      nocin_sigs <- as.data.frame(cbind(nocin_samps, NA, NA, NA, NA, 'Resistant'))
      colnames(nocin_sigs) <- c('bcr_patient_barcode', 'CX3', 'CX5', 'sCX3', 'sCX5', 'prediction')
      sigs_to_predict <- rbind(sigs_to_predict, nocin_sigs)
    }
    in_df <- left_join(in_df, sigs_to_predict)
  }
  if(chemo=='doxorubicin'){
    sigs_to_predict <- sigs_to_predict %>%
      mutate(prediction = ifelse(CX8 <= 0.01 & CX9 <= 0.009 & CX13 <= 0.009,
                                 'Sensitive', 'Resistant'))
    if(cohort=='hmf'){
      sigs_to_predict <- sigs_to_predict %>%
        select(sampleId=sample, CX8, CX9, CX13, prediction)
      nocin_samps <- readRDS(file.path(input_dir, 'hmf_all_segs_nocin_samps.rds')) #download copy number profiles from HMF & apply dCIN threshold
      nocin_samps <- unique(nocin_samps$sample[nocin_samps$sample %in% in_df$sampleId])
      if(length(nocin_samps)>0){
        nocin_sigs <- as.data.frame(cbind(nocin_samps, NA, NA, NA, 'Sensitive'))
        colnames(nocin_sigs) <- c('sampleId', 'CX8', 'CX9', 'CX13', 'prediction')
        sigs_to_predict <- rbind(sigs_to_predict, nocin_sigs)
      }
    }else if(grepl('tcga', cohort)){
      sigs_to_predict <- sigs_to_predict %>%
        select(bcr_patient_barcode, CX8, CX9, CX13, prediction)
      metadata_cin_tcga <- readRDS(file.path(input_dir, 'Metadata_ALL_TCGA.rds'))
      nocin_samps <- metadata_cin_tcga$sample[metadata_cin_tcga$QC=="YES" & metadata_cin_tcga$CIN=="NO"]
      nocin_samps <- nocin_samps[nocin_samps %in% in_df$bcr_patient_barcode]
      nocin_sigs <- as.data.frame(cbind(nocin_samps, NA, NA, NA, 'Sensitive'))
      colnames(nocin_sigs) <- c('bcr_patient_barcode', 'CX8', 'CX9', 'CX13', 'prediction')
      sigs_to_predict <- rbind(sigs_to_predict, nocin_sigs)
    }
    in_df <- left_join(in_df, sigs_to_predict)
  }
  return(in_df)
}

identify_sa_patients <- function(in_df, chemo, cohort){
  if(cohort=='hmf'){
    chemo_str <- c('platin', 'taxel', 'rubicin')
    names(chemo_str) <- c('platinum', 'taxane', 'doxorubicin')
    sa_df <- in_df %>% filter(grepl(chemo_str[chemo], name) & !grepl('/', name))
  }else if(grepl('tcga', cohort)){
    # Check that no prior treatment line includes chemo in combination
    if(cohort=='tcga-ov'){
      sa_df <- in_df %>% filter(grepl(chemo, drug_name) & !grepl('/', drug_name) & !grepl('platin', drug_name))
    }else{
      sa_df <- in_df %>% filter(grepl(chemo, drug_name) & !grepl('/', drug_name))
    }
    if(chemo=='taxane' & cohort=='tcga-ov'){
    }else{
      preceded_by_combo <- rep(FALSE, nrow(sa_df))
      for(i in 1:nrow(sa_df)){
        earlier_combo <- in_df %>%
          filter(bcr_patient_barcode==sa_df$bcr_patient_barcode[i] &
                   grepl(chemo, drug_name) &
                   as.numeric(substr(treatment_line, 1, 1))<as.numeric(substr(sa_df$treatment_line[i], 1, 1)))
        if(nrow(earlier_combo)>0){
          preceded_by_combo[i] <- TRUE
        }
      }
      sa_df <- sa_df[!preceded_by_combo,]
    }
  }
  return(sa_df)
}

identify_pre_biopsy_sa_chemo <- function(sa_df, chemo){
  chemo_str <- c('platin', 'taxel', 'rubicin')
  names(chemo_str) <- c('platinum', 'taxane', 'doxorubicin')
  pre_biopsy_chemo_patients <- unique(hmf_clin %>%
                                        filter(patientIdentifier %in% sa_df$patientIdentifier) %>%
                                        filter(biopsyState!='Post' & grepl(chemo_str[chemo], name)) %>%
                                        pull(patientIdentifier))
  sa_df$preceded_by_chemo <- FALSE
  sa_df$preceded_by_chemo[sa_df$patientIdentifier %in% pre_biopsy_chemo_patients] <- TRUE
  return(sa_df)
}

identify_soc_patients <- function(in_df, chemo, cohort, SA_control, control_pattern){
  if(cohort=='hmf'){
    chemo_str <- c('platin', 'taxel', 'rubicin')
    names(chemo_str) <- c('platinum', 'taxane', 'doxorubicin')
    soc_df <- in_df %>%
      filter(!grepl(chemo_str[chemo], name))
    if(!is.null(control_pattern)){
      soc_df <- soc_df %>% filter(grepl(control_pattern, name))
    }
    if(SA_control){
      soc_df <- soc_df %>% filter(!grepl('/', name))
    }
  }else if(grepl('tcga', cohort)){
    if(cohort=='tcga-ov'){
      soc_df <- in_df %>% filter(!grepl(chemo, drug_name) & !grepl('platin', drug_name))
    }else{
      soc_df <- in_df %>% filter(!grepl(chemo, drug_name))
    }
    if(!is.null(control_pattern)){
      soc_df <- soc_df %>% filter(grepl(control_pattern, drug_name))
    }
    if(SA_control){
      soc_df <- soc_df %>% filter(!grepl('/', drug_name))
    }
    if(chemo=='taxane' & cohort=='tcga-ov'){
    }else{
      preceded_by_combo <- rep(FALSE, nrow(soc_df))
      for(i in 1:nrow(soc_df)){
        earlier_combo <- in_df %>%
          filter(bcr_patient_barcode==soc_df$bcr_patient_barcode[i] &
                   grepl(chemo, drug_name) &
                   as.numeric(substr(treatment_line, 1, 1))<as.numeric(substr(soc_df$treatment_line[i], 1, 1)))
        if(nrow(earlier_combo)>0){
          preceded_by_combo[i] <- TRUE
        }
      }
      soc_df <- soc_df[!preceded_by_combo,]
    }
  }
  return(soc_df)
}

filter_short_cycles <- function(in_df, cohort){
  if(cohort=='tcga'){
    in_df <- in_df %>%
      mutate(t_length=days_to_end - days_to_start) %>%
      filter(cycles>=3 | t_length>=28)
  }else if(cohort=='hmf'){
    in_df <- in_df %>%
      mutate(t_length=endDate - startDate) %>%
      filter(t_length>28)
  }else if(cohort=='tcga-ov'){
    in_df <- in_df %>%
      filter(`Cycles of adjuvant therapy`>=3)
  }
  return(in_df)
}

identify_overlapping_patients <- function(sa_df, soc_df, cohort){
  if(cohort=='hmf'){
    overlapping_patients <- unique(intersect(sa_df$patientIdentifier,
                                             soc_df$patientIdentifier))
    soc_df <- soc_df %>% filter(!(patientIdentifier %in% overlapping_patients))
  }else if(grepl('tcga', cohort)){
    overlapping_patients <- unique(intersect(sa_df$bcr_patient_barcode,
                                             soc_df$bcr_patient_barcode))
    soc_df <- soc_df %>% filter(!(bcr_patient_barcode %in% overlapping_patients))
  }
  return(soc_df)
}

identify_overlapping_patients2 <- function(sa_df, soc_df, cohort){
  overlapping_patients <- unique(intersect(sa_df$bcr_patient_barcode,
                                           soc_df$bcr_patient_barcode))
  overlapping_patients_1 <- sample(overlapping_patients, floor(length(overlapping_patients)/2))
  overlapping_patients_2 <- overlapping_patients[!(overlapping_patients %in% overlapping_patients_1)]
  soc_df <- soc_df %>% filter(!(bcr_patient_barcode %in% overlapping_patients_1))
  sa_df <- sa_df %>% filter(!(bcr_patient_barcode %in% overlapping_patients_2))
  return(list(soc_df, sa_df))
}

top_N_SoC <- function(soc_res, soc_sen, cohort, n=5){
  cname <- ifelse(grepl('tcga', cohort), 'drug_name', 'name')
  top_n <- names(sort(table(unlist(strsplit(c(soc_res[[cname]], soc_sen[[cname]]), '/'))), decreasing=TRUE))
  top_n <- top_n[1:n]
  return(top_n)
}

identify_top_SoC <- function(soc_df, top_soc, cohort){
  if(cohort=='hmf'){
    top_soc_indexes <- which(unlist(lapply(strsplit(soc_df$name, '/'),
                                           FUN=function(x) all(x %in% top_soc)))) 
  }else if(grepl('tcga', cohort)){
    top_soc_indexes <- which(unlist(lapply(strsplit(soc_df$drug_name, '/'),
                                           FUN=function(x) all(x %in% top_soc))))
  }
  soc_df <- soc_df[top_soc_indexes,]
  return(soc_df)
}

preceded_by_chemo <- function(in_df, chemo){
  chemo_str <- c('platin', 'taxel', 'rubicin')
  names(chemo_str) <- c('platinum', 'taxane', 'doxorubicin')
  first_chemo_lines <- hmf_clin %>%
    filter(patientIdentifier %in% in_df$patientIdentifier & grepl(chemo_str[chemo], name)) %>%
    group_by(patientIdentifier) %>%
    summarise(first_chemo=min(treatment_line))
  in_df <- left_join(in_df, first_chemo_lines)
  in_df <- in_df %>% mutate(preceded_by_chemo=first_chemo<treatment_line)
  in_df$preceded_by_chemo[is.na(in_df$preceded_by_chemo)] <- FALSE
  in_df <- in_df %>% select(!first_chemo)
  return(in_df)
}

annotate_with_covariates <- function(in_df){
  metadata_df <- metadata %>%
    mutate(dup=duplicated(patientIdentifier),
           age=as.numeric(substr(biopsyDate, 1, 4)) - as.numeric(birthYear)) %>%
    filter(!dup) %>%
    dplyr::select(patientIdentifier, gender, primaryTumorSubType, biopsySite, age)
  return(left_join(in_df, metadata_df))
}

get_year_of_treatment <- function(bcr_barcode, tline){
  tline <- as.numeric(substr(tline, 1, 1))
  tcga_clin <- prepare_tcga_ov()
  years_treated <- floor(sum((tcga_clin %>% filter(bcr_patient_barcode==bcr_barcode & as.numeric(substr(treatment_line, 1, 1)<=tline)) %>% pull(TTF)), na.rm=T) / 12)
  year_of_diag <- unique(tcga_clin %>% filter(bcr_patient_barcode==bcr_barcode) %>% pull(year_of_initial_pathologic_diagnosis))
  return(years_treated + year_of_diag)
}

exclude_ttf_long_followup <- function(in_df, cohort){
  tcga_clin <- readRDS(file.path(input_dir, 'TCGA_clinical_data.RDS'))
  if(cohort == 'tcga-ov'){
    tcga_ov <- prepare_tcga_ov()
    tcga_ov <- tcga_ov %>% filter(!is.na(TTF))
    # Check which patients are using their last treatment line
    final_line <- merge(tcga_ov %>%
                          group_by(bcr_patient_barcode) %>%
                          summarise(max_line=max(as.numeric(substr(treatment_line, 1, 1)))),
                        in_df %>%
                          mutate(treatment_line=as.numeric(substr(treatment_line, 1, 1))) %>%
                          select(bcr_patient_barcode, treatment_line),
                        by='bcr_patient_barcode') %>%
      mutate(agree=treatment_line==max_line) %>%
      filter(agree)
    last_line_survival <- left_join(tcga_clin %>%
                                      filter(bcr_patient_barcode %in% final_line$bcr_patient_barcode) %>%
                                      group_by(bcr_patient_barcode) %>%
                                      filter(treatment_line==max(treatment_line)) %>%
                                      mutate(followup_time=overall_survival - days_to_end) %>%
                                      select(bcr_patient_barcode, followup_time),
                                    tcga_ov %>%
                                      filter(bcr_patient_barcode %in% final_line$bcr_patient_barcode) %>%
                                      group_by(bcr_patient_barcode) %>%
                                      filter(treatment_line==max(treatment_line)) %>%
                                      mutate(TTF=TTF*(365/12)) %>%
                                      select(bcr_patient_barcode, TTF, drug_name))
  }
  if(cohort == 'tcga'){
    # Check which patients are using their last treatment line
    final_line <- merge(tcga_clin %>%
                          group_by(bcr_patient_barcode) %>%
                          summarise(max_line=max(as.numeric(substr(treatment_line, 1, 1)))),
                        in_df %>%
                          mutate(treatment_line=as.numeric(substr(treatment_line, 1, 1))) %>%
                          select(bcr_patient_barcode, treatment_line),
                        by='bcr_patient_barcode') %>%
      mutate(agree=treatment_line==max_line) %>%
      filter(agree)
    last_line_survival <- tcga_clin %>%
                            filter(bcr_patient_barcode %in% final_line$bcr_patient_barcode) %>%
                            group_by(bcr_patient_barcode) %>%
                            filter(treatment_line==max(treatment_line)) %>%
                            mutate(followup_time=overall_survival - days_to_end) %>%
                            select(bcr_patient_barcode, followup_time, TTF)
  }
  patients_to_exclude <- last_line_survival %>% filter(followup_time>365 & TTF>730) %>% pull(bcr_patient_barcode)
  in_df <- in_df[!in_df$bcr_patient_barcode%in%patients_to_exclude,] #no followup info
  return(in_df)
}


phase_3_func <- function(cancer, chemo, cohort, SA_exp=TRUE, SA_control=FALSE, control_pattern=NULL){
  if(cohort=="hmf"){
    hmf_clin <- readRDS(file.path(input_dir, 'HMF_clinical_data.RDS')) #output of "Curating_HMF_Clinical_Data.R"
    metadata <- read_delim(file.path(input_dir, 'metadata.tsv'), #download from HMF
                           delim = "\t",
                           escape_double = FALSE,
                           trim_ws = TRUE)
    metadata$patientIdentifier <- gsub('TI*V*$', '', metadata$sampleId, perl=TRUE)
    cancer_df <- restrict_cancer(hmf_clin, cancer)
    cancer_df <- immediate_post_biopsy_samples_only(cancer_df, metadata)
  }else if(cohort=="tcga"){
    tcga_clin <- readRDS(file.path(input_dir, 'TCGA_clinical_data.RDS'))
    cancer_df <- restrict_cancer(tcga_clin, cancer)
    cancer_df <- remove_na_ttf(cancer_df)
    cancer_df <- remove_stage1(cancer_df)
  }else if(cohort=='tcga-ov'){
    tcga_clin <- prepare_tcga_ov()
    cancer_df <- remove_na_ttf(tcga_clin)
  }
  if(cancer=='Breast' | cancer=='BRCA'){
    cancer_df <- remove_herceptin(cancer_df, cohort)
  }
  if(cancer=='OV'){
    if(cohort=='tcga'){
      cancer_df <- cancer_df %>% filter(treatment_line!=1)
    }else if(cohort=='tcga-ov'){
      cancer_df <- cancer_df %>% filter(treatment_line!='1st regimen')
    }
  }
  cancer_df <- predict_chemo_resistance(cancer_df, chemo, cohort)
  # Resistant SA
  cancer_df_resistant <- cancer_df %>% filter(prediction=='Resistant')
  if(SA_exp){
    cancer_sa_resistant <- identify_sa_patients(cancer_df_resistant, chemo, cohort)
  }else{
    if(cohort=='hmf'){
      chemo_str <- c('platin', 'taxel', 'rubicin')
      names(chemo_str) <- c('platinum', 'taxane', 'doxorubicin')
      cancer_sa_resistant <- cancer_df_resistant %>% filter(grepl(chemo_str[[chemo]], name))
    }else if(cohort=="tcga-ov"){
      cancer_sa_resistant <- cancer_df_resistant %>% filter(grepl(chemo, drug_name) & !grepl('platin', drug_name))
    }else{
      cancer_sa_resistant <- cancer_df_resistant %>% filter(grepl(chemo, drug_name))
    }
  }
  # Sensitive SA
  cancer_df_sensitive <- cancer_df %>% filter(prediction=='Sensitive')
  if(SA_exp){
    cancer_sa_sensitive <- identify_sa_patients(cancer_df_sensitive, chemo, cohort)
  }else{
    if(cohort=='hmf'){
      chemo_str <- c('platin', 'taxel', 'rubicin')
      names(chemo_str) <- c('platinum', 'taxane', 'doxorubicin')
      cancer_sa_sensitive <- cancer_df_sensitive %>% filter(grepl(chemo_str[[chemo]], name))
    }else if(cohort=="tcga-ov"){
      cancer_sa_sensitive <- cancer_df_sensitive %>% filter(grepl(chemo, drug_name) & !grepl('platin', drug_name))
    }else{
      cancer_sa_sensitive <- cancer_df_sensitive %>% filter(grepl(chemo, drug_name))
    }
  }
  # Generate SoC
  cancer_soc_resistant <- identify_soc_patients(cancer_df_resistant, chemo, cohort, SA_control, control_pattern)
  cancer_soc_sensitive <- identify_soc_patients(cancer_df_sensitive, chemo, cohort, SA_control, control_pattern)
  # Filter out patients who received too little chemo
  cancer_soc_resistant <- filter_short_cycles(cancer_soc_resistant, cohort)
  cancer_soc_sensitive <- filter_short_cycles(cancer_soc_sensitive, cohort)
  cancer_sa_resistant <- filter_short_cycles(cancer_sa_resistant, cohort)
  cancer_sa_sensitive <- filter_short_cycles(cancer_sa_sensitive, cohort)
  # Remove patients from SoC that appear in SA
  cancer_soc_resistant <- identify_overlapping_patients(cancer_sa_resistant,
                                                        cancer_soc_resistant,
                                                        cohort)
  cancer_soc_sensitive <- identify_overlapping_patients(cancer_sa_sensitive,
                                                        cancer_soc_sensitive,
                                                        cohort)
  # Top 5 SoCs
  if(cancer=='BRCA'){
    cancer_soc_resistant$drug_name <- gsub('cytoxan', 'cyclophosphamide', cancer_soc_resistant$drug_name)
    cancer_soc_sensitive$drug_name <- gsub('cytoxan', 'cyclophosphamide', cancer_soc_sensitive$drug_name)
  }
  top_soc <- top_N_SoC(cancer_soc_resistant, cancer_soc_sensitive, cohort=cohort)
  if(cancer=='BRCA'){
    top_soc[top_soc=='tamoxifen'] <- 'methotrexate'
    top_soc[top_soc=='arimidex'] <- 'epirubicin'
  }
  cancer_soc_resistant <- identify_top_SoC(cancer_soc_resistant, top_soc, cohort)
  cancer_soc_sensitive <- identify_top_SoC(cancer_soc_sensitive, top_soc, cohort)
  # Keep only 1st available treatment line
  if(cohort=='hmf'){
    cancer_soc_resistant <- cancer_soc_resistant %>% arrange(patientIdentifier, treatment_line) %>% filter(!duplicated(patientIdentifier))
    cancer_soc_sensitive <- cancer_soc_sensitive %>% arrange(patientIdentifier, treatment_line) %>% filter(!duplicated(patientIdentifier))
    cancer_sa_resistant <- cancer_sa_resistant %>% arrange(patientIdentifier, treatment_line) %>% filter(!duplicated(patientIdentifier))
    cancer_sa_sensitive <- cancer_sa_sensitive %>% arrange(patientIdentifier, treatment_line) %>% filter(!duplicated(patientIdentifier))
  }else{
    cancer_soc_resistant <- cancer_soc_resistant %>% arrange(bcr_patient_barcode, treatment_line) %>% filter(!duplicated(bcr_patient_barcode))
    cancer_soc_sensitive <- cancer_soc_sensitive %>% arrange(bcr_patient_barcode, treatment_line) %>% filter(!duplicated(bcr_patient_barcode))
    cancer_sa_resistant <- cancer_sa_resistant %>% arrange(bcr_patient_barcode, treatment_line) %>% filter(!duplicated(bcr_patient_barcode))
    cancer_sa_sensitive <- cancer_sa_sensitive %>% arrange(bcr_patient_barcode, treatment_line) %>% filter(!duplicated(bcr_patient_barcode))
  }
  cancer_sa_resistant$Exp_treatment <- TRUE
  cancer_soc_resistant$Exp_treatment <- FALSE
  cancer_sa_sensitive$Exp_treatment <- TRUE
  cancer_soc_sensitive$Exp_treatment <- FALSE
  cancer_resistant <- rbind(cancer_sa_resistant, cancer_soc_resistant)
  cancer_sensitive <- rbind(cancer_sa_sensitive, cancer_soc_sensitive)
  if(cohort=='hmf'){
    # Annotate with covariates
    cancer_resistant <- annotate_with_covariates(cancer_resistant)
    cancer_sensitive <- annotate_with_covariates(cancer_sensitive)
    # Add censoring
    cancer_resistant$censoring <- ifelse(cancer_resistant$reason=='RECIST' | cancer_resistant$reason=='Length of treatment', 0, 1)
    cancer_sensitive$censoring <- ifelse(cancer_sensitive$reason=='RECIST' | cancer_sensitive$reason=='Length of treatment', 0, 1) 
  }else{
    # Exclude TCGA patients with long followup gaps with no clinical data for the last treatment line
    cancer_resistant <- exclude_ttf_long_followup(cancer_resistant, cohort)
    cancer_sensitive <- exclude_ttf_long_followup(cancer_sensitive, cohort)
    # Add censoring
    if(cohort=='tcga'){
      cancer_resistant$censoring <- ifelse(cancer_resistant$progression_reason %in% c('End of last treatment',
                                                                                      'Last followup',
                                                                                      'Last followup/Last followup/End of last treatment'), 0, 1)
      cancer_sensitive$censoring <- ifelse(cancer_sensitive$progression_reason %in% c('End of last treatment',
                                                                                      'Last followup',
                                                                                      'Last followup/Last followup/End of last treatment'), 0, 1)
    }
  }
  return(list(cancer_resistant,
              cancer_sensitive))
}

phase_2_func <- function(cancer, chemo, cohort, SA_exp=FALSE){
  if(cohort=="hmf"){
    hmf_clin <- readRDS(file.path(input_dir, 'HMF_clinical_data.RDS')) #output of "Curating_HMF_Clinical_Data.R"
    metadata <- read_delim(file.path(input_dir, 'metadata.tsv'), #downloaded from HMF
                           delim = "\t",
                           escape_double = FALSE,
                           trim_ws = TRUE)
    metadata$patientIdentifier <- gsub('TI*V*$', '', metadata$sampleId, perl=TRUE)
    cancer_df <- restrict_cancer(hmf_clin, cancer)
    cancer_df <- immediate_post_biopsy_samples_only(cancer_df, metadata)
  }else if(cohort=="tcga"){
    tcga_clin <- readRDS(file.path(input_dir, 'TCGA_clinical_data.RDS'))
    cancer_df <- restrict_cancer(tcga_clin, cancer)
    cancer_df <- remove_na_ttf(cancer_df)
    cancer_df <- remove_stage1(cancer_df)
  }else if(cohort=='tcga-ov'){
    tcga_clin <- prepare_tcga_ov()
    cancer_df <- remove_na_ttf(tcga_clin)
  }
  if(cancer=='Breast' | cancer=='BRCA'){
    cancer_df <- remove_herceptin(cancer_df, cohort)
  }
  # Predict resistance
  cancer_df <- predict_chemo_resistance(cancer_df, chemo, cohort)
  # Restrict to chemo of interest as SA or combo
  if((cancer=='OV' & chemo=='platinum') | (cancer=='SARC' & chemo=='doxorubicin')){
    cancer_df <- cancer_df %>% filter(grepl(chemo, drug_name))
  }else{
    if(SA_exp){
      cancer_df <- identify_sa_patients(cancer_df, chemo, cohort)
    }else{
      if(cohort=='hmf'){
        chemo_str <- c('platin', 'taxel', 'rubicin')
        names(chemo_str) <- c('platinum', 'taxane', 'doxorubicin')
        cancer_df <- cancer_df %>% filter(grepl(chemo_str[[chemo]], name))
      }else if(cohort=="tcga-ov"){
        cancer_df <- cancer_df %>% filter(grepl(chemo, drug_name) & !grepl('platin', drug_name))
      }else{
        cancer_df <- cancer_df %>% filter(grepl(chemo, drug_name))
      }
    }
  }
  
  cancer_res <- cancer_df %>% filter(prediction=='Resistant')
  cancer_sen <- cancer_df %>% filter(prediction=='Sensitive')
  # Filter out patients who received too little chemo
  cancer_res <- filter_short_cycles(cancer_res, cohort)
  cancer_sen <- filter_short_cycles(cancer_sen, cohort)
  if(cancer=='OV' & chemo!='platinum'){
    if(cohort=='tcga-ov'){
      cancer_res <- cancer_res %>% filter(treatment_line!='1st regimen')
      cancer_sen <- cancer_sen %>% filter(treatment_line!='1st regimen')
    }else{
      cancer_res <- cancer_res %>% filter(treatment_line!=1)
      cancer_sen <- cancer_sen %>% filter(treatment_line!=1)
    }
  }
  # Keep only 1st treatment line from sensitive & resistant patients
  if(grepl('tcga', cohort)){
    cancer_sen <- cancer_sen[!duplicated(cancer_sen$bcr_patient_barcode),]
    cancer_res <- cancer_res[!duplicated(cancer_res$bcr_patient_barcode),]
  }else if(cohort=='hmf'){
    cancer_sen <- cancer_sen[!duplicated(cancer_sen$patientIdentifier),]
    cancer_res <- cancer_res[!duplicated(cancer_res$patientIdentifier),]
  }
  cancer_df <- rbind(cancer_sen, cancer_res)
  # Annotate
  if(cohort=='hmf'){
    # Annotate with covariates
    cancer_df <- annotate_with_covariates(cancer_df)
    # Add censoring
    cancer_df$censoring <- ifelse(cancer_df$reason=='RECIST' | cancer_df$reason=='Length of treatment', 0, 1)
  }else if(cohort=='tcga'){
    # Add censoring
    cancer_df$censoring <- ifelse(cancer_df$progression_reason %in% c('End of last treatment',
                                                                      'Last followup',
                                                                      'Last followup/Last followup/End of last treatment'), 0, 1)
  }
  cancer_df$prediction <- factor(cancer_df$prediction, levels=c('Sensitive', 'Resistant'))
  return(cancer_df)
}

##### Filtering Flowchart Functions #####

make_str <- function(in_df, title='', type=TRUE, progression=TRUE, treatment=FALSE, line=FALSE, cohort){
  if(grepl('tcga', cohort)){
    n_patients <- length(unique(in_df$bcr_patient_barcode))
  }else if(cohort=='hmf'){
    n_patients <- length(unique(in_df$patientIdentifier))
  }
  str_list <- list()
  if(type){
    str_list[[length(str_list)+1]] <- make_type_str(in_df, cohort)
  }
  if(line){
    str_list[[length(str_list)+1]] <- make_line_str(in_df)
  }
  if(treatment){
    str_list[[length(str_list)+1]] <- make_treatment_str(in_df, cohort)
  }
  if(progression){
    if(cohort!='tcga-ov'){
      str_list[[length(str_list)+1]] <- make_progression_str(in_df, cohort)
    }
    str_list[[length(str_list)+1]] <- make_censoring_str(in_df, cohort)
  }
  feature_str <- paste(unlist(str_list), collapse='\\n\\n')
  out_str <- paste0(title,
                    '\\nN=',
                    n_patients,
                    '\\n# of treatment lines=',
                    nrow(in_df),
                    '\\n\\n',
                    feature_str)
  return(out_str)
}

make_line_str <- function(in_df){
  return(paste0(names(table(in_df$treatment_line)),
                ' - ',
                table(in_df$treatment_line),
                '/',
                nrow(in_df),
                ' (',
                round(table(in_df$treatment_line)/nrow(in_df)*100, digits=2),
                '%)', collapse='\\n'))
}

make_type_str <- function(in_df, cohort){
  if(grepl('tcga', cohort)){
    in_df <- in_df[!duplicated(in_df$bcr_patient_barcode),]
  }else if(cohort=='hmf'){
    in_df <- in_df[!duplicated(in_df$patientIdentifier),]
    metadata <- read_delim(file.path(input_dir, 'metadata.tsv'),
                           delim = "\t",
                           escape_double = FALSE,
                           trim_ws = TRUE)
    metadata$patientIdentifier <- gsub('TI*V*$', '', metadata$sampleId, perl=TRUE)
    in_df$type <- metadata$primaryTumorSubType[match(in_df$patientIdentifier, metadata$patientIdentifier)]
  }
  return(paste0(names(sort(table(in_df$type), decreasing=TRUE)),
                ' - ',
                sort(table(in_df$type), decreasing=TRUE),
                '/',
                nrow(in_df),
                ' (',
                round((sort(table(in_df$type), decreasing=TRUE)/nrow(in_df)*100), digits=2),
                '%)', collapse='\\n'))
}

make_treatment_str <- function(in_df, cohort){
  if(grepl('tcga', cohort)){
    drugs <- sort(table(in_df$drug_name), decreasing=TRUE)
    drug_percs <- round((sort(table(in_df$drug_name), decreasing=TRUE)/nrow(in_df)*100), digits=2)
    # Only include treatments that account for more than 1% of treatment lines
    drugs <- drugs[drug_percs>1]
    drug_percs <- drug_percs[drug_percs>1]
    return(paste0(names(drugs),
                  ' - ',
                  drugs,
                  '/',
                  nrow(in_df),
                  ' (',
                  drug_percs,
                  '%)', collapse='\\n'))
  }else if(cohort=='hmf'){
    drugs <- sort(table(in_df$name), decreasing=TRUE)
    drug_percs <- round((sort(table(in_df$name), decreasing=TRUE)/nrow(in_df)*100), digits=2)
    # Only include treatments that account for more than 1% of treatment lines
    drugs <- drugs[drug_percs>1]
    drug_percs <- drug_percs[drug_percs>1]
    return(paste0(names(drugs),
                  ' - ',
                  drugs,
                  '/',
                  nrow(in_df),
                  ' (',
                  drug_percs,
                  '%)', collapse='\\n'))
  }
}

make_progression_str <- function(in_df, cohort){
  if(grepl('tcga', cohort)){
    return(paste0(names(sort(table(in_df$progression_reason), decreasing=TRUE)),
                  ' - ',
                  sort(table(in_df$progression_reason), decreasing=TRUE),
                  '/',
                  nrow(in_df),
                  ' (',
                  round((sort(table(in_df$progression_reason), decreasing=TRUE)/nrow(in_df)*100), digits=2),
                  '%)', collapse='\\n'))
  }else if(cohort=='tcga-ov'){
    return('')
  }else if(cohort=='hmf'){
    return(paste0(names(sort(table(in_df$reason), decreasing=TRUE)),
                  ' - ',
                  sort(table(in_df$reason), decreasing=TRUE),
                  '/',
                  nrow(in_df),
                  ' (',
                  round((sort(table(in_df$reason), decreasing=TRUE)/nrow(in_df)*100), digits=2),
                  '%)', collapse='\\n'))
  }
}

make_censoring_str <- function(in_df, cohort){
  if(cohort=='tcga'){
    n_censored <- sum(grepl('Last followup', in_df$progression_reason) | grepl('End of last treatment', in_df$progression_reason))
  }else if(cohort=='tcga-ov'){
    n_censored <- nrow(in_df) - sum(in_df$censoring, na.rm = TRUE)
  }else if(cohort=='hmf'){
    n_censored <- sum(grepl('Length of treatment', in_df$reason) | grepl('RECIST', in_df$reason))
  }
  return(paste0('Overall censoring - ',
                n_censored,
                '/',
                nrow(in_df),
                ' (',
                round(n_censored/nrow(in_df) * 100, digits=2),
                '%)')) 
}

phase_3_plot_func <- function(cancer, chemo, cohort, SA_exp=TRUE, SA_control=FALSE, control_pattern=NULL){
  g <- create_graph(attr_theme = NULL) %>%
    add_global_graph_attrs(attr='layout', value='dot', attr_type='graph') %>%
    add_global_graph_attrs(attr='shape', value='box', attr_type='node')
  if(cancer=='Breast' | cancer=='BRCA'){
    str_type <- TRUE
  }else{
    str_type <- FALSE
  }
  if(cohort=="hmf"){
    hmf_clin <- readRDS(file.path(input_dir, 'HMF_clinical_data.RDS')) #output of "Curating_HMF_Clinical_Data.R"
    metadata <- read_delim(file.path(input_dir, 'metadata.tsv'), #downloaded from HMF
                           delim = "\t",
                           escape_double = FALSE,
                           trim_ws = TRUE)
    metadata$patientIdentifier <- gsub('TI*V*$', '', metadata$sampleId, perl=TRUE)
    g <- g %>% add_node(label=paste0('All HMF patients with post-biopsy treatment\nN=', length(unique(hmf_clin$patientIdentifier))))
    cancer_df <- restrict_cancer(hmf_clin, cancer)
    g <- g %>%
      add_node(label=paste0('Keep only ', cancer, ' patients'), from=1) %>%
      add_node(label=make_str(cancer_df, title=paste0('HMF ', cancer, ' patients'), cohort=cohort, type=str_type), from=2)
    cancer_df <- immediate_post_biopsy_samples_only(cancer_df, metadata)
    g <- g %>%
      add_node(label='Keep only immediate post-biopsy treatment lines', from=3) %>%
      add_node(label=make_str(cancer_df, cohort=cohort, type=str_type), from=4)
  }else if(cohort=="tcga"){
    tcga_clin <- readRDS(file.path(input_dir, 'TCGA_clinical_data.RDS'))
    g <- g %>% add_node(label=paste0('All TCGA patients\nN=', length(unique(tcga_clin$bcr_patient_barcode))))
    cancer_df <- restrict_cancer(tcga_clin, cancer)
    g <- g %>%
      add_node(label=paste0('Keep only ', cancer, ' patients'), from=1) %>%
      add_node(label=make_str(cancer_df, title=paste0('TCGA ', cancer, ' patients'), cohort='tcga', type=str_type), from=2)
    cancer_df <- remove_na_ttf(cancer_df)
    g <- g %>%
      add_node(label='Remove patients with missing TTF values', from=3) %>%
      add_node(label=make_str(cancer_df, cohort='tcga', type=str_type), from=4)
    cancer_df <- remove_stage1(cancer_df)
    g <- g %>%
      add_node(label='Remove patients with Stage I cancer', from=5) %>%
      add_node(label=make_str(cancer_df, cohort='tcga', type=str_type), from=6)
  }else if(cohort=='tcga-ov'){
    tcga_clin <- prepare_tcga_ov()
    g <- g %>%
      add_node(label=make_str(tcga_clin, title='All HGSOC TCGA patients', type=str_type, cohort=cohort))
    cancer_df <- remove_na_ttf(tcga_clin)
    g <- g %>%
      add_node(label='Remove treatments with missing TTF values', from=1) %>%
      add_node(label=make_str(cancer_df, cohort=cohort, type=str_type), from=2)
  }
  if(cancer=='Breast' | cancer=='BRCA'){
    cancer_df <- remove_herceptin(cancer_df, cohort)
    g <- g %>%
      add_node(label='Remove patients treated with Anti-HER2s', from=max(get_node_ids(g))) %>%
      add_node(label=make_str(cancer_df, cohort=cohort), from=max(get_node_ids(g))+1)
  }
  if(cancer=='OV'){
    if(cohort=='tcga'){
      cancer_df <- cancer_df %>% filter(treatment_line!=1)
    }else if(cohort=='tcga-ov'){
      cancer_df <- cancer_df %>% filter(treatment_line!='1st regimen')
    }
  g <- g %>%
    add_node(label='Remove 1st line treatments', from=max(get_node_ids(g))) %>%
    add_node(label=make_str(cancer_df, cohort=cohort), from=max(get_node_ids(g))+1)
  }
  cancer_df <- predict_chemo_resistance(cancer_df, chemo, cohort)
  cancer_df_resistant <- cancer_df %>% filter(prediction=='Resistant')
  cancer_df_sensitive <- cancer_df %>% filter(prediction=='Sensitive')
  if(grepl('tcga', cohort)){
    res_str <- paste0('Predicted Resistant\\n\\nN=',
                      length(unique(cancer_df_resistant$bcr_patient_barcode)),
                      ' (',
                      round(length(unique(cancer_df_resistant$bcr_patient_barcode))/length(unique(c(cancer_df_resistant$bcr_patient_barcode, cancer_df_sensitive$bcr_patient_barcode))) * 100, digits=2),
                      '%)\\n\\n# of treatment lines=',
                      nrow(cancer_df_resistant))
    sen_str <- paste0('Predicted Sensitive\\n\\nN=',
                      length(unique(cancer_df_sensitive$bcr_patient_barcode)),
                      ' (',
                      round(length(unique(cancer_df_sensitive$bcr_patient_barcode))/length(unique(c(cancer_df_resistant$bcr_patient_barcode, cancer_df_sensitive$bcr_patient_barcode))) * 100, digits=2),
                      '%)\\n\\n# of treatment lines=',
                      nrow(cancer_df_sensitive))
  }else if(cohort=='hmf'){
    res_str <- paste0('Predicted Resistant\\n\\nN=',
                      length(unique(cancer_df_resistant$patientIdentifier)),
                      ' (',
                      round(length(unique(cancer_df_resistant$patientIdentifier))/length(unique(c(cancer_df_resistant$patientIdentifier, cancer_df_sensitive$patientIdentifier))) * 100, digits=2),
                      '%)\\n\\n# of treatment lines=',
                      nrow(cancer_df_resistant))
    sen_str <- paste0('Predicted Sensitive\\n\\nN=',
                      length(unique(cancer_df_sensitive$patientIdentifier)),
                      ' (',
                      round(length(unique(cancer_df_sensitive$patientIdentifier))/length(unique(c(cancer_df_resistant$patientIdentifier, cancer_df_sensitive$patientIdentifier))) * 100, digits=2),
                      '%)\\n\\n# of treatment lines=',
                      nrow(cancer_df_sensitive))
  }
  g <- g %>%
    add_node(label='Predict Resistance', from=max(get_node_ids(g))) %>%
    add_node(label=res_str, from=max(get_node_ids(g)) + 1, type='pred_res') %>%
    add_node(label=sen_str, from=max(get_node_ids(g)) + 1, type='pred_sen')
  # Resistant Exp
  if(SA_exp){
    cancer_sa_resistant <- identify_sa_patients(cancer_df_resistant, chemo, cohort)
  }else{
    if(cohort=='hmf'){
      chemo_str <- c('platin', 'taxel', 'rubicin')
      names(chemo_str) <- c('platinum', 'taxane', 'doxorubicin')
      cancer_exp_resistant <- cancer_df_resistant %>% filter(grepl(chemo_str[[chemo]], name))
    }else if(cohort=="tcga-ov"){
      cancer_exp_resistant <- cancer_df_resistant %>% filter(grepl(chemo, drug_name) & !grepl('platin', drug_name))
    }else{
      cancer_exp_resistant <- cancer_df_resistant %>% filter(grepl(chemo, drug_name))
    }
  }
  g <- g %>%
    add_node(label=make_str(cancer_exp_resistant, title='Exp Resistant', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='pred_res'), type='exp_res')
  # Sensitive Exp
  if(SA_exp){
    cancer_exp_sensitive <- identify_sa_patients(cancer_df_sensitive, chemo, cohort)
  }else{
    if(cohort=='hmf'){
      chemo_str <- c('platin', 'taxel', 'rubicin')
      names(chemo_str) <- c('platinum', 'taxane', 'doxorubicin')
      cancer_exp_sensitive <- cancer_df_sensitive %>% filter(grepl(chemo_str[[chemo]], name))
    }else if(cohort=="tcga-ov"){
      cancer_exp_sensitive <- cancer_df_sensitive %>% filter(grepl(chemo, drug_name) & !grepl('platin', drug_name))
    }else{
      cancer_exp_sensitive <- cancer_df_sensitive %>% filter(grepl(chemo, drug_name))
    }
  }
  g <- g %>%
    add_node(label=make_str(cancer_exp_sensitive, title='Exp Sensitive', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='pred_sen'), type='exp_sen')
  # Generate SoC
  cancer_soc_resistant <- identify_soc_patients(cancer_df_resistant, chemo, cohort, SA_control, control_pattern)
  g <- g %>%
    add_node(label=make_str(cancer_soc_resistant, title='SoC Resistant', type=str_type, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='pred_res'), type='soc_res')
  cancer_soc_sensitive <- identify_soc_patients(cancer_df_sensitive, chemo, cohort, SA_control, control_pattern)
  g <- g %>%
    add_node(label=make_str(cancer_soc_sensitive, title='SoC Sensitive', type=str_type, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='pred_sen'), type='soc_sen')
  # Filter out patients who received too little chemo
  cancer_soc_resistant <- filter_short_cycles(cancer_soc_resistant, cohort)
  cancer_soc_sensitive <- filter_short_cycles(cancer_soc_sensitive, cohort)
  cancer_exp_resistant <- filter_short_cycles(cancer_exp_resistant, cohort)
  cancer_exp_sensitive <- filter_short_cycles(cancer_exp_sensitive, cohort)
  g <- g %>%
    add_node(label='Keep only treatments with 3+ cycles', from=get_node_ids(g, type=='soc_res'), type='soc_res') %>%
    add_node(label='Keep only treatments with 3+ cycles', from=get_node_ids(g, type=='soc_sen'), type='soc_sen') %>%
    add_node(label='Keep only treatments with 3+ cycles', from=get_node_ids(g, type=='exp_res'), type='exp_res') %>%
    add_node(label='Keep only treatments with 3+ cycles', from=get_node_ids(g, type=='exp_sen'), type='exp_sen')
  g <- g %>%
    add_node(label=make_str(cancer_soc_resistant, title='SoC Resistant', type=str_type, line=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='soc_res')), type='soc_res') %>%
    add_node(label=make_str(cancer_soc_sensitive, title='SoC Sensitive', type=str_type, line=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='soc_sen')), type='soc_sen') %>%
    add_node(label=make_str(cancer_exp_resistant, title='Exp Resistant', type=str_type, line=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='exp_res')), type='exp_res') %>%
    add_node(label=make_str(cancer_exp_sensitive, title='Exp Sensitive', type=str_type, line=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='exp_sen')), type='exp_sen')
  # Remove patients from SoC that appear in Exp
  cancer_soc_resistant <- identify_overlapping_patients(cancer_exp_resistant,
                                                        cancer_soc_resistant,
                                                        cohort)
  cancer_soc_sensitive <- identify_overlapping_patients(cancer_exp_sensitive,
                                                        cancer_soc_sensitive,
                                                        cohort)
  g <- g %>%
    add_node(label='Remove patients who received any Exp treatment', from=max(get_node_ids(g, type=='soc_res')), type='soc_res') %>%
    add_node(label='Remove patients who received any Exp treatment', from=max(get_node_ids(g, type=='soc_sen')), type='soc_sen')
  g <- g %>%
    add_node(label=make_str(cancer_soc_resistant, type=str_type, line=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='soc_res')), type='soc_res') %>%
    add_node(label=make_str(cancer_soc_sensitive, type=str_type, line=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='soc_sen')), type='soc_sen')
  # Top 5 SoCs
  if(cancer=='BRCA'){
    cancer_soc_resistant$drug_name <- gsub('cytoxan', 'cyclophosphamide', cancer_soc_resistant$drug_name)
    cancer_soc_sensitive$drug_name <- gsub('cytoxan', 'cyclophosphamide', cancer_soc_sensitive$drug_name)
  }
  top_soc <- top_N_SoC(cancer_soc_resistant, cancer_soc_sensitive, cohort=cohort)
  if(cancer=='BRCA'){
    top_soc[top_soc=='tamoxifen'] <- 'methotrexate'
    top_soc[top_soc=='arimidex'] <- 'epirubicin'
  }
  if(grepl('tcga', cohort)){
    top_soc_str <- paste0('Treatment lines with at least 5 instances\\n\\n',
                          paste(rbind(cancer_soc_resistant, cancer_soc_sensitive) %>%
                                  group_by(drug_name) %>%
                                  summarise(count=n()) %>%
                                  arrange(desc(count)) %>%
                                  mutate(perc=round(count/sum(count) * 100, digits=2)) %>%
                                  mutate(s=paste0(drug_name, ' - ', count, '/', sum(count), ' (', perc, '%)')) %>%
                                  filter(count>=5) %>%
                                  pull(s),
                                collapse='\\n'))
  }else if(cohort=='hmf'){
    top_soc_str <- paste0('Treatment lines with at least 5 instances\\n\\n',
                          paste(rbind(cancer_soc_resistant, cancer_soc_sensitive) %>%
                                  group_by(name) %>%
                                  summarise(count=n()) %>%
                                  arrange(desc(count)) %>%
                                  mutate(perc=round(count/sum(count) * 100, digits=2)) %>%
                                  mutate(s=paste0(name, ' - ', count, '/', sum(count), ' (', perc, '%)')) %>%
                                  filter(count>=5) %>%
                                  pull(s),
                                collapse='\\n'))
  }
  g <- g %>%
    add_node(label=top_soc_str)
  cancer_soc_resistant <- identify_top_SoC(cancer_soc_resistant, top_soc, cohort)
  cancer_soc_sensitive <- identify_top_SoC(cancer_soc_sensitive, top_soc, cohort)
  g <- g %>%
    add_node(label='Keep only most common SoCs', from=max(get_node_ids(g, type=='soc_res')), type='soc_res') %>%
    add_node(label='Keep only most common SoCs', from=max(get_node_ids(g, type=='soc_sen')), type='soc_sen')
  g <- g %>%
    add_node(label=make_str(cancer_soc_resistant, type=str_type, line=TRUE, treatment=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='soc_res')), type='soc_res') %>%
    add_node(label=make_str(cancer_soc_sensitive, type=str_type, line=TRUE, treatment=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='soc_sen')), type='soc_sen')
  # Keep only matching treatment lines
  if(cohort=='hmf'){
    cancer_soc_resistant <- cancer_soc_resistant %>% filter(!duplicated(patientIdentifier))
    cancer_soc_sensitive <- cancer_soc_sensitive %>% filter(!duplicated(patientIdentifier))
    cancer_exp_resistant <- cancer_exp_resistant %>% filter(!duplicated(patientIdentifier))
    cancer_exp_sensitive <- cancer_exp_sensitive %>% filter(!duplicated(patientIdentifier))
  }else{
    cancer_soc_resistant <- cancer_soc_resistant %>% filter(!duplicated(bcr_patient_barcode))
    cancer_soc_sensitive <- cancer_soc_sensitive %>% filter(!duplicated(bcr_patient_barcode))
    cancer_exp_resistant <- cancer_exp_resistant %>% filter(!duplicated(bcr_patient_barcode))
    cancer_exp_sensitive <- cancer_exp_sensitive %>% filter(!duplicated(bcr_patient_barcode))
  }
  cancer_exp_resistant$Exp_treatment <- TRUE
  cancer_soc_resistant$Exp_treatment <- FALSE
  cancer_exp_sensitive$Exp_treatment <- TRUE
  cancer_soc_sensitive$Exp_treatment <- FALSE
  cancer_resistant <- rbind(cancer_exp_resistant, cancer_soc_resistant)
  cancer_sensitive <- rbind(cancer_exp_sensitive, cancer_soc_sensitive)
  g <- g %>%
    add_node(label='Keep only 1st available line', from=c(max(get_node_ids(g, type=='soc_res')), max(get_node_ids(g, type=='exp_res'))), type='match_res') %>%
    add_node(label='Keep only 1st available line', from=c(max(get_node_ids(g, type=='soc_sen')), max(get_node_ids(g, type=='exp_sen'))), type='match_sen')
  g <- g %>%
    add_node(label=make_str(cancer_resistant[!cancer_resistant$Exp_treatment,], title='SoC Resistant', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='match_res'), type='soc_res') %>%
    add_node(label=make_str(cancer_resistant[cancer_resistant$Exp_treatment,], title='Exp Resistant', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='match_res'), type='exp_res') %>%
    add_node(label=make_str(cancer_sensitive[!cancer_sensitive$Exp_treatment,], title='SoC Sensitive', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='match_sen'), type='soc_sen') %>%
    add_node(label=make_str(cancer_sensitive[cancer_sensitive$Exp_treatment,], title='Exp Sensitive', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='match_sen'), type='exp_sen')
  # Exclude TCGA patients with long followup gaps with no clinical data for the last treatment line
  if(cohort!='hmf'){
    cancer_resistant <- exclude_ttf_long_followup(cancer_resistant, cohort)
    cancer_sensitive <- exclude_ttf_long_followup(cancer_sensitive, cohort)
    g <- g %>%
      add_node(label='Exclude patients with long followup periods', from=c(max(get_node_ids(g, type=='soc_res')), max(get_node_ids(g, type=='exp_res'))), type='followup_res') %>%
      add_node(label='Exclude patients with long followup periods', from=c(max(get_node_ids(g, type=='soc_sen')), max(get_node_ids(g, type=='exp_sen'))), type='followup_sen')
    g <- g %>%
      add_node(label=make_str(cancer_resistant[!cancer_resistant$Exp_treatment,], title='SoC Resistant', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='followup_res'), type='soc_res') %>%
      add_node(label=make_str(cancer_resistant[cancer_resistant$Exp_treatment,], title='Exp Resistant', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='followup_res'), type='exp_res') %>%
      add_node(label=make_str(cancer_sensitive[!cancer_sensitive$Exp_treatment,], title='SoC Sensitive', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='followup_sen'), type='soc_sen') %>%
      add_node(label=make_str(cancer_sensitive[cancer_sensitive$Exp_treatment,], title='Exp Sensitive', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='followup_sen'), type='exp_sen')
  }
  return(g)
}

phase_2_plot_func <- function(cancer, chemo, cohort){
  g <- create_graph(attr_theme = NULL) %>%
    add_global_graph_attrs(attr='layout', value='dot', attr_type='graph') %>%
    add_global_graph_attrs(attr='shape', value='box', attr_type='node')
  if(cancer=='Breast' | cancer=='BRCA'){
    str_type <- TRUE
  }else{
    str_type <- FALSE
  }
  if(cohort=="hmf"){
    hmf_clin <- readRDS(file.path(input_dir, 'HMF_clinical_data.RDS')) #output of "Curating_HMF_Clinical_Data.R"
    metadata <- read_delim(file.path(input_dir, 'metadata.tsv'), #downloaded from HMF
                           delim = "\t",
                           escape_double = FALSE,
                           trim_ws = TRUE)
    metadata$patientIdentifier <- gsub('TI*V*$', '', metadata$sampleId, perl=TRUE)
    g <- g %>% add_node(label=paste0('All HMF patients with post-biopsy treatment\nN=', length(unique(hmf_clin$patientIdentifier))))
    cancer_df <- restrict_cancer(hmf_clin, cancer)
    g <- g %>%
      add_node(label=paste0('Keep only ', cancer, ' patients'), from=1) %>%
      add_node(label=make_str(cancer_df, title=paste0('HMF ', cancer, ' patients'), cohort=cohort, type=str_type), from=2)
    cancer_df <- immediate_post_biopsy_samples_only(cancer_df, metadata)
    g <- g %>%
      add_node(label='Keep only immediate post-biopsy treatment lines', from=3) %>%
      add_node(label=make_str(cancer_df, cohort=cohort, type=str_type), from=4)
  }else if(cohort=="tcga"){
    tcga_clin <- readRDS(file.path(input_dir, 'TCGA_clinical_data.RDS'))
    g <- g %>% add_node(label=paste0('All TCGA patients\nN=', length(unique(tcga_clin$bcr_patient_barcode))))
    cancer_df <- restrict_cancer(tcga_clin, cancer)
    g <- g %>%
      add_node(label=paste0('Keep only ', cancer, ' patients'), from=1) %>%
      add_node(label=make_str(cancer_df, title=paste0('TCGA ', cancer, ' patients'), cohort='tcga', type=str_type), from=2)
    cancer_df <- remove_na_ttf(cancer_df)
    g <- g %>%
      add_node(label='Remove patients with missing TTF values', from=3) %>%
      add_node(label=make_str(cancer_df, cohort='tcga', type=str_type), from=4)
    cancer_df <- remove_stage1(cancer_df)
    g <- g %>%
      add_node(label='Remove patients with Stage I cancer', from=5) %>%
      add_node(label=make_str(cancer_df, cohort='tcga', type=str_type), from=6)
  }else if(cohort=='tcga-ov'){
    tcga_clin <- prepare_tcga_ov()
    g <- g %>%
      add_node(label=make_str(tcga_clin, title='All HGSOC TCGA patients', type=str_type, cohort=cohort))
    cancer_df <- remove_na_ttf(tcga_clin)
    g <- g %>%
      add_node(label='Remove treatments with missing TTF values', from=1) %>%
      add_node(label=make_str(cancer_df, cohort=cohort, type=str_type), from=2)
  }
  if(cancer=='Breast' | cancer=='BRCA'){
    cancer_df <- remove_herceptin(cancer_df, cohort)
    g <- g %>%
      add_node(label='Remove patients treated with Anti-HER2s', from=max(get_node_ids(g))) %>%
      add_node(label=make_str(cancer_df, cohort=cohort), from=max(get_node_ids(g))+1)
  }
  # Predict resistance
  cancer_df <- predict_chemo_resistance(cancer_df, chemo, cohort)
  if(cancer=='OV' & chemo=='platinum'| (cancer=='SARC' & chemo=='doxorubicin')){
    # Make single-chemo datasets
    cancer_df <- cancer_df %>% filter(grepl(chemo, drug_name))
    g <- g %>%
      add_node(label=paste0('Keep only ', chemo, '-containing treatment lines'), from=max(get_node_ids(g))) %>%
      add_node(label=make_str(cancer_df, cohort=cohort, type=str_type, treatment=TRUE), from=max(get_node_ids(g))+1)
  }else{
    if(SA_exp){
      # Restrict to SA chemo only
      cancer_df <- identify_sa_patients(cancer_df, chemo, cohort)
      g <- g %>%
        add_node(label=paste0('Keep only SA ', chemo, ' treatment lines'), from=max(get_node_ids(g))) %>%
        add_node(label=make_str(cancer_df, cohort=cohort, type=str_type, treatment=TRUE), from=max(get_node_ids(g))+1)
    }else{
      # Make single-chemo datasets
      if(cohort=='hmf'){
        chemo_str <- c('platin', 'taxel', 'rubicin')
        names(chemo_str) <- c('platinum', 'taxane', 'doxorubicin')
        cancer_df <- cancer_df %>% filter(grepl(chemo_str[[chemo]], name))
        g <- g %>%
          add_node(label=paste0('Keep only ', chemo_str[[chemo]], '-containing treatment lines'), from=max(get_node_ids(g))) %>%
          add_node(label=make_str(cancer_df, cohort=cohort, type=str_type, treatment=TRUE), from=max(get_node_ids(g))+1)
      }else if(cohort=="tcga-ov"){
        cancer_df <- cancer_df %>% filter(grepl(chemo, drug_name) & !grepl('platin', drug_name))
        g <- g %>%
          add_node(label=paste0('Keep only ', chemo, '-containing treatment lines'), from=max(get_node_ids(g))) %>%
          add_node(label=make_str(cancer_df, cohort=cohort, type=str_type, treatment=TRUE), from=max(get_node_ids(g))+1)
      }else{
        cancer_df <- cancer_df %>% filter(grepl(chemo, drug_name))
        g <- g %>%
          add_node(label=paste0('Keep only ', chemo, '-containing treatment lines'), from=max(get_node_ids(g))) %>%
          add_node(label=make_str(cancer_df, cohort=cohort, type=str_type, treatment=TRUE), from=max(get_node_ids(g))+1)
      }
    }
  }
  cancer_res <- cancer_df %>% filter(prediction=='Resistant')
  cancer_sen <- cancer_df %>% filter(prediction=='Sensitive')
  if(grepl('tcga', cohort)){
    res_str <- paste0('Predicted Resistant\\n\\nN=',
                      length(unique(cancer_res$bcr_patient_barcode)),
                      ' (',
                      round(length(unique(cancer_res$bcr_patient_barcode))/length(unique(c(cancer_res$bcr_patient_barcode, cancer_sen$bcr_patient_barcode))) * 100, digits=2),
                      '%)\\n\\n# of treatment lines=',
                      nrow(cancer_res))
    sen_str <- paste0('Predicted Sensitive\\n\\nN=',
                      length(unique(cancer_sen$bcr_patient_barcode)),
                      ' (',
                      round(length(unique(cancer_sen$bcr_patient_barcode))/length(unique(c(cancer_res$bcr_patient_barcode, cancer_sen$bcr_patient_barcode))) * 100, digits=2),
                      '%)\\n\\n# of treatment lines=',
                      nrow(cancer_sen))
  }else if(cohort=='hmf'){
    res_str <- paste0('Predicted Resistant\\n\\nN=',
                      length(unique(cancer_res$patientIdentifier)),
                      ' (',
                      round(length(unique(cancer_res$patientIdentifier))/length(unique(c(cancer_res$patientIdentifier, cancer_sen$patientIdentifier))) * 100, digits=2),
                      '%)\\n\\n# of treatment lines=',
                      nrow(cancer_res))
    sen_str <- paste0('Predicted Sensitive\\n\\nN=',
                      length(unique(cancer_sen$patientIdentifier)),
                      ' (',
                      round(length(unique(cancer_sen$patientIdentifier))/length(unique(c(cancer_res$patientIdentifier, cancer_sen$patientIdentifier))) * 100, digits=2),
                      '%)\\n\\n# of treatment lines=',
                      nrow(cancer_sen))
  }
  g <- g %>%
    add_node(label='Predict Resistance', from=max(get_node_ids(g))) %>%
    add_node(label=res_str, from=max(get_node_ids(g)) + 1, type='pred_res') %>%
    add_node(label=sen_str, from=max(get_node_ids(g)) + 1, type='pred_sen')
  # Filter out patients who received too little chemo
  cancer_res <- filter_short_cycles(cancer_res, cohort)
  cancer_sen <- filter_short_cycles(cancer_sen, cohort)
  g <- g %>%
    add_node(label='Keep only treatment lines with 3+ cycles', from=max(get_node_ids(g, type=='pred_res')), type='pred_res') %>%
    add_node(label='Keep only treatment lines with 3+ cycles', from=max(get_node_ids(g, type=='pred_sen')), type='pred_sen')
  g <- g %>%
    add_node(label=make_str(cancer_res, type=str_type, line=TRUE, treatment=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='pred_res')), type='pred_res') %>%
    add_node(label=make_str(cancer_sen, type=str_type, line=TRUE, treatment=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='pred_sen')), type='pred_sen')
  if(cancer=='OV' & chemo!='platinum'){
    if(cohort=='tcga-ov'){
      cancer_res <- cancer_res %>% filter(treatment_line!='1st regimen')
      cancer_sen <- cancer_sen %>% filter(treatment_line!='1st regimen')
    }else{
      cancer_res <- cancer_res %>% filter(treatment_line!=1)
      cancer_sen <- cancer_sen %>% filter(treatment_line!=1)
    }
    g <- g %>%
      add_node(label='Remove 1st line treatments', from=max(get_node_ids(g, type=='pred_res')), type='pred_res') %>%
      add_node(label='Remove 1st line treatments', from=max(get_node_ids(g, type=='pred_sen')), type='pred_sen')
    g <- g %>%
      add_node(label=make_str(cancer_res, type=str_type, line=TRUE, treatment=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='pred_res')), type='pred_res') %>%
      add_node(label=make_str(cancer_sen, type=str_type, line=TRUE, treatment=TRUE, cohort=cohort), from=max(get_node_ids(g, type=='pred_sen')), type='pred_sen')
  }
  # Keep only 1st available treatment line
  if(grepl('tcga', cohort)){
    cancer_sen <- cancer_sen[!duplicated(cancer_sen$bcr_patient_barcode),]
    cancer_res <- cancer_res[!duplicated(cancer_res$bcr_patient_barcode),]
  }else if(cohort=='hmf'){
    cancer_sen <- cancer_sen[!duplicated(cancer_sen$patientIdentifier),]
    cancer_res <- cancer_res[!duplicated(cancer_res$patientIdentifier),]
  }
  g <- g %>%
    add_node(label='Keep only 1st available line', from=c(max(get_node_ids(g, type=='pred_res')), max(get_node_ids(g, type=='pred_sen'))), type='match')
  g <- g %>%
    add_node(label=make_str(cancer_res, title='Resistant', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='match'), type='pred_res') %>%
    add_node(label=make_str(cancer_sen, title='Sensitive', type=str_type, treatment=TRUE, line=TRUE, cohort=cohort), from=get_node_ids(g, type=='match'), type='pred_sen')
  return(g)
}

#### Functions for processing segment tables ###################################

idSmoothingTargets = function(dfAllSegs, WIGGLE, colNameSegVal, colNameChr, IGNOREDELS = TRUE) {
  
  ### Check column name
  testSegVal = dfAllSegs[[colNameSegVal]][1]
  testChr = dfAllSegs[[colNameChr]][1]
  if(!is.numeric(testSegVal)){stop("Segment Value column has no numeric value in it. Supplied correct column name? Forgot conversion?")}
  if(is.null(testSegVal)){stop("Chromosome column has no numeric value in it. Supplied correct column name?")}
  
  # take differences to segment down below
  dfAllSegs$diffs = c(abs(dfAllSegs[[colNameSegVal]][1:(nrow(dfAllSegs)-1)] - dfAllSegs[[colNameSegVal]][2:nrow(dfAllSegs)]), WIGGLE+1)
  # set TRUE if difference to next segment is smaller than the user supplied cutoff
  dfAllSegs$smooth = dfAllSegs$diffs <= WIGGLE
  # set all segments which are last in a chromosome to FALSE. This also prevents leaking to other samples and cohorts.
  dfAllSegs$smooth[cumsum(rle(as.character(dfAllSegs[[colNameChr]]))$lengths)] = FALSE
  
  # Ignore deletions if wished
  if(IGNOREDELS) {dfAllSegs$smooth[dfAllSegs[[colNameSegVal]] == 0] = FALSE}
  
  return(dfAllSegs)
}


smoothSegments = function(lRaw, CORES, SMOOTHINGFACTOR, colNameMerge, colNameChr, colNameStart, colNameEnd,
                          IGNOREDELS = TRUE, asDf = FALSE) {
  
  ### Check column names
  test = lRaw[[1]]
  testMerge = test[[colNameMerge]][1]
  testChr = test[[colNameChr]][1]
  testStart = test[[colNameStart]][1]
  testEnd = test[[colNameEnd]][1]
  if(! is.numeric(testMerge)) { stop("Merge column has no numeric value in it. Supplied correct column name?")}
  if(is.null(testChr)) { stop("Chromosome column has no numeric value in it. Supplied correct column name?")}
  if(! is.numeric(testStart)) { stop("Start column has no numeric value in it. Supplied correct column name?")}
  if(! is.numeric(testEnd)) { stop("End column has no numeric value in it. Supplied correct column name?")}
  
  # Add diff column to names we want to keep when merging (comes from function "idSmoothingTargets").
  colNameMerge = c(colNameMerge, "diffs")
  
  registerDoMC(CORES)
  lSmooth = foreach(thisSample = lRaw, .final = function(x) setNames(x, names(lRaw)) ) %dopar% {
    
    thisOut = thisSample
    stillSmoothing = sum(thisOut$smooth)
    while( stillSmoothing > 0 ) {
      # For the while loop:
      # Read lines from thisSample and change in thisOut. Hence for a new iteration I need to sync the two.
      thisSample = thisOut
      
      rleRaw = rle(thisSample$smooth)
      # This takes the indeces of the FALSE chains and adds 1. This should give you the next segment which is TRUE.
      # Two challenges:
      # 1) Last segment always FALSE (see above), hence removal of the last number as this would indicate to a segment outside the df.
      # 2) If it starts with a TRUE segment, this would not be found when looking at the FALSE chains. Hence, adding index 1 manually if chain starts with TRUE.
      indRaw = cumsum(rleRaw$lengths)[ ! rleRaw$values ] + 1
      indRaw = indRaw[ -length(indRaw) ]
      if( rleRaw$values[1] ) { indRaw = c(1, indRaw) }
      
      # loop over start indices of TRUE chains.
      for(i in indRaw) {
        # detect length of segments to smooth. add 1 as the last segment has a FALSE value in it but still belongs to this chain.
        endOfStreak = i + rle(thisSample$smooth[i:nrow(thisSample)])$lengths[1]
        # extract reads
        dfMerge = thisSample[i:endOfStreak,]
        
        # too stupid to make this work with data.table
        newElement = as.data.frame( dfMerge[1,] )
        # Get new end and check first wether valid number.
        newEnd = dfMerge[nrow(dfMerge),][[colNameEnd]]
        if(! is.null(newEnd)) {
          newElement[[colNameEnd]] = newEnd
        } else {
          stop("New end coordinate is null. Supplied correct column name?")
        }
        ## Column "segVal" will be dealt with in a minute. Column "diffs" later when running again idSmoothingTargets.
        
        # Merge cn specifically by taking the length of the elements into consideration
        widthWeights = dfMerge[[colNameEnd]] - dfMerge[[colNameStart]]
        # save those cases where length of the segment is zero in both segs to merge (that happens when SV breaks are included in the release)
        if(widthWeights[1]==0 && widthWeights[2]==0){
          newElement[[colNameMerge[1]]] = mean(dfMerge[[colNameMerge[1]]])
        } else {
          newElement[[colNameMerge[1]]] = weighted.mean(dfMerge[[colNameMerge[1]]], widthWeights)
        }
        
        # Replace all to merge segments with the new merged segment. Later delete duplicated.
        thisOut[i:endOfStreak,] = newElement
      }
      
      # as we have replaced all segments with the new mean segment, we need to remove the duplicates
      thisOut = thisOut[ ! duplicated(thisOut), ]
      # again detect segments which needs smoothing
      thisOut = idSmoothingTargets(thisOut, SMOOTHINGFACTOR, colNameSegVal = colNameMerge[[1]], colNameChr = colNameChr,
                                   IGNOREDELS = IGNOREDELS)
      stillSmoothing = sum(thisOut$smooth)
    }
    
    # after smoothing is finished, change name of cohort
    thisOut$smooth = NULL
    thisOut$diffs = NULL
    return( thisOut )
  }
  
  if( isTRUE(asDf) ) {
    dfSmooth = setDT( rbindlist( lSmooth ) )
    return( dfSmooth )
  } else {
    return( lSmooth )
  }
  
}


#### Functions for swimmer plots ###############################################

prepare_swimmer<-function(chemo_history,clinical_history,chemo){
  samps = unique(chemo_history$study_subject_id)
  
  out=list()
  
  # main bars
  main.data = chemo_history[chemo_history$study_subject_id%in%samps,]
  
  # dots 
  point.data = as.data.frame(do.call(rbind, lapply(samps, function(s){
    clinical_history_patient <- clinical_history %>%
      filter(study_subject_id == s)
    
    if(chemo == "Platinum"){
      clinical_history_patient_line <- clinical_history_patient %>%
        filter(treat_line == 1)
    }else{
      clinical_history_patient_line <- clinical_history_patient %>%
        filter(treat_line == main.data$treatment_line[main.data$study_subject_id==s])
      new_zero=min(clinical_history_patient_line$days_since_first_treatment)
      clinical_history_patient_line$days_since_first_treatment=clinical_history_patient_line$days_since_first_treatment-new_zero
      clinical_history_patient_line$treat_line_duration=clinical_history_patient_line$treat_line_duration-new_zero
    }
    
    clinical_treatment_patient <- clinical_history_patient_line %>%
      filter(event != 'CA125 reading')
    
    if(chemo == "Platinum"){
      clinical_treatment_patient <- clinical_treatment_patient %>% 
        filter(event %in% c("CARBOPLATIN","CISPLATIN","PLATINUM"))
      
    }
    if(chemo == "Taxane"){
      clinical_treatment_patient <- clinical_treatment_patient %>% 
        filter(event == "PACLITAXEL")
    }
    if(chemo == "Doxorubicin"){
      clinical_treatment_patient <- clinical_treatment_patient %>% 
        filter(event %in% c("DOXORUBICIN","EPIRUBICIN"))
    }
    
    treat_start = unique(min(clinical_treatment_patient$days_since_first_treatment))
    treat_end = unique(max(clinical_treatment_patient$days_since_first_treatment))
    treat_line = unique(clinical_treatment_patient$treat_line)
    treat_line_duration = unique(clinical_treatment_patient$treat_line_duration)
    next_line_start = treat_line_duration
    
    clinical_history_patient <- clinical_history_patient_line %>%
      filter(event == 'CA125 reading')
    clinical_history_patient = rbind(as.data.frame(clinical_history_patient),c(s,"Treatment_start",treat_start,treat_line,treat_line_duration))
    clinical_history_patient = rbind(as.data.frame(clinical_history_patient),c(s,"Treatment_end",treat_end,treat_line,treat_line_duration))
    clinical_history_patient = rbind(as.data.frame(clinical_history_patient),c(s,"Next_treatment",next_line_start,treat_line,treat_line_duration))
    
    if(main.data$censoring[main.data$study_subject_id==s]==0){
      censored_time = main.data$PFS[main.data$study_subject_id==s]
      clinical_history_patient = rbind(as.data.frame(clinical_history_patient),c(s,"Censored",censored_time,treat_line,treat_line_duration))
      colnames(clinical_history_patient) = c("study_subject_id", "event", "days_since_first_treatment","treat_line","treat_line_duration")
      clinical_history_patient = clinical_history_patient[clinical_history_patient$event!="Next_treatment",]
    }
    colnames(clinical_history_patient) = c("study_subject_id", "event", "days_since_first_treatment","treat_line","treat_line_duration")
    return(clinical_history_patient)
  })))
  point.data[,3:5]=sapply(3:5,function(x)as.numeric(point.data[,x]))
  
  out[["main"]]=main.data
  out[["points"]]=point.data
  out[["follow"]]=point.data[point.data$event=="CA125 reading",]
  out[["treatment"]]=point.data[point.data$event%in%c("Treatment_start","Treatment_end"),]
  out[["next_treatment"]]=point.data[point.data$event%in%c("Next_treatment","Censored"),]
  out[["order"]]=main.data$study_subject_id[order(main.data$PFS_from_initial_treatment)]
  return(out)
}

