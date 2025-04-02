# This script is for generating the curated TCGA clinical history file 'TCGA_clinical_history.RDS'

# Clean environment
freshr::freshr()

# Libraries
library(this.path)
library(dplyr)

# Files
pancancer_dir <- dirname(this.path())
base_dir <- dirname(dirname(pancancer_dir))
input_dir <- file.path(base_dir, 'Input_Data')
dir.create(file.path(input_dir, 'TCGA_XMLs'))
xml_dirs <- file.path(input_dir, 'TCGA_XMLs')
helper_dir <- file.path(base_dir, 'Analysis/Helper_Scripts')
source(file.path(helper_dir, 'Helper_functions.R'))

# Download TCGA clinical XMLs
tcga_types <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
                "LUAD", "LUSC", "MESO", "PAAD", "PCPG", "PRAD", "READ", "SARC", "OV", # Note for donwstream analyses we will use the curated clinical data from Villalobos, Wang, & Sikic
                "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
for(i in tcga_types){
  download_gdc(i)
}

# Parse XMLs by TCGA cancer type
for(i in 1:length(tcga_types)){
  print(tcga_types[i])
  xml_paths <- list.files(file.path(xml_dirs, tcga_types[i]), full.names=TRUE)
  if(i==1){
    tcga_clin_df <- parse_xmls(xml_paths)
    tcga_clin_df$cancer_type <- tcga_types[i] 
    tcga_clin_df <- cbind(tcga_clin_df,
                          data.frame(Kalecky=NA,
                                     Thennavan=NA,
                                     TBLab=NA,
                                     count=NA,
                                     ER=NA,
                                     PR=NA,
                                     HER2=NA,
                                     type=NA))
  }else{
    temp_df <- parse_xmls(xml_paths)
    if(nrow(temp_df)>0){
      temp_df$cancer_type <- tcga_types[i]
      if(tcga_types[i]=='BRCA'){
        brca_types <- get_brca_types()
        temp_df <- merge(temp_df, brca_types, by='bcr_patient_barcode', all.x = TRUE)
      }else{
        temp_df <- cbind(temp_df,
                         data.frame(Kalecky=NA,
                                    Thennavan=NA,
                                    TBLab=NA,
                                    count=NA,
                                    ER=NA,
                                    PR=NA,
                                    HER2=NA,
                                    type=NA))
      }
      tcga_clin_df <- rbind(tcga_clin_df, temp_df)
    }
  }
}
tcga_clin_df <- tcga_clin_df %>% arrange(bcr_patient_barcode, treatment_line)

# Add BRCA mutant status
brca_status <- readRDS(file.path(input_dir, 'TCGA_BRCA_MutStatus.rds'))
tcga_clin_df <- merge(tcga_clin_df, brca_status[,c(1,3)], by.x='bcr_patient_barcode', by.y='Sample', all.x=TRUE)
colnames(tcga_clin_df)[which(colnames(tcga_clin_df)=='Status')] <- 'brca_mut_status'

# Add wGII
wgii <- get_wgii_tcga()
tcga_clin_df <- left_join(tcga_clin_df, wgii)

# Write to input dir
tcga_clin_df <- tcga_clin_df %>% arrange(bcr_patient_barcode, treatment_line)
saveRDS(tcga_clin_df, 'TCGA_clinical_history.RDS')
