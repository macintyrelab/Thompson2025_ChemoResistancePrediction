################################################################################
### This script is to identify the signatures that have predictive signal for
### taxanes. To this, we use cell lines from the DepMap database, which have
### been treated with different drugs
################################################################################


rm(list = ls(all = TRUE))

## PACKAGES
library(this.path)
library(dplyr)
library(ggplot2)
library(ggthemes)

## DIRECTORIES 
drug_dir <- dirname(this.path())
base_dir <- dirname(dirname(drug_dir))
input_dir <- file.path(base_dir, "Input_Data")
output_dir <- file.path(drug_dir, "outputs")
dir.create(output_dir)
helper_dir <- file.path(dirname(drug_dir), 'Helper_Scripts')
figs_dir <- file.path(base_dir, "Figures")

## PLOTTING THEME
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6),
             axis.line = element_line(linewidth = 0.5), 
             axis.ticks = element_line(linewidth = 0.5),
             axis.ticks.length = unit(.1, "cm"))


## LOAD FILES
cell_exp <- readRDS(file.path(input_dir, '4_Activities_CompendiumCINSigs_THRESH95_DepMap.rds'))
cell_id_mapping <- read.table(file.path(input_dir, 'cellLine_CEL_file_mapping.tsv'), sep = '\t', header = T, stringsAsFactors = F)
samp_info <- read.csv(file.path(input_dir, 'sample_info.csv'), header = T)
drug_screen <- read.csv(file.path(input_dir, 'secondary-screen-dose-response-curve-parameters.csv'), stringsAsFactors = F)
cell_line_info <- read.csv(file.path(input_dir, 'secondary-screen-cell-line-info.csv'))

## PROCESS SIGNATURE ACTIVITIES
#remove duplicate entries
rownames(cell_exp) <- cell_id_mapping[match(rownames(cell_exp),
                                            cell_id_mapping[cell_id_mapping$study == 'CCLE', 'fileid']),
                                      'cellid']
cell_exp <- cell_exp[!duplicated(rownames(cell_exp)), ]

rownames(cell_exp) <- samp_info[unlist(sapply(rownames(cell_exp),
                                              function(x) grep(x, samp_info$cell_line_name, ignore.case = T)[1])), 1]
cell_exp <- cell_exp[!is.na(rownames(cell_exp)), ]
cell_exp <- cell_exp[!duplicated(rownames(cell_exp)), ]

# Total number of cell lines with signature activities
nrow(cell_exp)


## PREPARE DRUG SCREEN & SAMPLE INFO DATA
drug_screen <- drug_screen[drug_screen$screen_id == 'HTS002', ]
cell_exp_filt_drug <- cell_exp[rownames(cell_exp) %in% drug_screen$depmap_id, ]

drug_lookup <- c()
unique_target <- unique(drug_screen[, c('name', 'target')])
for(i in 1:nrow(unique_target)){
  drug <- unique_target[i, 'name']
  targets <- unlist(strsplit(unique_target[i, 'target'], ', '))
  drug_lookup <- rbind(drug_lookup, cbind(drug, targets))
}

drug_lookup <- unique(drug_lookup)
drug_lookup <- data.frame(drug_lookup)
colnames(drug_lookup) <- c('drug', 'target')
drug_lookup <- drug_lookup[!is.na(drug_lookup$target), ]

# Number of cell lines with drug response data
nrow(cell_exp_filt_drug)


## SELECT CYTOTOXIC CHEMOTHERAPIES
platins <- unique(drug_lookup$drug[grep('platin', drug_lookup$drug)])
taxanes <- unique(drug_lookup$drug[grep('tax', drug_lookup$drug)])
anthracyclines <- unique(drug_lookup$drug[grep('rubicin', drug_lookup$drug)])
drugs <- c(platins, taxanes, anthracyclines)

# Number of cytotoxic chemotherapies included in the drug screen
length(drugs)


## FILTER CYTOTOXIC CHEMOTHREAPIS
ranges <- c()
for(d in drugs){
  screen <- drug_screen[drug_screen$name == d, ]
  screen <- screen[screen$auc <= 1, ]
  ranges <- rbind(ranges, c(d,
                            mean(screen$auc),
                            sd(screen$auc),
                            max(screen$auc),
                            min(screen$auc)))
  p <- ggplot(screen) +
    geom_density(aes(x = auc)) + 
    labs(title = 'AUC distribution in all cell lines') + 
    xlim(min(screen$auc) - 0.1, max(screen$auc) + 0.1) +
    scale_color_manual(values = alpha(c('gray40', 'red'), 1)) +  
    scale_fill_manual(values = alpha(c('gray40', 'red'), 0.6)) + 
    scale_x_continuous(breaks = seq(round(min(screen$auc)) - 0.1,
                                    round(max(screen$auc)) + 0.1,
                                    0.1)) + 
    xlab(paste0('Response to ', d, ' (AUC)')) +
    theme(text = element_text(size = 8), 
          axis.text = element_text(size = 6),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, 'cm'),
          legend.position = 'right')
  ggsave(file.path(figs_dir, paste0('Additional/Distribution_PanCancer_AUCs_', d, '.png')),
         width = 80, height = 80, units = 'mm', dpi = 300)
}
colnames(ranges) <- c('drug', 'mean', 'sd', 'max', 'min')
ranges <- as.data.frame(ranges)
ranges[, 2:5] <- sapply(2:5, function(x) as.numeric(ranges[, x]))

# Distribution of AUC values in all cells per drug
print(ranges)

# We cannot use *in vitro* response data for cisplatin since AUC values of all samples are 
# ~1 indicating that most cells are sensitive to this treatment at the dose tested. 
# On the other hand, the range of AUC values for anthracyclines is not wide enough 
# to perform correlations (AUC sd < 0.1). In addition, there is not enough number 
# of cells with AUC close to 0 and thus resistant to these treatments.
# Therefore, we have to limit the correlation analysis of signatures with drug response 
# data to taxanes. We also included nemorubicin (just in case we see something interesting) 
# since the AUC range is wider than the rest of anthracyclines included in the drug screening.

drugs <- c('docetaxel', 'cabazitaxel', 'paclitaxel', 'nemorubicin')

## CORRELATION ANALYSIS WITH DRUG RESPONSE PAN-CANCER
sigs <- colnames(cell_exp_filt_drug)

drug_response <- lapply(1:length(sigs), function(s){
  response <- c()
  for (d in drugs){
    pdat <- left_join(drug_screen[drug_screen$name == d, ],
                    data.frame(depmap_id = rownames(cell_exp_filt_drug),
                                    cell_exp_filt_drug),
                    by = 'depmap_id')
    pdat <- left_join(pdat, cell_line_info, by = 'depmap_id')
    pdat <- pdat[pdat$passed_str_profiling.y == T, ]
    pdat <- pdat[!is.na(pdat[, paste0('CX', s)]), ]
    moa <- unique(pdat$moa)
    #correlation
    res <- cor.test(pdat[, 'auc'], pdat[, paste0('CX', s)],
                    method = 'k', alternative = 'less')
    response <- rbind(response, c(s, d, moa, res$estimate, res$p.value))
    
    if(res$estimate < -0.07 & res$p.value < 0.05){
      dat <- as.data.frame(cbind(pdat[, 'auc'],
                                      pdat[, paste0('CX', s)],
                                      pdat$primary_tissue))
      colnames(dat) <- c('auc', 'signature', 'tissue')
      dat[, 1:2] <- sapply(1:2, function(x) as.numeric(dat[, x]))
      p <- ggplot(dat, aes(x = auc, y = signature, color = tissue)) + 
        geom_point(size = 1.5) +
        geom_smooth(method = 'lm', se = TRUE, color = 'darkred', fill = 'grey50') +
        xlim(min(dat$auc, na.rm = T) - 0.1, max(dat$auc, na.rm = T) + 0.1) +
        labs(title = 'Pancancer cell lines',
             subtitle = paste0('tau=',
                               round(res$estimate, 4),
                               '; p-value=',
                               round(res$p.value, 4))) +
        ylab(paste0('CX', s, ' activity levels')) +
        xlab(paste0('Response to ', d, ' (AUC)')) +
        theme(text = element_text(size = 8),
              axis.text = element_text(size = 6),
              axis.line = element_line(size = 0.5), 
              axis.ticks = element_line(size = 0.5),
              axis.ticks.length = unit(.1, 'cm'),
              legend.position = 'right')
      ggsave(file.path(figs_dir, paste0('Additional/PanCancer_CX', s, '_', d, '.png')),
             width = 200, height = 80, units = 'mm', dpi = 300)
    }
  }
  colnames(response) <- c('signature', 'drug', 'moa', 'tau', 'p.value')
  response <- cbind(response, q.value = p.adjust(response[, 5], method = 'BH'))
  response <- data.frame(response, stringsAsFactors = F)
  response[, 4:6] <- sapply(4:6, function(x) as.numeric(response[, x]))
  response <- response[order(response$signature, response$q.value), ]
  return(response)
  
})
drug_response_pan_sigs <- data.table::rbindlist(drug_response)

#Filtering
tau_thresh <- -0.07
drug_response_pan_filt <- drug_response_pan_sigs[drug_response_pan_sigs$tau < tau_thresh, ]
drug_response_pan_filt <- drug_response_pan_filt[drug_response_pan_filt$q.value < 0.05, ]

# Remove CX7 (seems to be an artefact)
drug_response_pan_filt <- drug_response_pan_filt[drug_response_pan_filt$signature != 7, ]

# Signature name
drug_response_pan_filt$signature <- paste0('CX', drug_response_pan_filt$signature)

#save
write.table(drug_response_pan_filt, file.path(output_dir, 'DrugResponse_PanCancer.txt'),
            sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)

# Number of significant hits for drug screen at pan-cancer level
print(table(drug_response_pan_filt$signature, drug_response_pan_filt$drug))


## PACLITAXEL PREDICTION: CX5 vs CX3
sigs <- colnames(cell_exp_filt_drug)
d <- 'paclitaxel'

pdat <- left_join(drug_screen[drug_screen$name == d, ],
                data.frame(depmap_id = rownames(cell_exp_filt_drug),
                                cell_exp_filt_drug),
                by = 'depmap_id')
pdat <- left_join(pdat, cell_line_info, by = 'depmap_id')
pdat <- pdat[pdat$passed_str_profiling.y == T, ]
pdat <- pdat[pdat$auc <= 1, ]

# linear model
lm <- lm(auc~CX3 + CX5 + CX3:CX5, pdat)
summary(lm)

# CX3 and CX5 present predictive signal as biomarkers for taxane response at pan-cancer level. 
# However, CX5 seems to have higher signal for predicting taxane response.
# We then aimed to create a prediction model and set a threhsold for classifying patients. 
# To do so, we first scaled signature activities of cell lines to the TCGA activities (reference cohort). 
# This aims to generalize the model for all clinical scenarios.

## TRAINING THE CLASSIFIER
# Load TCGA activities
TCGA.EXPS=readRDS(paste0(input_dir,"/4_Activities_CompendiumCINSigs_THRESH95_TCGA.rds"))
mTCGA = as.matrix(TCGA.EXPS[ , c("CX3", "CX5")])
# Generate a model for scaling
lModel = list(mean = attributes(scale(mTCGA))$`scaled:center`, 
              scale = attributes(scale(mTCGA))$`scaled:scale`)


# Scale Cell line activities
sc_cell_exp = sweep(sweep(cell_exp[,c("CX3", "CX5")], 2, lModel$mean, FUN = '-'), 2, lModel$scale, FUN = "/")
colnames(sc_cell_exp)=c("sCX3","sCX5")
sc_cell_exp=data.frame(depmap_id=rownames(sc_cell_exp),sc_cell_exp)

# Add scale values to pdat
pdat<-left_join(pdat,sc_cell_exp,by="depmap_id")
pdat<-pdat[!is.na(pdat$auc) & !is.na(pdat$CX5),]
pdat.filt<-pdat[pdat$auc<=quantile(pdat$auc,0.98),]

# Find optimal threshold for sCX5
thrs=seq(-1,1,0.1)

statistics=c()
for(thr in thrs){
  pdat.filt$pred=ifelse(pdat.filt$sCX5<thr, 'Resistant', 'Sensitive')
  tab=table(pdat.filt$pred)
  if(length(tab)!=1){
    # Compare AUC values
    n.res=tab[1]
    n.tot=sum(tab)
    mean.res=mean(pdat.filt$auc[pdat.filt$pred=="Resistant"])
    mean.sen=mean(pdat.filt$auc[pdat.filt$pred=="Sensitive"])
    test=t.test(pdat.filt$auc[pdat.filt$pred=="Resistant"],pdat.filt$auc[pdat.filt$pred=="Sensitive"])
    pval=test$p.value
    out=c(thr,n.tot,n.res,n.res/n.tot,mean.res,mean.sen,mean.res-mean.sen,pval)
    statistics=rbind(statistics,out)
  }
}
statistics=as.data.frame(statistics)
colnames(statistics)=c("thrs","n","n.res","f.res","mean.res","mean.sen","mean.diff","p.val")
statistics$signature="CX5"
statistics$significant=ifelse(statistics$p.val<0.1,"yes","no")


png(paste0(figs_dir,"/CellLines_Taxane_OptThresh.png"), width = 50, height = 50, units="mm", res=300)
p<-ggplot(statistics, aes(x=thrs, y=mean.res)) + 
  geom_vline(xintercept = -0.1, linetype = "dashed", colour = "grey70", size = 0.5) +
  geom_vline(xintercept = 0.7, linetype = "dashed", colour = "grey70", size = 0.5) +
  geom_point(size=0.8) +
  geom_line(size=0.5) +
  labs(title="", subtitle="") +
  scale_color_manual(values = c("black","red")) +
  ylab("Mean AUC values of cell lines\n predicted as resistant") +
  xlab("Activity threshold") +
  theme(text = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.line = element_line(size = 0.5), 
        axis.ticks = element_line(size = 0.5),
        axis.ticks.length = unit(.1, "cm"),
        legend.position="none")
print(p)
dev.off()
ggsave(file.path(figs_dir, "CellLines_Taxane_OptThresh.svg"),p,
       width = 50, height = 50, units = 'mm', dpi = 300)

# Classify cell lines with optimal threshold
pdat.filt$pred=ifelse(pdat.filt$sCX5<0, 'Resistant','Sensitive')

# Compare AUC values
test=t.test(pdat.filt$auc[pdat.filt$pred=="Resistant"],pdat.filt$auc[pdat.filt$pred=="Sensitive"])
pval=test$p.value


# Select resistant and sensitive cell lines
# Explore the distribution of AUC values & set threshold for classifying samples
p<-ggplot(pdat.filt) +
  geom_density(aes(x = auc, color=pred))+
  labs(title="IC50 distribution in all cell lines",
       subtitle=paste0("P-value t-test = ",round(pval,3)))+
  scale_x_continuous(breaks=seq(round(min(pdat$auc))-0.1,1,0.1), limits = c(0,1))+
  xlab("Response to paclitaxel (AUC)") +
  theme(text = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.5), 
        axis.ticks = element_line(size = 0.5),
        axis.ticks.length = unit(.1, "cm"),
        legend.position="right")
ggsave(file.path(figs_dir, paste0('Supplementary/CellLines_Taxane_AUC.svg')),
       width = 100, height = 80, units = 'mm', dpi = 300)
