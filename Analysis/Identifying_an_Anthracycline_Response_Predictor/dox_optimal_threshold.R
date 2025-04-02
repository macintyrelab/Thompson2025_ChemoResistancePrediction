##### OPTIMAL THRESHOLD FOR DOXORUBICIN RESPONSE PREDICTION #####

rm(list=ls(all=TRUE))

## Libraries
library(dplyr)
library(reshape2)
library(readr)
library(survival)
library(ggplot2)
library(ggthemes)
library(this.path)

## Plotting theme
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(text = element_text(size = 6),
             axis.text = element_text(size = 6))

## Paths
DIR <- dirname(this.path())
BASE <- dirname(dirname(DIR))
INDIR <- file.path(BASE, 'Input_Data')
FIGDIR <- file.path(BASE, 'Figures')


## Load experimental data
exps.org=as.data.frame(readRDS(paste0(INDIR, "/4_Activities_CompendiumCINSigs_THRESH95_OV_organoids.rds")))
exps.sph=as.data.frame(readRDS(paste0(INDIR, "/4_Activities_CompendiumCINSigs_THRESH95_OV_spheroids.rds")))


## Correct exposures => Apply thresholds computed for trusting in low activity values (already done)
# Organoids
corr.exps.org=exps.org
corr.exps.org$sample=row.names(corr.exps.org)

# Spheroids
corr.exps.sph = exps.sph
corr.exps.sph$`JBLAB ID`=substr(row.names(corr.exps.sph),1,nchar(row.names(corr.exps.sph))-7)

## Metadata for excluding samples from the same patient
meta=read_csv(paste0(INDIR,"/possible_ascites_list_V20180319.csv"))[,c(1,11)]
corr.exps.sph=left_join(corr.exps.sph,meta,by="JBLAB ID")
colnames(corr.exps.sph)[ncol(corr.exps.sph)]="sample"

map=read_csv(paste0(INDIR,"/spheroid_organoid_sample_overlap.csv"))
map$spheroid_id=gsub(".bwamem","",map$spheroid_id)
map=map[map$organoid_id%in%corr.exps.org$sample | map$spheroid_id%in%corr.exps.sph$`JBLAB ID`,]

## Add response data
# Organoids
response=cbind(sample=c("118976org","119127org","54327org","119127org","54288org","54276org","151723org","23868org",
                        "32077org","151773org"),
               response=c(rep("resistant",8),
                          rep("sensitive",2)))
response=as.data.frame(response)
data.org=left_join(response,corr.exps.org,by="sample")
colnames(map)[2]="sample"
map=left_join(map,data.org[,1:2],by="sample")

# Spheroids
response=cbind(sample=c("119120","119178","54289","118902","54327","80601","119136","119016","119025","118976","118947",
                        "54356","32072","80630","80720"),
               response=c(rep("resistant",11),
                          rep("sensitive",4)))
response=as.data.frame(response)
response$sample=as.numeric(response$sample)
data.sph=left_join(response,corr.exps.sph,by="sample")
colnames(map)[2:4]=c("organoid_id","JBLAB ID","organoid_resp")
map=left_join(map,data.sph[,c(20,2)],by="JBLAB ID")

## Threshold for signature response predictor
# For those patients with two samples, we selected spheroids since they are directly collected from patient
dupl=map$organoid_id[!is.na(map$organoid_resp) & !is.na(map$response)] #response is the same in organoids and spheroids
data.org=data.org[!data.org$sample%in%dupl,]
data.exp = rbind(data.org,data.sph[,1:19]) #23 samples

# searching space 
CX8.thrs=seq(0.001,0.02,0.001)
CX9.thrs=seq(0.001,0.02,0.001)
CX13.thrs=seq(0.001,0.02,0.001)

statistics=c()
for(thrs1 in 1:length(CX8.thrs)){
    for(thrs2 in 1:length(CX9.thrs)){
        for(thrs3 in 1:length(CX13.thrs)){
            data.exp$pan.prediction="pred.sensitive"
            data.exp$pan.prediction[data.exp$CX8>CX8.thrs[thrs1] | data.exp$CX9>CX9.thrs[thrs2] | data.exp$CX13>CX13.thrs[thrs3]]="pred.resistant"
            #statistics
            table=melt(table(data.exp$response,data.exp$pan.prediction))
            colnames(table)=c("response","pan.prediction","counts")
            sensitivity=table$counts[table$response=="sensitive" & table$pan.prediction=="pred.sensitive"]/
                (table$counts[table$response=="sensitive" & table$pan.prediction=="pred.sensitive"]+
                     table$counts[table$response=="sensitive" & table$pan.prediction=="pred.resistant"])
            specificity=table$counts[table$response=="resistant" & table$pan.prediction=="pred.resistant"]/
                (table$counts[table$response=="resistant" & table$pan.prediction=="pred.resistant"]+
                     table$counts[table$response=="resistant" & table$pan.prediction=="pred.sensitive"])
            youden = sensitivity + specificity - 1
            out=c(CX8.thrs[thrs1],CX9.thrs[thrs2],CX13.thrs[thrs3],sensitivity,specificity,youden)
            statistics=rbind(statistics,out)
        }
    }
}
statistics=as.data.frame(statistics)
colnames(statistics)=c("CX8.thrs","CX9.thrs","CX13.thrs","sensitivity","specificity","youden")

# best solutions => the aim is to maximize sensitivity!
optimal.thres.sens=statistics[statistics$sensitivity==max(statistics$sensitivity),]
optimal.thres.sens=optimal.thres.sens[optimal.thres.sens$specificity==max(optimal.thres.sens$specificity),]

optimal.thres.spec=statistics[statistics$specificity==max(statistics$specificity),]
optimal.thres.spec=optimal.thres.spec[optimal.thres.spec$sensitivity==max(optimal.thres.spec$sensitivity),]


# 17 CX sigs
#   max sensitivity => CX8 > 0.01 (middle number) & CX9/CX10 > 0.009 (middle number)
#   max specificity => CX8 = 0.005 & CX9/CX13 > 0.009

#### Test in tumours ####
meta=readRDS(paste0(INDIR,"/tailor_patient_plat+dox_data.rds"))
corr.exps=as.data.frame(readRDS(paste0(INDIR, "/4_Activities_CompendiumCINSigs_THRESH95_OV_tissue_old_50kb.rds"))) # Applied signature-specific thresholds 
corr.exps$primary_sample=substr(row.names(corr.exps),13,nchar(row.names(corr.exps)))
corr.exps$primary_sample=gsub("-Tseq-53","",corr.exps$primary_sample)
corr.exps$primary_sample=gsub("-Tseq-26","",corr.exps$primary_sample)
corr.exps$primary_sample=gsub("-Tseq-46","",corr.exps$primary_sample)
corr.exps$primary_sample=gsub("-Tseq-43","",corr.exps$primary_sample)
corr.exps$primary_sample=gsub("-Tseq-48","",corr.exps$primary_sample)
data.tum=left_join(meta,corr.exps,by="primary_sample")
data.tum=data.tum[!is.na(data.tum$PFS.response),]

# Which is the optimal in a independent set?
performance.sens=c()
for(thisOpt in 1:nrow(optimal.thres.sens)){
    data.tum$pan.prediction="pred.sensitive"
    data.tum$pan.prediction[data.tum$CX8>optimal.thres.sens[thisOpt,1] | data.tum$CX9>optimal.thres.sens[thisOpt,2] | data.tum$CX13>optimal.thres.sens[thisOpt,3]]="pred.resistant"
    #survival
    data.tum$pan.prediction=factor(data.tum$pan.prediction,levels = c("pred.sensitive","pred.resistant"))
    fit.coxph=coxph(Surv(PFS.response,progressed) ~ pan.prediction + pan.prediction:PFS_from_diagnosis, data = data.tum)
    out=c(optimal.thres.sens[thisOpt,1],optimal.thres.sens[thisOpt,2],optimal.thres.sens[thisOpt,3],summary(fit.coxph)$coefficients[1,])
    performance.sens=rbind(performance.sens,out)
}
performance.sens=as.data.frame(performance.sens)
colnames(performance.sens)=c("CX8.thrs","CX9.thrs","CX13.thrs","coef","exp(coef)","se(coef)","z","Pr(>|z|)")

performance.spec=c()
for(thisOpt in 1:nrow(optimal.thres.spec)){
    data.tum$pan.prediction="pred.sensitive"
    data.tum$pan.prediction[data.tum$CX8>optimal.thres.spec[thisOpt,1] | data.tum$CX9>optimal.thres.spec[thisOpt,2] | data.tum$CX13>optimal.thres.spec[thisOpt,3]]="pred.resistant"
    #survival
    data.tum$pan.prediction=factor(data.tum$pan.prediction,levels = c("pred.sensitive","pred.resistant"))
    fit.coxph=coxph(Surv(PFS.response,progressed) ~ pan.prediction + pan.prediction:PFS_from_diagnosis, data = data.tum)
    out=c(optimal.thres.spec[thisOpt,1],optimal.thres.spec[thisOpt,2],optimal.thres.spec[thisOpt,3],summary(fit.coxph)$coefficients[1,])
    performance.spec=rbind(performance.spec,out)
}
performance.spec=as.data.frame(performance.spec)
colnames(performance.spec)=c("CX8.thrs","CX9.thrs","CX13.thrs","coef","exp(coef)","se(coef)","z","Pr(>|z|)")

# 17 CX sigs
#   max sensitivity in organoids/spehroids is working fine => CX8 > 0.01 & CX9/CX10 > 0.009
#   max specificity in organoids/spehroids is working fine => CX8 = 0.005 & CX9/CX13 > 0.009


### CONCLUSION! 
# Optimal threshold is dependent of what performance metric we want to maximize
# In this case, we want to maximize sensitivity of the biomarker

STATS="Sensitivity"

#### Apply filters & show results ####
## Organoids
data.exp$pan.prediction="pred.sensitive"
if(STATS=="Sensitivity"){data.exp$pan.prediction[data.exp$CX8>0.01 | data.exp$CX9>0.009 | data.exp$CX13>0.009]="pred.resistant"}
if(STATS=="Specificity"){data.exp$pan.prediction[data.exp$CX8>0.005 | data.exp$CX9>0.009 | data.exp$CX13>0.009]="pred.resistant"}

table=melt(table(data.exp$response,data.exp$pan.prediction))
colnames(table)=c("response","pan.prediction","counts")

# Plot
p=ggplot(table %>% group_by(response) %>% 
             mutate(perc = round(counts/sum(counts),2)), 
         aes(x = response, y = perc, 
             fill = pan.prediction, cumulative = TRUE)) +
    geom_col() +
    labs(title="Organoids & Spheroids response to doxo", subtitle = "Observed vs Predicted") +
    geom_text(aes(label = counts), 
              position = position_stack(vjust = 0.5), size = 2) +
    theme(plot.title = element_text(size = 8),
          plot.subtitle = element_text(size = 7),
          axis.title = element_text(size = 7),
          axis.line = element_line(size = 0.5), 
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(.1, "cm"),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5, size = 6),
          axis.text.y = element_text(size = 6))
ggsave(paste0(FIGDIR,"/DoxResponse_Organoids&Spheroids_Obs_vs_Pred.svg"), p, width = 45, height = 45, units = "mm")

sensitivity=6/6
specificity=14/17

## Test this threshold in OV04
data.tum$pan.prediction="pred.sensitive"
data.tum$pan.prediction[data.tum$CX8>0.01 | data.tum$CX9>0.009 | data.tum$CX13>0.009]="pred.resistant"
data.tum$pan.prediction=factor(data.tum$pan.prediction,levels = c("pred.sensitive","pred.resistant"))
# survival
fit.coxph=coxph(Surv(PFS.response,progressed) ~ pan.prediction + pan.prediction*PFS_from_diagnosis, data = data.tum) #corrected by plat_PFS
summary(fit.coxph)

# We can see the prediction capacity even not including all confounding covariates!

