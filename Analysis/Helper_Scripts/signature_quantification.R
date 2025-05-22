## Workflow for quantifying our compendium of CIN signatures from copy number profiles


######################## MODIFY THIS AS NEEDED!!! ###############################

### INPUT DIRECTORY
INPDIR="Thompson2025_ChemoResistancePrediction/Input_Data" # Where your input file is located?

### INPUT FILE NAME
# Input file with segmented absolute copy number values (QDNAseq object or segTab)
# If the input is a segTab, the segments per sample should be ordered by chromosome, start and end 
# The input file should be a .rds file
ABS.CN="TSO500_OV04_absoluteCN.rds" 

### OUTPUT DIRECTORY
OUTDIR="Thompson2025_ChemoResistancePrediction/Input_Data" # Where do you want to save output files?

### ID NAME OF COHORT/EXPERIMENT/PROJECT
OUTNAME="OV_tissue_TSO500" # Name to identify the output files

### PARAMETERS
CORES=1 # You can add more if you have a linux local machine
SAVEINTERMEDIATE=FALSE # Add TRUE in case you want to save all intermediate files
REMOVEQUIET=TRUE # We cannot quantify activities of samples with <20 CNAs

################################################################################


### LIBRARIES
library(this.path)
library(data.table)
library(QDNAseq)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(tidyverse)
library(foreach)
library(doMC)
library(ComplexHeatmap)
library(YAPSA)


### PATHS
BASE=dirname(this.path())


### FUNCTIONS
source(file.path(BASE, "sig_scripts/main_functions.R"))


########################### WORKFLOW ###########################################

#### Step 1: Process input data ################################################
ABS.CN=readRDS(file.path(INPDIR,ABS.CN))

### Step 1.1: Get segTables with copy number profiles
CN.TABLE=c()
if(is(ABS.CN,"QDNAseqCopyNumbers")==TRUE){
    samples<-colnames(ABS.CN)
    #get copy number profiles
    for(sample in samples){
        tab<-as.data.frame(CNpare::getSegTable(ABS.CN[,sample]))
        CN.TABLE=rbind(CN.TABLE,cbind(tab,sample=rep(sample,nrow(tab))))
    }
    CN.TABLE[,2:4]=sapply(2:4,function(x)as.numeric(CN.TABLE[,x]))
    tab=NULL
}

if(is(ABS.CN,"data.frame")==TRUE){
    if(ncol(ABS.CN)!=5){
        print("Format of input data is not correct!")
        CN.TABLE=NULL
    }
    CN.TABLE=as.data.frame(ABS.CN)
    CN.TABLE[,2:4]=sapply(2:4,function(x)as.numeric(CN.TABLE[,x]))
}

### Step 1.2: Smoothing normal segments
WIGGLE=0.1
IGNOREDELS=FALSE

# Set everything very close to 2 to 2
CN.TABLE$segVal[CN.TABLE$segVal > (2-WIGGLE) & CN.TABLE$segVal < (2+WIGGLE)] = 2
# Merge segments only when two normal follow each other -> SMOOTHINGFACTOR = 0
CN.TABLE = idSmoothingTargets(CN.TABLE, WIGGLE = 0, colNameSegVal = "segVal", colNameChr = "chromosome", IGNOREDELS = IGNOREDELS)
# Split by sample name
lRaw = split(CN.TABLE, CN.TABLE$sample)
# Smooth segments by taking the weighted average of the segVal and their lengths
lCN.TABLE = smoothSegments(lRaw, CORES=1, SMOOTHINGFACTOR = 0, colNameMerge = "segVal", colNameChr = "chromosome",
                           colNameStart = "start", colNameEnd = "end", IGNOREDELS = IGNOREDELS, asDf = FALSE)

CN.TABLE = as.data.table(rbindlist(lCN.TABLE))


### Step 1.3: Remove quiet samples (<20 CNAs)
# These samples do not have enough copy number alterations to quantify signatures

# Identify CNAs per sample
dtCNAs = CN.TABLE[CN.TABLE$segVal != 2, ]
quietSamples = names(table(dtCNAs$sample))[table(dtCNAs$sample) < 20]

# Identify those samples with all segVal=2 => filtered out 
allSamples = unique(CN.TABLE$sample)
normalSamples = allSamples[!allSamples%in%dtCNAs$sample]
quietSamples = c(quietSamples,normalSamples)

# Saving info of detectable CIN
if(length(quietSamples)!=0){
    dfQuiets=as.data.frame(rbind(cbind(sample=quietSamples,CIN="NO"),cbind(sample=unique(CN.TABLE$sample),CIN="YES")))
    if(SAVEINTERMEDIATE==TRUE){saveRDS(dfQuiets, paste0(OUTDIR, "/1_CINGroup_",OUTNAME,".rds"))}
} else {
    dfQuiets=as.data.frame(cbind(sample=unique(CN.TABLE$sample),CIN="YES"))
    if(SAVEINTERMEDIATE==TRUE){saveRDS(dfQuiets, paste0(OUTDIR, "/1_CINGroup_",OUTNAME,".rds"))}
}
dfQuiets=NULL

if(REMOVEQUIET) {
    # Remove quiet samples
    CN.TABLE = CN.TABLE[!CN.TABLE$sample %in% quietSamples, ]
}

# Saving segTable
if(SAVEINTERMEDIATE==TRUE){saveRDS(CN.TABLE, paste0(OUTDIR, "/1_SegTables_",OUTNAME,".rds"))}


#### Step 2: Extract copy number features ######################################
PREPATH=paste0(BASE, "/refgenome/")
RMNORM=TRUE
lProfiles=split(CN.TABLE,factor(CN.TABLE$sample))

lECNF=extractCopynumberFeatures.2024(lProfiles, cores=CORES, prePath=PREPATH, rmNorm=RMNORM)
if(SAVEINTERMEDIATE==TRUE){saveRDS(lECNF, paste0(OUTDIR, "/2_ECNF_",OUTNAME,".rds"))}

#### Step 3: Compute the sum-of-posterior imput matrix #########################
# Mixture models
MODELS=readRDS(paste0(BASE, "/data/Mixmodels_merged_components.rds"))

# Probability of each segment to belong to each category & sum up per sample
UNINFPRIOR="TRUE"
allFeatures=names(MODELS)

lMats = lapply(allFeatures, function(thisFeature) {
    
    print(thisFeature)
    thisEcnf = lECNF[[thisFeature]]
    thisModel = as.data.frame(MODELS[[thisFeature]])
    
    dat = as.numeric(thisEcnf[,2] )
    # We want a posterior, hence likelihood (density) times prior (weight)
    if( ncol(thisModel) == 2 ) {
        # Poisson model
        print("Poisson-based posterior")
        if(UNINFPRIOR){
            postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) )
        } else {
            postDatUnscaled = sapply(1:nrow(thisModel), function(x) dpois(x = dat, lambda = thisModel[[x,"Mean"]]) * thisModel[[x, "Weight"]] )
        }
        
    } else {
        # Gaussian model
        print("Gauss-based posterior")
        if(UNINFPRIOR){
            postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) )
        } else {
            postDatUnscaled = sapply(1:nrow(thisModel), function(x) dnorm(x = dat, mean = thisModel[[x,"Mean"]], sd = thisModel[[x,"SD"]]) * thisModel[[x, "Weight"]] )
        }
    }
    
    # Normalise densities to probabilities
    postDatScaled = data.frame(postDatUnscaled / rowSums(postDatUnscaled))
    # Make sure that focal amplifications with huge number of copies have weight in the chp=10 
    id = which(rowSums(is.na(postDatScaled)) == ncol(postDatScaled))
    if (thisFeature == "changepoint"){
        if(length(id)!=0){postDatScaled[id,]=matrix(rep(c(rep(0,9),1), length(id)), ncol = 10, nrow =length(id), byrow = T)}
    }
    # Create component matrix per sample
    postDatScaled$Sample = thisEcnf[,1]
    matSxC = aggregate(. ~ Sample, postDatScaled, sum)
    rownames(matSxC) = matSxC$Sample
    matSxC$Sample = NULL
    matSxC = as.matrix(matSxC)
    
    # Should be sorted but just to be sure
    matSxC = matSxC[ ,order(thisModel[,"Mean"])]
    colnames(matSxC) = paste0(thisFeature, 1:ncol(matSxC))
    
    return(matSxC)
})
SxC = do.call(cbind, lMats)
if(SAVEINTERMEDIATE==TRUE){saveRDS(SxC, paste0(OUTDIR, "/3_SxC_",OUTNAME,".rds"))}


#### Step 4: Quantify copy number signatures ###################################
# Definitions of our compendium of signatures 
SIGS.DEFINITIONS=readRDS(paste0(BASE, "/data/Signature_Compendium_v5_Cosine-0.74_Signatures_NAMESAPRIL21.rds"))

# Quantify 17 CIN signatures
V = t(SxC)
W = SIGS.DEFINITIONS

# Sanity check signature matrix
if(nrow(W) < ncol(W)) W = t(W)

# Check order of components and fix if necessary
if(!identical(rownames(W), rownames(V))) {
    W = W[match(rownames(V), rownames(W)), ]
}

### YAPSA needs:
## Full matrix V        mutCatalogue        components (rows) by samples (cols)    <= HAVE
## Left matrix W        sigCatalogue        components (rows) by signature (cols)   <= HAVE
## Right matrix H       expCatalogue        signature (rows) by samples (cols)     <= WANT
Hraw = as.matrix(LCD(in_mutation_catalogue_df = V, in_signatures_df = W, in_per_sample_cutoff = 0))
H = t(apply(Hraw, 2, function(x) x/sum(x)))
if(SAVEINTERMEDIATE==TRUE){saveRDS(H, file = paste0(OUTDIR, "/4_Activities_CompendiumCINSigs_",OUTNAME,".rds"))}

# Apply signature-specific thresholds
thresh=readRDS(paste0(BASE,"/data/Jul2024_TCGA_Signature_Thresholds.rds"))
EXPS.CORRECTED = do.call(cbind,lapply(1:ncol(H), function(thisSig){
    acts=H[,thisSig]
    thrs=thresh[thisSig]
    acts[acts<thrs]=0
    return(acts)
}))
colnames(EXPS.CORRECTED)=colnames(H)
EXPS.CORRECTED=t(apply(EXPS.CORRECTED,1,function(x) x/sum(x)))

# Save
saveRDS(EXPS.CORRECTED, file = paste0(OUTDIR, "/4_Activities_CompendiumCINSigs_THRESH95_",OUTNAME,".rds"))
