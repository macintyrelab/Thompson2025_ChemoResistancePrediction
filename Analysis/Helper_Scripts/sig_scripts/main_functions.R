library("NMF")
library("quadprog")
#library("ComplexHeatmap")
#library("circlize")
# for the manual parallelisation and also used by NMF for their internal parallelisation
library("foreach")
# needed for faster nmf runs based on the developer's own recommendations:
# https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
# page 16 ff.
library("doParallel")
# Libraries "bigmemory" and "synchronicity" are recommended for non-Windows machines
# for remodelling data structures
library("reshape2")
# for plotting
library("ggplot2")
#library("Cairo")
# Correlation of signatures
library("nnls")
library("corrplot")


locationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
    this.file = NULL
    # This file may be 'sourced'
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }
    
    if (!is.null(this.file)) return(dirname(this.file))
    
    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
    
    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))
    
    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}
this_path<-locationOfThisScript()
source(paste(this_path,"helper_functions.R",sep="/"))


extractCopynumberFeatures.2024<-function(CN_data, cores = 1, prePath="data/", allowedError = 0.1, rmNorm = FALSE)
{
    #get chromosome lengths
    chrlen<-read.table(paste0(prePath, "hg19.chrom.sizes.txt"),sep="\t",stringsAsFactors = F)[1:24,]
    
    #get centromere locations
    gaps<-read.table(paste0(prePath, "gap_hg19.txt"),sep="\t",header=F,stringsAsFactors = F)
    centromeres<-gaps[gaps[,8]=="centromere",]
    
    if(cores > 1) {
        require(foreach)
        require(doMC)
        
        registerDoMC(cores)
        
        temp_list = foreach(i=1:6) %dopar% {
            if(i == 1){
                list(segsize = getSegsize(CN_data, rmNorm = rmNorm))
            } else if (i == 2) {
                list(bp10MB = getBPnum(CN_data,chrlen))
            } else if (i == 3) {
                list(osCN = getOscilation.2024(CN_data,chrlen))
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCounts(CN_data,centromeres,chrlen))
            } else if (i == 5) {
                list(changepoint = getChangepointCN(CN_data, allowedError, rmNorm = rmNorm))
            } else {
                list(copynumber = getCN(CN_data, rmNorm = rmNorm))
            }
            
        }
        
        # Another failsafe that the outcome is definitely numeric
        temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(temp_list, function(thisDF) { 
            thisDF[,2] = as.numeric(thisDF[,2])
            return(thisDF)
        })
        return( outList )
        
    } else {  
        
        segsize<-getSegsize(CN_data, rmNorm = rmNorm)
        bp10MB<-getBPnum(CN_data,chrlen)
        osCN<-getOscilation.2024(CN_data,chrlen)
        bpchrarm<-getCentromereDistCounts(CN_data,centromeres,chrlen)
        changepoint<-getChangepointCN(CN_data, allowedError, rmNorm = rmNorm)
        copynumber<-getCN(CN_data, rmNorm = rmNorm)
        
        temp_list = list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
        #temp_list = unlist( temp_list, recursive = FALSE )
        outList = lapply(names(temp_list), function(thisDF) { 
            DF=temp_list[[thisDF]]
            DF[,2] = as.numeric(DF[,2])
            return(DF)
        })
        names(outList)=names(temp_list)
        return(outList)
    }
}


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
