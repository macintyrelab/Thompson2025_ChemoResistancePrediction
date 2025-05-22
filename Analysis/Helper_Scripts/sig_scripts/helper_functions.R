
getSegsize<-function(abs_profiles, rmNorm = FALSE)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        # First tap.
        segTab$segVal = as.numeric(segTab$segVal)
        # If wished, don't consider normal segments
        if(rmNorm) { segTab = segTab[ segTab$segVal != 2, ] }
        # Avoiding potential artefact
        segTab$segVal[segTab$segVal<0]<-0
        seglen<-segTab$end-segTab$start
        seglen<-seglen[seglen>0]
        # Double tap.
        out<-rbind(out,cbind(ID=rep(i,length(seglen)),value=as.numeric(seglen)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}

getBPnum<-function(abs_profiles,chrlen, SIZE = 10000000)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        chrs<-unique(segTab$chromosome)
        allBPnum<-c()
        for(c in chrs)
        {
            currseg<-segTab[segTab$chromosome==c,]
            intervals<-seq(1,chrlen[chrlen[,1]==paste0("chr",c),2]+SIZE,SIZE)
            res <- hist(as.numeric(currseg$end[-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
            allBPnum<-c(allBPnum,res)
        }
        # Make sure it's really numeric
        out<-rbind(out,cbind(ID=rep(i,length(allBPnum)),value=as.numeric(allBPnum)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


getOscilation.2024<-function(abs_profiles,chrlen)
{
  out<-c()
  samps<-getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    chrs<-unique(segTab$chromosome)
    oscCounts<-c()
    for(c in chrs)
    {
      currseg<-as.numeric(segTab$segVal[segTab$chromosome==c])
      currseg<-round(as.numeric(currseg))
      if(length(currseg)>3)
      {
        prevval<-currseg[1]
        count=0
        for(j in 3:length(currseg))
        {
          if(currseg[j]==prevval&currseg[j]!=currseg[j-1])
          {
            count<-count+1
            #fixing issue #22
            if (j == length(currseg)) {
              oscCounts = c(oscCounts, count)
              count = 0
            }
          }else{
            oscCounts<-c(oscCounts,count)
            count=0
          }
          prevval<-currseg[j-1]
        }
      }
    }
    # Make sure it's really numeric
    out<-rbind(out,cbind(ID=rep(i,length(oscCounts)),value=as.numeric(oscCounts)))
    if(length(oscCounts)==0)
    {
      out<-rbind(out,cbind(ID=i,value=0))
    }
  }
  rownames(out)<-NULL
  data.frame(out,stringsAsFactors = F)
}

getCentromereDistCounts<-function(abs_profiles,centromeres,chrlen)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }else
        {
            segTab<-as.data.frame(abs_profiles[[i]])
            colnames(segTab)[4]<-"segVal"
        }
        chrs<-unique(segTab$chromosome)
        all_dists<-c()
        for(c in chrs)
        {
            if(nrow(segTab)>1)
            {
                starts<-as.numeric(segTab[segTab$chromosome==c,2])[-1]
                segstart<-as.numeric(segTab[segTab$chromosome==c,2])[1]
                ends<-as.numeric(segTab[segTab$chromosome==c,3])
                segend<-ends[length(ends)]
                ends<-ends[-length(ends)]
                centstart<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,3])
                centend<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,4])
                chrend<-chrlen[substr(chrlen[,1],4,5)==c,2]
                ndist<-cbind(rep(NA,length(starts)),rep(NA,length(starts)))
                ndist[starts<=centstart,1]<-(centstart-starts[starts<=centstart])/(centstart-segstart)*-1
                ndist[starts>=centend,1]<-(starts[starts>=centend]-centend)/(segend-centend)
                ndist[ends<=centstart,2]<-(centstart-ends[ends<=centstart])/(centstart-segstart)*-1
                ndist[ends>=centend,2]<-(ends[ends>=centend]-centend)/(segend-centend)
                ndist<-apply(ndist,1,min)
                
                all_dists<-rbind(all_dists,sum(ndist>0))
                all_dists<-rbind(all_dists,sum(ndist<=0))
            }
        }
        if(nrow(all_dists)>0)
        {
            # Make sure it's really numeric
            out<-rbind(out,cbind(ID=i,ct1=as.numeric(all_dists[,1])))
        }
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


getChangepointCN<-function(abs_profiles, allowedError = 0.1, rmNorm = FALSE)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers") {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        } else {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal = as.numeric(segTab$segVal)
        segTab$segVal[segTab$segVal<0]<-0
        chrs<-unique(segTab$chromosome)
        allcp<-c()
        for(c in chrs) {
            currseg<-as.numeric(segTab$segVal[segTab$chromosome==c])
            firstSeg = abs(2 - currseg[1] )
            # As we look only at the left end of a CNA, we might miss a changepoint at the beginning of the p-arm
            # That's why we check manually but only regard this value if it is higher than an allowed error rate.
            if(firstSeg <= allowedError) {
                theseChanges = abs(currseg[-1]-currseg[-length(currseg)])
                if(rmNorm) { theseChanges = theseChanges[ currseg[-1] != 2 ] }
                allcp<-c(allcp, theseChanges)
            } else {
                theseChanges = c( firstSeg, abs(currseg[-1]-currseg[-length(currseg)]) )
                if(rmNorm) { theseChanges = theseChanges[ currseg != 2 ] }
                allcp<-c(allcp, theseChanges)
            }
            
        }
        if(length(allcp)==0) {
            allcp<-0 #if there are no changepoints
        }
        # Make sure it's really numeric
        out<-rbind(out,cbind(ID=rep(i,length(allcp)),value=as.numeric(allcp)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


getCN<-function(abs_profiles, rmNorm = FALSE)
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal[as.numeric(segTab$segVal)<0]<-0
        # If wished, don't consider normal segments
        if(rmNorm) { segTab = segTab[ segTab$segVal != 2, ] }
        cn<-as.numeric(segTab$segVal)
        # Make sure it's really numeric.
        out<-rbind(out,cbind(ID=rep(i,length(cn)),value=as.numeric(cn)))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


getSampNames<-function(abs_profiles)
{
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
        samps<-colnames(abs_profiles)
    }
    else
    {
        samps<-names(abs_profiles)
    }
    samps
}


getSegTable<-function(x)
{
    dat<-x
    sn<-Biobase::assayDataElement(dat,"segmented")
    fd <- Biobase::fData(dat)
    fd$use -> use
    fdfiltfull<-fd[use,]
    sn<-sn[use,]
    segTable<-c()
    for(c in unique(fdfiltfull$chromosome))
    {
        snfilt<-sn[fdfiltfull$chromosome==c]
        fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
        sn.rle<-rle(snfilt)
        starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
        ends <- cumsum(sn.rle$lengths)
        lapply(1:length(sn.rle$lengths), function(s) {
            from <- fdfilt$start[starts[s]]
            to <- fdfilt$end[ends[s]]
            segValue <- sn.rle$value[s]
            c(fdfilt$chromosome[starts[s]], from, to, segValue)
        }) -> segtmp
        segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
        segTable<-rbind(segTable,segTableRaw)
    }
    colnames(segTable) <- c("chromosome", "start", "end", "segVal")
    segTable
}


getPloidy<-function(abs_profiles)
{
  out<-c()
  samps<-getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }
    else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    segLen<-(as.numeric(segTab$end)-as.numeric(segTab$start))
    ploidy<-sum((segLen/sum(segLen))*as.numeric(segTab$segVal))
    out<-c(out,ploidy)
  }
  data.frame(out,stringsAsFactors = F)
}

