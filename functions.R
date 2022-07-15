
logit<-function(xx)log(xx)-log(1-xx)
invLogit<-function(xx)1/(1+exp(-xx))
meanCrI<-function(xx,quants=c(.025,.975))c(mean(xx),quantile(xx,quants))

runCountyStan<-function(countTab,adjacencyTab,mod,regions=rep(1,nrow(countTab)),nChain=50,nIter=10000){
  if(any(colnames(adjacencyTab)!=rownames(adjacencyTab))||any(rownames(countTab)!=rownames(adjacencyTab)))stop('Mismatch between counts and adjacency')
  regionId<-structure(1:length(unique(regions)),.Names=sort(unique(regions)))
  dat<-list(
    nCounty=nrow(countTab),
    nRegion=length(unique(regions)),
    regions=regionId[regions],
    counts=apply(countTab,1,sum),
    positives=countTab[,'TRUE'],
    adjacency=adjacencyTab
  )
  stan_sample <- rstan::sampling(
    mod,
    data=dat,
    iter=nIter,
    chains=nChain,
    thin=2
  )
  return(list('stan'=stan_sample,'dat'=dat,'mod'=mod,'tab'=countTab,'region'=regionId))
}

insetScale<-function(breaks,col,insetPos=c(.025,.015,.04,.25),main='',offset=1e-3,at=NULL,labels=NULL,cex=1,labXOffset=0,labYOffset=0){
  if(length(breaks)!=length(col)+1)stop('Number of breaks must be one more than colors')
  insetPos<-c(graphics::grconvertY(insetPos[1],'nfc','user'),graphics::grconvertX(insetPos[2],'nfc','user'),graphics::grconvertY(insetPos[3],'nfc','user'),graphics::grconvertX(insetPos[4],'nfc','user'))
  breakPos<-((breaks)-(min(breaks)))/max((breaks)-(min(breaks)))*(insetPos[4]-insetPos[2])+insetPos[2]
  #add a bit of offset to avoid pdf viewers displaying breaks between exact rectangle border meeting
  offsetPos<-breakPos[-1]+c(rep(offset*diff(range(breakPos)),length(breakPos)-2),0)
  graphics::rect(breakPos[-length(breakPos)],insetPos[1],offsetPos,insetPos[3],col=col,xpd=NA,border=NA)
  graphics::rect(insetPos[2],insetPos[1],insetPos[4],insetPos[3],xpd=NA)
  if(is.null(at)){
    at<-pretty(breaks)
    at<-at[at<=max(breaks)&at>=min(breaks)]
  }
  if(is.null(labels))labels<-at
  convertPos<-(at-(min(breaks)))/((max(breaks))-(min(breaks)))*(insetPos[4]-insetPos[2])+insetPos[2]
  graphics::segments(convertPos,insetPos[1],convertPos,insetPos[1]-diff(insetPos[c(1,3)])*.1,xpd=NA)
  graphics::text(convertPos+labXOffset*diff(insetPos[c(2,4)]),insetPos[1]-diff(insetPos[c(1,3)])*.175+labYOffset*diff(insetPos[c(1,3)]),labels,xpd=NA,adj=c(.5,1),cex=.85*cex)
  graphics::text(mean(insetPos[c(2,4)]),insetPos[3]+diff(insetPos[c(1,3)])*.45,main,xpd=NA,adj=c(.5,0),cex=cex)
  invisible(NULL)
}

fillDown<-function(x,emptyStrings=c(NA,''),errorIfFirstEmpty=TRUE){
  #depending on %in% to catch NAs if necessary
  isEmpty<-x %in% emptyStrings
  if(isEmpty[1]&errorIfFirstEmpty)stop(simpleError('First value empty'))
  #if first is empty and we don't want errors then have to just fill down from it anyway
  isEmpty[1]<-FALSE
  ids<-1:length(x)
  ids[isEmpty]<-0
  ids<-cummax(ids)
  return(x[ids])
}


