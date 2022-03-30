options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)
logit<-function(xx)log(xx)-log(1-xx)
invLogit<-function(xx)1/(1+exp(-xx))
meanCrI<-function(xx,quants=c(.025,.975))c(mean(xx),quantile(xx,quants))


deer<-read.csv('WF swab samples 1-18-22 version final deer_fixedWithMissingCounties.csv')
deer<-deer[!is.na(deer$WF.ID..)&deer$WF.ID..!=''&deer$PCR.date!='duplicate',]
deer$pos<-deer$N1.Ct!='Undetermined'  | deer$N2.Ct!='Undetermined'
if(any(is.na(deer$County.Location)|deer$County.Location==''))stop('Deer with missing counties')


runCountyStan<-function(countTab,adjacencyTab,mod,nChain=50,nIter=10000){
  if(any(colnames(adjacencyTab)!=rownames(adjacencyTab))||(rownames(countTab)!=rownames(adjacencyTab)))stop('Mismatch between counts and adjacency')
  dat<-list(
    nCounty=nrow(countTab),
    counts=apply(countTab,1,sum),
    positives=countTab[,'TRUE'],
    adjacency=adjacencyTab
  )
  stan_sample <- rstan::sampling(
    mod,
    data=dat,
    iter=nIter,
    chains=nChain,
    thin=2,
    control=list(max_treedepth=15) #adapt_delta=.99
    #pars=c('means'),
    #include=FALSE
  )
  return(list('stan'=stan_sample,'dat'=dat,'mod'=mod,'tab'=countTab))
}

#https://www.census.gov/geographies/reference-files/2010/geo/county-adjacency.html
adj<-read.table('county_adjacency.txt',sep='\t')
colnames(adj)<-c('county','id','neighbor','neighborId')
adj$county<-dnar::fillDown(adj$county)
adj<-adj[grepl(', PA$',adj$county)&grepl(', PA$',adj$neighbor),]
adj$simple<-sub(' County, PA','',adj$county)
adj$neighborSimple<-sub(' County, PA','',adj$neighbor)
if(any(!deer$County.Location %in% adj$simple))stop('Unknown county')
adjacencyAll<-table(adj$simple,adj$neighborSimple)
diag(adjacencyAll)<-0
adj<-adj[adj$simple %in% deer$County.Location &adj$neighborSimple %in% deer$County.Location,]
adjacency<-table(adj$simple,adj$neighborSimple)
diag(adjacency)<-0

penn <- maps::map("county","Pennsylvania",plot=FALSE)
penn$prettyName<-sub('pennsylvania,','',penn$name)
substring(penn$prettyName,1,1)<-toupper(substring(penn$prettyName,1,1))
if(any(!colnames(props) %in% penn$names))stop('Unknown county name')
counts<-table(deer$County.Location,deer$pos)
allCounts<-rbind(counts,matrix(0,nrow=length(penn$names)-nrow(counts),ncol=2,dimnames=list(penn$prettyName[!penn$prettyName %in% rownames(counts)],c('FALSE','TRUE'))))
allCounts<-allCounts[penn$prettyName,]

mod <- rstan::stan_model("counties.stan")
stan<-runCountyStan(counts,adjacency,mod,nIter=3000)
matNoAdj<-as.matrix(stan$stan)
propsNoAdj<-apply(invLogit(matNoAdj[,'overallProp']+matNoAdj[,grep('^countyProp\\[',colnames(matNoAdj))]),2,meanCrI)
colnames(propsNoAdj)<-sprintf('pennsylvania,%s',tolower(rownames(adjacency)))

mod2 <- rstan::stan_model("counties_adjacent.stan")
stan2<-runCountyStan(counts,adjacency,mod2,nIter=20000)
stan3<-runCountyStan(allCounts,adjacencyAll,mod2,nIter=20000)
mat<-as.matrix(stan2$stan)
props<-apply(invLogit(mat[,'overallProp']+mat[,grep('^countyProp\\[',colnames(mat))]),2,meanCrI)
colnames(props)<-sprintf('pennsylvania,%s',tolower(rownames(adjacency)))

matAll<-as.matrix(stan3$stan)
propsAll<-apply(invLogit(matAll[,'overallProp']+matAll[,grep('^countyProp\\[',colnames(matAll))]),2,meanCrI)
colnames(propsAll)<-sprintf('pennsylvania,%s',tolower(rownames(adjacencyAll)))

cols<-1+penn$names %in% colnames(props)
breaks<-seq(0,ceiling(max(props[1,])*10)/10,.02)
cuts<-cut(props[1,],breaks)
cuts2<-cut(props[2,],breaks)
cuts3<-cut(propsAll[1,],breaks)
cutsNoAdj<-cut(propsNoAdj[1,],breaks)
propCol<-structure(dnar::rainbow.lab(length(levels(cuts))),.Names=levels(cuts))
cols<-structure(propCol[cuts],.Names=colnames(props))[penn$names]
cols2<-structure(propCol[cuts2],.Names=colnames(props))[penn$names]
cols3<-structure(propCol[cuts3],.Names=colnames(propsAll))[penn$names]
colsNoAdj<-structure(propCol[cutsNoAdj],.Names=colnames(props))[penn$names]
colsBase<-structure(ifelse(counts[,2]>0,propCol[length(propCol)],propCol[1]),.Names=colnames(props))[penn$names]
#cols[is.na(cols)]<-'grey'
#names(cols)<-penn$names
labels<-ifelse(rowSums(allCounts)>0,sprintf('%s\n%d/%d',rownames(allCounts),allCounts[,'TRUE'],rowSums(allCounts)),'')

pdf('test.pdf')
  #positive/negative coloring
  maps::map("county","Pennsylvania",col=colsBase,fill=TRUE)
  maps::map.text("county","Pennsylvania",add=TRUE,label=labels)
  #simple model
  maps::map("county","Pennsylvania",col=colsNoAdj,fill=TRUE)
  maps::map.text("county","Pennsylvania",add=TRUE,label=labels)
  dnar::insetScale(breaks,propCol,main='Estimated proportion positive',insetPos = c(0.025, 0.035, 0.04, 0.3))
  #adjacency model
  maps::map("county","Pennsylvania",col=cols,fill=TRUE)
  maps::map.text("county","Pennsylvania",add=TRUE,label=labels)
  dnar::insetScale(breaks,propCol,main='Estimated proportion positive',insetPos = c(0.025, 0.035, 0.04, 0.3))
  #lower credible interval
  maps::map("county","Pennsylvania",col=cols2,fill=TRUE)
  maps::map.text("county","Pennsylvania",add=TRUE,label=labels)
  dnar::insetScale(breaks,propCol,main='Estimated lower credible',insetPos = c(0.025, 0.035, 0.04, 0.3))
  #all estimates
  maps::map("county","Pennsylvania",col=cols3,fill=TRUE)
  maps::map.text("county","Pennsylvania",add=TRUE,label=labels)
  dnar::insetScale(breaks,propCol,main='Estimated proportion positive',insetPos = c(0.025, 0.035, 0.04, 0.3))
dev.off()



