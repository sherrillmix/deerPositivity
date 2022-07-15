options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)
source('functions.R')

deer<-read.csv('TableS2.csv',skip=1)
deer$pos<-deer$N1.Ct!='Undetermined' | deer$N2.Ct!='Undetermined'
deer$Region<-trimws(deer$Region)
if(any(is.na(deer$County.Location)|deer$County.Location==''))stop('Deer with missing counties')

#https://www.census.gov/geographies/reference-files/2010/geo/county-adjacency.html
adj<-read.table('county_adjacency.txt',sep='\t')
colnames(adj)<-c('county','id','neighbor','neighborId')
adj$county<-fillDown(adj$county)
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
counts<-table(deer$County.Location,deer$pos)

#https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
mod <- rstan::stan_model("counties_adjacent.stan")
stan<-runCountyStan(counts,adjacency,mod,nChain=8,nIter=3000)
mat<-as.matrix(stan$stan)
props<-apply(invLogit(mat[,'overallProp']+mat[,grep('^countyProp\\[',colnames(mat))]),2,meanCrI)
colnames(props)<-sprintf('pennsylvania,%s',tolower(rownames(adjacency)))
if(any(!colnames(props) %in% penn$names))stop('Unknown county name')

cols<-1+penn$names %in% colnames(props)
breaks<-seq(0,ceiling(max(props[1,])*10)/10,.005)
cuts<-cut(props[1,],breaks)
cuts2<-cut(props[2,],breaks)
nCut<-length(levels(cuts))
nCutTop<-0;nCutBottom<-150
propCol<-structure(tail(head(rev(colorRampPalette(c(viridis::rocket(30),'white'),space='Lab')(nCut+nCutTop+nCutBottom)),nCut+nCutTop),nCut),.Names=levels(cuts))
cols<-structure(propCol[cuts],.Names=colnames(props))[penn$names]
labels<-rep('',length(penn$prettyName))
select<-penn$prettyName %in% rownames(counts)
labels[select]<-sprintf('%s\n%d/%d',penn$prettyName[select],counts[penn$prettyName[select],'TRUE'],rowSums(counts[penn$prettyName[select],]))
pdf('countyPositivity.pdf',height=8,width=12)
  naCol<-'#EAEAEA'
  maps::map("county","Pennsylvania",col=ifelse(is.na(cols),naCol,cols),fill=TRUE)
  maps::map.text("county","Pennsylvania",add=TRUE,label=labels)
  insetScale(breaks+rep(c(0,.0001),c(length(breaks)-1,1)),propCol,main='Estimated proportion positive',insetPos = c(0.025, 0.035, 0.04, 0.3))
dev.off()



