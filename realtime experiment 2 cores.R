library("parallel")
library("foreach")
library("doParallel")

cl <- makeCluster(detectCores() )
registerDoParallel(cl, cores = detectCores() )
data = foreach(i = 1:length(filenames), .packages = c("ncdf","chron","stats"),
               .combine = rbind) %dopar% {
                       try({
                               
                               # make sure, the  right columns are addressed (same ordering in overview and df)
                               df.unrevised=df.unrevised[,as.character(overview.nr[,1])]
                               variables.tlag=overview.nr[overview.nr$Publication.Lag!=0,'Abbreviation']
                               for (var.tlag in variables.tlag){
                                       df.unrevised[,var.tlag]=lag.exact(df.unrevised[,var.tlag,drop=F]
                                                                         ,overview.nr[overview.nr$Abbr==var.tlag
                                                                                      ,'Publication.Lag'])
                               }
                               # getting vintage dates that correspond to dates in rownames
                               vint.names=names(sets)
                               vint.last19=grep('99M12',vint.names)# last vintage of the 20th century
                               vint.dates=vint.names
                               vint.dates[1:vint.last19]=paste('19',vint.dates[1:vint.last19],sep='')
                               vint.dates[(vint.last19+1):length(vint.dates)]=paste('20',vint.dates[(vint.last19+1):length(vint.dates)],sep='')
                               vint.missing.zeros=grep('M[1-9]$',vint.dates)
                               # vint.dates[vint.missing.zeros]=paste('0',vint.dates[vint.missing.zeros],sep='')
                               vint.dates[vint.missing.zeros]=gsub('M',':0',vint.dates[vint.missing.zeros])
                               vint.dates=gsub('M',':',vint.dates)
                               vintage=data.frame(name=vint.names,date=vint.dates)
                               
                               # dates in data matrices, both realtime and unrevised, must 
                               # at least contain vintages (observation 2015.1 not existing yet, but vintage 2015.1).
                               # If not a line needs to be added
                               aux.match=match(vintage$date,row.names(df.unrevised))
                               missing.dates=vintage$date[is.na(aux.match)]
                               missing.dates=as.character(missing.dates)
                               df.unrevised[missing.dates,]=NA
                               
                               # are sets to short, too?
                               set.rt=sets[[193]]
                               aux.match=match(vintage$date,row.names(set.rt))
                               missing.dates=vintage$date[is.na(aux.match)]
                               missing.dates=as.character(missing.dates)
                               # add missing lines
                               sets=lapply(sets,function(x){x[missing.dates,]=NA
                                                            return(x)})
                               
                               # sourcing necassary scripts for estimation and forecast (for a description see "olsbmalag.Rmd")
                               source(paste(DirCode,'/bmafo.R',sep=''))
                               library(BMA)
                               push.down=function(variable,set){
                                       # pushes all variables to the forecast origin 
                                       aux=set[,variable]
                                       naux=length(aux)
                                       values=aux[is.na(aux)==F]
                                       nvalues=length(values)
                                       aux2=rep(NA,naux)
                                       aux2[(naux-nvalues+1):naux]=values
                                       return(aux2)
                               }
                               forecast.all=vector('list',length(sets))
                               for (vint.num in 1:length(sets)){
                                       set.rt=sets[[vint.num]]
                                       set.last.date=vintage[
                                               grep(paste(vintage$name[vint.num],'$',sep=''),vintage$name)
                                               # $ is needed to get exactly the vintage required. 
                                               # 1998M1 and not 1998M10,1998M11, 1998M12
                                               ,'date']
                                       set.lst.obs=grep(set.last.date,row.names(set.rt))
                                       set.rt=set.rt[1:set.lst.obs,]
                                       set.unrevised=df.unrevised[row.names(set.rt),]
                                       set=cbind(set.rt,set.unrevised)
                                       
                                       # filtering those variables that cannot be transformed to change rates
                                       # due to zeros.#### attention: this might lead to an unbalanced panel!
                                       set.aux=set
                                       set.aux[is.na(set)]=1
                                       chg.ok=colSums(set.aux==0)==0
                                       
                                       # adding transformations
                                       set=cbind(set,diff(set,1,1)
                                                 ,diff(set,2,1)
                                                 ,chg(set[,chg.ok],1)
                                                 ,chg(set[,chg.ok],12)
                                       )
                                       
                                       # window of max.obs month back
                                       set=set[(nrow(set)-max.obs+1):nrow(set),]
                                       
                                       # pushing down all variables        
                                       variables=colnames(set)
                                       dates=row.names(set)
                                       set=data.frame(sapply(variables,push.down,set))
                                       row.names(set)=dates
                                       
                                       # spliting inflation and rest of set
                                       colnames(set)[grep(target,'_chg_12',colnames(set))]='infl'
                                       infl.col=grep('infl',colnames(set))
                                       df.y=set[,infl.col,drop=F]
                                       set.woinfl=set[,-infl.col]# woinfl=without inflation
                                       
                                       for (i in 1:ncol(set.woinfl)){
                                               df.x=cbind(df.y,set.woinfl[,i,drop=F])
                                               if (i==1){forecast=bmafo(df.y,df.x,horizon,max.lag)}else{
                                                       forecast=rbind(forecast,bmafo(df.y,df.x,horizon,max.lag))
                                               }
                                       }   #end model loop
                                       forecast.all[[vint.num]]=forecast
                                       
                               }#end vintage
                       })stopCluster(cl)
               }
stopCluster(cl)