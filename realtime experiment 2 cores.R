library("parallel")
library("foreach")
library("doParallel")

cl <- makeCluster(detectCores() )
registerDoParallel(cl, cores = detectCores() )
data = foreach(i = 1:length(filenames), .packages = c("ncdf","chron","stats"),
               .combine = rbind) %dopar% {
                       try({
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
                       })
               }
stopCluster(cl)