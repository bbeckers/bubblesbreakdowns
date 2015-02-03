bmafo <-function(df.y,df.x,horizon,max.lag){
        # building on lag.exact, lag.hormax and plugin.values, the functions estimates a 
        # ols model for predicting # of horizons ahead with a maximum lag of max.lag.
        # requires BMA package. df.y and df.x should be of same length and complete cases.
        # The resulting forecast is row name is the row name of the last observation and
        # "H"+horizon. 
        #         df.y=df.y 
        #         df.x=df.x 
        #         horizon=12 
        #         max.lag=12
        lag.exact<-function(df.x,lag.length){
                # Returns a dataframe of the lags of df.x with lag.length
                x.n=ncol(df.x)
                x.obs=nrow(df.x)
                df.x.lag=data.frame(matrix(NA,nrow(df.x),x.n))
                x.names=colnames(df.x)
                df.x.lag[(lag.length+1):nrow(df.x),]=df.x[1:(nrow(df.x)-lag.length),]
                rownames(df.x.lag)=rownames(df.x)
                colnames(df.x.lag)=paste(colnames(df.x),'L',lag.length,sep='')
                return(df.x.lag)
        }
        
        lag.hormax=function(df,horizon,maxlag){
                # returns a lag matrix of the lags of dataframe df from lag=horizon to lag=maxlag.
                x.n=ncol(df)
                x.obs=nrow(df)
                for (lag in horizon:(horizon+maxlag-1)){
                        if (lag==horizon){df.hormax=lag.exact(df,horizon)}else{
                                df.hormax=cbind(df.hormax,lag.exact(df,lag))
                        }
                }
                return(df.hormax)
        }
        plugin.values<-function(df,max.lag){
                #  creates the most recent values to be plugged in 
                # the regression results of lag=1 to max.lag
                x.obs=nrow(df)
                df.lag=lag.hormax(df,1,max.lag-1)
                plugin.values=cbind(df[x.obs,],df.lag[x.obs,])
                return(plugin.values)
        }
        df.x.lag=lag.hormax(df.x,horizon,max.lag)
        
        # making sure, that all lags employed have at least 20 observations.
        ind20=colSums(is.na(df.x.lag)==F)>=20
        df.x.lag=df.x.lag[,ind20]
        
        # restricting to complete cases
        xy=cbind(df.x.lag,df.y)
        xy.complete=complete.cases(xy)
        nobs=sum(xy.complete)
        df.x.lag=as.matrix(df.x.lag[xy.complete,])
        y=df.y[xy.complete,]
        y=as.vector(y)
        # estimation and prediction
        bma.res=bicreg(df.x.lag,y)
        plugin.values=plugin.values(df.x,max.lag)
        colnames(plugin.values)=colnames(df.x.lag)
        bma.fc=predict(bma.res,newdata=plugin.values,topmodels=1)$mean  
        names(bma.fc)=paste(row.names(df.y)[xy.complete][nobs],horizon,sep='H')
        bma.coefs=names(bma.res$ols[1,-1])[bma.res$ols[1,-1]!=0]
        coefs=paste(bma.coefs,collapse=',')
        useful=length(grep(colnames(df.y),bma.coefs))<length(bma.coefs)
        result=data.frame(bma.fc
                          ,nobs
                          ,residvar=bma.res$residvar[1]
                          ,names=paste(colnames(df.x)
                                       ,collapse=',')
                          ,useful=useful
                          ,coefficients=coefs
        )
        return(result)
}
