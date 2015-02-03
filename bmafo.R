bmafo <-
function(y.df,x.df,horizon,max.lag){
        # building on lag.exact, lag.hormax and plugin.values, the functions estimates a 
        # ols model for predicting # of horizons ahead with a maximum lag of max.lag.
        # requires BMA package. y.df and x.df should be of same length and complete cases.
        # The resulting forecast is row name is the row name of the last observation and
        # "H"+horizon. 
#         y.df=df.y 
#         x.df=df.x 
#         horizon=12 
#         max.lag=12
        lag.exact<-function(x.df,lag.length){
                # Returns a dataframe of the lags of x.df with lag.length
                x.n=ncol(x.df)
                x.obs=nrow(x.df)
                x.df.lag=data.frame(matrix(NA,nrow(x.df),x.n))
                x.names=colnames(x.df)
                x.df.lag[(lag.length+1):nrow(x.df),]=x.df[1:(nrow(x.df)-lag.length),]
                rownames(x.df.lag)=rownames(x.df)
                colnames(x.df.lag)=paste(colnames(x.df),'L',lag.length,sep='')
                return(x.df.lag)
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
        x.df.lag=lag.hormax(x.df,horizon,max.lag)
        
        # making sure, that all lags employed have at least 20 observations.
        ind20=colSums(is.na(x.df.lag)==F)>=20
        x.df.lag=x.df.lag[,ind20]
        
        # restricting to complete cases
        xy=cbind(x.df.lag,y.df)
                xy.complete=complete.cases(xy)
                nobs=sum(xy.complete)
                x.df.lag=as.matrix(x.df.lag[xy.complete,])
                y=y.df[xy.complete,]
        y=as.vector(y)
        bma.res=bicreg(x.df.lag,y)
        plugin.values=plugin.values(x.df,max.lag)
        colnames(plugin.values)=colnames(x.df.lag)
        bma.fc=predict(bma.res,newdata=plugin.values,topmodels=1)$mean  
        names(bma.fc)=paste(row.names(y.df)[xy.complete][nobs],horizon,sep='H')
        bma.coefs=names(bma.res$ols[1,-1])[bma.res$ols[1,-1]!=0]
        coefs=paste(bma.coefs,collapse=',')
        useful=length(grep(colnames(y.df),bma.coefs))<length(bma.coefs)
        result=data.frame(bma.fc
                          ,nobs
                          ,residvar=bma.res$residvar[1]
                          ,names=paste(colnames(x.df)
                                       ,collapse=',')
                          ,useful=useful
                          ,coefficients=coefs
                          )
        return(result)
        }
