# source('chg.R')

earliest.chg<-function(df.name){
        # df.name is the name of a dataframe of vintages in columns, rownames are observation 
        # dates, columnes vintages. This code is 
        # intended to check, if vintage data are needed at all.
        # Returns the first observation that is contained in all vintages 
        # and the first change to the vintage.
        df=eval(parse(text=df.name))
        
        result=list()
        
        # get observations that appear in every vintage
        df.compl=df[complete.cases(df),]
        # get the last unrevised vintage
        nvint.compl=ncol(df.compl)
        nobs.compl=ncol(df.compl)
        df.compl.uc=matrix(NA,nrow=nobs.compl,ncol=nvint.compl)
        
        for (obs in 1:nobs.compl){
                for (vint in 1:nvint.compl){
                        df.compl.uc[obs,vint]=df.compl[obs,1]==df.compl[obs,vint]
                }
        }
        rev.start=colSums(df.compl.uc)==nvint.compl
        unrevised.vint=colnames(df.compl)[rev.start]
        if (length(unrevised.vint)==0){result$unrevised='none'}else{
                result$unrevised=unrevised.vint[length(unrevised.vint)]       
        }
        
        # get the first complete observation appearing in all vintages and get
        # the first revised observation
        test=apply(df.compl,1,function(x) max(x)-min(x))
        result$first.complete=row.names(df.compl)[1]
        result$first.chg=row.names(df.compl)[test!=0][1]
        
        #         df12=chg(df,12)# creating 12 mth. change rates
        #         df12.compl=df12[complete.cases(df12),]
        #         test12=apply(df12.compl,1,function(x) max(x)-min(x))
        #         
        #         row.names(df12.compl)[test!=0][1]
        #         df1=chg(df,1)# creating 1 mth. change rates
        #         df1.compl=df1[complete.cases(df1),]
        #         test1=apply(df1.compl,1,function(x) max(x)-min(x))  
        return(result)
}



