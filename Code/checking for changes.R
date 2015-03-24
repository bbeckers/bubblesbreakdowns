source('chg.R')

earliest.chg<-function(df){
        # df is a dataframe of vintages in columns, rownames are observation 
        # dates, columnes vintages. This code is 
        # intended to check, if vintage data are needed at all.
        # Returns the first observation that is contained in all vintages 
        # and the first change to the vintage.
        result=list()
        df.compl=df[complete.cases(df),]
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

earliest.chg(PPPI)

