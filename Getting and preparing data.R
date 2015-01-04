#  ------------------------------------------------------------------------
# Starting values (do not comment) ---------------------------------------------------------
#  ------------------------------------------------------------------------

wd='C:/Users/Dirk/Documents/GitHub/bubblesbreakdowns'
setwd(wd)
dir.rt=paste(wd,'/data',sep='')

#  ------------------------------------------------------------------------
# Downloading and storing Data ---------------------------------------------------------
#  ------------------------------------------------------------------------

library("XLConnect")

dir.create(dir.rt)


# searching the index page of data base -----------------------------------

startaddress='http://www.philadelphiafed.org/research-and-data/real-time-center/real-time-data/data-files/'
fileaddress=paste(startaddress,'files/',sep='')
starthtml=readLines(startaddress)
# write.csv(starthtml,'starthtml.txt') # for inspection of starthtml only 

variable_index=regexec(paste('\\(','([A-Z]+)','\\)',sep=''),starthtml)
variable_vector=regmatches(starthtml,variable_index)#[[1]][2]

for (i in length(variable_vector):1){
        if (length(variable_vector[[i]])==0){variable_vector[[i]]=NULL}
}
rm(i)
variables=sapply(variable_vector,function(x) x[[2]])
rm(variable_index,variable_vector,starthtml)
variables=unique(variables)
nvar=length(variables)


# downloads all files -----------------------------------------------------

for (j in 1:nvar){
        dir.var=paste(dir.rt,'/',variables[j],sep='')
        dir.create(dir.var)
        
        varaddress=paste(startaddress,variables[j],sep='')
        varhtml=readLines(varaddress)
        # write.csv(varhtml,'varthtml.txt')
        
        xls_index=regexec(paste('/','(.*)','.xl',sep=''),varhtml)
        xls_index=regexec(paste('([A-Za-z0-9\\_]{1,}\\.xls[x]?)',sep=''),varhtml)
        
        xls_vector=regmatches(varhtml,xls_index)#[[1]][2]
        nuller=sapply(xls_vector, function(x)length(x)==0)
        xls=xls_vector[nuller!=T]
        xls=sapply(xls,function(x) x[[1]])
        #         grep(,xls)
        # write.csv(xls_vector,'varhtml.txt')
        for (i in 1:length(xls)){
                download.file(paste(fileaddress,xls[i],sep='')
                              , paste(dir.var,'/',xls[i],sep='')
                              ,mode='wb')                    
        }
}

# Failed download? --------------------------------------------------------

variables_downloaded=list.files(paste(dir.rt))
not_downl=match(variables,variables_downloaded)
variables[which(is.na(not_downl)==T)]

# Reading files, binding files of each variable together ------------------
no_mth_variables=character(0)
for (j in 1:nvar){
        #         if (variables[j]=="EMPLOY"){next}
        dir.var=paste(dir.rt,'/',variables[j],sep='')
        files=list.files(dir.var)
        mfiles=files[grep('MvMd',files)]
        if (length(mfiles)==0){
                no_mth_variables=c(no_mth_variables,variables[j])
                warning(paste('no monthly data available for ',variables[j],sep=''))
                next
        }
        cleanmfiles=gsub('[0-9]','',mfiles)
        varname=unique(cleanmfiles)
        if (length(mfiles)==1){
                #                 only one file
                df=read.xlsx2(paste(dir.var,'/',varname,sep=''),sheetIndex=1,row.names=1)
                eval(parse(text=paste(variables[j],'=df',sep='')))
        }else{
                #                         more than one file
                
                if(length(varname)>1){warning('not only one variable')}
                
                varname=gsub('\\.xls','',varname)
                numfiles=gsub(varname,'',mfiles)
                numfiles=gsub('\\.xls','',numfiles)
                numfiles=max(as.numeric(numfiles))
                
                for (i in 1:numfiles){
                        df_add=read.xlsx2(paste(dir.var,'/',varname,i,'.xls',sep=''),sheetIndex=1,row.names=1)
                        if (i==1){
                                df=df_add
                        }else{
                                if (nrow(df)!=nrow(df_add)){
                                        df_dates=row.names(df)
                                        df_add_dates=row.names(df_add)
                                        matchdates=match(df_add_dates,df_dates)
                                        add_dates=df_add_dates[which(is.na(matchdates)==T)]
                                        add_lines=df[1:length(add_dates),]
                                        add_lines[1:length(add_dates),]=NA
                                        row.names(add_lines)=add_dates
                                        df=rbind(df,add_lines)
                                }
                                df=cbind(df,df_add)
                        }
                        if (i==numfiles){
                                eval(parse(text=paste(variables[j],'=df',sep='')))
                        }
                }
        }
}
rm(no_mth_variables,variables,add_lines,df,df_add,cleanmfiles,i,j,mfiles,numfiles,add_dates,df_add_dates,df_dates,dir.rt,dir.rt,dir.var,fileaddress,files,matchdates,nvar,varname,startaddress)
# Needed repairs: EMPLOY 4 and CUM 2, zipped versions were all right.

list=ls()
list=list[-c(grep('wd',list),grep('list',list))]
for (variable in list){
        cat(paste('write.csv(',variable,',\'',wd,'/',variable,'.csv\')',sep=''),file='out')
        #         test=paste(wd,'/',variable,'.csv',sep='')
        eval(parse(file='out'))
}

save.image(paste(wd,"/realtime.RData",sep=''))



# #  ------------------------------------------------------------------------
# # Processing Data ---------------------------------------------------------
# #  ------------------------------------------------------------------------
# 
# load(paste(wd,"/realtime.RData",sep=''))
# # checking availability of each variable ----------------------------------
# # overview is a dataframe helping here
# variables=ls()
# variables=variables[grep('[A-Z]',variables)]
# overview=data.frame('nrows'=rep(NA,length(variables)),row.names=variables)
# startaddress='http://www.philadelphiafed.org/research-and-data/real-time-center/real-time-data/data-files/'
# 
# for (variable in variables){
#         aux=eval(parse(text=variable))
#         overview[variable,'nrows']=nrow(aux)
#         overview[variable,'ncols']=ncol(aux)   
#         overview[variable,'first obs']=row.names(aux)[1]    
#         overview[variable,'first obs year']=as.numeric(strsplit(overview[variable,'first obs'],':')[[1]][1])
#         overview[variable,'first obs month']=as.numeric(strsplit(overview[variable,'first obs'],':')[[1]][2])
#         
#         overview[variable,'last obs']=row.names(aux)[nrow(aux)]
#         overview[variable,'last obs year']=as.numeric(strsplit(overview[variable,'last obs'],':')[[1]][1])
#         overview[variable,'last obs month']=as.numeric(strsplit(overview[variable,'last obs'],':')[[1]][2])
#         
#         
#         overview[variable,'first vintage']=colnames(aux)[1]  
#         overview[variable,'first vintage']=gsub(variable,'',overview[variable,'first vintage'])
#         
#         overview[variable, 'first vintage year']=as.numeric(strsplit(overview[variable, 'first vintage'],'M')[[1]][1])
#         if (overview[variable, 'first vintage year']<=14){
#                 overview[variable, 'first vintage year']=overview[variable, 'first vintage year']+2000}else{
#                         overview[variable, 'first vintage year']=overview[variable, 'first vintage year']+1900}
#         overview[variable, 'first vintage month']=as.numeric(strsplit(overview[variable, 'first vintage'],'M')[[1]][2])
#         
#         overview[variable,'last vintage']=colnames(aux)[ncol(aux)]  
#         overview[variable,'last vintage']=gsub(variable,'',overview[variable,'last vintage'])
#         
#         overview[variable, 'last vintage year']=as.numeric(strsplit(overview[variable, 'last vintage'],'M')[[1]][1])
#         if (overview[variable, 'last vintage year']<=14){
#                 overview[variable, 'last vintage year']=overview[variable, 'last vintage year']+2000}else{
#                         overview[variable, 'last vintage year']=overview[variable, 'last vintage year']+1900}
#         overview[variable, 'last vintage month']=as.numeric(strsplit(overview[variable, 'last vintage'],'M')[[1]][2])
#         
#         
#         varaddress=paste(startaddress,variable,sep='')
#         varhtml=readLines(varaddress)
#         name=varhtml[grep('title',varhtml)]
#         overview[variable,'description']=gsub('<title>||</title>||- historical real-time data - Philadelphia Fed','',name)
#         rm(name,varaddress,varhtml,aux)
# }
# rm(startaddress)
# 
# # Getting rid of all variables that start in 2009-8 or later  -------------
# # (mainly household spending variables)
# overview.complete=overview
# overview=overview[overview$ncols>64,]
# variables=row.names(overview)
# variables.rm=row.names(overview.complete[overview.complete$ncols<=64,])
# rm(list=variables.rm)
# 
# overview=overview[order(overview$'first vintage year'),]
# 
# 
# # resizing variables to the smallest common (vintage) sample --------------
# vintages=colnames(POP)
# vintages=gsub('POP','',vintages)
# 
# for (variable in variables){
#         colselection=paste(variable,vintages,sep='')
#         comm=paste(variable,'=',variable,'[,colselection]',sep='')
#         eval(parse(text=comm))   
# }
# 
# # eliminating non-numeric elements ----------------------------------------
# 
# for (variable in variables){
#         aux=eval(parse(text=variable))
#         aux[aux=='#N/A']=NA
#         write.csv(aux,'test.csv')
#         aux=read.csv('test.csv',row.names=1)
#         eval(parse(text=paste(variable,'=aux',sep='')))
# }
# # creating specimen and resizing dataframes -------------------------------------------------------
# specimen.nrow=nrow(IPM)
# specimen.ncol=ncol(IPM)
# specimen=data.frame(matrix(NA,nrow=specimen.nrow,ncol=specimen.ncol))
# row.names(specimen)=row.names(IPM)
# 
# for (variable in variables){
#         aux=eval(parse(text=variable))
#         aux.specimen=specimen
#         colnames(aux.specimen)=colnames(aux)
#         aux.specimen[row.names(aux),colnames(aux)]=aux
#         aux=aux.specimen
#         eval(parse(text=paste(variable,'=aux',sep='')))
#         
# }
# rm(aux.specimen,aux,specimen)
# 
# # First observation in remaining variables and vintages, deleting useless rows -------------------
# 
# for (variable in variables){
#         overview[variable,'first obs subsample']=row.names(eval(parse(text=variable)))[complete.cases(eval(parse(text=variable)))][1]
#         overview[order(overview$'first obs subsample'),]
# }
# 
# fobs=overview[order(overview$'first obs subsample'),'first obs subsample'][1]
# fobs.row.number=grep(fobs,row.names(IPM))
# for (variable in variables){
#         aux=eval(parse(text=variable))
#         aux=aux[fobs.row.number:nrow(aux),]
#         eval(parse(text=paste(variable,'=aux',sep='')))
# }
# rm(aux)
# 
# # Creating sets of vintages -----------------------------------------------
# vintage=vintages[10]
# 
# vintage.extractor=function(variable,vintage){
#         # returns a vintage of a variable
#         vintage.colnum=grep(vintage,colnames(CUM))
#         aux=eval(parse(text=variable))[vintage.colnum]
#         colnames(aux)=variable
#         return(aux)        
# } 
# set.maker=function(vintage){
#         # puts corresponding vintages of all variables together in one set
#         aux=data.frame(matrix(NA,nrow=nrow(H),ncol=length(variables)))
#         colnames(aux)=variables
#         row.names(aux)=row.names(POP)
#         for (variable in variables){
#                 aux[,variable]=vintage.extractor(variable,vintage)
#         }
#         return(aux)
# }
# 
# 
# # creating a list containing all vintage sets -----------------------------
# 
# sets=vector('list',length=193)
# names(sets)=vintages
# 
# for (vintage in vintages){
#         sets[[vintage]]=set.maker(vintage)
# }
# 
# # pushing all values to forecast origin -----------------------------------
# # (eliminating NAs)
# # Here for one test-set
# 
# set=sets[[vintage]]
# push.up=function(variable,set){
#         aux=set[,variable]
#         naux=length(aux)
#         values=aux[is.na(aux)==F]
#         nvalues=length(values)
#         aux2=rep(NA,naux)
#         aux2[(naux-nvalues+1):naux]=values
#         return(aux2)
# }
# set.pushed=sapply(variables,push.up,set)
# row.names(set.pushed)=row.names(set)
