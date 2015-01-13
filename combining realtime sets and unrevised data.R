
# loading realtime data sets and unrevised data 
DirCode='C:/Users/Dirk/Documents/GitHub/bubblesbreakdowns'
source(paste(DirCode,'/lag.exact.R',sep=''))
source(paste(DirCode,'/diff.R',sep=''))
source(paste(DirCode,'/chg.R',sep=''))

load(paste(DirCode,'/data/sets.Rdata',sep=''))
df.unrevised=read.csv(paste(DirCode,'/data/unrevised data.csv'
                            ,sep='')
                      ,sep=';'
                      ,na.strings='NaN'
                      ,row.names=1
                      )
# overview=read.csv(paste(DirCode,'/overview.csv',sep=''),row.names=1)


                      


# transforming rownames (dates) of unrevised data accordingly
df.ur.rownames=row.names(df.unrevised)
dates=strsplit(df.ur.rownames,'\\.')
dates=sapply(dates,function(x) x[c(3,2)])
dates=apply(dates,2,function(x) paste(x,collapse=':'))
row.names(df.unrevised)=dates

# getting vintage dates that correspond to dates in rownames
vint.names=names(sets)
vint.last19=grep('99M12',vint.names)
vint.dates=vint.names
vint.dates[1:vint.last19]=paste('19',vint.dates[1:vint.last19],sep='')
vint.dates[(vint.last19+1):length(vint.dates)]=paste('20',vint.dates[(vint.last19+1):length(vint.dates)],sep='')
vint.missing.zeros=grep('M[1-9]$',vint.dates)
# vint.dates[vint.missing.zeros]=paste('0',vint.dates[vint.missing.zeros],sep='')
vint.dates[vint.missing.zeros]=gsub('M',':0',vint.dates[vint.missing.zeros])
vint.dates=gsub('M',':',vint.dates)
vintage=data.frame(name=vint.names,date=vint.dates)

# dates in data matrices, both realtime and unrevised, must at least contain vintages. If not
# a line needs to be added
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

# comparing cpi rt and in unrevised
test.cpi=df.unrevised[row.names(set.rt),'CPI',drop=F]
test.cpi$rt=NA
test.cpi[row.names(set.rt),'rt']=set.rt[row.names(set.rt),'PCPI']

# possible loop -----------------------------------------------------------


set.rt=sets[[1]]

set.last.date=vintage[grep(vintage$name[1],vintage$name),'date']
set.lst.obs=grep(set.last.date,row.names(set.rt))
set.rt=set.rt[1:set.lst.obs,]
set.unrevised=df.unrevised[row.names(set.rt),]

set=cbind(set.rt,set.unrevised)

# adding transformations

set=cbind(set,diff(set,1),diff(diff(set,1),1),chg(set,1),chg(set,12))
