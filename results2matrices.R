test=sapply(forecast.all,function(x) nrow(x))
which.min(test)
tt=forecast.all[[which.min(test)]]
ttt=forecast.all[[1]]
tt1=tt$names
ttt1=ttt$names
beast=ttt1[!ttt1%in%tt1]

tt=sapply(forecast.all,function(x){
  sel=!x$names%in%beast
  x=x[sel,]
  return(x)
} )
ttt=tt[4,]
tttt=sapply(ttt,function(x) x)
namesm=forecast.all[[185]]$names
row.names(tttt)=namesm
msr.arx=tttt