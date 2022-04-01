VARSEDIM<-function(data, variables, group, method="overlap", stepwise=TRUE, VARSEDIG=TRUE,
minimum=TRUE, kernel="gaussian", cor=TRUE, ellipse=FALSE, convex=TRUE, file="Plots VARSEDIG.pdf",
na="NA", dec=",", row.names=FALSE){

#Selection of variables
datosT<-data.frame(subset(data, select=group), subset(data, select=variables))
datos<-na.exclude(datosT)
groups<-as.character(unique(datos[,1]))
n<-length(groups)


pdf(file=file)

###Running VARSEDIG among all taxa
for(zz in 1:(n-1)){
for(tt in (zz+1):n){

datosT<-datos

if(dec=="."){
write.csv(x=datosT,file = "Temp.csv", fileEncoding = "", row.names=row.names,na=na)
datosT<-read.csv(file="Temp.csv" ,header=TRUE)
}
else{
write.csv2(x = datosT,file = "Temp.csv", fileEncoding = "", row.names=row.names,na=na)
datosT<-read.csv2(file="Temp.csv" ,header=TRUE)
}


VARSEDIG::VARSEDIG(data=datosT, variables=variables, group=group, group1=groups[zz], group2=groups[tt], ellipse=ellipse, convex=convex,
method=method, stepwise=stepwise, VARSEDIG=VARSEDIG, minimum=minimum, kernel=kernel, cor=cor, devnew=FALSE,
na=na, dec=dec)

}#End Loop tt
}#End Loop zz

if(file.exists("Temp.csv")){
file.remove("Temp.csv")
}

list<-dev.list()
borrar<-which(names(list)=="pdf")
dev.off(list[borrar])

}