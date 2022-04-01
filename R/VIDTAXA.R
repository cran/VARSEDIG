VIDTAXA<-function(data, var,  labels,  cat=NULL, analysis="PCA", por=80, k=NULL,  pthreshold=0.05, ellipse=FALSE, convex=FALSE, dim=c(1,2), size=c(1,5), showCluster=TRUE,
VIF=FALSE, VARSEDIG=TRUE, BUBBLE=TRUE,  threshold=10,  method="overlap", minimum=TRUE, ResetPAR=TRUE, PAR=NULL, PCA=NULL, SCATTERPLOT=NULL, HCLUST=NULL,
CLUSTER=NULL, BOXPLOT=NULL, mfrowBOXPLOT=NULL, LabelCat=NULL, COLOR=NULL, COLORC=NULL, COLORB=NULL, PCH=NULL,
XLIM=NULL, YLIM=NULL,  XLAB=NULL, YLAB=NULL, ylabBOXPLOT=NULL,
LEGEND=NULL, MTEXT= NULL, TEXTvar=NULL, TEXTlabels=NULL, arrows=TRUE, larrow=0.7, colArrows="black", quadratic=FALSE,
file1="Output.txt", file2="Cat loadings.csv", file3="Descriptive statistics of clusters.csv", 
file4="Original data and cluster number.csv", file5="Var loadings-Linear.csv", file6="Cat loadings-Linear.csv",
file7="Table cross-validation-Linear.csv", file8="Cases cross-validation-Linear.csv", 
file9="Table cross-validation-Quadratic.csv", file10="Cases cross-validation-Quadratic.csv", file11="Plots VARSEDIG.pdf", file12="U Mann-Whitney test.csv",
na="NA", dec=",", row.names=TRUE){

if(is.null(k)){
kmethod="automatically"
}
else{
kmethod="manually"
}

count<-function(x){ 
length(na.omit(x)) 
} 
######Final


#####Discriminat Function

DA<-function(data, cat, var, ellipse=TRUE, convex=FALSE, quadratic=FALSE, expand=TRUE, dimS=c(1,2), ResetPAR=TRUE, PAR=NULL,
CANDISC1=NULL, CANDISC2=NULL, CANPLOT=NULL, SCATTERPLOT=NULL, COLOR=NULL,
PCH=NULL,TEXT=NULL, LEGEND=NULL, AXIS=NULL, MTEXT= NULL, arrows=TRUE, larrow=0.95, colArrows="black",
 file1="Var loadings-Linear.csv", file2="Cat loadings-Linear.csv",
file3="Table cross-validation-Linear.csv", file4="Cases cross-validation-Linear.csv", 
file5="Table cross-validation-Quadratic.csv", file6="Cases cross-validation-Quadratic.csv", 
na="NA", dec=",", row.names=TRUE){


datos<-data

datosT<-data.frame(subset(datos, select=cat), subset(datos, select=var))

datos<-na.exclude(datosT)

datos[,1]<-as.character(datos[,1])


#LINEAR  DISCRIMINANT FUNCTION ANALYSIS

lenf1<-length(unique(datos[,cat]))

if(ResetPAR==TRUE){
#Resetear par() a las opciones por defecto
resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}
par(resetPar()) }
else{
}

if(!is.null(PAR)){
parexe<-paste("par(new,",toString(x=PAR), ")")
eval(parse(text=parexe))
}
else{
par(font.lab=2, mar=c(5,5,3,2),cex.lab=1.5)
}


zz<-dim(datos)
fo<-paste(names(datos)[2],",", names(datos)[3] )

for(i in 4:zz[2]){
fo<-paste(fo, ",", names(datos)[i])
}

fo<-paste("(cbind(",fo,")","~",names(datos)[1],")")

mod<-lm(fo, data=datos)

if(!is.null(CANDISC1)){
candiscexe<-paste("candisc::candisc(","mod=mod,", "term=cat,","ndim=1,", toString(x=CANDISC1), ")")
can1<-eval(parse(text=candiscexe))
}
else{
candiscexe<-paste("candisc::candisc(","mod=mod,", "term=cat,","ndim=1",  ")")
can1<-eval(parse(text=candiscexe))
}


Resultado<-list("Variance explanied", paste("Can", 1, round(can1$pct[1], digits=2), "%"))
tyu<-length(can1$pct)
for(h in 2:tyu){
Resultado<-append(Resultado,paste("Can", h, round(can1$pct[h], digits=2), "%"))
}

devi<-as.numeric(length(dev.list()))
if(devi<=1) asss<-1 else dev.new()

if(!is.null(CANPLOT)){
canplotexe<-paste("plot(","x=can1,", toString(x=CANPLOT), ")")
eval(parse(text=canplotexe))
}
else{
canplotexe<-paste("plot(","x=can1", ")")
eval(parse(text=canplotexe))
}


if(dec=="."){
write.csv(x= can1$scores,file = file2, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = can1$scores ,file = file2, fileEncoding = "", row.names=row.names,na=na)
}

if(dec=="."){
write.csv(x= can1$structure,file = file1, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = can1$structure ,file = file1, fileEncoding = "", row.names=row.names,na=na)
}



lenf<-length(unique(datos[,cat]))

if(lenf>2){


if(!is.null(CANDISC2)){
candiscexe<-paste("candisc::candisc(","mod=mod,", "term=cat,","ndim=(lenf-1),", toString(x=CANDISC2), ")")
can2<-eval(parse(text=candiscexe))
}
else{
candiscexe<-paste("candisc::candisc(","mod=mod,", "term=cat,","ndim=(lenf-1)",  ")")
can2<-eval(parse(text=candiscexe))
}

datos2<-data.frame(can2$scores[order(can2$scores[,1]), ])
dim2<-dim(datos2)

if(dim2[2]>dimS[2]){

Xmax<-max(datos2[,(dimS[1]+1)])
Xmin<-min(datos2[,(dimS[1]+1)])
Ymin<-min(datos2[,(dimS[2]+1)])
Ymax<-max(datos2[,(dimS[2]+1)])

if(!is.null(COLOR)){
color1<-COLOR
}
else{
if(convex==TRUE) color1<-rainbow(length(unique(datos2[,1])), alpha=0.4) else color1<-rainbow(length(unique(datos2[,1])))
}

catn<-length(unique(datos2[,1]))

pch1<-14+1
for (t in 2:catn){
pcht<-14+t
pch1<-append(pch1,pcht)
}




if(!is.null(PCH)) pcht<-PCH else pcht<-pch1


Resultado<-list("Variance explanied", paste("Can", 1, round(can1$pct[1], digits=2), "%"))
tyu<-length(can2$pct)
for(h in 2:tyu){
Resultado<-append(Resultado,paste("Can", h, round(can1$pct[h], digits=2), "%"))
}

can2$structure<-unique(can2$structure)

facX<-((Xmax-Xmin)/(max(can2$structure[,1])-min(can2$structure[,1])))/2

facY<-((Ymax-Ymin)/(max(can2$structure[,2])-min(can2$structure[,2])))/2

conX<-can2$structure[,1]*(facX)
conY<-can2$structure[,2]*(facY)

if(expand==TRUE) conX<-conX else conX<-can2$structure[,1] 
if(expand==TRUE) conY<-conY else conY<-can2$structure[,2] 

datos3<-data.frame(datos2[,1],datos2[,(dimS[1]+1)],datos2[,(dimS[2]+1)])

Xmax<-max(datos3[,2], max(conX))
Xmin<-min(datos3[,2], min(conX))
Ymin<-min(datos3[,3], min(conY))
Ymax<-max(datos3[,3],max(conY))


colnames(datos3)<-c("Categ", "CanX", "CanY")
Categ<-datos3[,1]
dev.new()

if(!is.null(PAR)){
parexe<-paste("par(new,",toString(x=PAR), ")")
eval(parse(text=parexe))
}
else{
par(font.lab=2, mar=c(5,5,3,2),cex.lab=1.5)
}

if(!is.null(SCATTERPLOT)){
scatterplotexe<-paste("car::scatterplot(","CanY~CanX| Categ,", "data=datos3,", toString(x=SCATTERPLOT), ")")
eval(parse(text=scatterplotexe))
}
else{
scatterplotexe<-paste("car::scatterplot(","CanY~CanX| Categ,", "data=datos3,", "regLine=FALSE,",
"smooth=FALSE,","grid=FALSE,", "xlim=c(Xmin,Xmax),", "ylim=c(Ymin,Ymax),",
"boxplots=FALSE,", "by.groups=TRUE,", "ellipse=ellipse,", "col=color1,", "pch=pcht,",
"xlab=paste('CAN',dimS[1],round(can2$pct[dimS[1]],digits=2),'%'),","ylab=paste('CAN',dimS[2],round(can2$pct[dimS[2]],digits=2),'%'),","legend=list(coords=c(x=500000000,y=500000000))", ")")
eval(parse(text=scatterplotexe))
}

if(!is.null(cat)){
if(convex==TRUE){
dati<-as.character(unique(datos3[,1]))
for(zz in 1:catn){
datis<-subset(datos3,(Categ == dati[zz]))
hpts <- chull(x=datis[,"CanX"], y=datis[,"CanY"])
hpts <- c(hpts, hpts[1])
datiss<-datis[hpts,]
polygon(datiss[,2],datiss[,3], col=color1[zz], border=NA)
}
}
}




if(!is.null(TEXT)){
textexe<-paste("text(","x=conX,", "y=conY," ,"labels = rownames(can2$structure),", toString(x=TEXT), ")")
eval(parse(text=textexe))
}
else{
textexe<-paste("text(","x=conX,", "y=conY," ,"labels = rownames(can2$structure),",
"xlim=c(Xmin,Xmax),", "ylim=c(Ymin,Ymax),", "col='black',","cex=1", ")")
eval(parse(text=textexe))
}


meanx<-(Xmin+Xmax)/2
meany<-(Ymin+Ymax)/2

minar<-min((Ymax-Ymin),(Xmax-Xmin))/4

dima<-length(conX)



if(arrows==TRUE){
for(aa in 1:dima){

x2<-conX[aa]*larrow
y2<-conY[aa]*larrow

IDPmisc::Arrows(x1=meanx, y1=meany, x2=x2, y2=y2, h.col=colArrows, sh.col=colArrows, h.col.bo=colArrows, open=FALSE)

}#END BUCLE arrows


}#END arrows TRUE





if(!is.null(LEGEND)){
legendexe<-paste("legend(",toString(x=LEGEND), ")")
eval(parse(text=legendexe))
}
else{
legendexe<-paste("legend(","x='topleft',","legend=unique(datos3[,1]),", "bty='n',", "col=color1,",
"pch=pcht", ")")
eval(parse(text=legendexe))
}

}

}

if(!is.null(AXIS)){
axisexe<-paste("axis(",toString(x=AXIS), ")")
eval(parse(text=axisexe))
}
else{
}

if(!is.null(MTEXT)){
mtextexe<-paste("mtext(",toString(x=MTEXT), ")")
eval(parse(text=mtextexe))
}
else{
}

if(lenf>2){
if(dec=="."){
write.csv(x= can2$scores,file = file2, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = can2$scores ,file = file2, fileEncoding = "", row.names=row.names,na=na)
}

if(dec=="."){
write.csv(x= can2$structure,file = file1, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = can2$structure ,file = file1, fileEncoding = "", row.names=row.names,na=na)
}

}
else{
if(dec=="."){
write.csv(x= can1$scores,file = file2, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = can1$scores ,file = file2, fileEncoding = "", row.names=row.names,na=na)
}

if(dec=="."){
write.csv(x= can1$structure,file = file1, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = can1$structure ,file = file1, fileEncoding = "", row.names=row.names,na=na)
}


}


#Porcentaje de aciertos
lenf<-lenf1
repli<-rep(1,lenf)
ad <- MASS::lda(datos[,var], grouping= datos[,cat], prior = repli/lenf)
t <- table(predict(ad, datos[,var])$class, datos[,cat])

acierto<-round(100-sum(predict(ad, datos[,var])$class!=datos[,cat])/ad$N*100, digits=2)


#Validacion cruzada
ad <- MASS::lda(datos[,var], grouping= datos[,cat], prior = repli/lenf, CV=TRUE)
t<- table(ad$class,datos[,cat])
#prop.table(t,2)*100
acierto2<-round(100-sum(ad$class!=datos[,cat])/sum(t)*100, digits=2)

indent<-data.frame(datos[,cat],ad$class)
colnames(indent)<-c("Real","Prediction")


if(dec=="."){
write.csv(x= t,file = file3, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = t ,file = file3, fileEncoding = "", row.names=row.names,na=na)
}

if(dec=="."){
write.csv(x= indent,file = file4, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = indent,file = file4, fileEncoding = "", row.names=row.names,na=na)
}


#QUADRATIC  DISCRIMINANT FUNCTION ANALYSIS

acierto3<-""
acierto4<-""

if(quadratic==TRUE){

counts<- data.frame(aggregate(datos[,var[1]],by=list(datos[,cat]),length)) 
minimo<-min(counts[,2] , na.rm = TRUE)

if(minimo>4){
#Porcentaje de aciertos


lenf<-lenf1
repli<-rep(1,lenf)
ad <- MASS::qda(datos[,var], grouping= datos[,cat], prior = repli/lenf)
t <- table(predict(ad, datos[,var])$class, datos[,cat])

acierto3<-round(100-sum(predict(ad, datos[,var])$class!=datos[,cat])/ad$N*100, digits=2)


#Validacion cruzada
ad <- MASS::qda(datos[,var], grouping= datos[,cat], prior = repli/lenf, CV=TRUE)
t<- table(ad$class,datos[,cat])
#prop.table(t,2)*100
acierto4<-round(100-sum(ad$class!=datos[,cat])/sum(t)*100, digits=2)

indent<-data.frame(datos[,cat],ad$class)
colnames(indent)<-c("Real","Prediction")


if(dec=="."){
write.csv(x= t,file = file5, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = t ,file = file5, fileEncoding = "", row.names=row.names,na=na)
}

if(dec=="."){
write.csv(x= indent,file = file6, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = indent,file = file6, fileEncoding = "", row.names=row.names,na=na)
}

}
}#end quadratic


Resultado<-c(Resultado, "LINEAR  DISCRIMINANT FUNCTION ANALYSIS", "Cases correctly identified",
acierto, "Cases correctly identified by cross-validation", acierto2,
"QUADRATIC  DISCRIMINANT FUNCTION ANALYSIS", "Cases correctly identified", acierto3,
"Cases correctly identified by cross-validation",acierto4)

print(Resultado)

}

#####End Discriminant function





####Function Bubble
Bubble<-function(data, varY, varX, varSize=NULL, size=c(1,5), legSpos = "bottomright" , orientation = "horizontal",
digitsS=1, k=NULL, XLAB=NULL, YLAB=NULL, XLIM=NULL, YLIM=NULL){


datos<-data

datosT<-data.frame(subset(datos, select=varX), subset(datos, select=varY), subset(datos, select=varSize))

datos<-na.exclude(datosT)



if(!is.null(XLAB)) xlab<-XLAB else xlab<-varX

if(!is.null(YLAB)) ylab<-YLAB else ylab<-varY



if(!is.null(XLIM)){
minsx<-XLIM[1]
maxsx<-XLIM[2]
}
else{
minsx<-min(datos[,varX])
maxsx<-max(datos[,varX])
XLIM<-c(minsx,maxsx)
}


if(!is.null(YLIM)){
minsy<-YLIM[1]
maxsy<-YLIM[2]
}
else{
minsy<-min(datos[,varY])
maxsy<-max(datos[,varY])
YLIM<-c(minsy,maxsy)
}


maxS<-max(datos[,varSize])
minS<-min(datos[,varSize])
matriz<-matrix(c(minS,maxS, size[1], size[2]),nrow = 2 , ncol = 2)
regS<-lm(matriz[,2]~matriz[,1])



if(!is.null(varSize)) cex<-regS$coefficients[1]+regS$coefficients[2]*datos[1,varSize] else cex<-1



plot(x=datos[1,varX], y=datos[1,varY], cex=cex, col='black', xlim=XLIM, yaxt="n", xaxt="n",  ylim=YLIM, xlab=xlab, ylab=ylab, pch=16)

axis(1, seq(1, k, by =1))
axis(2, seq(1, k, by = 1))

uni<-length(unique(datos[,3]))

dimS<-dim(datos)
for(zz in 2:dimS[1]){
if(!is.null(varSize)) cex<-regS$coefficients[1]+regS$coefficients[2]*datos[zz,varSize] else cex=1
if(uni==1) cex<-1 else cex<-cex
points(x=datos[zz,varX], y=datos[zz,varY], cex=cex, col='black', xlim=XLIM, ylim=YLIM, pch=16)

}


ranX<-abs(XLIM[2]-XLIM[1])
ranY<-abs(YLIM[2]-YLIM[1])



if(!is.null(varSize)){

if(uni==1){
val<-min(datos[,3], na.rm=TRUE)
size[2]<-size[1]
}
else{
val<-format((size[1]-regS$coefficients[1])/regS$coefficients[2],digits=digitsS)
for(zz in (size[1]+1):size[2]){
val<-append(val,format((zz-regS$coefficients[1])/regS$coefficients[2],digits=digitsS))
}
}

hh<-0

xx<-XLIM[2]
yy<-YLIM[1]
gy<-(1)
gx<-(-1)


for(zz in size[1]:size[2]){
hh<-hh+1
if(as.numeric(val[size[2]])>1000) mul<-10 else mul<-8
if(uni>1){
if(XLIM[1]==XLIM[2]){
points(xx+(zz+xx)/15, yy+abs(yy*1/100)*(gy), cex=zz, xlim=XLIM, ylim=YLIM)
text(xx+(zz+xx)/15, yy+ranY*9/100*(gy),label=val[hh])
}
else{
points(xx+ranX*((-size[1]+zz)*mul)/100*(gx), yy+abs(yy*1/100)*(gy), cex=zz, xlim=XLIM, ylim=YLIM)
text(xx+ranX*((-size[1]+zz)*mul)/100*(gx), yy+ranY*9/100*(gy),label=val[hh])
}
}
else{
text(xx+ranX*((-size[1]+zz)*mul)/100*(gx), yy+ranY*9/100*(gy),label=val[hh], pos=3)
}
}
}



}####end Bubble




###Descriptive statistics
stat<-function(data, var, grupo=NULL, file=NULL){

if(is.null(grupo)){
datos<-data.frame(subset(data, select=var))
}
else{
datos<-data.frame(subset(data, select=var), subset(data, select=grupo))
}

if(is.null(grupo)){
num<-as.data.frame(apply(X =datos[,var] , MARGIN = 2 , FUN = count))
min<-as.data.frame(apply(X =datos[,var] , MARGIN = 2 , FUN = min, na.rm=TRUE))
max<-as.data.frame(apply(X =datos[,var] , MARGIN = 2 , FUN = max, na.rm=TRUE))
media<-as.data.frame(apply(X =datos[,var] , MARGIN = 2 , FUN = mean, na.rm=TRUE))
mediana<-as.data.frame(apply(X =datos[,var] , MARGIN = 2 , FUN = median, na.rm=TRUE))
suma<-as.data.frame(apply(X =datos[,var] , MARGIN = 2 , FUN = sum, na.rm=TRUE))
sd<-as.data.frame(apply(X =datos[,var] , MARGIN = 2 , FUN = sd, na.rm=TRUE))
cv<-sd/media
}
else{

num<-aggregate(x=datos[,var], by = list(datos[ , grupo]), FUN=count)
min<-aggregate(x=datos[,var], by = list(datos[ , grupo]), FUN=min, na.rm=TRUE)
max<-aggregate(x=datos[,var], by = list(datos[ , grupo]), FUN=max, na.rm=TRUE)
media<-aggregate(x=datos[,var], by = list(datos[ , grupo]), FUN=mean, na.rm=TRUE)
mediana<-aggregate(x=datos[,var], by = list(datos[ , grupo]), FUN=median, na.rm=TRUE)
suma<-aggregate(x=datos[,var], by = list(datos[ , grupo]), FUN=sum, na.rm=TRUE)
sd<-aggregate(x=datos[,var], by = list(datos[ , grupo]), FUN=sd, na.rm=TRUE)
cv<-sd/media
}


if(is.null(grupo)){
final<-cbind(num,min,max,media,mediana, suma, sd,cv )
names(final)<-names(final)<-c("N\u00FAmero","M\u00EDnimos","M\u00E1ximo","Media","Mediana","Suma", "SD", "CV")

if(dec=="."){
write.csv(x=final,file = file, fileEncoding = "", row.names=FALSE,na=na)
}
else{
write.csv2(x = final ,file = file, fileEncoding = "", row.names=FALSE,na=na)
}

}
else{


###Numero
dim<-dim(num)
names<-names(num)
rep<-rep(names[2],dim[1])
vnum<-data.frame(num[,1],rep,num[,2])
names(vnum)<-c("Grupo","Variable","Numero")
for(zz in 3:dim[2]){
rep<-rep(names[zz],dim[1])
v1<-data.frame(media[,1],rep,num[,zz])
names(v1)<-c("Grupo","Variable","Numero")
vnum<-rbind(vnum,v1)
}

##Minimo
dim<-dim(min)
names<-names(min)
rep<-rep(names[2],dim[1])
vmin<-data.frame(min[,1],rep,min[,2])
names(vmin)<-c("Grupo","Variable","Minimo")
for(zz in 3:dim[2]){
rep<-rep(names[zz],dim[1])
v1<-data.frame(min[,1],rep,min[,zz])
names(v1)<-c("Grupo","Variable","Minimo")
vmin<-rbind(vmin,v1)
}


##Maximo
dim<-dim(max)
names<-names(max)
rep<-rep(names[2],dim[1])
vmax<-data.frame(max[,1],rep,max[,2])
names(vmax)<-c("Grupo","Variable","Maximo")
for(zz in 3:dim[2]){
rep<-rep(names[zz],dim[1])
v1<-data.frame(max[,1],rep,max[,zz])
names(v1)<-c("Grupo","Variable","Maximo")
vmax<-rbind(vmax,v1)
}

##Media
dim<-dim(media)
names<-names(media)
rep<-rep(names[2],dim[1])
vmedia<-data.frame(media[,1],rep,media[,2])
names(vmedia)<-c("Grupo","Variable","Media")
for(zz in 3:dim[2]){
rep<-rep(names[zz],dim[1])
v1<-data.frame(media[,1],rep,media[,zz])
names(v1)<-c("Grupo","Variable","Media")
vmedia<-rbind(vmedia,v1)
}


##Mediana
dim<-dim(mediana)
names<-names(mediana)
rep<-rep(names[2],dim[1])
vmedian<-data.frame(mediana[,1],rep,mediana[,2])
names(vmedian)<-c("Grupo","Variable","Mediana")
for(zz in 3:dim[2]){
rep<-rep(names[zz],dim[1])
v1<-data.frame(mediana[,1],rep,mediana[,zz])
names(v1)<-c("Grupo","Variable","Mediana")
vmedian<-rbind(vmedian,v1)
}


##Suma
dim<-dim(suma)
names<-names(suma)
rep<-rep(names[2],dim[1])
vsum<-data.frame(suma[,1],rep,suma[,2])
names(vsum)<-c("Grupo","Variable","Suma")
for(zz in 3:dim[2]){
rep<-rep(names[zz],dim[1])
v1<-data.frame(suma[,1],rep,suma[,zz])
names(v1)<-c("Grupo","Variable","Suma")
vsum<-rbind(vsum,v1)
}


##SD
dim<-dim(sd)
names<-names(sd)
rep<-rep(names[2],dim[1])
vsd<-data.frame(sd[,1],rep,sd[,2])
names(vsd)<-c("Grupo","Variable","SD")
for(zz in 3:dim[2]){
rep<-rep(names[zz],dim[1])
v1<-data.frame(sd[,1],rep,sd[,zz])
names(v1)<-c("Grupo","Variable","SD")
vsd<-rbind(vsd,v1)
}

##CV
dim<-dim(cv)
names<-names(cv)
rep<-rep(names[2],dim[1])
vcv<-data.frame(cv[,1],rep,cv[,2])
names(vcv)<-c("Grupo","Variable","SD")
for(zz in 3:dim[2]){
rep<-rep(names[zz],dim[1])
v1<-data.frame(cv[,1],rep,cv[,zz])
names(v1)<-c("Grupo","Variable","SD")
vcv<-rbind(vcv,v1)
}

final<-cbind(vnum,vmin[,3], vmax[,3],vmedia[,3],vmedian[,3], vsum[,3],vsd[,3], vcv[,3])

names(final)<-names(final)<-c("Cluster","Variable","Number","Minimum","Maximum","Mean","Median","Sum", "SD", "CV")

if(dec=="."){
write.csv(x=final,file = file, fileEncoding = "", row.names=FALSE,na=na)
}
else{
write.csv2(x = final,file = file, fileEncoding = "", row.names=FALSE,na=na)
}


}

}

###End function stat



if(analysis=="CA"){

####Convert to numeric the variables

conv<-NULL
len<-length(var)
for(yy in 1:len){
if(inherits(data[,var[yy]],"factor") | inherits(data[,var[yy]], "character")){
data[,var[yy]]<-as.character(data[,var[yy]])
cate<-unique(data[,var[yy]])
lt<-length(unique(data[,var[yy]]))
se<-seq(1,lt,1)
conv1<-data.frame(var[yy],cate,se)
names(conv1)<-c("Variable", "Category","Numeric value")
if(is.null(conv)) conv<-conv1 else conv<-rbind(conv,conv1)
for(ff in 1:lt){
data[,var[yy]]<-replace(data[,var[yy]], data[,var[yy]]==cate[ff], se[ff])
}
data[,var[yy]]<-as.numeric(data[,var[yy]])
}###End if
}###End for
}


###Selection of variables

if(!is.null(cat)){
if(cat=="Cluster"){
noms<-names(data)
noms<-replace(noms, noms=="Cluster", "ClusterC")
names(data)<-noms
cat<-"ClusterC"
}
}

if(ResetPAR==TRUE){
#Resetear par() a las opciones por defecto
resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}
par(resetPar()) }
else{
}

if(!is.null(PAR)){
parexe<-paste("par(new,",toString(x=PAR), ")")
eval(parse(text=parexe))
}
else{
par(font.lab=2,cex.lab=1.5)
}

bc<-COLORC

datos1<-data.frame(subset(data, select=var))

if(analysis=="PCA"){

###VIF
if(VIF==TRUE){
VIF<-usdm::vif(datos1)
VIF1<-subset(VIF, VIF>threshold)
var<-as.character(unique(VIF1[,1]))
}
else{
VIF<-"No estimated"
}

datos<-data

datosT<-data.frame(subset(datos, select=var))
num<-0

if(!is.null(labels)){
datosT<-data.frame(subset(datos, select=labels),datosT)
num<-num+1
}

if(!is.null(cat)){
datosT<-data.frame(subset(datos, select=cat), datosT)
num<-num+1
}

datos<-na.exclude(datosT)
if(!is.null(cat)){
datos[,1]<-as.character(datos[,1])
}


datos1<-data.frame(subset(datos, select=var))


###kmo test

kmo<-psych::KMO(datos1)


###Bartlett test
bar<-REdaS::bart_spher(datos1)

###Multivariate analysis

if(!is.null(PCA)){
pcaexe<-paste("prcomp(","x=datos1,",toString(x=PCA), ")")
pca<-eval(parse(text=pcaexe))
}
else{
pcaexe<-paste("prcomp(","x=datos1,","scale.=TRUE,", ")")
pca<-eval(parse(text=pcaexe))
}


sum<-summary(pca)

corr<-ltm::rcor.test(datos1)


###File cat loadings

d<-dim(pca$x)

seq<-seq(1,d[1],1)

rownames(pca$x)<-seq

if(!is.null(cat)){
zz<-data.frame(datos[, cat],datos[, labels],pca$x)
colnames(zz)<-c(cat,labels, colnames(pca$x))
catn<-length(unique(zz[,1]))
zz<-zz[order(zz[,1]), ]
dati<-as.character(unique(zz[,1]))
}
else{
zz<-data.frame(datos[, labels],pca$x)
colnames(zz)<-c(labels, colnames(pca$x))
}



}####End Principal Components analysis


if(analysis=="CA"){


datos<-data

datosT<-data.frame(subset(datos, select=var))
num<-0

if(!is.null(labels)){
datosT<-data.frame(subset(datos, select=labels),datosT)
datosZ<-data.frame(subset(datos, select=labels),datosT)
num<-num+1
}

if(!is.null(cat)){
datosT<-data.frame(subset(datos, select=cat), datosT)
num<-num+1
}

datos<-na.exclude(datosT)
datosZ<-na.exclude(datosZ)
if(!is.null(cat)){
datos[,1]<-as.character(datos[,1])
}


datos1<-data.frame(subset(datosZ, select=var))



###Correspondence analysis

zz<-as.character(datosZ[,labels])
datos2<-as.data.frame(datos1,row.names=zz)
corre<-ca::ca(datos2)

sum<-corre


###File cat loadings
zz<-data.frame(datos[, cat],datos[, labels],corre$rowcoord)
colnames(zz)<-c(cat,labels, colnames(corre$rowcoord))
dimz<-dim(zz)
seqz<-seq(1,dimz[1],1)
rownames(zz)<-seqz

if(!is.null(cat)){
catn<-length(unique(zz[,1]))
zz<-zz[order(zz[,1]), ]
dati<-as.character(unique(zz[,1]))
}



}####End Correspondence analysis


###Plot 

if(!is.null(COLOR)){
color1<-COLOR
}
else{
if(!is.null(cat)){
if(convex==TRUE) color1<-rainbow(length(unique(zz[,1])), alpha=0.4) else color1<-rainbow(length(unique(zz[,1])))
}
else{
color1<-"tomato"
}
}


if(!is.null(cat)){
pch1<-14+1
for (t in 2:catn){
pcht<-14+t
pch1<-append(pch1,pcht)
}
}
else{
pch1<-15
}

names<-names(zz)

if(!is.null(PCH)) pcht<-PCH else pcht<-pch1

if(!is.null(XLAB)){
xlab<-XLAB
}
else{
xlab<-names[dim[1]+num]
}


if(!is.null(YLAB)){
ylab<-YLAB
}
else{
ylab<-names[dim[2]+num]
}

if(!is.null(XLIM)){
xlim<-XLIM
}
else{
xlim<-c(min(zz[,dim[1]+num],na.rm=TRUE),max(zz[,dim[1]+num],na.rm=TRUE))
}

if(!is.null(YLIM)){
ylim<-YLIM
}
else{
ylim<-c(min(zz[,dim[2]+num],na.rm=TRUE),max(zz[,dim[2]+num],na.rm=TRUE))
}

dev.new()

if(!is.null(cat)){
by=TRUE
if(!is.null(SCATTERPLOT)){
scatterplotexe<-paste("car::scatterplot(","x=zz[, dim[1]+num],", "y=zz[, dim[2]+num],", "group=zz[,cat],", toString(x=SCATTERPLOT), ")")
eval(parse(text=scatterplotexe))
}
else{
scatterplotexe<-paste("car::scatterplot(","x=zz[, dim[1]+num],", "y=zz[, dim[2]+num],", "group=zz[,cat],", "regLine=FALSE,",
"smooth=FALSE,","grid=FALSE,", "xlim=xlim,", "ylim=ylim,",
"boxplots=FALSE,", "by.groups=by,", "ellipse=ellipse,", "col=color1,", "pch=pcht,",
"xlab=xlab,","ylab=ylab,","legend=list(coords=c(x=500000000,y=500000000))", ")")
eval(parse(text=scatterplotexe))
}

}
else{
by=FALSE
if(!is.null(SCATTERPLOT)){
scatterplotexe<-paste("car::scatterplot(","x=zz[, dim[1]+num],", "y=zz[, dim[2]+num],", toString(x=SCATTERPLOT), ")")
eval(parse(text=scatterplotexe))
}
else{
scatterplotexe<-paste("car::scatterplot(","x=zz[, dim[1]+num],", "y=zz[, dim[2]+num],", "regLine=FALSE,",
"smooth=FALSE,", "grid=FALSE,","xlim=xlim,", "ylim=ylim,",
"boxplots=FALSE,", "by.groups=by,", "ellipse=ellipse,", "col=color1,", "pch=pcht,",
"xlab=xlab,","ylab=ylab,","legend=list(coords=c(x=500000000,y=500000000))", ")")
eval(parse(text=scatterplotexe))
}

}

if(!is.null(cat)){
if(convex==TRUE){
for(zh in 1:catn){
datis<-subset(zz,(zz[,1] == dati[zh]))
hpts <- chull(x=datis[,dim[1]+num], y=datis[,dim[2]+num])
hpts <- c(hpts, hpts[1])
datiss<-datis[hpts,]
polygon(datiss[,dim[1]+num],datiss[,dim[2]+num], col=color1[zh], border=NA)
}
}
}


if(!is.null(cat)){
if(!is.null(LEGEND)){
legendexe<-paste("legend(",toString(x=LEGEND), ")")
eval(parse(text=legendexe))
}
else{
legendexe<-paste("legend(","x='topleft',","legend=dati,", "bty='n',", "col=color1,",
"pch=pcht", ")")
eval(parse(text=legendexe))
}
}


meanx<-(xlim[1]+xlim[2])/2
meany<-(ylim[1]+ylim[2])/2

minar<-min((ylim[2]-ylim[1]),(xlim[2]-xlim[1]))/4

if(analysis=="PCA"){
ma<-as.data.frame(pca$rotation)
}
else{
ma<-as.data.frame(corre$colcoord)
}

dima<-length(rownames(ma))


rx1<-abs(xlim[2]-xlim[1])
ry1<-abs(ylim[2]-ylim[1])

rx2<-abs(max(ma[,dim[1]],na.rm=TRUE)-min(ma[,dim[1]],na.rm=TRUE))
ry2<-abs(max(ma[,dim[2]],na.rm=TRUE)-min(ma[,dim[2]],na.rm=TRUE))

px<-rx1/rx2
ma[,dim[1]]<-ma[,dim[1]]*px*larrow

py<-ry1/ry2
ma[,dim[2]]<-ma[,dim[2]]*py*larrow


if(arrows==TRUE){
for(aa in 1:dima){

x2<-ma[aa, dim[1]]
y2<-ma[aa, dim[2]]

IDPmisc::Arrows(x1=meanx, y1=meany, x2=x2, y2=y2, h.col=colArrows, sh.col=colArrows, h.col.bo=colArrows, open=FALSE)

}#END BUCLE arrows

}#END arrows TRUE



if(!is.null(MTEXT)){
mtextexe<-paste("mtext(",toString(x=MTEXT), ")")
eval(parse(text=mtextexe))
}
else{
}

if(!is.null(TEXTvar)){
textexe<-paste("text(","x=ma[,dim[1]],", "y=ma[,dim[2]]," ,"labels = rownames(ma),", toString(x=TEXTvar), ")")
eval(parse(text=textexe))
}
else{
textexe<-paste("text(","x=ma[,dim[1]],", "y=ma[,dim[2]]," ,"labels = rownames(ma),",
"xlim=xlim,", "ylim=ylim,", "col='black',","cex=1,","pos=3,", "offset=0.25", ")")
eval(parse(text=textexe))
}

if(!is.null(labels)){
if(!is.null(TEXTlabels)){
textexe<-paste("text(","x=zz[,dim[1]+num],", "y=zz[,dim[2]+num]," ,"labels = zz[,labels],", toString(x=TEXTlabels), ")")
eval(parse(text=textexe))
}
else{
textexe<-paste("text(","x=zz[,dim[1]+num],", "y=zz[,dim[2]+num]," ,"labels = zz[,labels],",
"xlim=xlim,", "ylim=ylim,", "col='black',","cex=1,","pos=3,", "offset=0.25", ")")
eval(parse(text=textexe))
}
}


fin<-0

if(kmethod=="automatically"){
k=2
}


while(fin<1){

###Cluster

if(analysis=="PCA"){
if(por==100){
di<-dim(zz)
if(is.null(cat)) corte<-(di[2]-2) else corte<-(di[2]-3)
}
else{
corte<-min(which(sum$importance[3,] > por/100))-1
}

}
else{
if(por==100){
di<-dim(zz)
if(is.null(cat)) corte<-(di[2]-2) else corte<-(di[2]-3)
}
else{
porce<-corre$sv^2*100/sum(corre$sv^2)
lenp<-length(porce)
for(gj in 1:lenp){
if(gj==1) porce[gj]<-porce[gj] else porce[gj]<-porce[gj]+porce[gj-1]
}
corte<-min(which(porce > por))-1
}
}


if(is.null(cat)) ini<-2 else ini<-3
seq<-seq(ini,corte+ini,1)
datos1<-zz[,seq]


if(!is.null(HCLUST)){
hclustexe<-paste("hclust(","d=dist(datos1, method='euclidea'),",toString(x=HCLUST), ")")
hc<-eval(parse(text=hclustexe))
}
else{
hclustexe<-paste("hclust(","d=dist(datos1, method='euclidea'),", ")")
hc<-eval(parse(text=hclustexe))
}


if(kmethod=="manually" | fin>0){
dev.new()


if(!is.null(PAR)){
parexe<-paste("par(new,",toString(x=PAR), ")")
eval(parse(text=parexe))
}
else{
par(font.lab=2,cex.lab=1.5)
}

if(!is.null(CLUSTER)){
hclustexe<-paste("plot(","x=hc,","label=zz[,labels],", toString(x=CLUSTER), ")")
eval(parse(text=hclustexe))
}
else{
hclustexe<-paste("plot(","x=hc,","label=zz[,labels],", "xlab='',", "sub='',", ")")
eval(parse(text=hclustexe))
}
}###End 1 if kmethod


if(!is.null(k)){

if(is.null(COLORC)){
COLORC<-rep("red",k)
}


Cluster <- cutree(hc, k=k)

if(kmethod=="manually" | fin>0){
rect.hclust(hc, k=k, border=COLORC)

if(showCluster==TRUE){
for(rr in 1:k){
ors<-mean(which(Cluster[hc$order]==rr))
text(x=ors, y=max(hc$height)/2.9, label=paste("C", rr, sep=" "), font=2, col=COLORC[rr])
}
}

}###End 2 if kmethod

COLORC<-bc


m<-data.frame(Nueva=as.numeric(rownames(zz)),Cluster) 
m<-m[order(m[,1]),]

nombres<-names(datos)
temp<-data.frame(datos,m[,2])
names(temp)<-c(nombres,"Cluster")


if(kmethod=="manually" | fin>0){


stat(data=temp, var, grupo="Cluster", file=file3)


if(dec=="."){
write.csv(x=temp,file = file4, fileEncoding = "", row.names=FALSE,na=na)
}
else{
write.csv2(x = temp ,file = file4, fileEncoding = "", row.names=FALSE,na=na)
}

}###End 3 if kmethod


}###End null k

len<-length(var)

if(kmethod=="manually" | fin>0){

###Results
if(analysis=="PCA"){
Resultados<-list("VIF values", VIF,  kmo,  bar, "Correlation matrix", corr,  "Summary Multivariate Analysis", sum)
}
else{
Resultados<-list("Numerical labels assigned to variables with factors or characters", conv, "Summary Multivariate Analysis", sum)
}

print(Resultados)

###Save the files

if(!is.null(file)){
sink(file1)
print(Resultados)
sink()
}

if(dec=="."){
write.csv(x=zz,file = file2, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = zz ,file = file2, fileEncoding = "", row.names=row.names,na=na)
}


if(length(var)>25){
print("It is no possible to depict the boxplot due to the high number of variables")
}
else{


###Combined boxplots

if(!is.null(k)){
temp<-temp[order(temp[,"Cluster"]),]
temp[,"Cluster"]<-factor(temp[,"Cluster"], levels = unique(temp[,"Cluster"]), labels = unique(temp[,"Cluster"]))

if(!is.null(LabelCat)) temp[,"Cluster"]<-factor(temp[,"Cluster"], levels = unique(temp[,"Cluster"]), labels = LabelCat) else temp<-temp


dev.new()

if(!is.null(PAR)){
parexe<-paste("par(new,",toString(x=PAR), ")")
eval(parse(text=parexe))
}
else{
par(font.lab=2,cex.lab=1.5)
}

if(is.null(mfrowBOXPLOT)){
div<-ceiling(length(var)/2)
if(div>5){
div<-5
}
ini<-ceiling(len/div)
par(mfrow=c(ini,div))
}
else{
par(mfrow=mfrowBOXPLOT)
}

for(hh in 1:len){

if(!is.null(COLORB)){
color<-COLORB
}
else{
color<-terrain.colors(length(unique(temp[,"Cluster"])))
}

if(is.null(ylabBOXPLOT)){
ylab<-var[hh]
}
else{
ylab<-ylabBOXPLOT[hh]
}


if(!is.null(BOXPLOT)){
boxplotexe<-paste("boxplot(","temp[,var[hh]]~temp[,'Cluster'],",toString(x=BOXPLOT), ")")
eval(parse(text=boxplotexe))
}
else{
boxplotexe<-paste("boxplot(","temp[,var[hh]]~temp[,'Cluster'],",
"xlab='Cluster',", "ylab=ylab,","col=color",")")
eval(parse(text=boxplotexe))
}


}###En loop combined boxplots

}##End null k

}###End var higher than 25

####Discriminant Analysis

DA(data=temp, cat="Cluster", var=var, ellipse=ellipse, convex=convex, quadratic=quadratic,
COLOR=COLORB, PCH=PCH, arrows=arrows, larrow=larrow, colArrows=colArrows, file1=file5,
file2=file6, file3=file7, file4=file8, file5=file9, file6=file10, na=na, dec=dec, row.names=row.names)

}###End 4 if kmethod


#Test U Mann-Whitney

filas<-sum(seq(k-1,1,-1))*len
sum(seq(k,1,-1))
matriz<-matrix(0,filas,4)
matriz<-as.data.frame(matriz)
names(matriz)<-c("Variable","ClusterA","ClusterB","Probability")
ro<-0

if(is.null(LabelCat)){
for(dd in 1:len){
for(gg in 1:(k-1)){
data1<-subset(temp,Cluster == gg)
for(jj in (gg+1):k){
ro<-ro+1
data2<-subset(temp,Cluster == jj)
W<-wilcox.test(data1[,var[dd]],data2[,var[dd]],  exact = FALSE)
matriz[ro,1]<-var[dd]
matriz[ro,2]<-gg
matriz[ro,3]<-jj
matriz[ro,4]<-W$p.value
}
}
}
}
else{
for(dd in 1:len){
for(gg in 1:(k-1)){
data1<-subset(temp,Cluster == LabelCat[gg])
for(jj in (gg+1):k){
ro<-ro+1
data2<-subset(temp,Cluster == LabelCat[jj])
W<-wilcox.test(data1[,var[dd]],data2[,var[dd]])
matriz[ro,1]<-var[dd]
matriz[ro,2]<-gg
matriz[ro,3]<-jj
matriz[ro,4]<-W$p.value
}
}
}
}


if(dec=="."){
write.csv(x=matriz,file = file12, fileEncoding = "", row.names=FALSE,na=na)
}
else{
write.csv2(x = matriz ,file = file12, fileEncoding = "", row.names=FALSE,na=na)
}

matriz[is.na(matriz)]<-1

if(kmethod=="manually"){
fin<-1
}
else{
Probability<-matriz[,4]
matrizP<-subset(matriz,(Probability < pthreshold))
dimP<-dim(matrizP)
if(dimP[1]>0){
datA<-aggregate(x=matrizP[,c("Variable")], by = list(matrizP[ , c("Variable")],matrizP[ , c("ClusterA")]), FUN=count)
names(datA)<-c("Variable","Cluster","Number")
datB<-aggregate(x=matrizP[,c("Variable")], by = list(matrizP[ , c("Variable")],matrizP[ , c("ClusterB")]), FUN=count)
names(datB)<-c("Variable","Cluster","Number")
dta<-rbind(datA, datB)
dta<-aggregate(x=dta[,c("Number")], by = list(dta[ , c("Variable")],dta[ , c("Cluster")]), FUN=sum)
names(dta)<-c("Variable","Cluster","Number")
write.csv2(dta,"dta.csv")
print(k)
jj<-k
for(ff in 1:jj){
dat1<-subset(dta,(Cluster==ff))
texy<-any(dat1[,3]==(k-1))
if(texy==FALSE){
kmethod<-"manually"
}
}

}
else{
fin<-1
k<-k-1
}

if(kmethod=="manually"){
k<-k-1
}
else{
k<-k+1
}



}

}###End while


Probability<-matriz[,4]
matriz<-subset(matriz,(Probability <= pthreshold))
dats<-aggregate(x=matriz[,c("Variable")], by = list(matriz[ , c("ClusterA")],matriz[ , c("ClusterB")]), FUN=count)
names(dats)<-c("ClusterA","ClusterB","Number")


if(BUBBLE==TRUE){

dev.new()

if(!is.null(PAR)){
parexe<-paste("par(new,",toString(x=PAR), ")")
eval(parse(text=parexe))
}
else{
par(font.lab=2,cex.lab=1.5)
}

Bubble(data = dats , varY = "ClusterB" , varX = "ClusterA" , varSize = "Number" , size=size, XLAB="Cluster", YLAB="Cluster", k=k)

}###End BUBBLE TRUE

####VARSEDIG

if(VARSEDIG==TRUE){

temp[,"Cluster"]<-paste("Cluster", temp[,"Cluster"])

temp1<-aggregate(x=temp[,1], by = list(temp[ , "Cluster"]), FUN=count)

if(min(temp1[,2], na.rm=TRUE)>4){
VARSEDIG::VARSEDIM(data=temp, variables=var, group="Cluster", method=method, minimum=minimum, ellipse=ellipse, convex=convex, file=file11,
na=na, dec=dec)
}
else{
print("Due to low number of data in some clusters, VARSEDIG algorithm was not run")
}

}###END VARSEDIG TRUE


}
