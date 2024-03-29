VARSEDIG<-function(data, variables, group, group1, group2, method="overlap", stepwise=TRUE, VARSEDIG=TRUE,
minimum=TRUE,kernel="gaussian", cor=TRUE, ellipse=FALSE, convex=TRUE, DPLOT=NULL, SCATTERPLOT=NULL,
BIVTEST12=NULL, BIVTEST21=NULL, Pcol="red", colbiv="lightblue", br=20, sub="", lty=1, lwd=2.5, 
ResetPAR=TRUE,  PAR=NULL,  XLABd=NULL, YLABd=NULL, XLIMd=NULL,
YLIMd=NULL, COLORd=NULL, COLORB=NULL,  LEGENDd=NULL, AXISd=NULL, MTEXTd= NULL, TEXTd=NULL,
XLABs=NULL, YLABs=NULL, XLIMs=NULL, YLIMs=NULL, PCHs=NULL, COLORs=NULL, 
LEGENDs=NULL, MTEXTs= NULL, TEXTs=NULL, LEGENDr=NULL, MTEXTr= NULL, TEXTr=NULL,
arrows=TRUE, larrow=1, ARROWS=NULL, TEXTa=NULL, devnew=TRUE, model="Model.rda", 
file1="Overlap.csv",  file2="Coefficients.csv", file3="Predictions.csv", 
file4="Polar coordinates.csv", file="Output.txt", na="NA", dec=",", row.names=FALSE, removefiles=FALSE){


if(is.null(group) & (ellipse==TRUE | convex==TRUE)) {
stop("It is necessary to specify the group")
}

codl<-length(variables)

valorX12<-"NO"
valorY12<-"NO"
valorX21<-"NO"
valorY21<-"NO"

corte<-"SI"

X12<-0
Y12<-0
TX12<-1
TY12<-1

X21<-0
Y21<-0
TX21<-1
TY21<-1


n<-length(variables)
dist<-0


for(zz in 1:n){

datos<-data

var<-variables[zz]

datosT<-data.frame(subset(datos, select=var), subset(datos, select=group))

datosT<-subset(datosT,(datosT[, group] == group1) |  (datosT[, group] == group2))


datos<-na.exclude(datosT)


dati<-as.character(unique(datos[,2]))

tot<-0

catn<-length(unique(dati))

llty<-length(lty)

if(llty==1){
pch1<-lty
for (t in 2:catn){
pcht<-lty
pch1<-append(pch1,pcht)
}
lty<-pch1
}
else{
}

overlap<-TRUE

if(overlap==TRUE){

Variable1<-"M"
Overlap<-1
Variable2<-"M"


minimo<-min(datos[,1])
maximo<-max(datos[,1])


for(h in 1:catn){
datos3<-subset(datos, datos[,group]==dati[h])
d1<-density(datos3[,1], from=minimo, to=maximo, kernel=kernel)
for(i in 1:catn){
datos3<-subset(datos, datos[,group]==dati[i])
d2<-density(datos3[,1], from=minimo, to=maximo, kernel=kernel)
p<-d1$y-d2$y
p[which(p<0)]=0
p3<-d1$y-p
p3[which(p3<0)]=0
auc1<-kulife::auc(d1$x,d1$y)
auc3<-kulife::auc(d1$x,p3)
P1<-auc3*100/auc1
if(dati[h]==dati[i]){
}
else{
Variable1<-append(Variable1,dati[h])
Variable2<-append(Variable2,dati[i])
Overlap<-append(Overlap,P1)
}
}
}
}

salida<-data.frame(Variable1[-1],Overlap[-1],Variable2[-1],var)

colnames(salida)<-c("Group1","Overlap","Group2","Variable")

if(zz==1) salidaF<-salida else salidaF<-rbind(salidaF,salida)

}


salidaF<-salidaF[order(salidaF[,2]),]


salidaF<-salidaF[!duplicated(salidaF[,4]),]


if(dec=="."){
write.csv(x=salidaF,file = file1, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = salidaF,file = file1, fileEncoding = "", row.names=row.names,na=na)
}



#Density plot

datos<-data

var<-as.character(salidaF[1,4])

datosT<-data.frame(subset(datos, select=var), subset(datos, select=group))

datosT<-subset(datosT,(datosT[, group] == group1) |  (datosT[, group] == group2))

datos<-na.exclude(datosT)

if(devnew==TRUE){
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

}

dati<-as.character(unique(datos[,2]))

if(!is.null(COLORB)){
}
else{
COLORB<-rainbow(length(dati))
}


if(!is.null(COLORd)){
}
else{
COLORd<-rainbow(length(dati),alpha=0.4)
}


if(!is.null(XLABd)) xlab<-XLABd else xlab<-var


if(!is.null(YLABd)){
ylab<-YLABd
}
else{
ylab<-"Density"
}


if(!is.null(PAR)){
parexe<-paste("par(new,",toString(x=PAR), ")")
eval(parse(text=parexe))
}
else{
par(font.lab=2, mar=c(5,5,3,2),cex.lab=1.5)
}


if(!is.null(XLIMd)){
}
else{
maxx<-max(datos[,1])
minx<-min(datos[,1])
XLIMd<-c(minx,maxx)
}

catn<-length(unique(dati))

maxden<-0

for(h in 1:catn){
datos3<-subset(datos, datos[,2]==dati[h])
den<-density(datos3[,1], kernel=kernel)
maxd<-max(den$y)
if(maxd>maxden) maxden<-maxd else maxden<-maxden
}


if(!is.null(YLIMd)){
}
else{

YLIMd<-c(0,maxden)
}

if(devnew==TRUE){
dev.new()
}

if(!is.null(DPLOT)){
tx<-paste("plot.default(","x=0,", "y=0,", "type='n',", toString(x=DPLOT), ")")
eval(parse(text=tx))
}
else{
tx<-paste("plot.default(","x=0,", "y=0,", "type='n',",
"xlim=XLIMd,", "ylim=YLIMd,", "xlab=xlab,", "ylab=ylab,", "main=''", ")")
eval(parse(text=tx))
}


llty<-length(lty)

if(llty==1){
pch1<-lty
for (t in 2:catn){
pcht<-lty
pch1<-append(pch1,pcht)
}
lty<-pch1
}
else{
}


for(h in 1:catn){
datos3<-subset(datos, datos[,group]==dati[h])
den<-density(datos3[,1], kernel=kernel)
polygon(x=den$x, y=den$y, col=COLORd[h],border=COLORB[h], lwd=lwd,lty=lty[h])
}

if(!is.null(LEGENDd)){
legendexe<-paste("legend(",toString(x=LEGENDd), ")")
ifelse(LEGENDd==FALSE, paso<-"NO", eval(parse(text=legendexe)))
}
else{
legendexe<-paste("legend(","x='topleft',","legend=dati,","bty='n',", "col=COLORB,",
"lty=lty", ")")
eval(parse(text=legendexe))
}


if(!is.null(AXISd)){
axisexe<-paste("axis(",toString(x=AXISd), ")")
eval(parse(text=axisexe))
}
else{
}

if(!is.null(MTEXTd)){
mtextexe<-paste("mtext(",toString(x=MTEXTd), ")")
eval(parse(text=mtextexe))
}
else{
}

if(!is.null(TEXTd)){
textexe<-paste("text(",toString(x=TEXTd), ")")
eval(parse(text=textexe))
}
else{
}


if(method=="logistic regression"){

#Logistic regression

vector<-datos[,group]
vector<-vector[!vector == group1] 
cod<-as.character(unique(vector))


#Crear variable binaria

datosT<-data.frame(subset(data, select=group),subset(data, select=variables))
datosT<-subset(datosT,(datosT[, group] == group1) |  (datosT[, group] == group2))
datosT<-na.exclude(datosT)

varbin<-datosT[,group]
varbinaria<-ifelse(varbin==group1,1,0)

#nuevo conjunto de datos

datos4<-data.frame(subset(datosT, select=variables))


datos3<-data.frame(varrespuesta=varbinaria,datos4)

colnames(datos3)<-c(group,colnames(datos4))

#Modelo de regresi�n log�stica.


leng<-length(variables)
fo<-variables[1]
for(i in 2:leng){
fo<-paste(fo,variables[i],sep="+")
}


formula<-as.formula(paste(group,"~",fo))
modelo<-glm(formula,family=binomial("logit"), data=datos3)

if(stepwise==TRUE){
modelo.2<-MASS::stepAIC(modelo, data=datos3)
}
else {
modelo.2<-modelo
}


ModelB<-modelo.2
save(ModelB,file=model)

Resumenmodelo<-(summary(modelo.2))

#Aciertos totales
pred1<-round(predict(modelo.2, type="response"))
prob2<-data.frame(datos3,pred1,pred1)
colnames(prob2)<-c(colnames(datos3),"Probability","Predictions")

dimty<-dim(prob2)
for(i in 1:dimty[1]){
if(prob2[i,"Probability"]>0.4999){
prob2[i, "Predictions"]<-1
}
else{
prob2[i, "Predictions"]<-0
}
}
prctacierto<-sum(prob2[,group]==prob2[,"Predictions"])/ dim(prob2)[1]*100


#Aciertos grupo referencia
datos1<-subset(datos3, datos3[,group]==1)
pred1<-round(predict(modelo.2, type="response", newdata=datos1))
prob2<-data.frame(datos1,pred1,pred1)

colnames(prob2)<-c(colnames(datos1),"Probability","Predictions")
dimty<-dim(prob2)
for(i in 1:dimty[1]){
if(prob2[i,"Probability"]>0.4999){
prob2[i, "Predictions"]<-1
}
else{
prob2[i, "Predictions"]<-0
}
}
prctacierto1<-sum(prob2[,group]==prob2[,"Predictions"])/ dim(prob2)[1]*100


#Aciertos grupo no referencia
datos1<-subset(datos3, datos3[,group]==0)
pred1<-round(predict(modelo.2, type="response", newdata=datos1))

prob2<-data.frame(datos1,pred1,pred1)

colnames(prob2)<-c(colnames(datos1),"Probability","Predictions")
dimty<-dim(prob2)
for(i in 1:dimty[1]){
if(prob2[i,"Probability"]>0.4999){
prob2[i, "Predictions"]<-1
}
else{
prob2[i, "Predictions"]<-0
}
}

prctacierto0<-sum(prob2[,group]==prob2[,"Predictions"])/ dim(prob2)[1]*100

resultadocls<-c("Percentage of cases correctly identified: All cases",prctacierto)

resultadocls1<-c(paste("Percentage of cases correctly identified:",group,group1,collapse=""),prctacierto1)

resultadocls0<-c(paste("Percentage of cases correctly identified:",group, cod,collapse=""),prctacierto0)

### Secci�n que muestra los resultados

coef<-as.data.frame(modelo.2$coefficients)

coef<-data.frame(Variables=rownames(coef), coef)

colnames(coef)<-c("Variables","Coefficients")

prob<-(predict(modelo.2, type="response"))

prob2<-data.frame(datosT,prob,prob)

colnames(prob2)<-c(colnames(datosT),"Probability","Predictions")

dimty<-dim(prob2)


for(i in 1:dimty[1]){
if(prob2[i,"Probability"]>0.4999){
prob2[i, "Predictions"]<-as.character(group1)
}
else{
prob2[i, "Predictions"]<-cod
}


}


if(dec=="."){
write.csv(x=coef,file = file2, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = coef,file = file2, fileEncoding = "", row.names=row.names,na=na)
}


if(dec=="."){
write.csv(x=prob2,file = file3, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = prob2,file = file3, fileEncoding = "", row.names=row.names,na=na)
}


}#End logistic regression


#Buscar la probabilidad menor con randomization test

lenvar<-length(variables)

for(rr in 1:lenvar){

datos<-data

var<-variables[rr]

datosT<-data.frame(subset(datos, select=var), subset(datos, select=group))

datosT<-subset(datosT,(datosT[, group] == group1) |  (datosT[, group] == group2))

datosP<-na.exclude(datosT)


grupo1<-subset(datosP,(datosP[, group] == group1))

grupo2<-subset(datosP,(datosP[, group] == group2))


X1<-mean(grupo1[,1])
X2<-mean(grupo2[,1])



#Comparing all values of group 2 with group 1

wz<-dim(grupo2)

tp1<-0

for(ww in 1:wz[1]){#bucle para grupo 2

if(X1>X2) randX<-ade4::as.randtest(as.numeric(grupo1[,1]), as.numeric(grupo2[ww,1]), alter="less") else randX<-ade4::as.randtest(as.numeric(grupo1[,1]), as.numeric(grupo2[ww,1]), alter="greater")


tp1<-tp1+randX$pvalue


}#fin bucle grupo 2


tp1<-tp1/wz

#Comparing all values of group 2 with group 1

wz<-dim(grupo1)

tp2<-0

for(ww in 1:wz[1]){#bucle para grupo 1


if(X2>X1) randX<-ade4::as.randtest(as.numeric(grupo2[,1]), as.numeric(grupo1[ww,1]), alter="less") else randX<-ade4::as.randtest(as.numeric(grupo2[,1]), as.numeric(grupo1[ww,1]), alter="greater")

tp2<-tp2+randX$pvalue


}#fin bucle grupo 1

tp2<-tp2/wz

Prob<-tp1*tp2

salida<-data.frame(variables[rr],Prob)

colnames(salida)<-c("Variable","Probability")

if(rr==1) salidaR<-salida else salidaR<-rbind(salidaR,salida)

}#fin blucle variables

salidaR<-salidaR[order(salidaR[,2]),]

salidaR<-salidaR[!duplicated(salidaR[,1]),]

#Polar coordinates


if(method=="logistic regression") Variables<-as.character(coef[-1,1]) else Variables<-as.character(salidaF[,4])

if(method=="Monte-Carlo") Variables<-as.character(salidaR[,1]) else Variables<-Variables


if(method=="logistic regression"){
}
else{
Resumenmodelo<-""
resultadocls<-""
resultadocls1<-""
resultadocls0<-""
}


n<-length(Variables)

datos<-data



for(zz in 1:n){#Bucle coordenadas polares


datos<-data

if(zz==1) var<-Variables[1] else var<-append(var,Variables[zz])

datosT<-data.frame(subset(datos, select=var), subset(datos, select=group))

datosT<-subset(datosT,(datosT[, group] == group1) |  (datosT[, group] == group2))




datos<-na.exclude(datosT)

#Organize by correlation

selection<-subset(datos, select=var)


dimdt<-dim(selection)

if(cor==TRUE){
if(dimdt[2]>3){
datosT2<-selection
datosT3<-selection[,1]

names<-names(selection)[1]

for(zz in 1:(dimdt[2]-2)){
cor<-cor(datosT2)
corT1<-cor[,-1]
datosT2<-datosT2[,-1]

df<-which.max(cor[!cor[,1]==1,1])

if(zz==(dimdt[2]-2)){
dimdim<-dim(datosT3)
cor1<-cor(datosT3[,dimdim[2]], datosT2[,1])
cor2<-cor(datosT3[,dimdim[2]], datosT2[,2])
if(cor1>cor2){
datosT3<-data.frame(datosT3, datosT2[,1],datosT2[,2])
names(datosT3)<-c(names,names(datosT2))
}
else{
datosT3<-data.frame(datosT3, datosT2[,2],datosT2[,1])
kk<-datosT2
kk<-data.frame(datosT2[,2],datosT2[,1])
names(kk)<-c(names(datosT2)[2],names(datosT2)[1])
names(datosT3)<-c(names,names(kk))
} 

}
else{
datosT3<-data.frame(datosT3, datosT2[,df])
names(datosT3)<-c(names,names(datosT2)[df])
}

names<-names(datosT3)


datosP<-datosT2[,-df]


datosZZZ<-datosT2
datosT2<-data.frame(datosT2[,names(datosT2)[df]], datosP)
names(datosT2)<-c(names(datosZZZ)[df], colnames(datosP))
}

selection<-datosT3
}
var<-colnames(selection)
}#end cor





#Standardization 1 to -1

a<-dim(selection)

datosE<-selection

for (z in 1:a[2]){
matrixE<-matrix(c(-1, 1, min(selection[,z],na.rm=TRUE),max(selection[,z],na.rm=TRUE)), nrow = 2 , ncol = 2)
reg<-lm(matrixE[,1]~matrixE[,2])
datosC<-reg$coefficients[1]+selection[,z]*reg$coefficients[2]
datosE<-cbind(datosE,datosC)
}


if(zz==1) datosE<-as.data.frame(datosE[,-1]) else datosE<-datosE[,-c(1:a[2])]


if(zz==1) colnames(datosE)<-var else colnames(datosE)<-colnames(selection)

selection<-datosE



ggg<-0

ggg<-ggg+1

datosT<-data.frame(subset(datos, select=group),selection)

datos<-datosT


#Estimation of polar coordinates

angle<-pi/a[2]

if(zz==1){
datosX<-datos[,1]
datosY<-datos[,1]
}
else{
datosX<-datos[,1:ggg]
datosY<-datos[,1:ggg]
}

h<-0
for (z in (1+ggg):(a[2]+ggg)){
h<-h+1
datosC<- ifelse(datos[,z] <=0, abs(datos[,z])*cos(angle*h+pi), abs(datos[,z])*cos(angle*h)) 
datosX<-cbind(datosX,datosC)
datosC<- ifelse(datos[,z] <=0, abs(datos[,z])*sin(angle*h+pi), abs(datos[,z])*sin(angle*h)) 
datosY<-cbind(datosY,datosC)
}

if(zz==1){
XX<-datosC
YY<-datosC 
}
else{
XX<-apply(datosX[,(ggg+1):(a[2]+ggg)],1,sum)
YY<-apply(datosY[,(ggg+1):(a[2]+ggg)],1,sum)
}
RX<-(max(XX)-min(XX))

RY<-(max(YY)-min(YY))


if(zz==1) datosF<-data.frame(datos[,1],XX,YY) else datosF<-data.frame(datos[,1:ggg],XX,YY)


colnames(datosF)<-c("Group","X","Y")


grupo1<-subset(datosF,(datosF[, "Group"] == group1))

grupo2<-subset(datosF,(datosF[, "Group"] == group2))


X1<-mean(grupo1[,"X"])
X2<-mean(grupo2[,"X"])
Y1<-mean(grupo1[,"Y"])
Y2<-mean(grupo2[,"Y"])

distA<-sqrt((X1-X2)^2+(Y1-Y2)^2)

if(!is.na(distA)) distA<-distA else distA<-0

#Comparing all values of group 2 with group 1

wz<-dim(grupo2)


if(distA>dist){


for(ww in 1:wz[1]){#bucle para grupo 2

if(X1>X2) randX<-ade4::as.randtest(as.numeric(grupo1[,2]), as.numeric(grupo2[ww,2]), alter="less") else randX<-ade4::as.randtest(as.numeric(grupo1[,2]), as.numeric(grupo2[ww,2]), alter="greater")
if(Y1>Y2) randY<-ade4::as.randtest(as.numeric(grupo1[,3]), as.numeric(grupo2[ww,3]), alter="less") else randY<-ade4::as.randtest(as.numeric(grupo1[,3]), as.numeric(grupo2[ww,3]), alter="greater")


if(X12<=randX$pvalue){
X12<-randX$pvalue
wx<-ww
}

if(Y12<=randY$pvalue){
Y12<-randY$pvalue
wy<-ww
}


}#fin bucle grupo 2

if(X1>X2) randX<-ade4::as.randtest(as.numeric(grupo1[,2]), as.numeric(grupo2[wx,2]), alter="less") else randX<-ade4::as.randtest(as.numeric(grupo1[,2]), as.numeric(grupo2[wx,2]), alter="greater")
if(Y1>Y2) randY<-ade4::as.randtest(as.numeric(grupo1[,3]), as.numeric(grupo2[wx,3]), alter="less") else randY<-ade4::as.randtest(as.numeric(grupo1[,3]), as.numeric(grupo2[wx,3]), alter="greater")

X12<-min(randX$pvalue,randY$pvalue)

if(X1>X2) randX<-ade4::as.randtest(as.numeric(grupo1[,2]), as.numeric(grupo2[wy,2]), alter="less") else randX<-ade4::as.randtest(as.numeric(grupo1[,2]), as.numeric(grupo2[wy,2]), alter="greater")
if(Y1>Y2) randY<-ade4::as.randtest(as.numeric(grupo1[,3]), as.numeric(grupo2[wy,3]), alter="less") else randY<-ade4::as.randtest(as.numeric(grupo1[,3]), as.numeric(grupo2[wy,3]), alter="greater")

Y12<-min(randX$pvalue,randY$pvalue)



if(minimum==FALSE){
if(TX12>=X12){
valorX12<-"NO"
TX12<-X12
}
}
else{
if(TX12>X12){
valorX12<-"NO"
TX12<-X12
}
}




#Comparing all values of group 1 with group 2

wz<-dim(grupo1)


for(ww in 1:wz[1]){#bucle para grupo 1

if(X2>X1) randX<-ade4::as.randtest(as.numeric(grupo2[,2]), as.numeric(grupo1[ww,2]), alter="less") else randX<-ade4::as.randtest(as.numeric(grupo2[,2]), as.numeric(grupo1[ww,2]), alter="greater")
if(Y2>Y1) randY<-ade4::as.randtest(as.numeric(grupo2[,3]), as.numeric(grupo1[ww,3]), alter="less") else randY<-ade4::as.randtest(as.numeric(grupo2[,3]), as.numeric(grupo1[ww,3]), alter="greater")

if(X21<=randX$pvalue){
X21<-randX$pvalue
wx<-ww
}

if(Y21<=randY$pvalue){
Y21<-randY$pvalue
wy<-ww
}


}#fin bucle grupo 1


if(X2>X1) randX<-ade4::as.randtest(as.numeric(grupo2[,2]), as.numeric(grupo1[wx,2]), alter="less") else randX<-ade4::as.randtest(as.numeric(grupo2[,2]), as.numeric(grupo1[wx,2]), alter="greater")
if(Y2>Y1) randY<-ade4::as.randtest(as.numeric(grupo2[,3]), as.numeric(grupo1[wx,3]), alter="less") else randY<-ade4::as.randtest(as.numeric(grupo2[,3]), as.numeric(grupo1[wx,3]), alter="greater")

X21<-min(randX$pvalue,randY$pvalue)

if(X2>X1) randX<-ade4::as.randtest(as.numeric(grupo2[,2]), as.numeric(grupo1[wy,2]), alter="less") else randX<-ade4::as.randtest(as.numeric(grupo2[,2]), as.numeric(grupo1[wy,2]), alter="greater")
if(Y2>Y1) randY<-ade4::as.randtest(as.numeric(grupo2[,3]), as.numeric(grupo1[wy,3]), alter="less") else randY<-ade4::as.randtest(as.numeric(grupo2[,3]), as.numeric(grupo1[wy,3]), alter="greater")

Y21<-min(randX$pvalue,randY$pvalue)


corte<-"SI"

if(TX12==X12 & TY12==Y12 & TX21==X21 & TY21==Y21){
corte<-"NO"
}

if(corte=="SI"){
if(TX12>=X12 & TY12>=Y12 & TX21>=X21 & TY21>=Y21){
valorX12<-"NO"
TX12<-X12
valorY12<-"NO"
TY12<-Y12

valorX21<-"NO"
TX21<-X21
valorY21<-"NO"
TY21<-Y21
}
}


if(minimum==FALSE){
if(X12<0.07 | Y12<0.07){
valorX12<-"NO"
TX12<-X12
valorY12<-"NO"
TY12<-Y12
}
}




if(minimum==FALSE){
if(X21<0.07 | Y21<0.07){
valorX21<-"NO"
TX21<-X21
valorY21<-"NO"
TY21<-Y21
}
}


}#final distA


if(VARSEDIG==TRUE){
if(valorX12=="NO" & valorY12=="NO" & valorX21=="NO" & valorY21=="NO"){
X12<-0
Y12<-0
X21<-0
Y21<-0


valorX12<-"SI"
valorY12<-"SI"
valorX21<-"SI"
valorY21<-"SI"

dist<-distA

if(dec=="."){
write.csv(x=datosF,file = file4, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = datosF,file = file4, fileEncoding = "", row.names=row.names,na=na)
}
}

else{
X12<-0
Y12<-0
X21<-0
Y21<-0


valorX12<-"SI"
valorY12<-"SI"
valorX21<-"SI"
valorY21<-"SI"

var<-var[-length(var)]
}
}#TRUE VARSEDIG
else{

if(dec=="."){
write.csv(x=datosF,file = file4, fileEncoding = "", row.names=row.names,na=na)
}
else{
write.csv2(x = datosF,file = file4, fileEncoding = "", row.names=row.names,na=na)
}

dist<-distA

X12<-0
Y12<-0
X21<-0
Y21<-0


valorX12<-"SI"
valorY12<-"SI"
valorX21<-"SI"
valorY21<-"SI"

}






}#Fin bucle coordenadas polares




if(method=="logistic regression"){
if(devnew==TRUE){
print(Resumenmodelo)
print(resultadocls)
print(resultadocls1)
print(resultadocls0)
}
}
else{
Resumenmodelo<-""
resultadocls<-""
resultadocls1<-""
resultadocls0<-""
}

if(dec=="."){
datosF<-read.csv(file=file4 ,header=TRUE)
}
else{
datosF<-read.csv2(file=file4 ,header=TRUE)
}

datosF<-datosF[order(datosF[,1]), ]


if(!is.null(COLORs)){
color1<-COLORs
}
else{
if(convex==TRUE) color1<-rainbow(length(unique(datosF[,1])), alpha=0.4) else color1<-rainbow(length(unique(datosF[,1])))
}


catn<-length(unique(datosF[,1]))

pch1<-14+1
for (t in 2:catn){
pcht<-14+t
pch1<-append(pch1,pcht)
}

if(!is.null(PCHs)) pcht<-PCHs else pcht<-pch1





if(devnew==TRUE){
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
}



if(!is.null(XLABs)) xlab<-XLABs else xlab<-"POLAR COORDINATES X"


if(!is.null(YLABs)){
ylab<-YLABs
}
else{
ylab<-"POLAR COORDINATES Y"
}

if(devnew==TRUE){
dev.new()
}

if(!is.null(PAR)){
parexe<-paste("par(new,",toString(x=PAR), ")")
eval(parse(text=parexe))
}
else{
par(font.lab=2, mar=c(5,5,3,2),cex.lab=1.5)
}

datosF[,1]<-as.factor(datosF[,1])

if(!is.null(SCATTERPLOT)){
scatterplotexe<-paste("car::scatterplot(","Y~X| Group,", "data=datosF,", toString(x=SCATTERPLOT), ")")
eval(parse(text=scatterplotexe))
}
else{
scatterplotexe<-paste("car::scatterplot(","Y~X| Group,", "data=datosF,", "regLine=FALSE,",
"smooth=FALSE,",  "grid=FALSE,", "xlab=xlab,","ylab=ylab,","xlim=XLIMs,","ylim=YLIMs,",
"boxplots=FALSE,", "by.groups=TRUE,", "ellipse=ellipse,", "col=color1,", "pch=pcht,","legend=list(coords=c(x=500000000,y=500000000))", ")")
eval(parse(text=scatterplotexe))
}

dati<-as.character(unique(datosF[,1]))

if(!is.null(group)){
if(convex==TRUE){
for(zz in 1:catn){
datis<-subset(datosF,(datosF[,"Group"] == dati[zz]))
hpts <- chull(x=datis[,"X"], y=datis[,"Y"])
hpts <- c(hpts, hpts[1])
datiss<-datis[hpts,]
polygon(datiss[,2],datiss[,3], col=color1[zz], border=NA)
}
}
}



if(!is.null(TEXTs)){
textexe<-paste("text(",toString(x=TEXTs), ")")
eval(parse(text=textexe))
}
else{
}


if(!is.null(LEGENDs)){
legendexe<-paste("legend(",toString(x=LEGENDs), ")")
eval(parse(text=legendexe))
}
else{
legendexe<-paste("legend(","x='topleft',","legend=unique(datosF[,'Group']),", "bty='n',", "col=color1,",
"pch=pcht", ")")
eval(parse(text=legendexe))
}




if(!is.null(MTEXTs)){
mtextexe<-paste("mtext(",toString(x=MTEXTs), ")")
eval(parse(text=mtextexe))
}
else{
}



#ARROWS

if(!is.null(XLIMs)){
minsx<-XLIMs[1]
maxsx<-XLIMs[2]
}
else{
minsx<-min(datosF[,2])
maxsx<-max(datosF[,2])
}

if(!is.null(YLIMs)){
minsy<-YLIMs[1]
maxsy<-YLIMs[2]
}
else{
minsy<-min(datosF[,3])
maxsy<-max(datosF[,3])
}

meanx<-(minsx+maxsx)/2
meany<-(minsy+maxsy)/2

minar<-min((maxsy-minsy),(maxsx-minsx))/4

dima<-length(var)

angle<-pi/dima


if(arrows==TRUE){
for(aa in 1:dima){
ag<-angle*aa
x2<-minar*cos(ag)*larrow+meanx
y2<-minar*sin(ag)*larrow+meany
if(!is.null(ARROWS)){
arrowsexe<-paste("IDPmisc::Arrows(","x1=meanx,", "y1=meany,", "x2=-x2,", "y2=y2,", toString(x=ARROWS), ")")
eval(parse(text=arrowsexe))
}
else{
arrowsexe<-paste("IDPmisc::Arrows(","x1=meanx,", "y1=meany,", "x2=x2,", "y2=y2,", "open=FALSE", ")")
eval(parse(text=arrowsexe))
}

if(!is.null(TEXTa)){
textexe<-paste("text(","x=x2,", "y=y2,", "labels = paste('+',var[aa], sep=''),", "pos=3,", toString(x=TEXTa), ")")
eval(parse(text=textexe))
}
else{
textexe<-paste("text(","x=x2,", "y=y2,", "labels = paste('+',var[aa], sep=''),", "pos=3", ")")
eval(parse(text=textexe))
}

}#END BUCLE arrows

for(aa in 1:dima){
ag<-angle*aa
x2<-minar*cos(ag+pi)*larrow+meanx
y2<-minar*sin(ag+pi)*larrow+meany
if(!is.null(ARROWS)){
arrowsexe<-paste("IDPmisc::Arrows(","x1=meanx,", "y1=meany,", "x2=-x2,", "y2=y2,", toString(x=ARROWS), ")")
eval(parse(text=arrowsexe))
}
else{
arrowsexe<-paste("IDPmisc::Arrows(","x1=meanx,", "y1=meany,", "x2=x2,", "y2=y2,", "open=FALSE", ")")
eval(parse(text=arrowsexe))
}

if(!is.null(TEXTa)){
textexe<-paste("text(","x=x2,", "y=y2,", "labels = paste('-',var[aa], sep=''),", "pos=1,", toString(x=TEXTa), ")")
eval(parse(text=textexe))
}
else{
textexe<-paste("text(","x=x2,", "y=y2,", "labels = paste('-',var[aa], sep=''),", "pos=1", ")")
eval(parse(text=textexe))
}

}#END BUCLE arrows

}#END arrows TRUE





#Random test


if(dec=="."){
datosF<-read.csv(file=file4 ,header=TRUE)
}
else{
datosF<-read.csv2(file=file4 ,header=TRUE)
}



grupo1<-subset(datosF,(datosF[, "Group"] == group1))

grupo2<-subset(datosF,(datosF[, "Group"] == group2))


X1<-mean(grupo1[,"X"])
X2<-mean(grupo2[,"X"])
Y1<-mean(grupo1[,"Y"])
Y2<-mean(grupo2[,"Y"])

matrix1<-data.frame(0,0)
colnames(matrix1)<-c("po","ED")

wz<-dim(grupo2)

for(ww in 1:wz[1]){

distP<-sqrt((X1-grupo2[ww,"X"])^2+(Y1-grupo2[ww,"Y"])^2)

matrix1<-rbind(matrix1,c(ww,distP))

}

matrix1<-matrix1[-1,]


matrix1<-na.exclude(matrix1[order(matrix1[,2], decreasing = FALSE),])




po1<-matrix1[1,1]

point1 <- c(grupo2[po1,2], grupo2[po1,3])



if(devnew==TRUE){
dev.new()
}

if(!is.null(BIVTEST12)){
tx<-paste("adehabitatHS::biv.test(","dfxy=grupo1[,2:3],", "point=point1,",  toString(x=BIVTEST12), ")")
eval(parse(text=tx))
}
else{
tx<-paste("adehabitatHS::biv.test(","dfxy=grupo1[,2:3],", "point=point1,", "br=br,",
"Pcol=Pcol,", "col=colbiv,", "sub=sub", ")")
eval(parse(text=tx))
}




dati<-c(group2,group1)
col<-c(Pcol, colbiv)

if(!is.null(LEGENDr)){
legendexe<-paste("legend(",toString(x=LEGENDr), ")")
ifelse(LEGENDr==FALSE, paso<-"NO", eval(parse(text=legendexe)))
}
else{
legendexe<-paste("legend(","x='topright',","legend=dati,","bty='n',", "col=col,","pch=c(16,16)",  ")")
eval(parse(text=legendexe))
}


if(!is.null(MTEXTr)){
mtextexe<-paste("mtext(",toString(x=MTEXTr), ")")
eval(parse(text=mtextexe))
}
else{
}

if(!is.null(TEXTr)){
textexe<-paste("text(",toString(x=TEXTr), ")")
eval(parse(text=textexe))
}
else{
}



matrix2<-data.frame(0,0)
colnames(matrix2)<-c("po","ED")


wz<-dim(grupo1)

for(ww in 1:wz[1]){

distP<-sqrt((X2-grupo1[ww,"X"])^2+(Y2-grupo1[ww,"Y"])^2)

matrix2<-rbind(matrix2,c(ww,distP))

}

matrix2<-matrix2[-1,]

matrix2<-na.exclude(matrix2[order(matrix2[,2], decreasing = FALSE),])



po2<-matrix2[1,1]

point2 <- c(grupo1[po2,2], grupo1[po2,3])


if(devnew==TRUE){
dev.new()
}

if(!is.null(BIVTEST21)){
tx<-paste("adehabitatHS::biv.test(","dfxy=grupo2[,2:3],", "point=point2,",  toString(x=BIVTEST21), ")")
eval(parse(text=tx))
}
else{
tx<-paste("adehabitatHS::biv.test(","dfxy=grupo2[,2:3],", "point=point2,", "br=br,",
"Pcol=Pcol,", "col=colbiv,", "sub=sub,", ")")
eval(parse(text=tx))
}



dati<-c(group1,group2)
col<-c(Pcol, colbiv)

if(!is.null(LEGENDr)){
legendexe<-paste("legend(",toString(x=LEGENDr), ")")
ifelse(LEGENDr==FALSE, paso<-"NO", eval(parse(text=legendexe)))
}
else{
legendexe<-paste("legend(","x='topright',","legend=dati,","bty='n',", "col=col,","pch=c(16,16)",  ")")
eval(parse(text=legendexe))
}


if(!is.null(MTEXTr)){
mtextexe<-paste("mtext(",toString(x=MTEXTr), ")")
eval(parse(text=mtextexe))
}
else{
}

if(!is.null(TEXTr)){
textexe<-paste("text(",toString(x=TEXTr), ")")
eval(parse(text=textexe))
}
else{
}


if(devnew==TRUE){
cat(c("Value of ", group2,"with highest probability to belong to", group1), "\n")
print("")
cat(c("Row:", po1, "Values: ", point1), "\n")
print("")
cat(c("Value of ", group1,"with highest probability to belong to", group2), "\n")
print("")
cat(c("Row:", po2, "Values: ", point2), "\n")
print("")
print("VARIABLES SELECTED")
print(var)
print("")
print("EUCLIDEAN DISTANCE BETWEEN GROUPS")
print(dist)
}


if(method=="logistic regression"){

if(devnew==TRUE){
if(!is.null(file)){
sink(file)
 print(Resumenmodelo)
 print(resultadocls)
 print(resultadocls1)
 print(resultadocls0)
cat(c("Value of ", group2,"with highest probability to belong to", group1), "\n")
print("")
cat(c("Row:", po1, "Values: ", point1), "\n")
print("")
cat(c("Value of ", group1,"with highest probability to belong to", group2), "\n")
print("")
cat(c("Row:", po2, "Values: ", point2), "\n")
print("")
print("VARIABLES SELECTED")
print(var)
print("")
print("EUCLIDEAN DISTANCE BETWEEN GROUPS")
print(dist)
sink()
}
}
}
else{
if(devnew==TRUE){
if(!is.null(file)){
sink(file)
cat(c("Value of ", group2,"with highest probability to belong to", group1), "\n")
print("")
cat(c("Row:", po1, "Values: ", point1), "\n")
print("")
cat(c("Value of ", group1,"with highest probability to belong to", group2), "\n")
print("")
cat(c("Row:", po2, "Values: ", point2), "\n")
print("")
print("VARIABLES SELECTED")
print(var)
print("")
print("EUCLIDEAN DISTANCE BETWEEN GROUPS")
print(dist)
sink()
}
}
}

if(removefiles==TRUE){

if(file.exists(file1)==TRUE){
file.remove(file1)
}

if(file.exists(file2)==TRUE){
file.remove(file2)
}

if(file.exists(file3)==TRUE){
file.remove(file3)
}

if(file.exists(file4)==TRUE){
file.remove(file4)
}

if(file.exists(file)==TRUE){
file.remove(file)
}

}


}
