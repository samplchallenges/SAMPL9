
########
#SAMPL9#
########

#06.02.23
#William Zamora

library("Metrics")

#Submissions 
#logP_MD.csv, logP_Mixed.csv,and logP_QM.csv 
#reported transfer energies with positive values when they should be negatives and vice versa
#here these values were multiplied by -1

data=read.csv("sampl9_submissions.csv")

#Root mean square deviation (RSMD)
rsmd=NULL
for (i in c(5:length(data))) {
y=data[,4]
r0=round(rmse(y,data[,i]),2)
rsmd=append(rsmd,r0)}

statdata=data.frame(colnames(data)[5:length(data)],rsmd)
colnames(statdata)[1]="Method"

#coefficient of determination (R2),
r2=NULL
for (i in c(5:length(data))) {
y=data[,4]
r0=round((cor(y,data[,i]))^2,2)
r2=append(r2,r0)}

statdata$r2=r2

#mean signed error (MSE)
mse=NULL
for (i in c(5:length(data))) {
y=data[,4]
r0=round(sum((y-data[,i])/length(data[,i])),2)
mse=append(mse,r0)}

statdata$mse=mse

#mean unsigned error (MUE)
mue=NULL
for (i in c(5:length(data))) {
y=data[,4]
r0<-round(mae(data[,i],y),2)
mue=append(mue,r0)}

statdata$mue=mue

#data ordered by RSMD
statdata=statdata[order(statdata$rsmd),]

write.table(statdata,"sampl9_stat.csv",sep=",",row.names=FALSE)



