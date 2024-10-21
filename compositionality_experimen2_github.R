#test of compositionality corrections
setwd("C:/Users/brandonh/Documents/analysis/compositionality")
library(ggpubr)
library(gridExtra)
library(edgeR)
library(compositions)
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(vegan)

#save.image("compos.Rdata")
#load("compos.Rdata")

#save.image("compos_ex2.Rdata")
load("compos_ex2.Rdata")

meta = data.frame(timepoint = rep(c("1","2"),each=20),
                  treatment = rep(c("c","tr","c","tr"),each=10))

absdata <- data.frame(s1=round(exp(rnorm(100,50,10)*0.25)))
#absdata <- data.frame(s1=round(exp(rnorm(100,10,5)*0.5)))
for(i in 2:20){
absdata[,paste("s",i,sep="")] <- absdata$s1 * exp(rnorm(100,0,0.5))
}
absdata[absdata<0]<-0
colnames(absdata) <- paste(colnames(absdata),rep(c("c","tr"),each=10),sep="_")

windows();barplot(as.matrix(absdata[rev(order(rowMeans(absdata))),]),border=F,cex.names = 0.8,col=1:100,las=2)


test <- data.frame(s1=round(exp(rnorm(100,50,10)*0.25)))
#test <- data.frame(s1=round(exp(rnorm(100,10,5)*0.5)))
for(i in 2:20){
  test[,paste("s",i,sep="")] <- test$s1 * exp(rnorm(100,0,0.5)) * exp(rnorm(100,0,1))
}
test[test<0]<-0
colnames(test) <- paste(colnames(test),rep(c("c","tr"),each=10),sep="_")

windows();barplot(as.matrix(test[rev(order(rowMeans(test))),]),border=F,cex.names = 0.8,col=1:100,las=2)

windows();barplot(round(exp(rnorm(100,100,10)*0.25)),border=F,cex.names = 0.8,col=1:100,las=2)

absdata <-test

#windows();plot(hclust(as.dist(1-cor(log(absdata+1)))))

palette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F",
          
          brewer.pal(9,"Greens")[-1],#"plum4","plum","plum1","violet","purple1","purple3","purple4",
          brewer.pal(9,"Purples")[-1],#brewer.pal(9,"Reds"),brewer.pal(9,"YlOrBr"),
         # brewer.pal(9,"Greys"),brewer.pal(9,"PuBuGn"),
          brewer.pal(9,"YlGnBu")[-1],rep(c("gray90","gray40"),80)))
windows();barplot(as.matrix(absdata[rev(order(rowMeans(absdata))),]),border=F,cex.names = 0.5,col=1:100)
windows();barplot(as.matrix(t(t(absdata[rev(order(rowMeans(absdata))),])/colSums(absdata[rev(order(rowMeans(absdata))),]))),border=F,cex.names = 0.5,col=1:100)

treatmentdata <- cbind(absdata,absdata)
for(i in 21:40){
  treatmentdata[,i]<-round(treatmentdata[,i]*exp(rnorm(100,0,sd=0.25)))#0.1
}
colnames(treatmentdata) <- paste(colnames(treatmentdata),rep(c("1","2"),each=20),sep="_")

windows();barplot(as.matrix(treatmentdata[rev(order(rowMeans(treatmentdata))),]),border=F,cex.names = 0.5,col=1:100,las=2)
windows();barplot(as.matrix(t(t(treatmentdata[rev(order(rowMeans(treatmentdata))),])/colSums(treatmentdata[rev(order(rowMeans(treatmentdata))),]))),border=F,cex.names = 0.5,col=1:100)

palette(c("#053061", "#2166AC", "#4393C3", "#F4A582", "#B2182B", "plum4","violet","purple1",brewer.pal(3,"Greens")[-1]))

#palette(c(brewer.pal(9,"Set1"),"#053061"))

windows(width=3,height = 4.5);
plot(log(t(treatmentdata)[meta$treatment=="c",4])~meta$timepoint[meta$treatment=="c"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Control group",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines(log(t(treatmentdata)[meta$treatment=="c",4][c(i,i+10)])~meta$timepoint[meta$treatment=="c"][c(i,i+10)],
        col=i,lwd=2)
}

windows(width=3,height = 4.5);
plot(log(t(treatmentdata)[meta$treatment=="tr",3])~meta$timepoint[meta$treatment=="tr"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Control group",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines(log(t(treatmentdata)[meta$treatment=="tr",3][c(i,i+10)])~meta$timepoint[meta$treatment=="tr"][c(i,i+10)],
        col=i,lwd=2)
}


windows();plot(hclust(as.dist(1-cor(log(treatmentdata+1)))))
windows();plot(hclust(as.dist(1-cor((treatmentdata)))))
#0.5,1,2,5,10

treatmentdatalist <- list()
for(k in c(1:100)){ 
for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){    
for(i in c(0.1,0.5)){
    treatmentdatalist[[paste(k,i,j)]] <- treatmentdata
    treatmentdatalist[[paste(k,i,j)]][1:round(100*(j/2)),31:40] <- round(treatmentdatalist[[paste(k,i,j)]][1:round(100*(j/2)),31:40]*(1+rnorm(n=10,i,i/2)))
    treatmentdatalist[[paste(k,i,j)]][floor(100*(j/2)+1):(j*100),31:40] <- round(treatmentdatalist[[paste(k,i,j)]][floor(100*(j/2)+1):(j*100),31:40]*(1-rnorm(n=10,i,i/2)))
    treatmentdatalist[[paste(k,i,j)]][treatmentdatalist[[paste(k,i,j)]]<0] <- 1
  }  
  for(i in c(1,5,10)){
    treatmentdatalist[[paste(k,i,j)]] <- treatmentdata
    treatmentdatalist[[paste(k,i,j)]][1:round(100*(j/2)),31:40] <- round(treatmentdatalist[[paste(k,i,j)]][1:round(100*(j/2)),31:40]*(1+rnorm(n=10,i,i/2)))
    treatmentdatalist[[paste(k,i,j)]][floor(100*(j/2)+1):(j*100),31:40] <- round(treatmentdatalist[[paste(k,i,j)]][floor(100*(j/2)+1):(j*100),31:40]/(1+rnorm(n=10,i,i/2)))
    treatmentdatalist[[paste(k,i,j)]][treatmentdatalist[[paste(k,i,j)]]<0] <- 1
     }

}
}

for(k in c(1:100)){ 
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){    
    for(i in c(0.1,0.5,1,5,10)){
      treatmentdatalist[[paste(k,i,j)]][treatmentdatalist[[paste(k,i,j)]]<0]<-0   
    }}}


windows(height=4);barplot(as.matrix(treatmentdata[rev(order(rowMeans(treatmentdata))),21:40]),border=F,cex.names = 0.5,col=1:100,las=2)

windows(height=4);barplot(as.matrix(treatmentdatalist[['3 0.1 0.01']][rev(order(rowMeans(treatmentdatalist[['3 0.1 0.01']]))),21:40]),border=F,cex.names = 0.5,col=1:100,las=2)
windows(height=4);barplot(as.matrix(treatmentdatalist[['3 10 0.75']][rev(order(rowMeans(treatmentdatalist[['3 10 0.75']]))),21:40]),border=F,cex.names = 0.5,col=1:100,las=2)


windows();plot(hclust(as.dist(1-cor(log(treatmentdatalist[[1]])))))
windows();plot(hclust(as.dist(1-cor(log(treatmentdatalist[[25]])))))


windows();boxplot(log(t(treatmentdatalist[[25]])[,1])~paste(meta$treatment,meta$timepoint))
windows();boxplot(log(t(treatmentdatalist[[1]])[,1])~paste(meta$treatment,meta$timepoint))


windows(width=3,height = 4.5);
plot(log(t(treatmentdatalist[['3 0.1 0.1']])[meta$treatment=="tr",3])~meta$timepoint[meta$treatment=="tr"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Treatment group 10%",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines(log(t(treatmentdatalist[['3 0.1 0.1']])[meta$treatment=="tr",3][c(i,i+10)])~meta$timepoint[meta$treatment=="tr"][c(i,i+10)],
        col=i,lwd=2)
}

windows(width=3,height = 4.5);
plot(log(t(treatmentdatalist[['1 10 0.75']])[meta$treatment=="tr",1])~meta$timepoint[meta$treatment=="tr"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Treatment group 10%",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines(log(t(treatmentdatalist[['1 10 0.75']])[meta$treatment=="tr",1][c(i,i+10)])~meta$timepoint[meta$treatment=="tr"][c(i,i+10)],
        col=i,lwd=2)
}




windows(width=3,height = 4.5);
plot(log(t(treatmentdatalist[['0.5 0.05']])[meta$treatment=="tr",3])~meta$timepoint[meta$treatment=="tr"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Treatment group 10%",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines(log(t(treatmentdatalist[['0.5 0.05']])[meta$treatment=="tr",3][c(i,i+10)])~meta$timepoint[meta$treatment=="tr"][c(i,i+10)],
        col=i,lwd=2)
}



windows(width=3,height = 4.5);
plot(log(t(treatmentdatalist[['2 100 0.5']])[meta$treatment=="tr",4])~meta$timepoint[meta$treatment=="tr"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Treatment group",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines(log(t(treatmentdatalist[['2 100 0.5']])[meta$treatment=="tr",4][c(i,i+10)])~meta$timepoint[meta$treatment=="tr"][c(i,i+10)],
        col=i,lwd=2)
}


windows(width=3,height = 4.5);
plot(log(t(treatmentdatalist[['10 0.25']])[meta$treatment=="tr",1])~meta$timepoint[meta$treatment=="tr"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Treatment group 10x",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines(log(t(treatmentdatalist[['10 0.25']])[meta$treatment=="tr",1][c(i,i+10)])~meta$timepoint[meta$treatment=="tr"][c(i,i+10)],
        col=i,lwd=2)
}

windows();hist(log(t(treatmentdatalist[[1]])[,1]))
windows();hist(log(t(treatmentdatalist[[1]])[,100]))
windows();hist(log(t(treatmentdatalist[[25]])[,1]))

FC <- list()
for(k in c(1:100)){ 
for(i in c(0.1,0.5,1,5,10)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(treatmentdatalist[[paste(k,i,j)]]))
    tr1[tr1<1]<-1
    tr1$timepoint <- rep(c("1","2"),each=20)
    tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
    
FC[[paste(k,i,j)]] <- data.frame(c=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]),
                               tr=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100]),
                               bw=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100])/
                                 colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]))
 }
}
}

gt_bw <- list()

for(k in c(1:100)){ 
for(i in c(0.1,0.5,1,5,10)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
tr1 <- log(as.data.frame(t(treatmentdatalist[[paste(k,i,j)]])))
tr1$timepoint <- rep(c("1","2"),each=20)
tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
tr1[tr1<0]<-0
delta <- (exp(tr1[tr1$timepoint==2,1:100])-exp(tr1[tr1$timepoint==1,1:100]))/exp(tr1[tr1$timepoint==1,1:100])
delta$treatment <- tr1$treatment[1:20]
gt_bw[[paste(k,i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
for(h in 1:100){
gt_bw[[paste(k,i,j)]][gt_bw[[paste(k,i,j)]]$taxon==h,"p"]<-summary(lm(as.formula(paste("V",h,"~treatment",sep="")),data=delta))$coef[2,4]
}
}
}
}

#i=0.5
#j=0.1


peardist_log <- function(x,distance){
  as.dist(1-cor(t(log(x+min(x[x>0])))))
}

peardist <- function(x,distance){
  as.dist(1-cor(t((x))))
}



palette(scales::alpha(c("plum1","purple","skyblue","royalblue"),0.750))


  ad_abs<-list()
  ad_abs_bc<-list()
  for(k in 1:100){
ad_abs[[k]]<- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
ad_abs_bc[[k]]<- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
  }

    for(k in 1:100){
  for(i in c(0.1,0.5,1,5,10)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
tr1 <- as.data.frame(t(treatmentdatalist[[paste(k,i,j)]]))
ad_abs[[k]][ad_abs[[k]]$effsize==i&ad_abs[[k]]$Nchange==j,"R2"] <- adonis2(peardist_log(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$R2[1]
ad_abs_bc[[k]][ad_abs_bc[[k]]$effsize==i&ad_abs_bc[[k]]$Nchange==j,"R2"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "bray")$R2[1]
}}}

ad_abs_mean <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
abss <- cbind(ad_abs[[1]][,3],ad_abs[[2]][,3])
for(k in 3:50){
  abss <- cbind(abss,ad_abs[[k]][,3])
}
ad_abs_mean$R2 <- rowMeans(abss)


windows(width=4,height = 3);
ggplot(ad_abs_mean, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="firebrick",limits=c(-0.1,1)) 

ad_abs_bc_mean <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
abss <- cbind(ad_abs_bc[[1]][,3],ad_abs_bc[[2]][,3])
for(k in 3:30){
  abss <- cbind(abss,ad_abs_bc[[k]][,3])
}
ad_abs_bc_mean$R2 <- rowMeans(abss)

windows(width=4,height = 3);
ggplot(ad_abs_bc_mean, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="turquoise4",limits=c(-0.1,1)) 

result_abs <-list()
for(k in 1:100){
result_abs[[k]] <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
for(i in c(0.1,0.5,1,5,10)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    pvals <- as.data.frame(gt_bw[[paste(k,i,j)]])
    pvals$adj <- p.adjust(pvals[,2],"fdr")
    folds <- FC[[paste(k,i,j)]]
    absfolds <- abs(log(folds[,"bw"]))
    fccut <- quantile(abs(relfolds[pvals$adj<0.1]),0.1)
    result_abs[[k]]$true_pos[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j]<- j*100
    result_abs[[k]]$true_neg[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j]<- 100-(j*100)
    result_abs[[k]]$discovered_pos_fc[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[absfolds>fccut&pvals$adj<0.1,])
    result_abs[[k]]$discovered_neg_fc[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[absfolds<fccut|pvals$adj>0.1,])
    result_abs[[k]]$discovered_neg[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[pvals[,"adj"]>0.1,])
    result_abs[[k]]$discovered_pos[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[pvals[,"adj"]<0.1,])
    
    result_abs[[k]]$correct_pos[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.1,])
    result_abs[[k]]$false_pos[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.1,])
    result_abs[[k]]$false_neg[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.1,])
    result_abs[[k]]$correct_neg[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.1,])
    
    result_abs[[k]]$correct_pos_fc[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.1&absfolds[1:(100*j)]>fccut,])
    result_abs[[k]]$false_pos_fc[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.1&absfolds[(100*j+1):100]>fccut,])
    result_abs[[k]]$false_neg_fc[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.1|absfolds[1:(100*j)]<fccut,])
    result_abs[[k]]$correct_neg_fc[result_abs[[k]]$effsize==i&result_abs[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.1|absfolds[(100*j+1):100]<fccut,])
  }}

result_abs[[k]]$false_pos_prop <- round(result_abs[[k]]$false_pos/(result_abs[[k]]$false_pos+result_abs[[k]]$true_neg),digits=2)
result_abs[[k]]$false_neg_prop <- round(result_abs[[k]]$false_neg/(result_abs[[k]]$true_pos+result_abs[[k]]$false_neg),digits=2)

result_abs[[k]]$false_pos_prop_fc <- round(result_abs[[k]]$false_pos_fc/(result_abs[[k]]$false_pos_fc+result_abs[[k]]$true_neg),digits=2)
result_abs[[k]]$false_neg_prop_fc <- round(result_abs[[k]]$false_neg_fc/(result_abs[[k]]$true_pos+result_abs[[k]]$false_neg_fc),digits=2)
}

result_abs_mean <-data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
abss <- cbind(result_abs[[1]][,c("false_pos_prop","false_neg_prop","false_pos_prop_fc","false_neg_prop_fc")],
              result_abs[[2]][,c("false_pos_prop","false_neg_prop","false_pos_prop_fc","false_neg_prop_fc")])
for(k in 3:100){
  abss <- cbind(abss,result_abs[[k]][,c("false_pos_prop","false_neg_prop","false_pos_prop_fc","false_neg_prop_fc")])
}
result_abs_mean$false_pos_prop <- rowMeans(abss[,grep(names(abss),pattern="pos_prop$",value=T)])
result_abs_mean$false_pos_prop_fc <- rowMeans(abss[,grep(names(abss),pattern="pos_prop_fc$",value=T)])
result_abs_mean$false_neg_prop <- rowMeans(abss[,grep(names(abss),pattern="neg_prop$",value=T)])
result_abs_mean$false_neg_prop_fc <- rowMeans(abss[,grep(names(abss),pattern="neg_prop_fc$",value=T)])


windows(width=4,height = 3);
ggplot(result_abs_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False pos"))+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4,height =3);
ggplot(result_abs_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False neg"))+
  scale_fill_gradient(low="white", high="royalblue",limits=c(0,1)) 


windows(width=4,height = 3);
ggplot(result_abs_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop_fc)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False pos"))+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4,height =3);
ggplot(result_abs_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop_fc)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False neg"))+
  scale_fill_gradient(low="white", high="royalblue",limits=c(0,1)) 


#----------------------------------------------------------------------------------

#relative data--------------------------------------------------------------



palette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F",
          brewer.pal(9,"Greens")[-1],
          brewer.pal(9,"Purples")[-1],
          brewer.pal(9,"YlGnBu")[-1],rep(c("gray90","gray40"),80)))

seqdatalist <- treatmentdatalist
for(k in 1:100){
for(i in c(0.1,0.5,1,5,10)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
    seqdata <- t(t(treatmentdatalist[[paste(k,i,j)]])/rowSums(t(treatmentdatalist[[paste(k,i,j)]])))
    for(h in c(1:40)) {
      seqdatalist[[paste(k,i,j)]][,h]<-round(seqdata[,h]*rnorm(n=1,mean=50000,sd=10000))
    }
    seqdatalist[[paste(k,i,j)]][seqdatalist[[paste(k,i,j)]]<0]<-0
    }}}

windows();plot((log(t(seqdatalist[[paste(k,i,j)]][1,])+1))~t(log((treatmentdatalist[[paste(k,i,j)]][1,]))), 
              ylim=c(0,12),xlim=c(8,22))
for(h in c(2:100)){
  points((log(t(seqdatalist[[paste(k,i,j)]][h,])+1))~t(log((treatmentdatalist[[paste(k,i,j)]][h,]))))
}


windows(width=7,height = 4);
par(mgp=c(3.5,1,0),mar=c(4,5,2,2),cex.axis=0.7)
barplot(as.matrix(treatmentdatalist[[paste(k,i,j)]][rev(order(rowMeans(treatmentdatalist[[paste(k,i,j)]]))),]),
                                   border=F,cex.names = 0.5,col=1:100,las=2,ylab="Genomes")

windows(width=7,height = 4);
par(mgp=c(3.5,1,0),mar=c(4,5,2,2))
barplot(as.matrix(seqdatalist[[paste(k,i,j)]][rev(order(rowMeans(treatmentdatalist[[paste(k,i,j)]]))),]),border=F,
        cex.names = 0.5,col=1:100,las=2,ylab="Reads")


windows(width=5,height = 4);barplot(as.matrix(seqdata[rev(order(rowMeans(seqdata))),]),border=F,cex.names = 0.5,col=1:100,las=2)

reldatalist <- list()
for(k in c(1:100)){ 
for(i in c(0.1,0.5,1,5,10)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
    reldatalist[[paste(k,i,j)]] <- t(t(seqdatalist[[paste(k,i,j)]])/rowSums(t(seqdatalist[[paste(k,i,j)]])))
    }}
}



windows(width=4,height = 4);plot(log(t(reldatalist[[paste(k,i,j)]])[,1]+0.0001)~log(t(treatmentdatalist[[paste(k,i,j)]])[,1]),
              ylim=c(-10,0),xlim=c(4,22),col=alpha("red",0.5),pch=20,
              ylab="log(relative)",xlab="log(absolute)")
for(h in c(2:100)){
  points(log(t(reldatalist[[paste(k,i,j)]])[,h]+0.0001)~log(t(treatmentdatalist[[paste(k,i,j)]])[,h]),
         col=alpha("red",0.5),pch=20)
}

windows(width=4,height = 4);plot((((reldatalist[[paste(k,i,j)]][1,]+0.0001)))~t(((treatmentdatalist[[paste(k,i,j)]][1,]))),
             ylim=c(0,0.8),xlim=c(0,exp(21)),col=alpha("red",0.5),pch=20,
             ylab="Relative)",xlab="Absolute")
for(h in c(2:100)){
  points((((reldatalist[[paste(k,i,j)]][h,]+0.0001)))~t(((treatmentdatalist[[paste(k,i,j)]][h,]))),col=alpha("red",0.5),pch=20)
}

cor_rel_mean <- 0.5
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
cors <- cbind(cor(t(reldatalist[[paste(k,i,j)]])[,1],t(treatmentdatalist[[paste(k,i,j)]])[,1]),
              cor(t(reldatalist[[paste(k,i,j)]])[,2],t(treatmentdatalist[[paste(k,i,j)]])[,2]))
for(k in 3:100){
  cors <- cbind(cors,cor(t(reldatalist[[paste(k,i,j)]])[,h],t(treatmentdatalist[[paste(k,i,j)]])[,h]))
}
cor_rel_mean <- cbind(cor_rel_mean,mean(cors))   
}}}

cor_rel_mean_log <- 0.5
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
      cors_log <- cbind(cor(log(t(reldatalist[[paste(k,i,j)]])[,1]+0.0001),log(t(treatmentdatalist[[paste(k,i,j)]])[,1])),
                    cor(t(reldatalist[[paste(k,i,j)]])[,2],t(treatmentdatalist[[paste(k,i,j)]])[,2]))
      for(k in 3:100){
        cors_log <- cbind(cors_log,cor(log(t(reldatalist[[paste(k,i,j)]])[,h]+0.0001),log(t(treatmentdatalist[[paste(k,i,j)]])[,h])))
      }
      cor_rel_mean_log <- cbind(cor_rel_mean_log,mean(cors_log))      
    }}}



windows();plot(log(reldatalist[[paste(k,i,j)]][,1])~log(treatmentdatalist[[paste(k,i,j)]][,1]))
for(h in c(2:100)){
  points(log(reldatalist[[paste(k,i,j)]][,h])~log(treatmentdatalist[[paste(k,i,j)]][,h]))
}

windows(width=7,height = 4);
par(mgp=c(3.5,1,0),mar=c(4,5,2,2))
barplot(as.matrix(reldatalist[[paste(k,i,j)]][rev(order(rowMeans(treatmentdatalist[[paste(k,i,j)]]))),]),border=F,
        cex.names = 0.5,col=1:100,las=2,ylab="Reads")

windows(height=4);
barplot(as.matrix(reldatalist[[1]][rev(order(rowMeans(reldatalist[[1]]))),])[,c(1,21)],border=F,cex.names = 0.5,col=1:100,las=2)

FC_rel <- list()
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
      tr1 <- as.data.frame(t(reldatalist[[paste(k,i,j)]]))
      tr1 <- tr1+0.0001
      tr1$timepoint <- rep(c("1","2"),each=20)
      tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
      
      FC_rel[[paste(k,i,j)]] <- data.frame(c=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]),
                                       tr=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100]),
                                       bw=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100])/
                                         colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]))
    }
  }
}

gt_bw_rel <- list()
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      tr1 <- (as.data.frame(t(reldatalist[[paste(k,i,j)]]+0.00001)))
      tr1$timepoint <- rep(c("1","2"),each=20)
      tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
      delta <- ((tr1[tr1$timepoint==2,1:100])-(tr1[tr1$timepoint==1,1:100]))/(tr1[tr1$timepoint==1,1:100])
      delta$treatment <- tr1$treatment[1:20]
     gt_bw_rel[[paste(k,i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
      for(h in 1:100){
      gt_bw_rel[[paste(k,i,j)]][gt_bw_rel[[paste(k,i,j)]]$taxon==h,"p"]<-kruskal.test(g=delta[,"treatment"],x=delta[,paste("V",h,sep="")])$p.value
      }
    }
  }
}


ad_rel<-list()
ad_rel_bc<-list()
for(k in 1:100){
  ad_rel[[k]]<- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
  ad_rel_bc[[k]]<- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
}

for(k in 1:100){
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      tr1 <- as.data.frame(t(reldatalist[[paste(k,i,j)]]))
      ad_rel[[k]][ad_rel[[k]]$effsize==i&ad_rel[[k]]$Nchange==j,"R2"] <- adonis2(peardist_log(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$R2[1]
      ad_rel_bc[[k]][ad_rel_bc[[k]]$effsize==i&ad_rel_bc[[k]]$Nchange==j,"R2"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "bray")$R2[1]
    }}}

ad_rel_mean <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
rels <- cbind(ad_rel[[1]][,3],ad_rel[[2]][,3])
for(k in 3:50){
  rels <- cbind(rels,ad_rel[[k]][,3])
}
ad_rel_mean$R2 <- rowMeans(rels)

ad_rel_bc_mean <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
rels <- cbind(ad_rel_bc[[1]][,3],ad_rel_bc[[2]][,3])
for(k in 3:30){
  rels <- cbind(rels,ad_rel_bc[[k]][,3])
}
ad_rel_bc_mean$R2 <- rowMeans(rels)

windows(width=4,height = 3);
ggplot(ad_rel_mean, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="firebrick",limits=c(-0.1,1)) 

windows(width=4,height = 3);
ggplot(ad_rel_bc_mean, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="turquoise4",limits=c(-0.1,1)) 

result_rel <-list()
for(k in 1:100){
  result_rel[[k]] <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      pvals <- as.data.frame(gt_bw_rel[[paste(k,i,j)]])
      pvals$p[is.na(pvals$p)]<-1
      pvals$adj <- p.adjust(pvals[,2],"fdr")
      folds <- FC_rel[[paste(k,i,j)]]
      relfolds <- rel(log(folds[,"bw"]))
      fccut <- quantile(abs(relfolds[pvals$adj<0.1]),0.1)
      result_rel[[k]]$true_pos[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j]<- j*100
      result_rel[[k]]$true_neg[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j]<- 100-(j*100)
      result_rel[[k]]$discovered_pos_fc[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[relfolds>fccut&pvals$adj<0.1,])
      result_rel[[k]]$discovered_neg_fc[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[relfolds<fccut|pvals$adj>0.1,])
      result_rel[[k]]$discovered_neg[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[pvals[,"adj"]>0.1,])
      result_rel[[k]]$discovered_pos[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[pvals[,"adj"]<0.1,])
      
      result_rel[[k]]$correct_pos[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.1,])
      result_rel[[k]]$false_pos[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.1,])
      result_rel[[k]]$false_neg[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.1,])
      result_rel[[k]]$correct_neg[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.1,])
      
      result_rel[[k]]$correct_pos_fc[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.1&relfolds[1:(100*j)]>fccut,])
      result_rel[[k]]$false_pos_fc[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.1&relfolds[(100*j+1):100]>fccut,])
      result_rel[[k]]$false_neg_fc[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.1|relfolds[1:(100*j)]<fccut,])
      result_rel[[k]]$correct_neg_fc[result_rel[[k]]$effsize==i&result_rel[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.1|relfolds[(100*j+1):100]<fccut,])
    }}
  
  result_rel[[k]]$false_pos_prop <- round(result_rel[[k]]$false_pos/(result_rel[[k]]$false_pos+result_rel[[k]]$true_neg),digits=2)
  result_rel[[k]]$false_neg_prop <- round(result_rel[[k]]$false_neg/(result_rel[[k]]$true_pos+result_rel[[k]]$false_neg),digits=2)
  
  result_rel[[k]]$false_pos_prop_fc <- round(result_rel[[k]]$false_pos_fc/(result_rel[[k]]$false_pos_fc+result_rel[[k]]$true_neg),digits=2)
  result_rel[[k]]$false_neg_prop_fc <- round(result_rel[[k]]$false_neg_fc/(result_rel[[k]]$true_pos+result_rel[[k]]$false_neg_fc),digits=2)
}

result_rel_mean <-data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
rels <- cbind(result_rel[[1]][,c("false_pos_prop","false_neg_prop","false_pos_prop_fc","false_neg_prop_fc")],
              result_rel[[2]][,c("false_pos_prop","false_neg_prop","false_pos_prop_fc","false_neg_prop_fc")])
for(k in 3:100){
  rels <- cbind(rels,result_rel[[k]][,c("false_pos_prop","false_neg_prop","false_pos_prop_fc","false_neg_prop_fc")])
}
result_rel_mean$false_pos_prop <- rowMeans(rels[,grep(names(rels),pattern="pos_prop$",value=T)])
result_rel_mean$false_pos_prop_fc <- rowMeans(rels[,grep(names(rels),pattern="pos_prop_fc$",value=T)])
result_rel_mean$false_neg_prop <- rowMeans(rels[,grep(names(rels),pattern="neg_prop$",value=T)])
result_rel_mean$false_neg_prop_fc <- rowMeans(rels[,grep(names(rels),pattern="neg_prop_fc$",value=T)])


windows(width=4,height = 3);
ggplot(result_rel_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False pos"))+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4,height = 3);
ggplot(result_rel_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop_fc)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False pos"))+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4,height =3);
ggplot(result_rel_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False neg"))+
  scale_fill_gradient(low="white", high="royalblue",limits=c(0,1)) 

windows(width=4,height =3);
ggplot(result_rel_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop_fc)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_neg_prop_fc), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False neg"))+
  scale_fill_gradient(low="white", high="royalblue",limits=c(0,1)) 

fc_comp <- list()
for(k in 1:100){
  fc_comp[[k]] <- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=500),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),500))
for(i in c(0.5,1,5,10)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    fc_comp[[k]]$FC[fc_comp[[k]]$effsize==i&fc_comp[[k]]$Nchange==j] <- FC[[paste(k,i,j)]]$bw
    fc_comp[[k]]$FC_rel[fc_comp[[k]]$effsize==i&fc_comp[[k]]$Nchange==j]<- FC_rel[[paste(k,i,j)]]$bw
  }}}


cols = brewer.pal(9, "GnBu")[-c(1,2)]
pal = colorRampPalette(cols)

ord = findInterval(fc_comp[[k]]$effsize, sort(fc_comp[[k]]$effsize))
windows();plot(log(FC_rel)~log(FC),data=fc_comp[[k]],pch=20,col=scales::alpha(pal(nrow(fc_comp[[k]]))[ord],0.75))
legend("bottomright",legend=unique(fc_comp[[k]]$effsize),pch = 21, pt.bg = scales::alpha(pal(6),0.75),cex=2,col="white")

ord = findInterval(fc_comp[[k]]$Nchange, sort(fc_comp[[k]]$Nchange))
windows();plot(log(FC_rel)~log(FC),data=fc_comp[[k]],pch=20,col=scales::alpha(pal(nrow(fc_comp[[k]]))[ord],0.75))
legend("bottomright",legend=unique(fc_comp[[k]]$Nchange),pch = 21, pt.bg = scales::alpha(pal(6),0.75),cex=2,col="white")



#clr data--------------------------------------------------------------

clrdatalist <- list()
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
      clrdatalist[[paste(k,i,j)]] <- t(clrt(seqdatalist[[paste(k,i,j)]]))
  }}
}

windows(width=4,height = 4);plot((t(clrdatalist[[paste(k,i,j)]])[,1]+0.0001)~log(t(treatmentdatalist[[paste(k,i,j)]])[,1]),
                                ylim=c(-5,7),xlim=c(4,22),col=alpha("red",0.5),pch=20,
                                ylab="CLR",xlab="log(absolute)")
for(h in c(2:100)){
  points((t(clrdatalist[[paste(k,i,j)]])[,h]+0.0001)~log(t(treatmentdatalist[[paste(k,i,j)]])[,h]),
         col=alpha("red",0.5),pch=20)
}



FC_clr <- list()
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
      tr1 <- as.data.frame(exp(t(clrdatalist[[paste(k,i,j)]])))
      tr1$timepoint <- rep(c("1","2"),each=20)
      tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
      FC_clr[[paste(k,i,j)]] <- data.frame(c=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]),
                                           tr=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100]),
                                           bw=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100])/
                                             colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]))
    }
  }
}


gt_bw_clr <- list()
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      tr1 <- as.data.frame(t(clrdatalist[[paste(k,i,j)]]))
      tr1$timepoint <- rep(c("1","2"),each=20)
      tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
      delta <- (exp(tr1[tr1$timepoint==2,1:100])-exp(tr1[tr1$timepoint==1,1:100]))/exp(tr1[tr1$timepoint==1,1:100])
      delta$treatment <- tr1$treatment[1:20]
   gt_bw_clr[[paste(k,i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
      for(h in 1:100){
       gt_bw_clr[[paste(k,i,j)]][gt_bw_clr[[paste(k,i,j)]]$taxon==h,"p"]<-kruskal.test(g=delta[,"treatment"],x=delta[,paste("V",h,sep="")])$p.value
      }
      
    }
  }
}


ad_clr<-list()
ad_clr_bc<-list()
for(k in 1:100){
  ad_clr[[k]]<- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
  ad_clr_bc[[k]]<- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
}

for(k in 51:100){
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      tr1 <- as.data.frame(t(clrdatalist[[paste(k,i,j)]]))
      ad_clr[[k]][ad_clr[[k]]$effsize==i&ad_clr[[k]]$Nchange==j,"R2"] <- adonis2(peardist(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$R2[1]
      ad_clr_bc[[k]][ad_clr_bc[[k]]$effsize==i&ad_clr_bc[[k]]$Nchange==j,"R2"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "euc")$R2[1]
    }}}

ad_clr_mean <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
clrs <- cbind(ad_clr[[1]][,3],ad_clr[[2]][,3])
for(k in 3:50){
  clrs <- cbind(clrs,ad_clr[[k]][,3])
}
ad_clr_mean$R2 <- rowMeans(clrs,na.rm=T)

ad_clr_bc_mean <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
clrs <- cbind(ad_clr_bc[[1]][,3],ad_clr_bc[[2]][,3])
for(k in 3:30){
  clrs <- cbind(clrs,ad_clr_bc[[k]][,3])
}
ad_clr_bc_mean$R2 <- rowMeans(clrs,na.rm=T)

windows(width=4,height = 3);
ggplot(ad_clr_mean, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="firebrick",limits=c(-0.1,1)) 

windows(width=4,height = 3);
ggplot(ad_clr_bc_mean, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="orange3",limits=c(-0.1,1)) 


result_clr <-list()
for(k in 1:100){
  result_clr[[k]] <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      pvals <- as.data.frame(gt_bw_clr[[paste(k,i,j)]])
      pvals$adj <- p.adjust(pvals[,2],"fdr")
      result_clr[[k]]$true_pos[result_clr[[k]]$effsize==i&result_clr[[k]]$Nchange==j]<- j*100
      result_clr[[k]]$true_neg[result_clr[[k]]$effsize==i&result_clr[[k]]$Nchange==j]<- 100-(j*100)
      result_clr[[k]]$discovered_pos[result_clr[[k]]$effsize==i&result_clr[[k]]$Nchange==j] <- nrow(pvals[pvals[,"adj"]<0.2,])
      result_clr[[k]]$discovered_neg[result_clr[[k]]$effsize==i&result_clr[[k]]$Nchange==j] <- nrow(pvals[pvals[,"adj"]>0.2,])
      result_clr[[k]]$correct_pos[result_clr[[k]]$effsize==i&result_clr[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.2,])
      result_clr[[k]]$false_pos[result_clr[[k]]$effsize==i&result_clr[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.2,])
      result_clr[[k]]$false_neg[result_clr[[k]]$effsize==i&result_clr[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.2,])
      result_clr[[k]]$correct_neg[result_clr[[k]]$effsize==i&result_clr[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.2,])
    }}
  result_clr[[k]]$false_pos_prop <- round(result_clr[[k]]$false_pos/(result_clr[[k]]$false_pos+result_clr[[k]]$true_neg),digits=2)
  result_clr[[k]]$false_neg_prop <- round(result_clr[[k]]$false_neg/(result_clr[[k]]$true_pos+result_clr[[k]]$false_neg),digits=2)
}

result_clr_mean <-data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
clrs <- cbind(result_clr[[1]][,c(11,12)],result_clr[[2]][,c(11,12)])
for(k in 3:50){
  clrs <- cbind(clrs,result_clr[[k]][,c(11,12)])
}
result_clr_mean$false_pos_prop <- rowMeans(clrs[,grep(names(clrs),pattern="pos",value=T)])
result_clr_mean$false_neg_prop <- rowMeans(clrs[,grep(names(clrs),pattern="neg",value=T)])


windows(width=4,height = 3);
ggplot(result_clr_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False pos"))+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4,height =3);
ggplot(result_clr_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False neg"))+
  scale_fill_gradient(low="white", high="royalblue",limits=c(0,1)) 


#deseq data--------------------------------------------------------------

deseqdatalist <- list()
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
      dds <- DESeqDataSetFromMatrix(countData = seqdatalist[[paste(k,i,j)]], 
                                    colData =(meta),
                                    design = ~ timepoint + treatment + timepoint:treatment)
      dds <- estimateSizeFactors(dds)
      deseqdatalist[[paste(k,i,j)]]  <- counts(dds, normalized=TRUE)
    }}
}

windows(width=4,height = 4);plot(log(t(deseqdatalist[[paste(k,i,j)]])[,1]+0.0001)~log(t(treatmentdatalist[[paste(k,i,j)]])[,1]),
                                ylim=c(-1,11),xlim=c(4,22),
                                col=alpha("red",0.5),pch=20,
                                ylab="DESeq",xlab="log(absolute)")
for(h in c(2:100)){
  points(log(t(deseqdatalist[[paste(k,i,j)]])[,h]+0.0001)~log(t(treatmentdatalist[[paste(k,i,j)]])[,h]),
         col=alpha("red",0.5),pch=20)
}


windows(width=4,height = 4);plot((log((deseqdatalist[[paste(k,i,j)]][1,]+0.0001)))~t(log((treatmentdatalist[[paste(k,i,j)]][1,]))),
                                ylim=c(-10,0),xlim=c(4,22),col=alpha("red",0.5),pch=20,
                                ylab="log(deseqative)",xlab="log(absolute)")
for(h in c(2:100)){
  points((log((deseqdatalist[[paste(k,i,j)]][h,]+0.0001)))~t(log((treatmentdatalist[[paste(k,i,j)]][h,]))),
         col=alpha("red",0.5),pch=20)
}

FC_deseq <- list()
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
      tr1 <- as.data.frame(t(deseqdatalist[[paste(k,i,j)]]))
      tr1[tr1<1]<-1
      tr1$timepoint <- rep(c("1","2"),each=20)
      tr1$deseq <- rep(c("c","tr","c","tr"),each=10)
      
      FC_deseq[[paste(k,i,j)]] <- data.frame(c=colMeans(tr1[tr1$timepoint==2&tr1$deseq=="c",1:100]/tr1[tr1$timepoint==1&tr1$deseq=="c",1:100]),
                                           tr=colMeans(tr1[tr1$timepoint==2&tr1$deseq=="tr",1:100]/tr1[tr1$timepoint==1&tr1$deseq=="tr",1:100]),
                                           bw=colMeans(tr1[tr1$timepoint==2&tr1$deseq=="tr",1:100]/tr1[tr1$timepoint==1&tr1$deseq=="tr",1:100])/
                                             colMeans(tr1[tr1$timepoint==2&tr1$deseq=="c",1:100]/tr1[tr1$timepoint==1&tr1$deseq=="c",1:100]))
    }
  }
}

gt_bw_deseq <- list()

for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      tr1 <- as.data.frame(t(deseqdatalist[[paste(k,i,j)]])+0.0001)
      tr1$timepoint <- rep(c("1","2"),each=20)
      tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
      delta <- (tr1[tr1$timepoint==2,1:100]-tr1[tr1$timepoint==1,1:100])/(tr1[tr1$timepoint==1,1:100])
      delta$treatment <- tr1$treatment[1:20]
      gt_bw_deseq[[paste(k,i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
      
      for(h in 1:100){
        gt_bw_deseq[[paste(k,i,j)]][gt_bw_deseq[[paste(k,i,j)]]$taxon==h,"p"]<-kruskal.test(g=delta[,"treatment"],x=delta[,paste("V",h,sep="")])$p.value
      }
    }
  }
}


ad_deseq<-list()
ad_deseq_bc<-list()
for(k in 1:100){
  ad_deseq[[k]]<- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
  ad_deseq_bc[[k]]<- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
}

for(k in 1:100){
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      tr1 <- as.data.frame(t(deseqdatalist[[paste(k,i,j)]]))
      ad_deseq[[k]][ad_deseq[[k]]$effsize==i&ad_deseq[[k]]$Nchange==j,"R2"] <- adonis2(peardist_log(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$R2[1]
      ad_deseq_bc[[k]][ad_deseq_bc[[k]]$effsize==i&ad_deseq_bc[[k]]$Nchange==j,"R2"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "bray")$R2[1]
    }}}

ad_deseq_mean <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
deseqs <- cbind(ad_deseq[[1]][,3],ad_deseq[[2]][,3])
for(k in 3:100){
  deseqs <- cbind(deseqs,ad_deseq[[k]][,3])
}
ad_deseq_mean$R2 <- rowMeans(deseqs)

ad_deseq_bc_mean <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
deseqs <- cbind(ad_deseq_bc[[1]][,3],ad_deseq_bc[[2]][,3])
for(k in 3:100){
  deseqs <- cbind(deseqs,ad_deseq_bc[[k]][,3])
}
ad_deseq_bc_mean$R2 <- rowMeans(deseqs)

windows(width=4,height = 3);
ggplot(ad_deseq_mean, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="firebrick",limits=c(-0.1,1)) 

windows(width=4,height = 3);
ggplot(ad_deseq_bc_mean, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="turquoise4",limits=c(-0.1,1)) 

result_deseq <-list()
for(k in 1:100){
  result_deseq[[k]] <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      pvals <- as.data.frame(gt_bw_deseq[[paste(k,i,j)]])
      pvals$adj <- p.adjust(pvals[,2],"fdr")
      result_deseq[[k]]$true_pos[result_deseq[[k]]$effsize==i&result_deseq[[k]]$Nchange==j]<- j*100
      result_deseq[[k]]$true_neg[result_deseq[[k]]$effsize==i&result_deseq[[k]]$Nchange==j]<- 100-(j*100)
      result_deseq[[k]]$discovered_pos[result_deseq[[k]]$effsize==i&result_deseq[[k]]$Nchange==j] <- nrow(pvals[pvals[,"adj"]<0.1,])
      result_deseq[[k]]$discovered_neg[result_deseq[[k]]$effsize==i&result_deseq[[k]]$Nchange==j] <- nrow(pvals[pvals[,"adj"]>0.1,])
      result_deseq[[k]]$correct_pos[result_deseq[[k]]$effsize==i&result_deseq[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.1,])
      result_deseq[[k]]$false_pos[result_deseq[[k]]$effsize==i&result_deseq[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.1,])
      result_deseq[[k]]$false_neg[result_deseq[[k]]$effsize==i&result_deseq[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.1,])
      result_deseq[[k]]$correct_neg[result_deseq[[k]]$effsize==i&result_deseq[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.1,])
    }}
  result_deseq[[k]]$false_pos_prop <- round(result_deseq[[k]]$false_pos/(result_deseq[[k]]$false_pos+result_deseq[[k]]$true_neg),digits=2)
  result_deseq[[k]]$false_neg_prop <- round(result_deseq[[k]]$false_neg/(result_deseq[[k]]$true_pos+result_deseq[[k]]$false_neg),digits=2)
}

result_deseq_mean <-data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
deseqs <- cbind(result_deseq[[1]][,c(11,12)],result_deseq[[2]][,c(11,12)])
for(k in 3:100){
  deseqs <- cbind(deseqs,result_deseq[[k]][,c(11,12)])
}
result_deseq_mean$false_pos_prop <- rowMeans(deseqs[,grep(names(deseqs),pattern="pos",value=T)])
result_deseq_mean$false_neg_prop <- rowMeans(deseqs[,grep(names(deseqs),pattern="neg",value=T)])


windows(width=4,height = 3);
ggplot(result_deseq_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False pos"))+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4,height =3);
ggplot(result_deseq_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False neg"))+
  scale_fill_gradient(low="white", high="royalblue",limits=c(0,1)) 


#edger data--------------------------------------------------------------

edgerdatalist <- list()
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
      y <- DGEList(counts=round(seqdatalist[[paste(k,i,j)]]),group=meta$treatment)
      edgerdatalist[[paste(k,i,j)]]   <- seqdatalist[[paste(k,i,j)]]*normLibSizes(y)$samples$norm
    }}
}


windows(width=4,height = 4);plot(log(t(edgerdatalist[[paste(k,i,j)]])[,1]+0.0001)~log(t(treatmentdatalist[[paste(k,i,j)]])[,1]),
                                ylim=c(-1,11),xlim=c(4,22),
                                col=alpha("red",0.5),pch=20,
                                ylab="DESeq",xlab="log(absolute)")
for(h in c(2:100)){
  points(log(t(edgerdatalist[[paste(k,i,j)]])[,h]+0.0001)~log(t(treatmentdatalist[[paste(k,i,j)]])[,h]),
         col=alpha("red",0.5),pch=20)
}

FC_edger <- list()
for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
      tr1 <- as.data.frame(t(edgerdatalist[[paste(k,i,j)]]))
      tr1[tr1<1]<-1
      tr1$timepoint <- rep(c("1","2"),each=20)
      tr1$edger <- rep(c("c","tr","c","tr"),each=10)
      
      FC_edger[[paste(k,i,j)]] <- data.frame(c=colMeans(tr1[tr1$timepoint==2&tr1$edger=="c",1:100]/tr1[tr1$timepoint==1&tr1$edger=="c",1:100]),
                                             tr=colMeans(tr1[tr1$timepoint==2&tr1$edger=="tr",1:100]/tr1[tr1$timepoint==1&tr1$edger=="tr",1:100]),
                                             bw=colMeans(tr1[tr1$timepoint==2&tr1$edger=="tr",1:100]/tr1[tr1$timepoint==1&tr1$edger=="tr",1:100])/
                                               colMeans(tr1[tr1$timepoint==2&tr1$edger=="c",1:100]/tr1[tr1$timepoint==1&tr1$edger=="c",1:100]))
    }
  }
}

gt_bw_edger <- list()

for(k in c(1:100)){ 
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      tr1 <- as.data.frame(t(edgerdatalist[[paste(k,i,j)]])+0.0001)
      tr1$timepoint <- rep(c("1","2"),each=20)
      tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
      delta <- (tr1[tr1$timepoint==2,1:100]-tr1[tr1$timepoint==1,1:100])/(tr1[tr1$timepoint==1,1:100])
      delta$treatment <- tr1$treatment[1:20]
      gt_bw_edger[[paste(k,i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
      for(h in 1:100){
        gt_bw_edger[[paste(k,i,j)]][gt_bw_edger[[paste(k,i,j)]]$taxon==h,"p"]<-kruskal.test(g=delta[,"treatment"],x=delta[,paste("V",h,sep="")])$p.value
      }
    }
  }
}


ad_edger<-list()
ad_edger_bc<-list()
for(k in 1:100){
  ad_edger[[k]]<- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
  ad_edger_bc[[k]]<- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
}

for(k in 1:100){
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      tr1 <- as.data.frame(t(edgerdatalist[[paste(k,i,j)]]))
      ad_edger[[k]][ad_edger[[k]]$effsize==i&ad_edger[[k]]$Nchange==j,"R2"] <- adonis2(peardist_log(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$R2[1]
      ad_edger_bc[[k]][ad_edger_bc[[k]]$effsize==i&ad_edger_bc[[k]]$Nchange==j,"R2"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "bray")$R2[1]
    }}}

ad_edger_mean <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
edgers <- cbind(ad_edger[[1]][,3],ad_edger[[2]][,3])
for(k in 3:100){
  edgers <- cbind(edgers,ad_edger[[k]][,3])
}
ad_edger_mean$R2 <- rowMeans(edgers)

ad_edger_bc_mean <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
edgers <- cbind(ad_edger_bc[[1]][,3],ad_edger_bc[[2]][,3])
for(k in 3:100){
  edgers <- cbind(edgers,ad_edger_bc[[k]][,3])
}
ad_edger_bc_mean$R2 <- rowMeans(edgers)

windows(width=4,height = 3);
ggplot(ad_edger_mean, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="firebrick",limits=c(-0.1,1)) 

windows(width=4,height = 3);
ggplot(ad_edger_bc_mean, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="turquoise4",limits=c(-0.1,1)) 

result_edger <-list()
for(k in 1:100){
  result_edger[[k]] <- data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
  for(i in c(0.1,0.5,1,5,10)){
    for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
      pvals <- as.data.frame(gt_bw_edger[[paste(k,i,j)]])
      pvals$adj <- p.adjust(pvals[,2],"fdr")
      result_edger[[k]]$true_pos[result_edger[[k]]$effsize==i&result_edger[[k]]$Nchange==j]<- j*100
      result_edger[[k]]$true_neg[result_edger[[k]]$effsize==i&result_edger[[k]]$Nchange==j]<- 100-(j*100)
      result_edger[[k]]$discovered_pos[result_edger[[k]]$effsize==i&result_edger[[k]]$Nchange==j] <- nrow(pvals[pvals[,"adj"]<0.2,])
      result_edger[[k]]$discovered_neg[result_edger[[k]]$effsize==i&result_edger[[k]]$Nchange==j] <- nrow(pvals[pvals[,"adj"]>0.2,])
      result_edger[[k]]$correct_pos[result_edger[[k]]$effsize==i&result_edger[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.2,])
      result_edger[[k]]$false_pos[result_edger[[k]]$effsize==i&result_edger[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.2,])
      result_edger[[k]]$false_neg[result_edger[[k]]$effsize==i&result_edger[[k]]$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.2,])
      result_edger[[k]]$correct_neg[result_edger[[k]]$effsize==i&result_edger[[k]]$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.2,])
    }}
  result_edger[[k]]$false_pos_prop <- round(result_edger[[k]]$false_pos/(result_edger[[k]]$false_pos+result_edger[[k]]$true_neg),digits=2)
  result_edger[[k]]$false_neg_prop <- round(result_edger[[k]]$false_neg/(result_edger[[k]]$true_pos+result_edger[[k]]$false_neg),digits=2)
}

result_edger_mean <-data.frame(effsize=rep(c(0.1,0.5,1,5,10),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),5))
edgers <- cbind(result_edger[[1]][,c(11,12)],result_edger[[2]][,c(11,12)])
for(k in 3:100){
  edgers <- cbind(edgers,result_edger[[k]][,c(11,12)])
}
result_edger_mean$false_pos_prop <- rowMeans(edgers[,grep(names(edgers),pattern="pos",value=T)])
result_edger_mean$false_neg_prop <- rowMeans(edgers[,grep(names(edgers),pattern="neg",value=T)])


windows(width=4,height = 3);
ggplot(result_edger_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False pos"))+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4,height =3);
ggplot(result_edger_mean, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False neg"))+
  scale_fill_gradient(low="white", high="royalblue",limits=c(0,1)) 




windows(width=4,height = 4);
plot(log(FC_rel)~log(FC),data=fc_comp[[k]],pch=20,col=scales::alpha(pal(nrow(fc_comp[[k]]))[ord],0.75))
legend("topleft",legend=unique(fc_comp[[k]]$Nchange),pch = 21, pt.bg = scales::alpha(pal(6),0.75),cex=1,col="white",bty="n")


#Difference plots for the R2 and different distances----
ad_rel_diff<-ad_rel_mean
ad_rel_diff$R2<-ad_rel_bc_mean$R2-ad_rel_mean$R2

d1<-ggplot(ad_rel_diff, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  # geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  # scale_fill_gradient(low="white", high="turquoise4",limits=c(-0.1,.1)) 
  geom_text(aes(label = round(R2,digits = 2)), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "FP diff")) +
  scale_fill_distiller(palette = 'RdBu',limits=c(-0.5,0.5))+
  theme(legend.position="none")  


ad_clr_mean_diff<-ad_clr_mean
ad_clr_mean_diff$R2<-ad_clr_bc_mean$R2-ad_clr_mean$R2

d2<-ggplot(ad_clr_mean_diff, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  # geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  # scale_fill_gradient(low="white", high="turquoise4",limits=c(-0.1,.1)) 
  geom_text(aes(label = round(R2,digits = 2)), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "FP diff")) +
  scale_fill_distiller(palette = 'RdBu',limits=c(-0.5,0.5))+
  theme(legend.position="none")  

ad_deseq_mean_diff<-ad_deseq_mean
ad_deseq_mean_diff$R2<-ad_deseq_bc_mean$R2-ad_deseq_mean$R2
#FIG experiment 1 RElAbundance R2 Bray----
#windows(width=4,height = 3);
d3<-ggplot(ad_deseq_mean_diff, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  # geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  # scale_fill_gradient(low="white", high="turquoise4",limits=c(-0.1,.1)) 
  geom_text(aes(label = round(R2,digits = 2)), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "FP diff")) +
  scale_fill_distiller(palette = 'RdBu',limits=c(-0.5,0.5))+
  theme(legend.position="none")  

ad_edger_mean_diff<-ad_edger_mean
ad_edger_mean_diff$R2<-ad_edger_bc_mean$R2-ad_edger_mean$R2
#FIG experiment 1 RElAbundance R2 Bray----
#windows(width=4,height = 3);
d4<-ggplot(ad_edger_mean_diff, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  # geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  # scale_fill_gradient(low="white", high="turquoise4",limits=c(-0.1,.1)) 
  geom_text(aes(label = round(R2,digits = 2)), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "FP diff")) +
  scale_fill_distiller(palette = 'RdBu',limits=c(-0.5,0.5))+ 
  theme(legend.position="right")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(d4)

graphics.off()
pdf('FIG_EXP2_R2_differences_.pdf',width=18,height = 4)
tiff("FIG_EXP2_R2_differences.tiff", units="in", width=18,height = 4, res=300)
gridExtra::grid.arrange(d1,d2,d3,d4,nrow=1);dev.off()



pdf('FIG_EXP2_R2_differences_.pdf',width=18,height = 4)
ggarrange(d1, d2, d3, d4, ncol=1, nrow=1, common.legend = TRUE, legend="right")
dev.off()


