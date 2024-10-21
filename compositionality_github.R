#test of compositionality corrections

library(edgeR)
library(compositions)
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(vegan)



meta = data.frame(timepoint = rep(c("1","2"),each=20),
                  treatment = rep(c("c","tr","c","tr"),each=10))

absdata <- data.frame(s1=round(exp(rnorm(100,50,10)*0.25)))

for(i in 2:20){
absdata[,paste("s",i,sep="")] <- absdata$s1 * exp(rnorm(100,0,0.5))
}
absdata[absdata<0]<-0
colnames(absdata) <- paste(colnames(absdata),rep(c("c","tr"),each=10),sep="_")

windows();barplot(as.matrix(absdata[rev(order(rowMeans(absdata))),]),border=F,cex.names = 0.8,col=1:100,las=2)



palette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F",
          brewer.pal(9,"Greens")[-1],
          brewer.pal(9,"Purples")[-1],
          brewer.pal(9,"YlGnBu")[-1],
          rep(c("gray90","gray40"),80)))
pdf('Abundance.pdf',width=7,height = 10)
par(mar = c(7,7,1,1),mgp=c(4.5,1,1)) 
  barplot(as.matrix(absdata[rev(order(rowMeans(absdata))),]),
                  border=F,cex.names = 1,col=1:100, las=2,
                  ylab="Absolute abundance")
dev.off()
 
pdf('RelativeAbundance.pdf',width=7,height = 10)
par(mar = c(7,7,1,1),mgp=c(4.5,1,1)) 
barplot(as.matrix(absdata[rev(order(rowMeans(absdata))),]),
        border=F,cex.names = 1,col=1:100, las=2,
        ylab="Absolute abundance")
dev.off()

 windows();barplot(as.matrix(t(t(absdata[rev(order(rowMeans(absdata))),])/colSums(absdata[rev(order(rowMeans(absdata))),]))),border=F,cex.names = 0.5,col=1:100)

treatmentdata <- cbind(absdata,absdata)
for(i in 21:40){
  treatmentdata[,i]<-round(treatmentdata[,i]*exp(rnorm(100,0,sd=0.25)))#0.1
}
colnames(treatmentdata) <- paste(colnames(treatmentdata),rep(c("1","2"),each=20),sep="_")

windows();barplot(as.matrix(treatmentdata[rev(order(rowMeans(treatmentdata))),]),border=F,cex.names = 0.5,col=1:100,las=2)
windows();barplot(as.matrix(t(t(treatmentdata[rev(order(rowMeans(treatmentdata))),])/colSums(treatmentdata[rev(order(rowMeans(treatmentdata))),]))),border=F,cex.names = 0.5,col=1:100)

palette(c("#053061", "#2166AC", "#4393C3", "#F4A582", "#B2182B", "plum4","violet","purple1",brewer.pal(3,"Greens")[-1]))


windows(width=3,height = 4.5);
plot(log(t(treatmentdata)[meta$treatment=="c",1])~meta$timepoint[meta$treatment=="c"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Control group",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines(log(t(treatmentdata)[meta$treatment=="c",1][c(i,i+10)])~meta$timepoint[meta$treatment=="c"][c(i,i+10)],
        col=i,lwd=2)
}

windows(width=4,height = 4.5);
plot((t(treatmentdata)[meta$treatment=="c",1])~meta$timepoint[meta$treatment=="c"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Control group",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines((t(treatmentdata)[meta$treatment=="c",1][c(i,i+10)])~meta$timepoint[meta$treatment=="c"][c(i,i+10)],
        col=i,lwd=2)
}


windows();plot(hclust(as.dist(1-cor(log(treatmentdata+1)))))
windows();plot(hclust(as.dist(1-cor((treatmentdata)))))
#0.5,1,2,5,10

treatmentdatalist <- list()
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
   treatmentdatalist[[paste(i,j)]] <- treatmentdata
    treatmentdatalist[[paste(i,j)]][1:(100*j),31:40] <- round(treatmentdatalist[[paste(i,j)]][1:(100*j),31:40]*(1+rnorm(n=10,i,i/2)))
    treatmentdatalist[[paste(i,j)]][treatmentdatalist[[paste(i,j)]]<0] <- 1
     }
}


windows();plot(hclust(as.dist(1-cor(log(treatmentdatalist[[1]])))))
windows();plot(hclust(as.dist(1-cor(log(treatmentdatalist[[25]])))))


windows();boxplot(log(t(treatmentdatalist[[25]])[,1])~paste(meta$treatment,meta$timepoint))
windows();boxplot(log(t(treatmentdatalist[[1]])[,1])~paste(meta$treatment,meta$timepoint))


windows(width=3,height = 4.5);
plot(log(t(treatmentdatalist[['1 0.01']])[meta$treatment=="tr",1])~meta$timepoint[meta$treatment=="tr"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Treatment group 10%",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines(log(t(treatmentdatalist[['1 0.01']])[meta$treatment=="tr",1][c(i,i+10)])~meta$timepoint[meta$treatment=="tr"][c(i,i+10)],
        col=i,lwd=2)
}

windows(width=4,height = 4.5);
plot((t(treatmentdatalist[['100 0.25']])[meta$treatment=="tr",1])~meta$timepoint[meta$treatment=="tr"],
     col=1:10,pch=20,cex=2,xlim=c(0.5,2.5),xlab="Timepoint",ylab="Log abundance of feature 1",
     main="Treatment group",xaxt="n")
axis(side=1,at=c(1,2))
for(i in 1:10){
  lines((t(treatmentdatalist[['100 0.25']])[meta$treatment=="tr",1][c(i,i+10)])~meta$timepoint[meta$treatment=="tr"][c(i,i+10)],
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
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(treatmentdatalist[[paste(i,j)]]))
    tr1[tr1<1]<-1
    tr1$timepoint <- rep(c("1","2"),each=20)
    tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
    
FC[[paste(i,j)]] <- data.frame(c=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]),
                               tr=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100]),
                               bw=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100])/
                                 colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]))
 }
}

gt_tr <- list()
gt_c <- list()
gt_bw <- list()
gt_p <- list()
gt_fc <- list()

for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
tr1 <- log(as.data.frame(t(treatmentdatalist[[paste(i,j)]])))
tr1$timepoint <- rep(c("1","2"),each=20)
tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
delta <- (exp(tr1[tr1$timepoint==2,1:100])-exp(tr1[tr1$timepoint==1,1:100]))/exp(tr1[tr1$timepoint==1,1:100])
delta$treatment <- tr1$treatment[1:20]
gt_tr[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
gt_c[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
gt_bw[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
gt_p[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
gt_fc[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))

for(h in 1:100){
gt_c[[paste(i,j)]][gt_c[[paste(i,j)]]$taxon==h,"p"]<-summary(lm(as.formula(paste("V",h,"~timepoint",sep="")),data=tr1[tr1$treatment=="c",]))$coef[2,4]
gt_tr[[paste(i,j)]][gt_tr[[paste(i,j)]]$taxon==h,"p"]<-summary(lm(as.formula(paste("V",h,"~timepoint",sep="")),data=tr1[tr1$treatment=="tr",]))$coef[2,4]
gt_bw[[paste(i,j)]][gt_bw[[paste(i,j)]]$taxon==h,"p"]<-summary(lm(as.formula(paste("V",h,"~treatment",sep="")),data=delta))$coef[2,4]
gt_p[[paste(i,j)]] <- cbind(c=gt_c[[paste(i,j)]]$p,tr=gt_tr[[paste(i,j)]]$p,bw=gt_bw[[paste(i,j)]]$p) 
}

}
}




peardist_log <- function(x,distance){
  as.dist(1-cor(t(log(x+min(x[x>0])))))
}

peardist <- function(x,distance){
  as.dist(1-cor(t((x))))
}


i=0.5
j=0.25

correl.sym <- gt_p[[paste(i,j)]][,c(1:2)]
correl.sym[p.adjust(gt_p[[paste(i,j)]][,c(1)],"fdr")<0.05,1]<-"#"
correl.sym[p.adjust(gt_p[[paste(i,j)]][,c(1)],"fdr")>0.05,1]<-""
correl.sym[p.adjust(gt_p[[paste(i,j)]][,c(2)],"fdr")<0.05,2]<-"#"
correl.sym[p.adjust(gt_p[[paste(i,j)]][,c(2)],"fdr")>0.05,2]<-""
correl.sym1 <- correl.sym

correl.sym <- correl.sym
correl.sym[p.adjust(gt_p[[paste(i,j)]][,c(3)],"fdr")<0.05,2]<-"*"
correl.sym[p.adjust(gt_p[[paste(i,j)]][,c(3)],"fdr")>0.05,2]<-""

symb <- gt_p[[paste(i,j)]]
for(h in 1:nrow(gt_p[[paste(i,j)]])){
  for(k in 1:2){
    symb[h,k] <- paste(correl.sym1[h,k] ,correl.sym[h,k])
  }}

windows();
gplots::heatmap.2(as.matrix(FC[[paste(i,j)]][,c(1:2)]),
                  col=brewer.pal(n=9,"BuPu")[1:7],
                           density.info = "none",trace="none",dendrogram = "none",
                           Colv=F,cellnote=symb,notecol = "red",cexRow = 0.3,notecex = 0.8)

palette(scales::alpha(c("plum1","purple","skyblue","royalblue"),0.750))


ad_abs<- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))
ad_abs_bc<- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))

pdf("abs_pcoat.pdf")
par(mfcol=c(5,5),mar=c(0.1,0.1,0.1,0.1),mgp=c(0.1,0,0),tck=-0.01,oma=c(2,2,2,2))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
tr1 <- as.data.frame(t(treatmentdatalist[[paste(i,j)]]))
ad_abs[ad_abs$effsize==i&ad_abs$Nchange==j,"R2"] <- adonis2(peardist_log(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$R2[1]
ad_abs[ad_abs$effsize==i&ad_abs$Nchange==j,"p"] <- adonis2(peardist_log(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$Pr[1]
abs_pcoa <- capscale(peardist_log(tr1)~1)
plot(summary(abs_pcoa)$sites[,2]~summary(abs_pcoa)$sites[,1],
              col=as.factor(paste(meta$treatment,meta$timepoint)),pch=20,cex=2,
              ylab="Component 2",xlab="Component 1",axes=F)
axis(side=1,labels=F)
axis(side=2,labels=F)
mtext(side=3,line=-2,text=ad_abs[ad_abs$effsize==i&ad_abs$Nchange==j,"p"])
}}
mtext(side=2,adj=-5, text="Proportion responding",line=40)
mtext(side=1, adj=-5,text="Effect size",line=1)
  
dev.off()


windows(width=5,height = 4.5);
ggplot(ad_abs, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 


pdf("abs_pcoa_bc.pdf")
par(mfcol=c(5,5),mar=c(0.1,0.1,0.1,0.1),mgp=c(0.1,0,0),tck=-0.01,oma=c(2,2,2,2))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(treatmentdatalist[[paste(i,j)]]))
    ad_abs_bc[ad_abs_bc$effsize==i&ad_abs_bc$Nchange==j,"R2"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "bray")$R2[1]
    ad_abs_bc[ad_abs_bc$effsize==i&ad_abs_bc$Nchange==j,"p"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "bray")$Pr[1]
    abs_pcoa <- capscale((tr1)~1,method = "bray")
    plot(summary(abs_pcoa)$sites[,2]~summary(abs_pcoa)$sites[,1],
         col=as.factor(paste(meta$treatment,meta$timepoint)),pch=20,cex=2,
         ylab="Component 2",xlab="Component 1",axes=F)
    axis(side=1,labels=F)
    axis(side=2,labels=F)
    mtext(side=3,line=-2,text=ad_abs[ad_abs$effsize==i&ad_abs$Nchange==j,"p"])
  }}
mtext(side=2,adj=-5, text="Proportion responding",line=40)
mtext(side=1, adj=-5,text="Effect size",line=1)

dev.off()


windows(width=5,height = 4.5);
ggplot(ad_abs_bc, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 


result <- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    pvals <- as.data.frame(gt_p[[paste(i,j)]])
    pvals$adj <- p.adjust(pvals[,3],"fdr")
    result$true_pos[result$effsize==i&result$Nchange==j]<- j*100
    result$true_neg[result$effsize==i&result$Nchange==j]<- 100-(j*100)
    result$discovered_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[pvals[,"adj"]<0.05,])
    result$discovered_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[pvals[,"adj"]>0.05,])
    result$correct_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.05,])
    result$false_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.05,])
    result$false_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.05,])
    result$correct_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.05,])
  }}
result$false_pos_prop <- round(result$false_pos/result$discovered_pos,digits=2)
result$false_pos_prop[result$false_pos==0]<-0
result$false_neg_prop <- round(result$false_neg/result$discovered_neg,digits=2)
result$false_neg_prop[result$false_neg==0]<-0

result$false_pos2 <- result$discovered_pos2-result$true_pos
result$false_neg2 <- result$discovered_neg2-result$true_neg
result$false_pos2[result$false_pos2<0]<-0
result$false_neg2[result$false_neg2<0]<-0
result$false_pos2_prop <- round(result$false_pos2/result$discovered_pos2,digits=2)
result$false_pos2_prop[result$false_pos2==0]<-0
result$false_neg2_prop <- round(result$false_neg2/result$discovered_neg2,digits=2)
result$false_neg2_prop[result$false_neg2==0]<-0


windows(width=5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False pos"))+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False neg"))+
  scale_fill_gradient(low="white", high="blue",limits=c(0,1)) 


windows(width=5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False pos"))+
  scale_fill_gradient(low="white", high="red",limits=c(0,100)) 

windows(width=5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_neg), color = "white", size = 4) +
  guides(fill = guide_colourbar(title = "False neg"))+
  scale_fill_gradient(low="white", high="blue",limits=c(0,100)) 



#relative data--------------------------------------------------------------
palette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F",
          brewer.pal(9,"Greens")[-1],
          brewer.pal(9,"Purples")[-1],
          brewer.pal(9,"YlGnBu")[-1],rep(c("gray90","gray40"),80)))

seqdatalist <- treatmentdatalist
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
    seqdata <- t(t(treatmentdatalist[[paste(i,j)]])/rowSums(t(treatmentdatalist[[paste(i,j)]])))
    for(h in c(1:40)) {
      seqdatalist[[paste(i,j)]][,h]<-round(seqdata[,h]*rnorm(n=1,mean=50000,sd=10000))
    }
    }}

windows();barplot(as.matrix(treatmentdatalist[[1]][rev(order(rowMeans(treatmentdatalist[[1]]))),]),border=F,cex.names = 0.5,col=1:100)
windows();barplot(as.matrix(seqdata[rev(order(rowMeans(seqdata))),]),border=F,cex.names = 0.5,col=1:100)
windows();barplot(as.matrix(seqdatalist[[1]][rev(order(rowMeans(seqdatalist[[1]]))),]),border=F,cex.names = 0.5,col=1:100,las=2)

reldatalist <- list()
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
    reldatalist[[paste(i,j)]] <- t(t(seqdatalist[[paste(i,j)]])/rowSums(t(seqdatalist[[paste(i,j)]])))
    }}

windows(height=4);
barplot(as.matrix(reldatalist[[1]][rev(order(rowMeans(reldatalist[[1]]))),]),border=F,cex.names = 0.5,col=1:100,las=2)

windows(height=4);
barplot(as.matrix(reldatalist[[1]][rev(order(rowMeans(reldatalist[[1]]))),])[,c(1,21)],border=F,cex.names = 0.5,col=1:100,las=2)

FC_rel <- list()
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(reldatalist[[paste(i,j)]]))
    tr1[tr1==0]<-0.0000001
    tr1$timepoint <- rep(c("1","2"),each=20)
    tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
 FC_rel[[paste(i,j)]] <- data.frame(c=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]),
                                   tr=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100]),
                                   bw=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100])/
                                     colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]))
  }
}

gt_tr_rel <- list()
gt_c_rel <- list()
gt_bw_rel <- list()
gt_p_rel <- list()
gt_fc_rel <- list()

for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(reldatalist[[paste(i,j)]]))
    tr1[tr1==0]<-0.0000001
    tr1$timepoint <- rep(c("1","2"),each=20)
    tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
    delta <- tr1[tr1$timepoint==2,1:100]-tr1[tr1$timepoint==1,1:100]
    delta$treatment <- tr1$treatment[1:20]
    
    
    gt_tr_rel[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_c_rel[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_bw_rel[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_p_rel[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_fc_rel[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    
    for(h in 1:100){
      gt_c_rel[[paste(i,j)]][gt_c_rel[[paste(i,j)]]$taxon==h,"p"]<- kruskal.test(g=tr1[tr1$treatment=="c","timepoint"],x=tr1[tr1$treatment=="c",paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~timepoint",sep="")),data=tr1[tr1$treatment=="c",]))$coef[2,4]
      gt_tr_rel[[paste(i,j)]][gt_tr_rel[[paste(i,j)]]$taxon==h,"p"]<-kruskal.test(g=tr1[tr1$treatment=="tr","timepoint"],x=tr1[tr1$treatment=="tr",paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~timepoint",sep="")),data=tr1[tr1$treatment=="tr",]))$coef[2,4]
      gt_bw_rel[[paste(i,j)]][gt_bw_rel[[paste(i,j)]]$taxon==h,"p"]<-kruskal.test(g=delta[,"treatment"],x=delta[,paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~treatment",sep="")),data=delta))$coef[2,4]
      gt_p_rel[[paste(i,j)]] <- cbind(c=gt_c_rel[[paste(i,j)]]$p,tr=gt_tr_rel[[paste(i,j)]]$p,bw=gt_bw_rel[[paste(i,j)]]$p) 
    }
    
  
  }
}


i=2
j=0.1

correl.sym <- gt_p_rel[[paste(i,j)]][,c(1:2)]
correl.sym[p.adjust(gt_p_rel[[paste(i,j)]][,c(1)],"fdr")<0.05,1]<-"#"
correl.sym[p.adjust(gt_p_rel[[paste(i,j)]][,c(1)],"fdr")>0.05,1]<-""
correl.sym[p.adjust(gt_p_rel[[paste(i,j)]][,c(2)],"fdr")<0.05,2]<-"#"
correl.sym[p.adjust(gt_p_rel[[paste(i,j)]][,c(2)],"fdr")>0.05,2]<-""
correl.sym1 <- correl.sym

correl.sym <- correl.sym
correl.sym[p.adjust(gt_p_rel[[paste(i,j)]][,c(3)],"fdr")<0.05,2]<-"*"
correl.sym[p.adjust(gt_p_rel[[paste(i,j)]][,c(3)],"fdr")>0.05,2]<-""

symb <- gt_p_rel[[paste(i,j)]]
for(h in 1:nrow(gt_p_rel[[paste(i,j)]])){
  for(k in 1:2){
    symb[h,k] <- paste(correl.sym1[h,k] ,correl.sym[h,k])
  }} 
windows();
gplots::heatmap.2(as.matrix(FC_rel[[paste(i,j)]][,c(1:2)]),
                  col=brewer.pal(n=9,"BuPu")[1:7],
                  #col=c( "aliceblue", "#F7F7F7", "#FDDBC7","#F4A582","#D6604D", "#B2182B", "#67001F"),#rainbow(256, start=0,end=0.34),
                  density.info = "none",trace="none",dendrogram = "none",
                  Colv=F,cellnote=symb,notecol = "red",cexRow = 0.3,notecex = 0.8)


result <- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
  pvals <- as.data.frame(gt_p_rel[[paste(i,j)]])
  pvals$adj <- p.adjust(pvals[,3],"fdr")
  pvals[pvals=="NaN"]<-1
  result$true_pos[result$effsize==i&result$Nchange==j]<- j*100
  result$true_neg[result$effsize==i&result$Nchange==j]<- 100-(j*100)
  result$discovered_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[pvals[,"adj"]<0.05,])
  result$discovered_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[pvals[,"adj"]>0.05,])
  result$correct_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.05,])
  result$false_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.05,])
  result$false_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.05,])
  result$correct_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.05,])
  }}

result$false_pos_prop <- round(result$false_pos/result$discovered_pos,digits=2)
result$false_pos_prop[result$false_pos==0]<-0
result$false_neg_prop <- round(result$false_neg/result$discovered_neg,digits=2)
result$false_neg_prop[result$false_neg==0]<-0


windows(width=4.5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = "none")+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4.5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="blue",limits=c(0,1)) 

fc_comp <- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=500),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),500))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
fc_comp$FC[fc_comp$effsize==i&fc_comp$Nchange==j]<- FC[[paste(i,j)]]$bw
fc_comp$FC_rel[fc_comp$effsize==i&fc_comp$Nchange==j]<- FC_rel[[paste(i,j)]]$bw
  }}

X11();plot(log(FC_rel)~log(FC),data=fc_comp,pch=20)

ad_rel<- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))

pdf("rel_pcoa.pdf")
par(mfrow=c(5,5),mar=c(0.1,0.1,0.1,0.1),mgp=c(0.1,0,0),tck=-0.01,oma=c(2,2,2,2))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(reldatalist[[paste(i,j)]]))
    ad_rel[ad_rel$effsize==i&ad_rel$Nchange==j,"R2"] <- adonis2(peardist_log(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$R2[1]
    ad_rel[ad_rel$effsize==i&ad_rel$Nchange==j,"p"] <- adonis2(peardist_log(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$Pr[1]
    rel_pcoa <- capscale(peardist_log(tr1)~1)
    plot(summary(rel_pcoa)$sites[,2]~summary(rel_pcoa)$sites[,1],
         col=as.factor(paste(meta$treatment,meta$timepoint)),pch=20,cex=2,
         ylab="Component 2",xlab="Component 1",axes=F)
    axis(side=1,labels=F)
    axis(side=2,labels=F)
    mtext(side=3,line=-2,text=ad_rel[ad_rel$effsize==i&ad_rel$Nchange==j,"p"])
  }}
mtext(side=2,adj=-5, text="Proportion responding",line=40)
mtext(side=1, adj=-5,text="Effect size",line=1)

dev.off()


windows(width=4.5,height = 4.5);
ggplot(ad_rel, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

ad_rel_bc<- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))

pdf("rel_pcoa_bc.pdf")
par(mfrow=c(5,5),mar=c(0.1,0.1,0.1,0.1),mgp=c(0.1,0,0),tck=-0.01,oma=c(2,2,2,2))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(reldatalist[[paste(i,j)]]))
    ad_rel_bc[ad_rel_bc$effsize==i&ad_rel_bc$Nchange==j,"R2"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method="bray")$R2[1]
    ad_rel_bc[ad_rel_bc$effsize==i&ad_rel_bc$Nchange==j,"p"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method="bray")$Pr[1]
    rel_pcoa <- capscale((tr1)~1,method="bray")
  plot(summary(rel_pcoa)$sites[,2]~summary(rel_pcoa)$sites[,1],
         col=as.factor(paste(meta$treatment,meta$timepoint)),pch=20,cex=2,
         ylab="Component 2",xlab="Component 1",axes=F)
    axis(side=1,labels=F)
    axis(side=2,labels=F)
    mtext(side=3,line=-2,text=ad_rel_bc[ad_rel_bc$effsize==i&ad_rel$Nchange==j,"p"])
  }}
mtext(side=2,adj=-5, text="Proportion responding",line=40)
mtext(side=1, adj=-5,text="Effect size",line=1)

dev.off()


windows(width=4.5,height = 4.5);
ggplot(ad_rel_bc, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 


#clr data--------------------------------------------------------------


clrdatalist <- list()
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
    clrdatalist[[paste(i,j)]] <- t(seqdatalist[[paste(i,j)]])
    clrdatalist[[paste(i,j)]] <- t(clr(clrdatalist[[paste(i,j)]]))
      }}


FC_clr <- list()
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(clrdatalist[[paste(i,j)]]))
    tr1[tr1==0]<-0.0000001
    tr1$timepoint <- rep(c("1","2"),each=20)
    tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
    FC_clr[[paste(i,j)]] <- data.frame(c=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]),
                                       tr=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100]),
                                       bw=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100])/
                                         colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]))
  }
}

gt_tr_clr <- list()
gt_c_clr <- list()
gt_bw_clr <- list()
gt_p_clr <- list()
gt_fc_clr <- list()

for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(clrdatalist[[paste(i,j)]]))
    tr1[tr1==0]<-0.0000001
    tr1$timepoint <- rep(c("1","2"),each=20)
    tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
    delta <- tr1[tr1$timepoint==2,1:100]-tr1[tr1$timepoint==1,1:100]
    delta$treatment <- tr1$treatment[1:20]
    
    
    #windows();boxplot(V2~paste(treatment,timepoint),data=tr1)
    #windows();boxplot(V1~paste(treatment),data=delta)
    #windows();hist(tr1$V1)
    
    gt_tr_clr[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_c_clr[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_bw_clr[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_p_clr[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_fc_clr[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    
    for(h in 1:100){
      gt_c_clr[[paste(i,j)]][gt_c_clr[[paste(i,j)]]$taxon==h,"p"]<- kruskal.test(g=tr1[tr1$treatment=="c","timepoint"],x=tr1[tr1$treatment=="c",paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~timepoint",sep="")),data=tr1[tr1$treatment=="c",]))$coef[2,4]
      gt_tr_clr[[paste(i,j)]][gt_tr_clr[[paste(i,j)]]$taxon==h,"p"]<-kruskal.test(g=tr1[tr1$treatment=="tr","timepoint"],x=tr1[tr1$treatment=="tr",paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~timepoint",sep="")),data=tr1[tr1$treatment=="tr",]))$coef[2,4]
      gt_bw_clr[[paste(i,j)]][gt_bw_clr[[paste(i,j)]]$taxon==h,"p"]<-kruskal.test(g=delta[,"treatment"],x=delta[,paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~treatment",sep="")),data=delta))$coef[2,4]
      gt_p_clr[[paste(i,j)]] <- cbind(c=gt_c_clr[[paste(i,j)]]$p,tr=gt_tr_clr[[paste(i,j)]]$p,bw=gt_bw_clr[[paste(i,j)]]$p) 
    }
    

  }
}
   
 corclr.sym <- gt_p_clr[[paste(i,j)]][,c(1:2)]
    corclr.sym[p.adjust(gt_p_clr[[paste(i,j)]][,c(1)],"fdr")<0.05,1]<-"#"
    corclr.sym[p.adjust(gt_p_clr[[paste(i,j)]][,c(1)],"fdr")>0.05,1]<-""
    corclr.sym[p.adjust(gt_p_clr[[paste(i,j)]][,c(2)],"fdr")<0.05,2]<-"#"
    corclr.sym[p.adjust(gt_p_clr[[paste(i,j)]][,c(2)],"fdr")>0.05,2]<-""
    corclr.sym1 <- corclr.sym
    
    corclr.sym <- corclr.sym
    corclr.sym[p.adjust(gt_p_clr[[paste(i,j)]][,c(3)],"fdr")<0.05,2]<-"*"
    corclr.sym[p.adjust(gt_p_clr[[paste(i,j)]][,c(3)],"fdr")>0.05,2]<-""
    
    symb <- gt_p_clr[[paste(i,j)]]
    for(h in 1:nrow(gt_p_clr[[paste(i,j)]])){
      for(k in 1:2){
        symb[h,k] <- paste(corclr.sym1[h,k] ,corclr.sym[h,k])
      }}
    
windows();
gplots::heatmap.2(as.matrix(FC_clr[[paste(i,j)]][,c(1:2)]),
                  col=brewer.pal(n=9,"BuPu")[1:7],
                  #col=c( "aliceblue", "#F7F7F7", "#FDDBC7","#F4A582","#D6604D", "#B2182B", "#67001F"),#rainbow(256, start=0,end=0.34),
                  density.info = "none",trace="none",dendrogram = "none",
                  Colv=F,cellnote=symb,notecol = "red",cexRow = 0.3,notecex = 0.8)


result <- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    pvals <- as.data.frame(gt_p_clr[[paste(i,j)]])
    pvals$adj <- p.adjust(pvals[,3],"fdr")
    pvals[pvals=="NaN"]<-1
    result$true_pos[result$effsize==i&result$Nchange==j]<- j*100
    result$true_neg[result$effsize==i&result$Nchange==j]<- 100-(j*100)
    result$discovered_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[pvals[,"adj"]<0.05,])
    result$discovered_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[pvals[,"adj"]>0.05,])
    result$correct_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.05,])
    result$false_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.05,])
    result$false_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.05,])
    result$correct_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.05,])
  }}

result$false_pos_prop <- round(result$false_pos/result$discovered_pos,digits=2)
result$false_pos_prop[result$false_pos==0]<-0
result$false_neg_prop <- round(result$false_neg/result$discovered_neg,digits=2)
result$false_neg_prop[result$false_neg==0]<-0



windows(width=4.5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = guide_colourbar(title = "False positives"))+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = "none")+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4.5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
 # guides(fill = guide_colourbar(title = "False negatives"))+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="blue",limits=c(0,1)) 



ad_clr<- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))
ad_clr_pear<- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))

pdf("clr_pcoa.pdf")
par(mfrow=c(5,5),mar=c(0.1,0.1,0.1,0.1),mgp=c(0.1,0,0),tck=-0.01,oma=c(2,2,2,2))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(clrdatalist[[paste(i,j)]]))
    ad_clr[ad_clr$effsize==i&ad_clr$Nchange==j,"R2"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "euc")$R2[1]
    ad_clr[ad_clr$effsize==i&ad_clr$Nchange==j,"p"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "euc")$Pr[1]
    clr_pcoa <- capscale((tr1)~1,method = "euc")
    plot(summary(clr_pcoa)$sites[,2]~summary(clr_pcoa)$sites[,1],
         col=as.factor(paste(meta$treatment,meta$timepoint)),pch=20,cex=2,
         ylab="Component 2",xlab="Component 1",axes=F)
    axis(side=1,labels=F)
    axis(side=2,labels=F)
    #mtext(side=3,line=-2,text=ad_clr[ad_clr$effsize==i&ad_clr$Nchange==j,"p"])
  }}
#mtext(side=2,adj=-5, text="Proportion responding",line=40)
#mtext(side=1, adj=-5,text="Effect size",line=1)

dev.off()


windows(width=4.5,height = 4.5);
ggplot(ad_clr, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 


pdf("clr_pcoa_pear.pdf")
par(mfrow=c(5,5),mar=c(0.1,0.1,0.1,0.1),mgp=c(0.1,0,0),tck=-0.01,oma=c(2,2,2,2))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(clrdatalist[[paste(i,j)]]))
    ad_clr_pear[ad_clr_pear$effsize==i&ad_clr_pear$Nchange==j,"R2"] <- adonis2(peardist(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$R2[1]
    ad_clr_pear[ad_clr_pear$effsize==i&ad_clr_pear$Nchange==j,"p"] <- adonis2(peardist(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$Pr[1]
    clr_pcoa_pear <- capscale(peardist(tr1)~1)
    plot(summary(clr_pcoa_pear)$sites[,2]~summary(clr_pcoa_pear)$sites[,1],
         col=as.factor(paste(meta$treatment,meta$timepoint)),pch=20,cex=2,
         ylab="Component 2",xlab="Component 1",axes=F)
    axis(side=1,labels=F)
    axis(side=2,labels=F)
    #mtext(side=3,line=-2,text=ad_clr[ad_clr$effsize==i&ad_clr$Nchange==j,"p"])
  }}
#mtext(side=2,adj=-5, text="Proportion responding",line=40)
#mtext(side=1, adj=-5,text="Effect size",line=1)

dev.off()


windows(width=4.5,height = 4.5);
ggplot(ad_clr_pear, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 


#deseq data--------------------------------------------------------------


deseqdatalist <- list()
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
   dds <- DESeqDataSetFromMatrix(countData = (seqdatalist[[paste(i,j)]]), 
                                  colData =(meta),
                                  design = ~ timepoint + treatment + timepoint:treatment)
  dds <- estimateSizeFactors(dds)
   deseqdatalist[[paste(i,j)]]  <- counts(dds, normalized=TRUE)
   }}


FC_deseq <- list()
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(deseqdatalist[[paste(i,j)]]))
    tr1[tr1==0]<-0.0000001
    tr1$timepoint <- rep(c("1","2"),each=20)
    tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
    FC_deseq[[paste(i,j)]] <- data.frame(c=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]),
                                       tr=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100]),
                                       bw=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100])/
                                         colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]))
  }
}

gt_tr_deseq <- list()
gt_c_deseq <- list()
gt_bw_deseq <- list()
gt_p_deseq <- list()
gt_fc_deseq <- list()

for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(deseqdatalist[[paste(i,j)]]))
    tr1[tr1==0]<-0.0000001
    tr1$timepoint <- rep(c("1","2"),each=20)
    tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
    delta <- tr1[tr1$timepoint==2,1:100]-tr1[tr1$timepoint==1,1:100]
    delta$treatment <- tr1$treatment[1:20]
    
    
    #windows();boxplot(V2~paste(treatment,timepoint),data=tr1)
    #windows();boxplot(V1~paste(treatment),data=delta)
    #windows();hist(tr1$V1)
    
    gt_tr_deseq[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_c_deseq[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_bw_deseq[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_p_deseq[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_fc_deseq[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    
    for(h in 1:100){
      gt_c_deseq[[paste(i,j)]][gt_c_deseq[[paste(i,j)]]$taxon==h,"p"]<- kruskal.test(g=tr1[tr1$treatment=="c","timepoint"],x=tr1[tr1$treatment=="c",paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~timepoint",sep="")),data=tr1[tr1$treatment=="c",]))$coef[2,4]
      gt_tr_deseq[[paste(i,j)]][gt_tr_deseq[[paste(i,j)]]$taxon==h,"p"]<-kruskal.test(g=tr1[tr1$treatment=="tr","timepoint"],x=tr1[tr1$treatment=="tr",paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~timepoint",sep="")),data=tr1[tr1$treatment=="tr",]))$coef[2,4]
      gt_bw_deseq[[paste(i,j)]][gt_bw_deseq[[paste(i,j)]]$taxon==h,"p"]<-kruskal.test(g=delta[,"treatment"],x=delta[,paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~treatment",sep="")),data=delta))$coef[2,4]
      gt_p_deseq[[paste(i,j)]] <- cbind(c=gt_c_deseq[[paste(i,j)]]$p,tr=gt_tr_deseq[[paste(i,j)]]$p,bw=gt_bw_deseq[[paste(i,j)]]$p) 
    }
    
  }
}




cordeseq.sym <- gt_p_deseq[[paste(i,j)]][,c(1:2)]
cordeseq.sym[p.adjust(gt_p_deseq[[paste(i,j)]][,c(1)],"fdr")<0.05,1]<-"#"
cordeseq.sym[p.adjust(gt_p_deseq[[paste(i,j)]][,c(1)],"fdr")>0.05,1]<-""
cordeseq.sym[p.adjust(gt_p_deseq[[paste(i,j)]][,c(2)],"fdr")<0.05,2]<-"#"
cordeseq.sym[p.adjust(gt_p_deseq[[paste(i,j)]][,c(2)],"fdr")>0.05,2]<-""
cordeseq.sym1 <- cordeseq.sym

cordeseq.sym <- cordeseq.sym
cordeseq.sym[p.adjust(gt_p_deseq[[paste(i,j)]][,c(3)],"fdr")<0.05,2]<-"*"
cordeseq.sym[p.adjust(gt_p_deseq[[paste(i,j)]][,c(3)],"fdr")>0.05,2]<-""

symb <- gt_p_deseq[[paste(i,j)]]
for(h in 1:nrow(gt_p_deseq[[paste(i,j)]])){
  for(k in 1:2){
    symb[h,k] <- paste(cordeseq.sym1[h,k] ,cordeseq.sym[h,k])
  }}
windows();
gplots::heatmap.2(as.matrix(FC_deseq[[paste(i,j)]][,c(1:2)]),
                  col=brewer.pal(n=9,"BuPu")[1:7],
                  #col=c( "aliceblue", "#F7F7F7", "#FDDBC7","#F4A582","#D6604D", "#B2182B", "#67001F"),#rainbow(256, start=0,end=0.34),
                  density.info = "none",trace="none",dendrogram = "none",
                  Colv=F,cellnote=symb,notecol = "red",cexRow = 0.3,notecex = 0.8)


result <- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    pvals <- as.data.frame(gt_p_deseq[[paste(i,j)]])
    pvals$adj <- p.adjust(pvals[,3],"fdr")
    pvals[pvals=="NaN"]<-1
    result$true_pos[result$effsize==i&result$Nchange==j]<- j*100
    result$true_neg[result$effsize==i&result$Nchange==j]<- 100-(j*100)
    result$discovered_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[pvals[,"adj"]<0.05,])
    result$discovered_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[pvals[,"adj"]>0.05,])
    result$correct_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.05,])
    result$false_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.05,])
    result$false_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.05,])
    result$correct_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.05,])
  }}

result$false_pos_prop <- round(result$false_pos/result$discovered_pos,digits=2)
result$false_pos_prop[result$false_pos==0]<-0
result$false_neg_prop <- round(result$false_neg/result$discovered_neg,digits=2)
result$false_neg_prop[result$false_neg==0]<-0


windows(width=4.5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = "none")+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4.5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  # guides(fill = guide_colourbar(title = "False negatives"))+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="blue",limits=c(0,1)) 


fc_comp <- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=500),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),500))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    fc_comp$FC[fc_comp$effsize==i&fc_comp$Nchange==j]<- FC[[paste(i,j)]]$bw
    fc_comp$FC_deseq[fc_comp$effsize==i&fc_comp$Nchange==j]<- FC_deseq[[paste(i,j)]]$bw
  }}


windows();plot(log(FC_deseq)~log(FC),data=fc_comp,pch=20)

ad_deseq<- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))
ad_deseq_pear<- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))

pdf("deseq_pcoa.pdf")
par(mfrow=c(5,5),mar=c(0.1,0.1,0.1,0.1),mgp=c(0.1,0,0),tck=-0.01,oma=c(2,2,2,2))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(deseqdatalist[[paste(i,j)]]))
    ad_deseq[ad_deseq$effsize==i&ad_deseq$Nchange==j,"R2"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "bray")$R2[1]
    ad_deseq[ad_deseq$effsize==i&ad_deseq$Nchange==j,"p"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "bray")$Pr[1]
    deseq_pcoa <- capscale((tr1)~1,method = "bray")
    plot(summary(deseq_pcoa)$sites[,2]~summary(deseq_pcoa)$sites[,1],
         col=as.factor(paste(meta$treatment,meta$timepoint)),pch=20,cex=2,
         ylab="Component 2",xlab="Component 1",axes=F)
    axis(side=1,labels=F)
    axis(side=2,labels=F)
  }}

dev.off()


windows(width=4.5,height = 4.5);
ggplot(ad_deseq, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 


pdf("deseq_pcoa_pear.pdf")
par(mfrow=c(5,5),mar=c(0.1,0.1,0.1,0.1),mgp=c(0.1,0,0),tck=-0.01,oma=c(2,2,2,2))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(deseqdatalist[[paste(i,j)]]))
    ad_deseq_pear[ad_deseq_pear$effsize==i&ad_deseq_pear$Nchange==j,"R2"] <- adonis2(peardist(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$R2[1]
    ad_deseq_pear[ad_deseq_pear$effsize==i&ad_deseq_pear$Nchange==j,"p"] <- adonis2(peardist(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$Pr[1]
    deseq_pcoa_pear <- capscale(peardist(tr1)~1)
    plot(summary(deseq_pcoa_pear)$sites[,2]~summary(deseq_pcoa_pear)$sites[,1],
         col=as.factor(paste(meta$treatment,meta$timepoint)),pch=20,cex=2,
         ylab="Component 2",xlab="Component 1",axes=F)
    axis(side=1,labels=F)
    axis(side=2,labels=F)
    #mtext(side=3,line=-2,text=ad_deseq[ad_deseq$effsize==i&ad_deseq$Nchange==j,"p"])
  }}
#mtext(side=2,adj=-5, text="Proportion responding",line=40)
#mtext(side=1, adj=-5,text="Effect size",line=1)

dev.off()


windows(width=4.5,height = 4.5);
ggplot(ad_deseq_pear, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 



#edger data--------------------------------------------------------------

meta = data.frame(timepoint = rep(c("1","2"),each=20),
                  treatment = rep(c("c","tr","c","tr"),each=10))

edgerdatalist <- list()
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){  
    y <- DGEList(counts=round(seqdatalist[[paste(i,j)]]),group=meta$treatment)
    edgerdatalist[[paste(i,j)]]   <- seqdatalist[[paste(i,j)]]*normLibSizes(y)$samples$norm
 }}


FC_edger <- list()
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(edgerdatalist[[paste(i,j)]]))
    tr1[tr1==0]<-0.0000001
    tr1$timepoint <- rep(c("1","2"),each=20)
    tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
    FC_edger[[paste(i,j)]] <- data.frame(c=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]),
                                         tr=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100]),
                                         bw=colMeans(tr1[tr1$timepoint==2&tr1$treatment=="tr",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="tr",1:100])/
                                           colMeans(tr1[tr1$timepoint==2&tr1$treatment=="c",1:100]/tr1[tr1$timepoint==1&tr1$treatment=="c",1:100]))
  }
}

gt_tr_edger <- list()
gt_c_edger <- list()
gt_bw_edger <- list()
gt_p_edger <- list()
gt_fc_edger <- list()

for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(edgerdatalist[[paste(i,j)]]))
    tr1[tr1==0]<-0.0000001
    tr1$timepoint <- rep(c("1","2"),each=20)
    tr1$treatment <- rep(c("c","tr","c","tr"),each=10)
    delta <- tr1[tr1$timepoint==2,1:100]-tr1[tr1$timepoint==1,1:100]
    delta$treatment <- tr1$treatment[1:20]
    
    
    #windows();boxplot(V2~paste(treatment,timepoint),data=tr1)
    #windows();boxplot(V1~paste(treatment),data=delta)
    #windows();hist(tr1$V1)
    
    gt_tr_edger[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_c_edger[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_bw_edger[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_p_edger[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    gt_fc_edger[[paste(i,j)]] <- data.frame(taxon=c(1:100),p=rep(NA,100))
    
    for(h in 1:100){
      gt_c_edger[[paste(i,j)]][gt_c_edger[[paste(i,j)]]$taxon==h,"p"]<- kruskal.test(g=tr1[tr1$treatment=="c","timepoint"],x=tr1[tr1$treatment=="c",paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~timepoint",sep="")),data=tr1[tr1$treatment=="c",]))$coef[2,4]
      gt_tr_edger[[paste(i,j)]][gt_tr_edger[[paste(i,j)]]$taxon==h,"p"]<-kruskal.test(g=tr1[tr1$treatment=="tr","timepoint"],x=tr1[tr1$treatment=="tr",paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~timepoint",sep="")),data=tr1[tr1$treatment=="tr",]))$coef[2,4]
      gt_bw_edger[[paste(i,j)]][gt_bw_edger[[paste(i,j)]]$taxon==h,"p"]<-kruskal.test(g=delta[,"treatment"],x=delta[,paste("V",h,sep="")])$p.value#summary(lm(as.formula(paste("V",h,"~treatment",sep="")),data=delta))$coef[2,4]
      gt_p_edger[[paste(i,j)]] <- cbind(c=gt_c_edger[[paste(i,j)]]$p,tr=gt_tr_edger[[paste(i,j)]]$p,bw=gt_bw_edger[[paste(i,j)]]$p) 
    }
    
 
  }
}

coredger.sym <- gt_p_edger[[paste(i,j)]][,c(1:2)]
coredger.sym[p.adjust(gt_p_edger[[paste(i,j)]][,c(1)],"fdr")<0.05,1]<-"#"
coredger.sym[p.adjust(gt_p_edger[[paste(i,j)]][,c(1)],"fdr")>0.05,1]<-""
coredger.sym[p.adjust(gt_p_edger[[paste(i,j)]][,c(2)],"fdr")<0.05,2]<-"#"
coredger.sym[p.adjust(gt_p_edger[[paste(i,j)]][,c(2)],"fdr")>0.05,2]<-""
coredger.sym1 <- coredger.sym

coredger.sym <- coredger.sym
coredger.sym[p.adjust(gt_p_edger[[paste(i,j)]][,c(3)],"fdr")<0.05,2]<-"*"
coredger.sym[p.adjust(gt_p_edger[[paste(i,j)]][,c(3)],"fdr")>0.05,2]<-""

symb <- gt_p_edger[[paste(i,j)]]
for(h in 1:nrow(gt_p_edger[[paste(i,j)]])){
  for(k in 1:2){
    symb[h,k] <- paste(coredger.sym1[h,k] ,coredger.sym[h,k])
  }}  
windows();
gplots::heatmap.2(as.matrix(FC_edger[[paste(i,j)]][,c(1:2)]),
                  col=brewer.pal(n=9,"BuPu")[1:7],
                  #col=c( "aliceblue", "#F7F7F7", "#FDDBC7","#F4A582","#D6604D", "#B2182B", "#67001F"),#rainbow(256, start=0,end=0.34),
                  density.info = "none",trace="none",dendrogram = "none",
                  Colv=F,cellnote=symb,notecol = "red",cexRow = 0.3,notecex = 0.8)


result <- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    pvals <- as.data.frame(gt_p_edger[[paste(i,j)]])
    pvals$adj <- p.adjust(pvals[,3],"fdr")
    pvals[pvals=="NaN"]<-1
    result$true_pos[result$effsize==i&result$Nchange==j]<- j*100
    result$true_neg[result$effsize==i&result$Nchange==j]<- 100-(j*100)
    result$discovered_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[pvals[,"adj"]<0.05,])
    result$discovered_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[pvals[,"adj"]>0.05,])
    result$correct_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]<0.05,])
    result$false_pos[result$effsize==i&result$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]<0.05,])
    result$false_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[1:(100*j),][pvals[1:(100*j),"adj"]>0.05,])
    result$correct_neg[result$effsize==i&result$Nchange==j] <- nrow(pvals[(100*j+1):100,][pvals[(100*j+1):100,"adj"]>0.05,])
  }}

result$false_pos_prop <- round(result$false_pos/result$discovered_pos,digits=2)
result$false_pos_prop[result$false_pos==0]<-0
result$false_neg_prop <- round(result$false_neg/result$discovered_neg,digits=2)
result$false_neg_prop[result$false_neg==0]<-0


windows(width=4.5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_pos_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  geom_text(aes(label = false_pos_prop), color = "white", size = 4) +
  guides(fill = "none")+
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

windows(width=4.5,height = 4.5);
ggplot(result, aes(as.factor(effsize), as.factor(Nchange), fill= false_neg_prop)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  # guides(fill = guide_colourbar(title = "False negatives"))+
  geom_text(aes(label = false_neg_prop), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="blue",limits=c(0,1)) 


fc_comp <- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=500),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),500))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    fc_comp$FC[fc_comp$effsize==i&fc_comp$Nchange==j]<- FC[[paste(i,j)]]$bw
    fc_comp$FC_edger[fc_comp$effsize==i&fc_comp$Nchange==j]<- FC_edger[[paste(i,j)]]$bw
  }}


windows();plot(log(FC_edger)~log(FC),data=fc_comp,pch=20)

ad_edger<- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))
ad_edger_pear<- data.frame(effsize=rep(c(0.5,1,5,10, 50,100),each=6),Nchange=rep(c(0.01,0.05,0.1,0.25,0.5,0.75),6))

pdf("edger_pcoa.pdf")
par(mfrow=c(5,5),mar=c(0.1,0.1,0.1,0.1),mgp=c(0.1,0,0),tck=-0.01,oma=c(2,2,2,2))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(edgerdatalist[[paste(i,j)]]))
    ad_edger[ad_edger$effsize==i&ad_edger$Nchange==j,"R2"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "bray")$R2[1]
    ad_edger[ad_edger$effsize==i&ad_edger$Nchange==j,"p"] <- adonis2((tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2],method = "bray")$Pr[1]
    edger_pcoa <- capscale((tr1)~1,method = "bray")
    plot(summary(edger_pcoa)$sites[,2]~summary(edger_pcoa)$sites[,1],
         col=as.factor(paste(meta$treatment,meta$timepoint)),pch=20,cex=2,
         ylab="Component 2",xlab="Component 1",axes=F)
    axis(side=1,labels=F)
    axis(side=2,labels=F)
  }}

dev.off()


windows(width=4.5,height = 4.5);
ggplot(ad_edger, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 


pdf("edger_pcoa_pear.pdf")
par(mfrow=c(5,5),mar=c(0.1,0.1,0.1,0.1),mgp=c(0.1,0,0),tck=-0.01,oma=c(2,2,2,2))
for(i in c(0.5,1,5,10, 50,100)){
  for(j in c(0.01,0.05,0.1,0.25,0.5,0.75)){ 
    tr1 <- as.data.frame(t(edgerdatalist[[paste(i,j)]]))
    ad_edger_pear[ad_edger_pear$effsize==i&ad_edger_pear$Nchange==j,"R2"] <- adonis2(peardist_log(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$R2[1]
    ad_edger_pear[ad_edger_pear$effsize==i&ad_edger_pear$Nchange==j,"p"] <- adonis2(peardist_log(tr1[meta$timepoint==2,])~meta$treatment[meta$timepoint==2])$Pr[1]
    edger_pcoa_pear <- capscale(peardist_log(tr1)~1)
    plot(summary(edger_pcoa_pear)$sites[,2]~summary(edger_pcoa_pear)$sites[,1],
         col=as.factor(paste(meta$treatment,meta$timepoint)),pch=20,cex=2,
         ylab="Component 2",xlab="Component 1",axes=F)
    axis(side=1,labels=F)
    axis(side=2,labels=F)
    #mtext(side=3,line=-2,text=ad_edger[ad_edger$effsize==i&ad_edger$Nchange==j,"p"])
  }}
#mtext(side=2,adj=-5, text="Proportion responding",line=40)
#mtext(side=1, adj=-5,text="Effect size",line=1)

dev.off()


windows(width=4.5,height = 4.5);
ggplot(ad_edger_pear, aes(as.factor(effsize), as.factor(Nchange), fill= R2)) + 
  geom_tile(color = "white",
            lwd = 1.5,linetype = 1) +
  ylab("Proportion responding")+
  xlab("Effect size")+
  guides(fill = "none")+
  geom_text(aes(label = round(R2,digits=2)), color = "white", size = 4) +
  scale_fill_gradient(low="white", high="red",limits=c(0,1)) 

