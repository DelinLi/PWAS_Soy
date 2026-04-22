setwd("example") # please adjust the path

## PWAS example of pubescence color and soybean cyst nematode PWAS 
## By Delin Li 20260402

### loading the GAPIT function, this will automatically install GAPIT if you have not installed it.
source("http://zzlab.net/GAPIT/gapit_functions.txt")

### A manhattan plot function, not necessary if you only want to conduct PWAS
source("Plot.Function.R")

myGD<-read.delim("data/Pro.geno.txt.gz",check.names =F)
myGM<-read.delim("data/Pro.map.txt.gz")

Pheno <- read.delim("data/Sample.Information.Phenotype.TableS1.txt")

##########################
### pubescence color
##########################

## convert the qualitative trait category in to numeric values for PWAS
pub.color  <-Pheno[   ,c("ReSeq_ID","Pubescence.Color")]
pub.color$PubColor <- NA
pub.color$PubColor[pub.color$Pubescence.Color %in% "tawny" ] <- 1
pub.color$PubColor[pub.color$Pubescence.Color %in% "gray" ] <- 2

myY <- pub.color[ pub.color$ReSeq_ID  %in% myGD$taxa  & !is.na(pub.color$PubColor),]
myY <- myY[order(myY$ReSeq_ID),]
myGD.tem<-myGD[myGD$taxa %in% myY$ReSeq_ID,]
myGD.tem<-myGD.tem[order(myGD.tem$taxa),]

myGM.tem <- myGM[myGM[,1] %in% colnames(myGD.tem),]
myGD.tem<- myGD.tem[,colnames(myGD.tem) %in% c("taxa",as.character(myGM[,1]))]

myGAPIT <- GAPIT(Y=myY[,-2] ,
                 GD=myGD.tem, GM=myGM.tem,
                 PCA.total=3,
                 model= "CMLM", # can try with different model provided by GAPIT
                 SNP.MAF=0.05,
                 file.output=F
)

write.csv(myGAPIT$GWAS,"Soybean.PWAS_Seedling_PubescenceColor.CMLM.Results.csv",row.names = F)

cutOff.pubcolor <- 0.1 / nrow(myGAPIT$GWAS)

png("Soybean.PWAS_Seedling_PubescenceColor.CMLM.png", width = 18,height=5,units = "in",res = 200)
layout(matrix(1:1, 1, 1, byrow = TRUE))
par(mar = (c(4.2,5,2,2)+ 0.5), mgp=c(3,1.8,0))    
par(bty="l", lwd=1.5)  ## bty=l  the plot is coordinate instead of box

GAPIT.Manhattan.plot(myGAPIT$GWAS[,-1] ,name.of.trait= "Pubescence Color", cutOff= cutOff.pubcolor ,cex.sig = 2,
                     highliht.sig=T,color1="lightgrey", color2="darkgrey",  
                     cex.lab=2,color1.sig="#FFA500", color2.sig="#FFA500") 

dev.off()

print(myGAPIT$GWAS[which.min(myGAPIT$GWAS$P.value),] )
# should be the known protein T (ID  Glyma.06G202300.1.p) with a P.value of "3.166518e-10"


##########################
### soybean cyst nematode
##########################
SCN  <-Pheno[ ! is.na(Pheno$Grade_SCN3 ) ,c("ReSeq_ID","Grade_SCN3")]

myY <- SCN[ SCN$ReSeq_ID  %in% myGD$taxa ,]
myY <- myY[order(myY$ReSeq_ID),]
myGD.tem<-myGD[myGD$taxa %in% myY$ReSeq_ID,]
myGD.tem<-myGD.tem[order(myGD.tem$taxa),]

myGM.tem <- myGM[myGM[,1] %in% colnames(myGD.tem),]
myGD.tem<- myGD.tem[,colnames(myGD.tem) %in% c("taxa",as.character(myGM[,1]))]

myGAPIT <- GAPIT(Y=myY ,
                 GD=myGD.tem, GM=myGM.tem,
                 PCA.total=3,
                 model= "CMLM",
                 SNP.MAF=0.05,
                 file.output=F
)

write.csv(myGAPIT$GWAS,"Soybean.PWAS_Seedling_SCN.CMLM.Results.csv",row.names = F)

cutOff.SCN <- 0.001 # a fixed lose cutoff was used for the quantative trait, SCN

png("Soybean.PWAS_Seedling_SCN.CMLM.png", width = 18,height=5,units = "in",res = 200)
layout(matrix(1:1, 1, 1, byrow = TRUE))
par(mar = (c(4.2,5,2,2)+ 0.5), mgp=c(3,1.8,0))    
par(bty="l", lwd=1.5)  

GAPIT.Manhattan.plot(myGAPIT$GWAS[,-1] ,name.of.trait= "SCN", cutOff= cutOff.SCN,cex.sig = 2,
                     highliht.sig=T,color1="lightgrey", color2="darkgrey", 
                     cex.lab=2,color1.sig="#FFA500", color2.sig="#FFA500")  

dev.off()

print(myGAPIT$GWAS[myGAPIT$GWAS$P.value <= 0.001 ,] )
#Should include the known SCN gene: 
#Glyma.07G195900.1.p (NSF07) P.value 6.85e-04
#SoyL02_18G021000.m1 (Rhg1/Glyma.18G022500) P.value 1.38e-04
