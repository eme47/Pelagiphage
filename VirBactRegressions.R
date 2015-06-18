#import file and check structure, get general summary
DD=read.table("VirBact_DeepDOMtransect2.csv",sep=",",header=TRUE)
str(DD)
summary(DD)

#plot log_Vir vs log_bact by water mass color, run regression and check assumptions
plot(VLA~BLA,data=DD,xlab="Viral Abundance", ylab="Bacterial Abundance")
lm=lm(VLA~BLA,data=DD)
summary(lm)

#Set "Surface" as referent dummy variable for Water Mass
RLWM=relevel(DD$WaterMass, ref="Surface") 

#test other variables in bivariate linear models
#for further investigation include reset level water mass (RLWM),temp, sal_avg and lat_start
summary(lm(VLA~RLWM,data=DD))
summary(lm(VLA~temp,data=DD))
summary(lm(VLA~sal_avg,data=DD))
summary(lm(VLA~O2_mL_L,data=DD))
summary(lm(VLA~lat_start,data=DD))
summary(lm(VLA~long_start,data=DD)) 

summary(lm(VLA~BLA+temp,data=DD)) 
summary(lm(VLA~BLA+RLWM,data=DD)) #important
summary(lm(VLA~BLA+sal_avg,data=DD)) 
summary(lm(VLA~BLA+lat_start,data=DD)) 
plot(VLA~RLWM, data=DD)

##Multiple linear regression for viral abundance
#Set "Surface" as referent dummy variable for Water Mass
RLWM=relevel(DD$WaterMass, ref="Surface") 

summary(lm) #from above
cor.test(DD$VLA,DD$BLA)

mlr=lm(VLA~BLA+RLWM, data=DD)
summary(mlr)

mlr2=lm(VLA~BLA+RLWM, data=DD)
summary(mlr2)

mlr3=lm(VLA~BLA+RLWM+temp, data=DD)
summary(mlr3)

mlr3B=lm(VLA~BLA+RLWM+sal_avg, data=DD)
summary(mlr3B)

mlr4=lm(VLA~BLA+RLWM+temp+sal_avg,data=DD)
summary(mlr4)

mlr5=lm(VLA~BLA+temp+sal_avg,data=DD)
summary(mlr5)

mlr5B=lm(VLA~BLA+temp+sal_avg+O2_mL_L,data=DD)
summary(mlr5B)

mlr6=lm(VLA~BLA+lat_start,data=DD)
summary(mlr6)

mlr7=lm(VLA~BLA+RLWM+lat_start, data=DD)
summary(mlr7)

summary(lm)

#compare models for most parsimonious and best fit using Akaike information criterion (AIC)
AIC(lm,mlr,mlr4)
AIC(lm, mlr,mlr2,mlr3,mlr4,mlr5,mlr6,mlr7)
AIC(lm,mlr,mlr3, mlr3B,mlr4)
AIC(mlr,mlr3B,mlr4,mlr5B)

resm4=rstandard(mlr4)
predm4=predict(mlr4)
hist(resm4)
qqnorm(resm4)
qqline(resm4)
plot(predm4,resm4)

str(DD)


#color plot

mycols=c("purple4","blue","darkorange1","yellow3","green3","red")
nonalphcols=c("red","darkorange1","yellow3", "green3", "blue", "purple4")
plot(DD$VLA~DD$BLA, xlab="log Bacterial Abundance", ylab="log Viral Abundance",xlim=c(3,7),ylim=c(4,9),pch=16,col=mycols[DD$WaterMass])
#legend("topleft", c("AABW","AAIW","DCM","Mesopelagic","NADW","Surface"), pch=16,col=mycols)
legend("topleft", c("Surface","DCM","Mesopelagic","AAIW","NADW","AABW"), pch=16,col=nonalphcols)
abline(lm, col="black")
#black and white plot
mycolsbw=c("black","gray20","gray55","gray45","gray10","gray70")
nonalphcolsbw=c("gray10","gray55","gray45","gray20","gray10","black")
tiff(file="VLABLA_BW.tif", width=700,height=600,type=c("windows"))
plot(DD$VLA~DD$BLA, xlab="log Bacterial Abundance", ylab="log Viral Abundance",xlim=c(3,7),ylim=c(4,9),pch=16,col=mycolsbw[DD$WaterMass])
#legend("topleft", c("AABW","AAIW","DCM","Mesopelagic","NADW","Surface"), pch=16,col=mycolsbw)
legend("topleft", c("Surface","DCM","Mesopelagic","AAIW","NADW","AABW"), pch=16,col=nonalphcolsbw)
abline(lm, col="black")
dev.off()
#different pch VLA~BLA
mypch=c(0,3,16,15,17,8)
nonalphcolsbw=c("gray10","gray55","gray45","gray20","gray10","black")
tiff(file="VLABLA_pch.tif", width=700,height=600,type=c("windows"))
plot(DD$VLA~DD$BLA, xlab="log Bacterial Abundance", ylab="log Viral Abundance",xlim=c(3,7),ylim=c(4,9),pch=mypch[DD$WaterMass],col=nonalphcolsbw[DD$WaterMass])
legend("topleft", c("Surface","DCM","Mesopelagic","AAIW","NADW","AABW"),pch=c(8,16,15,3,17,0),col=c("black","gray45","gray20","gray55","gray10","gray10"))
abline(lm, col="black")
dev.off()

res=rstandard(lm)
pred=predict(lm)
hist(res)
qqnorm(res)
qqline(res)
plot(pred,res)
summary(lm)
#assess temp with relation to Vir_abun, Bact_abund and VBR
plot(DD$temp,DD$VLA)
lm2=lm(temp~VLA,data=DD)
res2=rstandard(lm2)
pred2=predict(lm2)
hist(res2)
qqnorm(res2)
qqline(res2)
plot(pred2,res2)
summary(lm2)

plot(DD$temp,DD$BLA)
lm3=lm(temp~BLA,data=DD)
res3=rstandard(lm3)
pred3=predict(lm3)
hist(res3)
qqnorm(res3)
qqline(res3)
plot(pred3,res3)
summary(lm3)

plot(DD$temp,DD$VBR)
lm4=lm(temp~VBR,data=DD)
res4=rstandard(lm4)
pred4=predict(lm4)
hist(res4)
qqnorm(res4)
qqline(res4)
plot(pred4,res4)
summary(lm4)

#one way ANOVA/Kruskall depth by Vir_abundance, Bact_abundance or VBR
plot(DD$WaterMass,DD$logVir)
KW1=kruskal.test(WaterMass~VLA, data=DD)
KW1
KW2=kruskal.test(WaterMass~BLA, data=DD)
KW2
KW3=kruskal.test(WaterMass~VBR, data=DD)
KW3


plot(DD$WaterMass,DD$VLA,xlab="Water Mass",ylab="log Viral Abundance")
plot(DD$WaterMass,DD$BLA, xlab="Water Mass",ylab="log Bacterial Abundance")
plot(DD$WaterMass,DD$VBR, xlab="Water Mass",ylab="Virus:Bacteria Ratio")

bartlett.test(VLA~RLWM, data=DD)
bartlett.test(BLA~RLWM, data=DD)
bartlett.test(VBR~RLWM, data=DD)

library(lsmeans)

A1=lm(VLA~RLWM, data=DD)
anova(A1)
summary(A1)
lsmeans(A1,~RLWM,data=DD)
lsmeans(A1,pairwise~RLWM,data=DD)
lsmeans(A1,pairwise~RLWM,data=DD, adjust="bonferroni")
AA=lm(Vir_abundance~RLWM,data=DD)
anova(AA)
summary(AA)
lsmeans(AA,~RLWM,data=DD)
lsmeans(AA,pairwise~RLWM,data=DD)
lsmeans(AA,pairwise~RLWM,data=DD, adjust="bonferroni")


A2=lm(BLA~RLWM, data=DD)
anova(A2)
summary(A2)
lsmeans(A2,~RLWM,data=DD)
lsmeans(A2,pairwise~RLWM,data=DD)
lsmeans(A2,pairwise~RLWM,data=DD, adjust="bonferroni")

plot(VBR~RLWM, data=DD)
A3=lm(VBR~RLWM, data=DD)
anova(A3)
summary(A3)
lsmeans(A3,~RLWM,data=DD)
lsmeans(A3,pairwise~RLWM,data=DD)
lsmeans(A3,pairwise~RLWM,data=DD, adjust="bonferroni")

summary(DD)
str(DD)


#plots of VLA, BLA and VBR by lattitude, need to fix legends
l=VPNE$lat_start
m=VPNE$VLA
n=VPNE$BLA
o=VPNE$VBR

SDVLA=VPNE$Vir_SE
SDBLA=VPNE$Bact_SE

VLA_lat=plot(VLA~lat_start,data=DD, xlab="Lattitude", ylab="log Viral Abundance",col=mycols[DD$WaterMass],pch=16)
legend("bottomleft", c("AABW","AAIW","DCM","Mesopelagic","NADW","Surface"), pch=16,col=mycols)
BLA_lat=plot(BLA~lat_start,data=DD, ylim=c(3,7), xlab="Lattitude", ylab="log Bacterial Abundance",col=mycols[DD$WaterMass],pch=16)
legend("bottomleft", c("AABW","AAIW","DCM","Mesopelagic","NADW","Surface"), pch=16,col=mycols)
VBR_lat=plot(VBR~lat_start,data=DD, xlab="Lattitude", ylab="Virus to Bacteria Ratio",col=mycols[DD$WaterMass],pch=16, ylim=c(0,80))
legend("topright", c("AABW","AAIW","DCM","Mesopelagic","NADW","Surface"), pch=16,col=mycols)

#temp and sal plots by lattitude
plot(DD$lat_start,DD$temp,col=mycols[DD$WaterMass],pch=16)
legend("bottomleft", c("AABW","AAIW","DCM","Mesopelagic","NADW","Surface"), pch=16,col=mycols)
str(DD)
plot(DD$lat_start,DD$sal_avg,col=mycols[DD$WaterMass],pch=16)
legend("bottomleft", c("AABW","AAIW","DCM","Mesopelagic","NADW","Surface"), pch=16,col=mycols)
