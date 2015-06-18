PC=read.table("VirBact_DeepDOMtransect2.csv",sep=",",header=TRUE)
str(PC)
summary(PC)

plot(temp~WaterMass, data=PC)
bartlett.test(temp~WaterMass, data=PC)
lsmeans
T1=lm(temp~RLWM, data=PC)
anova(T1)
summary(T1)
lsmeans(T1,~RLWM,data=PC)
lsmeans(T1,pairwise~RLWM,data=PC)
lsmeans(T1,pairwise~RLWM,data=PC, adjust="bonferroni")

kwTW=kruskal.test(PC$temp~PC$WaterMass)
kwTW
plot(sal_avg~WaterMass, data=PC)
KWSW=kruskal.test(PC$sal_avg~PC$WaterMass)
KWSW
plot(O2_mL_L~WaterMass,data=PC)
KWOW=kruskal.test(PC$O2_mL_L~PC$WaterMass)
KWOW
bartlett.test(O2_mL_L~WaterMass, data=PC)
plot(temp~lat_start,data=PC,pch=16, col=mycols[DD$WaterMass])
abline(lm(temp~lat_start,data=PC, DD$WaterMass=="Surface"))
abline(lm(temp~lat_start,data=PC, DD$WaterMass=="DCM"))
abline(lm(temp~lat_start,data=PC, DD$WaterMass=="Mesopelagic"))
abline(lm(temp~lat_start,data=PC, DD$WaterMass=="AAIW"))
abline(lm(temp~lat_start,data=PC, DD$WaterMass=="NADW"))
abline(lm(temp~lat_start,data=PC, DD$WaterMass=="AABW"))

plot(VLA~lat_start,data=PC,pch=16, col=mycols[DD$WaterMass])
plot(sal_avg~lat_start,data=PC,pch=16, col=mycols[DD$WaterMass])
abline(lm(sal_avg~lat_start,data=PC, DD$WaterMass=="Surface"))
abline(lm(sal_avg~lat_start,data=PC, DD$WaterMass=="DCM"))
abline(lm(sal_avg~lat_start,data=PC, DD$WaterMass=="Mesopelagic"))
abline(lm(sal_avg~lat_start,data=PC, DD$WaterMass=="AAIW"))
abline(lm(sal_avg~lat_start,data=PC, DD$WaterMass=="NADW"))
abline(lm(sal_avg~lat_start,data=PC, DD$WaterMass=="AABW"))
