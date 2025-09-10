#fine-gray regression model
library(prodlim)
library(survival)
library(cmprsk)
library (riskRegression) 
library(lava)


CRE3 <- read.csv("CRE3.csv")
View(CRE3)

fgr1<-FGR(Hist(time, status)~ atg + reint + mv + arf + kpc + crepre_60 + multipre + crepost_60+ multipost, data=CRE3, cause=1)

pred30<-predictRisk(fgr1, newdata=CRE3, time=30, cause=1)

#Draw histograms of predicted probabilities at day 30 and day 60

a <- ggplot(CRE3, aes(x = pred30))
a + geom_histogram(bins = 50, color = "black", fill = "grey") + geom_vline(aes(xintercept = mean(pred30)), color = "#FC4E07", linetype = "dashed", size = 1) + theme_minimal()
a + geom_violin(trim = FALSE)+ geom_jitter(position=position_jitter(0.2))

pred60<-predictRisk(fgr1, newdata=CRE3, time=60, cause=1)

a2 <- ggplot(CRE3, aes(x = pred60))
a2 + geom_histogram(bins = 50, color = "black", fill = "grey") + geom_vline(aes(xintercept = mean(pred60)), color = "#FC4E07", linetype = "dashed", size = 1) + theme_minimal()

#draw calibration curves

library(survival)

score1<-Score(list("Fine-Gray"=fgr1),formula=Hist(time,status)~1,data=CRE3, times=seq(60, 90, 180), plots="Calibration", summary="risks")
score2<-Score(list("Fine-Gray"=fgr1),formula=Hist(time,status)~1,data=CRE3, times=seq(60, 90, 180), plots="AUC", summary="risks")
score3<-Score(list("Fine-Gray"=fgr1),formula=Hist(time,status)~1,data=CRE3, plots="box", times=seq (60, 90, 180), null.model=FALSE, summary="risks")
dev.new(width=5,height=4)
dev.new(width=5,height=4)

plotCalibration(score1,times=60, pseudo = TRUE)

plotPredictRisk(fgr1, time=1:60, cause=1)
dev.new(width=5,height=4)
boxplot(score3, time=30)

#Model cross-validation
fgr1.cv<-FGR(Hist(time, status)~ atg + reint + mv + arf + mv + kpc + crepre_60 + multipre + crepost_60+ multipost, data=CRE3, cause=1)
score.cv2<-Score(list("Fine-Gray CV"=[fgr1.cv],(http://fgr1.cv/)),formula=Hist(time,status)~1,data=CRE3, times=seq(60, 90, 180), split.method="bootcv", B=100, plots="calibration")
dev.new(width=5,height=4)
plotCalibration(score.cv2,times=60, pseudo = TRUE)
#plot predicted risk
pred.fgr1<-predict(fgr1, newdata = CRE3[267,],time=1:60, cause=1)
dev.new(width=5,height=4)
plotPredictRisk(fgr1, newdata = CRE3[22,],time=1:60, cause=1)
#construct nomogram
library(survival)

# treat time to CRE infection and death and competing risk

etime <-with(CRE3, ifelse(status==0, futime, cretime))
event<-with(CRE3, ifelse(status==0, 2*death, 1))
event<-factor(event, 0:2, labels=c("censor", "pcm", "death"))
#FG model-Weighed multistate survival dataset
adata<-finegray(Surv(etime, event)~atg + reint + mv + arf + mv + kpc + crepre_60 + multipre + crepost_60+ multipost, data=CRE3, na.action=na.pass)
fgfit<-coxph(Surv(fgstart, fgstop, fgstatus) ~ atg + reint + mv + arf + mv + kpc + crepre_60 + multipre + crepost_60+ multipost, weight=fgwt, data=adata)
f.csc <-coxph(Surv(cretime,status==1) ~ atg + reint + arf + mv + kpc + crepre_60 + multipre + crepost_60+ multipost, data=CRE3)

# Draw nomogram for 60-day risk of CRE infection

library(regplot)
f.csc <-coxph(Surv(time,status==1) ~ atg + reint + arf + mv + kpc + crepre_60 + multipre + crepost_60+ multipost, data=CRE3)
dev.new(width=5,height=4)
regplot(f.csc,observation=CRE3[CRE3$id==22,],failtime=c(30, 60), center=TRUE, points=TRUE, rank = "range", subticks = TRUE, interval="confidence", title = "Risk of CRE Infection Post OLT", prfail=TRUE, droplines=T)
dev.new(width=5,height=4)
plotPredictRisk(fgr1, newdata = CRE3[267,],time=1:60, cause=1)
#weighted dataset for competing risk analysis to build nomogram
f.csc <-coxph(Surv(time,status==1) ~ reint + mv + arf + crepre_60 + crepost_60+ multipost, data=CRE3)
library(DynNom)
library(shiny)
DynNom(f.csc, data=CRE3, m.summary = c("formatted"), covariate = "slider", DNtitle = "Move slider to (0) when risk factor is absent; Move slider to (1) when risk factor is present",DNxlab = "Predicted CRE infection probability +/-95% CI", DNlimits = c(0, 180), DNylab = "Group",KMtitle="CRE infection probability", KMxlab="Days post transplant", KMylab="Probability of CRE infection", ptype=c("1-st"))
DNbuilder(f.csc, data=CRE3, m.summary = c("formatted"), covariate = "slider", DNtitle = "Move slider to (0) when risk factor is absent; Move slider to (1) when risk factor is present",DNxlab = "Predicted CRE infection probability +/-95% CI", DNlimits = c(0, 180), DNylab = "Group",KMtitle="CRE infection probability", KMxlab="Days post transplant", KMylab="Probability of CRE infection", ptype=c("1-st"))
library(rsconnect)
rsconnect::setAccountInfo(name='idbologna',token='DA368A4CBFC7435F67911388C6A3E9B4',secret='N4GLNTKq0RMD4HqdByj+CB6MD05VL1XqFnkRx8cp')
runApp("DynNomapp")
deployApp()
