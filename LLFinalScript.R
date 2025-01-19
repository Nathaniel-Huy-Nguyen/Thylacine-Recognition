library(patchwork)
install.packages("generalhoslem")
library(generalhoslem)
install.packages("sure")
library(sure)
install.packages("glmm")
library(glmm)
install.packages("pscl")
library(pscl)

LL<-read.csv("LLTime3.csv")
LLW<-LL[which(LL$Species == "Wallaby"),]
LLW
nrow(LLW)/nrow(LL)

LL1<-LL[1:917, ]
LL2<-LL[918:2187, ]

aggregate(LLW$Model, by=list(LLW$Model), FUN=length)
aggregate(LL$Species, by=list(LL$Species), FUN=length)

mean(LLWT$Investigating)

pairwise.t.test((LLW$Investigating), LLW$Model, NA.rm=TRUE)
sd(LLW$Investigating)
mean(LLW$Grazing)
sd(LLW$Grazing)

mean(LLW$T_loc)
sd(LLW$T_loc)


LLPie<-c(941, 488, 350, 141, 8, 58, 3, 3, 88, 10, 5, 92)
names(LLPie)<-c("Bennett's Wallaby","Tasmanian Pademelon","Common Brushtail Possum Possum","Common Wombat", "Short-Beaked Echidna", "Tasmanian Devil", "Spotted-Tail Quoll", "Eastern Spotted Quoll", "Feral Rabbit", "Feral Sambar Deer", "Feral Cat", "Unknown")
LLSpeciesPlot<-barchart(LLPie, col=c("palegreen", "palegreen3", "paleturquoise", "paleturquoise3", "paleturquoise4", "lavender", "plum1", "plum3", "papayawhip", 'peachpuff', "palevioletred", "white"), ylab="Species", xlab="Number of Observations", main="A) London Lakes")
names(LLPie)<-c("Bennett's Wallaby","Tasmanian Pademelon", "Common Brushtail Possum", "Other")
pie(LLPie, col=c("palegreen", "palegreen3", "paleturquoise", "paleturquoise3", "paleturquoise4", "lavender", "plum1", "plum3", "papayawhip", 'peachpuff', "palevioletred", "white"))

table(LLW$Model, LLW$Pref_Dist)
chisq.test(table(LLW$Model, LLW$Pref_Dist))

table(LLW$Model, LLW$Site)
chisq.test(table(LLW$Model, LLW$Site))

#Poisson GLM for total visits
LLn_by_model<-aggregate(LLW$Model, by=list(Model=LLW$Model, Site=LLW$Site, Moon_Phase=factor(LLW$Moon_Phase)), FUN=length, drop=FALSE)
LLn_by_model
LLn_by_model$x <- ifelse(is.na(LLn_by_model$x),0,LLn_by_model$x)

LLn_by_model$Moon_Phase <- as.numeric(as.character(LLn_by_model$Moon_Phase))
str(LLn_by_model$Moon_Phase)
LLVisitModel<-glmmTMB(x ~ Model + (1|Site) + Moon_Phase, data = LLn_by_model, family = "nbinom1")
check_model(LLVisitModel)
summary(LLVisitModel)
emmeans(LLVisitModel, ~ Model, type = "response")
pairs(emmeans(LLVisitModel, ~ Model, type = "response"))
plot(emmeans(LLVisitModel, ~ Model, type = "response"))

library(performance)
library(glmmTMB)
library(tidyverse)
install.packages("DHARMa")
library(DHARMa)

load('LLn_by_model.RData')

LLVisitModel <- glmmTMB(x ~ Model + (1|Site), data = LLn_by_model, family = "nbinom1")
## in hindsight, don't use this
## check_model(LLVisitModel)

## use DHARMa functions instead:
simres <- simulateResiduals(LLVisitModel)

## standard DHARMa diagnostic plots
plot(simulateResiduals(LLVisitModel))

## simulated residuals stratified by the predictor variable
plotResiduals(simres, form = LLn_by_model$Model)

## ## some tests, if you really want them (they are reported on plot(simres) anyway)
#testUniformity() - tests if the overall distribution conforms to expectations.
#testOutliers() - tests if there are more simulation outliers than expected.
#testDispersion() - tests if the simulated dispersion is equal to the observed dispersion.
#testQuantiles() - fits a quantile regression or residuals against a predictor (default predicted value),
##                   and tests of this conforms to the expected quantile.
#testCategorical(simulationOutput, catPred = testData$group) tests residuals against a categorical predictor.

LLn_by_model<-aggregate(LLW$Model, by=list(Model=LLW$Model, Site=LLW$Site), FUN=length, drop=FALSE)
LLn_by_model
LLn_by_model$x <- ifelse(is.na(LLn_by_model$x),0,LLn_by_model$x)
LLVisitModel<-glmmTMB(x ~ Model + (1|Site), data = LLn_by_model, family = "nbinom1")
check_model(LLVisitModel)
summary(LLVisitModel)
emmeans(LLVisitModel, ~ Model, type = "response")
pairs(emmeans(LLVisitModel, ~ Model, type = "response"))
LLTotalPlot<-plot(emmeans(LLVisitModel, ~ Model, type = "response"), ylab="Treatment Type", xlab="") + ggtitle("A) London Lakes")
LLTotalPlot
Anova(LLVisitModel)
summary(aov(x ~ Model + Site, data=LLn_by_model))


summary(LLn_by_model)

LLn_by_modelT<-LLn_by_model[which(LLn_by_model$Model=="Thylacine"), ]
mean(LLn_by_modelT$x)
var(LLn_by_modelT$x)
sqrt(var(LLn_by_modelT$x))
sd(LLn_by_modelT$x)




#Save this and send the line to Jeremy
save(LLn_by_model, file = 'LLn_by_model.RData')
LLVisitModel<-glmmTMB(x ~ Model + (1|Site), data = LLn_by_model, family = "nbinom1")

#Converting times to POSIX and creating Duration
LLTime <- mutate(LL,
                  ## create POSIXct objects for the first & last times
                  First_datetime = strptime(paste(Date,First),format = '%d/%m/%Y %H:%M:%S'),
                  Last_datetime = strptime(paste(Date,Last),format = '%d/%m/%Y %H:%M:%S'),
                  ## adjust if they span the date boundary
                  Last_datetime = Last_datetime + ifelse(Last_datetime >= First_datetime, 0, 60*60*24),
                  Duration = as.numeric(difftime(Last_datetime, First_datetime, units = 'secs')))

#GROUP SIZE STUFF
install.packages("ivs")
library(ivs)


Sites <- unique(LLTime$Site)
site <- Sites[1]
## create a new variable for the number of overlaps
LLTime$n_overlaps <- NA
for(site in Sites){
  irows <- which(LLTime$Site == site)
  ## create an "interval vector" (which contains the start/end times) for this site
  IV <- with(LLTime[irows,],
             ## if the start-time and end-time are identical, add 1 second to the end-time
             iv(First_datetime, Last_datetime + ifelse(First_datetime == Last_datetime, 1, 0)))
  ## Number of overlapping overservations *for this site*
  ##
  ## The '-1' excludes overlap with "self", so it is the count of
  ## the number of other observations that overlap with it
  LLTime$n_overlaps[irows] <- iv_count_overlaps(IV,IV) - 1
  ## to look at what the overlaps are (if you want to check):
  ## subset(iv_locate_overlaps(IV,IV), needles != haystack)
  ## the "subset( ..., needles != haystack)" code filters out self-overlaps
}
LLTimeW<-LLTime[which(LLTime$Species == "Wallaby"),]

mean(LLTimeW$n_overlaps)
pairwise.t.test((LLTimeW$n_overlaps), LLTimeW$Model, na.rm=TRUE)

LLOverlapModel<-glmmTMB(n_overlaps ~ Model + (1|Site) + Position + Moon_Phase, ziformula =~1, data=LLTimeW, family="poisson")
summary(LLOverlapModel)
plot(LLOverlapModel)
anova(LLOverlapModel)
plot(emmeans(LLOverlapModel, ~Model, typw="response"))
emmeans(LLOverlapModel, ~Model, type="response")
pairs(emmeans(LLOverlapModel, ~Model))
plot(simulateResiduals(LLOverlapModel))
check_model(LLOverlapModel)


#Zero inflated model for overlaps
LLOverZeroInf<-glmmTMB(n_overlaps ~ Model + (1|Site) + Moon_Phase, data = LLTimeW, ziformula = ~., family="nbinom1")

summary(LLOverZeroInf)
check_model(LLOverZeroInf, n_bins=10)
emmeans(LLOverZeroInf, ~Model, type = "response")
emmeans(LLOverZeroInf, ~Model, type = "conditional")
LLOverPlot<-plot(emmeans(LLOverZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="") + ggtitle("A) London Lakes")
LLOverPlot
pairs(emmeans(LLOverZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(LLOverZeroInf, ~Model, type = "response")))
plot(simulateResiduals(LLOverZeroInf))

emmip(LLOverZeroInf, ~Model, type = "response", CIs=TRUE)

LLDurationModel<-lmer(Duration ~ Model + (1|Site) + Position + Moon_Phase, data=LLTimeW)
check_model(LLDurationModel)
summary(LLDurationModel)
anova(LLDurationModel)
LLDurationPlot<-plot(emmeans(LLDurationModel, ~Model), ylab="Treatment Type", xlab="") + ggtitle("A) London Lakes")
emmeans(LLDurationModel, ~Model)
pairs(emmeans(LLDurationModel, ~Model))

#Zero inflated weighted binomial for grazing
LLGrazeZeroInf<-glmmTMB(Grazing ~ Model + (1|Site) + Position + Moon_Phase+n_overlaps, data = LLTimeW, ziformula = ~., weights = Total, family = binomial)
summary(LLGrazeZeroInf)
Anova(LLGrazeZeroInf)
check_model(LLGrazeZeroInf, n_bins=10)
emmeans(LLGrazeZeroInf, ~Model, type = "response")
emmeans(LLGrazeZeroInf, ~Model, type = "conditional")
LLGrazePlot<-plot(emmeans(LLGrazeZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="") + ggtitle("A) London Lakes")
pairs(emmeans(LLGrazeZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(LLGrazeZeroInf, ~Model, type = "response")))
plot(simulateResiduals(LLGrazeZeroInf))

LLGrazeWBinomial<-glmmTMB(Grazing ~ Model + (1|Site) + Position + Moon_Phase+n_overlaps, data = LLTimeW, weights = Total, family = binomial)
summary(LLGrazeWBinomial)
check_model(LLGrazeWBinomial, n_bins=10)
emmeans(LLGrazeWBinomial, ~Model, type = "response")

LLGrazeHurdle<-hurdle(Grazing ~ Model + (1|Site) + Position + Moon_Phase+n_overlaps, data = LLTimeW, weights = Total, family = binomial)


#Zero-Inflated weighted binomial for investigating
LLInvZeroInf<-glmmTMB(Investigating ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = LLTimeW, ziformula = ~(1|Site) + Position + Moon_Phase+n_overlaps, weights = Total, family = binomial)
summary(LLInvZeroInf)
check_model(LLInvZeroInf, n_bins=10)
emmeans(LLInvZeroInf, ~Model, type = "response")
emmeans(LLInvZeroInf, ~Model, type = "conditional")
LLInvPlot<-plot(emmeans(LLInvZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="") + ggtitle("A) London Lakes")
pairs(emmeans(LLInvZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(LLInvZeroInf, ~Model, type = "response")))
Anova(LLInvZeroInf)

LLInvWBinomial<-glmmTMB(Investigating ~ Model + (1|Site) + Position + Moon_Phase+n_overlaps, data = LLTimeW, weights = Total, family = binomial)
summary(LLInvWBinomial)
check_model(LLInvWBinomial, n_bins=10)
emmeans(LLInvWBinomial, ~Model, type = "response")
pairs(emmeans(LLInvWBinomial, ~Model, type = "response"))

#Zero-Inflated weighted binomial for total vigilance
LLTimeW$T_vig<-(LLTimeW$V_low+LLTimeW$V_mid+LLTimeW$V_high)
LLVigZeroInf<-glmmTMB(T_vig ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = LLTimeW, ziformula = ~., weights = Total, family = binomial)
summary(LLVigZeroInf)
check_model(LLVigZeroInf, n_bins=10)
emmeans(LLVigZeroInf, ~Model, type = "response")
emmeans(LLVigZeroInf, ~Model, type = "conditional")
LLVigPlot<-plot(emmeans(LLVigZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="") + ggtitle("A) London Lakes")
pairs(emmeans(LLVigZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(LLVigZeroInf, ~Model, type = "response")))
Anova(LLVigZeroInf)

LLTimeW$V_direction<-((LLTimeW$V_toward+1)/(LLTimeW$V_away+1))
LLVigDir<-glmmTMB(V_direction ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = LLTimeW)
summary(LLVigDir)
check_model(LLVigDir, n_bins=10)
emmeans(LLVigDir, ~Model, type="response")
pairs(emmeans(LLVigDir, ~Model, type="response"))

#Zero-Inflated weighted binomial for total locomotion
LLTimeW$T_loc<-(LLTimeW$L_towards+LLTimeW$L_away+LLTimeW$L_across)
LLLocZeroInf<-glmmTMB(T_loc ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = LLTimeW, ziformula = ~., weights = Total, family = binomial)
summary(LLLocZeroInf)
check_model(LLLocZeroInf, n_bins=10)
emmeans(LLLocZeroInf, ~Model, type = "response")
emmeans(LLLocZeroInf, ~Model, type = "conditional")
LLLocPlot<-plot(emmeans(LLLocZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="") + ggtitle("A) London Lakes")
pairs(emmeans(LLLocZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(LLLocZeroInf, ~Model, type = "response")))
Anova(LLLocZeroInf)
LLLocPlot


(LLWT$Grazing)
(LLWT$Grazing)
mean(LLWF$Grazing)
mean(LLWG$Grazing)
mean(LLWC$Grazing)
LLGrazePlot

#Cutting data back into photos and applying vigilance into an ordered factor
LLTime_by_frame <-
  LLTime %>%
  ## head %>%
  ## calculate numbers of frames by vigileace roup
  mutate(n_V_low = round(Total*V_low),
         n_V_mid = round(Total*V_mid),
         n_V_high = round(Total*V_high)) %>%
  group_by(Individual) %>%
  reframe(Site = Site, Camera = Camera, Position = Position, Date = Date, Model = Model, 
          Weather = Weather, Min_Temp = Min_Temp, Wind_Speed = Wind_Speed, 
          Rainfall = Rainfall, Moon_Phase = Moon_Phase, ## Individual = Individual, 
          Species = Species, Size = Size, Sex = Sex, Joey = Joey, n_overlaps = n_overlaps,
          ## for the three vigilence categories, represent these as 0's or 1's 
          V_low = rep(0:1,c(Total-n_V_low,n_V_low)),
          V_mid = c(rep(0:1,c((Total-n_V_low) - n_V_mid,n_V_mid)),rep(0,n_V_low)),
          V_high = c(rep(0:1,c((Total-n_V_low-n_V_mid) - n_V_high,n_V_high)),rep(0,n_V_low + n_V_mid))
  ) %>%
  ## create a new "vigilance" variable, which is an ordered factor
  mutate(vigilance_level = case_when(
    V_low == 1 ~ "low",
    V_mid == 1 ~ "mid",
    V_high == 1 ~ "high",
    .default = NA),
    vigilance_level = factor(vigilance_level, levels = c('low','mid','high'), ordered = TRUE))

#Create an ordinal logistical model for vigilance level
LLTimeW_by_frame<-LLTime_by_frame[which(LLTime_by_frame$Species == "Wallaby"),]
LLVigModelW<-clmm(vigilance_level~Model+(1|Site)+Moon_Phase+n_overlaps, data=LLTimeW_by_frame)
LLVigModelW
summary(LLVigModelW)
check_model(LLVigModelW)
Anova(LLVigModelW)

#Outputting modeled effect of statue type on vigilance with P-values
emmeans(LLVigModelW,~Model, type='response')
summary(emmeans(LLVigModelW, ~Model), type = 'response')
pairs(emmeans(LLVigModelW, ~Model))
LLVigLevelPlot<-plot(emmeans(LLVigModelW, ~Model, type="response", mode="exc.prob"), xlab="", ylab="Treatment Type") + ggtitle("A) London Lakes")
plot(pairs(emmeans(LLVigModelW, ~Model)), comparisons = TRUE)

LLVigLevelFullPlot<-plot(emmeans(LLVigModelW,~Model | vigilance_level, mode="prob", type='response'), horizontal=TRUE, xlab="", ylab="Treatment Type") + ggtitle("A) London Lakes")

#Defining preferred distance as a vector
LLTime$Pref_Dist<-factor(LLTime$Pref_Dist, levels=c("<5m", "5-10m", "10-15m", ">15m"), ordered=TRUE)
LLTime$Pref_Dist
LLTimeW<-LLTime[which(LLTime$Species=="Wallaby"),]


#Creating a model for preferred distance
LLDistModelW<-clmm(Pref_Dist~Model+(1|Site)+Moon_Phase+n_overlaps, data=LLTimeW)
LLDistModelW
emmeans(LLDistModelW, ~Model)
summary(LLDistModelW)
Anova(LLDistModelW)
check_model(LLDistModelW)
model_performance(LLDistModelW)

#Outputting and displaying modeled effect of statue on preferred distance
emmeans(LLDistModelW,~Model, mode="exc.prob", type='response')
summary(emmeans(LLDistModelW, ~Model, mode="exc.prob", type = 'response'))
pairs(emmeans(LLDistModelW, ~Model, mode="exc.prob", type = 'response'))
plot(emmeans(LLDistModelW, ~Model, mode="exc.prob", type = 'response'))
plot(pairs(emmeans(LLDistModelW, ~Model, mode="exc.prob", type = 'response')))
plot(pairs(emmeans(LLDistModelW, ~Model, mode="exc.prob", type = 'response')), comparisons = TRUE)

pairs(emmeans(LLDistModelW,~Model | Pref_Dist, mode="prob", type='response'))
plot(emmeans(LLDistModelW,~Model | Pref_Dist, mode="prob", type='response'), horizontal=TRUE, xlab="EMM of probability of selection", ylab="Treatment Type", main="A) London Lakes")
LLFullDistPlot<-plot(emmeans(LLDistModelW,~Model | Pref_Dist, mode="prob", type='response'), horizontal=TRUE, xlab="", ylab="Treatment Type", main="A) London Lakes") + ggtitle("A) London Lakes")
LLDistPlot<-plot(emmeans(LLDistModelW, ~Model, mode="exc.prob", type = 'response'), ylab="Treatment Type", xlab="") + ggtitle("A) London Lakes")
plot(emmeans(LLDistModelW, ~Model, mode="exc.prob", type = 'response'), ylab="Treatment Type", xlab="") + ggtitle("A) London Lakes")

Anova(DunLocZeroInf)


pulkrob.chisq(LLDistModelW, c("Model", "Site"))
lipsitz.test(LLDistModelW,g=10)
autoplot.clm(LLDistModelW, what = c("qq", "fitted", "covariate"))

hist(table(LLW$Model, LLW$Pref_Dist))

LLDistBarplot<-LLTimeW %>%
  drop_na(Pref_Dist) %>%
  ggplot(aes(x=Pref_Dist, fill=Model))+
  geom_bar(position = position_dodge())+
  labs(x = "Count",
       y = "Treatment")
LLDistBarplot


emmeans(LLDistModelW,~Position, type='response')
summary(emmeans(LLDistModelW, ~Position), type = 'response')
pairs(emmeans(LLDistModelW, ~Position))
plot(pairs(emmeans(LLDistModelW, ~Position)))
plot(pairs(emmeans(LLDistModelW, ~Position)), comparisons = TRUE)

#NEW idea for binary regression of distance
LLBinModelW<-glmmTMB(Pref_Dist2 ~Model+(1|Site)+Moon_Phase+Position, family="binomial", data = LLTimeW)
LLBinModelW
summary(LLBinModelW)
emmeans(LLBinModelW, ~Model, type="response")
pairs(emmeans(LLBinModelW, ~Model, type="response"))
check_model(LLBinModelW, n_bins=10)
model_performance(LLBinModelW)
Anova(LLBinModelW)
plot(emmeans(LLBinModelW, ~Model, type="response"))

