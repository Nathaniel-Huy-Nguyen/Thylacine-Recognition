Q<-read.csv("QTime3.csv")
QW<-Q[which(Q$Species == "Wallaby"),]
QW
nrow(QW)/nrow(Q)

Q1<-Q[1:2265, ]
Q2<-Q[2266:3600, ]

mean(QTimeW$Duration)
mean(QTimeW$Total)

aggregate(QW$Model, by=list(QW$Model), FUN=length)
aggregate(Q$Species, by=list(Q$Species), FUN=length)

QPie<-c(2526, 698, 164, 74, 2, 16, 2, 14, 26, 78)
names(QPie)<-c("Bennett's Wallaby","Tasmanian Pademelon", "Eastern Grey Kangaroo", "Common Brushtail Possum Possum","Common Wombat", "Tasmanian Devil", "Spotted-Tail Quoll", "Feral Rabbit", "Feral Sambar Deer", "Unknown")
QSpeciesPlot<-barchart(QPie, col=c("palegreen", "palegreen3", "palegreen4", "paleturquoise", "paleturquoise3", "lavender", "plum2", "papayawhip", 'peachpuff', "white"), ylab="Species", xlab="Number of Observations", main="B) The Quoin")
names(QPie)<-c("Bennett's Wallaby","Tasmanian Pademelon", "Eastern Grey Kangaroo")
pie(QPie, col=c("palegreen", "palegreen3", "palegreen4", "paleturquoise", "paleturquoise3", "lavender", "plum2", "papayawhip", 'peachpuff', "white"))


table(QW$Model, QW$Pref_Dist)
chisq.test(table(QW$Model, QW$Pref_Dist))

#Poisson GLM for total visits
Qn_by_model<-aggregate(QW$Model, by=list(Model=QW$Model, Site=QW$Site), FUN=length)
QVisitModel<-glmer(x ~ Model + (1|Site), data = Qn_by_model, family = "poisson")
check_model(QVisitModel)
summary(QVisitModel)
emmeans(QVisitModel, ~ Model, type = "response")
pairs(emmeans(QVisitModel, ~ Model, type = "response"))
plot(emmeans(QVisitModel, ~ Model, type = "response"))

## use DHARMa functions instead:
Qsimres <- simulateResiduals(QVisitModel)

## standard DHARMa diagnostic plots
plot(Qsimres)

Qn_by_model<-aggregate(QW$Model, by=list(Model=QW$Model, Site=QW$Site), FUN=length, drop=FALSE)
Qn_by_model
Qn_by_model$x <- ifelse(is.na(Qn_by_model$x),0,Qn_by_model$x)
QVisitModel<-glmmTMB(x ~ Model + (1|Site), data = Qn_by_model, family = "nbinom1")
check_model(QVisitModel)
summary(QVisitModel)
emmeans(QVisitModel, ~ Model, type = "response")
pairs(emmeans(QVisitModel, ~ Model, type = "response"))
QTotalPlot<-plot(emmeans(QVisitModel, ~ Model, type = "response"), comparisons=TRUE, ylab="Treatment Type", xlab="") + ggtitle("B) The Quoin")
QTotalPlot
Anova(QVisitModel)

#Converting times to POSIX and creating Duration
QTime <- mutate(Q,
                  ## create POSIXct objects for the first & last times
                  First_datetime = strptime(paste(Date,First),format = '%d/%m/%Y %H:%M:%S'),
                  Last_datetime = strptime(paste(Date,Last),format = '%d/%m/%Y %H:%M:%S'),
                  ## adjust if they span the date boundary
                  Last_datetime = Last_datetime + ifelse(Last_datetime >= First_datetime, 0, 60*60*24),
                  Duration = as.numeric(difftime(Last_datetime, First_datetime, units = 'secs')))

#GROUP SIZE STUFF
install.packages("ivs")
library(ivs)


Sites <- unique(QTime$Site)
site <- Sites[1]
## create a new variable for the number of overlaps
QTime$n_overlaps <- NA
for(site in Sites){
  irows <- which(QTime$Site == site)
  ## create an "interval vector" (which contains the start/end times) for this site
  IV <- with(QTime[irows,],
             ## if the start-time and end-time are identical, add 1 second to the end-time
             iv(First_datetime, Last_datetime + ifelse(First_datetime == Last_datetime, 1, 0)))
  ## Number of overlapping overservations *for this site*
  ##
  ## The '-1' excludes overlap with "self", so it is the count of
  ## the number of other observations that overlap with it
  QTime$n_overlaps[irows] <- iv_count_overlaps(IV,IV) - 1
  ## to look at what the overlaps are (if you want to check):
  ## subset(iv_locate_overlaps(IV,IV), needles != haystack)
  ## the "subset( ..., needles != haystack)" code filters out self-overlaps
}
QTimeW<-QTime[which(QTime$Species == "Wallaby"),]

mean(QTimeW$n_overlaps)
pairwise.t.test((QTimeW$n_overlaps), QTimeW$Model, na.rm=TRUE)

QOverlapModel<-glm(n_overlaps ~ Model + Site + Position + Moon_Phase + Min_Temp + Rainfall + Wind_Speed, data=QTimeW)
summary(QOverlapModel)
anova(QOverlapModel)
plot(emmeans(QOverlapModel, ~Model))
emmeans(QOverlapModel, ~Model)
pairs(emmeans(QOverlapModel, ~Model))

QOverZeroInf<-glmmTMB(n_overlaps ~ Model + (1|Site) + Moon_Phase, data = QTimeW, ziformula = ~., family=poisson)
QOverZeroInf
summary(QOverZeroInf)
check_model(QOverZeroInf, n_bins=10)
emmeans(QOverZeroInf, ~Model, type = "response")
emmeans(QOverZeroInf, ~Model, type = "conditional")
QOverPlot<-plot(emmeans(QOverZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="")  + ggtitle("B) The Quoin")
pairs(emmeans(QOverZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(QOverZeroInf, ~Model, type = "response")))

QDurationModel<-lmer(Duration ~ Model + (1|Site) + Position + Moon_Phase, data=QTimeW)
summary(QDurationModel)
QDurationModel
check_model(QDurationModel)
anova(QDurationModel)
QDurationPlot<-plot(emmeans(QDurationModel, ~Model), ylab="Treatment Type", xlab="")  + ggtitle("B) The Quoin")
emmeans(QDurationModel, ~Model)
pairs(emmeans(QDurationModel, ~Model))

#Zero inflated weighted binomial for grazing
QGrazeZeroInf<-glmmTMB(Grazing ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = QTimeW, ziformula = ~., weights = Total, family = binomial)
summary(QGrazeZeroInf)
check_model(QGrazeZeroInf, n_bins=10)
emmeans(QGrazeZeroInf, ~Model, type = "response")
emmeans(QGrazeZeroInf, ~Model, type = "conditional")
QGrazePlot<-plot(emmeans(QGrazeZeroInf, ~Model, type = "response"), comparisons=TRUE, ylab="Treatment Type", xlab="")  + ggtitle("B) The Quoin")
QGrazePlot
pairs(emmeans(QGrazeZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(QGrazeZeroInf, ~Model, type = "response")))
Anova(QGrazeZeroInf)

#Zero-Inflated weighted binomial for investigating
QInvZeroInf<-glmmTMB(Investigating ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = QTimeW, ziformula = ~., weights = Total, family = binomial)
summary(QInvZeroInf)
check_model(QInvZeroInf, n_bins=10)
emmeans(QInvZeroInf, ~Model, type = "response")
emmeans(QInvZeroInf, ~Model, type = "conditional")
QInvPlot<-plot(emmeans(QInvZeroInf, ~Model, type = "response"), comparisons=TRUE, ylab="Treatment Type", xlab="")  + ggtitle("B) The Quoin")
QInvPlot
pairs(emmeans(QInvZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(QInvZeroInf, ~Model, type = "response")))
Anova(QInvZeroInf)

#Zero-Inflated weighted binomial for total vigilance
QTimeW$T_vig<-(QTimeW$V_low+QTimeW$V_mid+QTimeW$V_high)
QVigZeroInf<-glmmTMB(T_vig ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = QTimeW, ziformula = ~., weights = Total, family = binomial)
summary(QVigZeroInf)
check_model(QVigZeroInf, n_bins=10)
emmeans(QVigZeroInf, ~Model, type = "response")
emmeans(QVigZeroInf, ~Model, type = "conditional")
QVigPlot<-plot(emmeans(QVigZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="", xlim=0.55)  + ggtitle("B) The Quoin")
pairs(emmeans(QVigZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(QVigZeroInf, ~Model, type = "response")))
Anova(QVigZeroInf)

QTimeW$V_direction<-((QTimeW$V_toward+1)/(QTimeW$V_away+1))
QVigDir<-glmmTMB(V_direction ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = QTimeW)
summary(QVigDir)
check_model(QVigDir, n_bins=10)
emmeans(QVigDir, ~Model, type="response")
pairs(emmeans(QVigDir, ~Model, type="response"))

#Zero-Inflated weighted binomial for total locomotion
QTimeW$T_loc<-(QTimeW$L_towards+QTimeW$L_away+QTimeW$L_across)
QLocZeroInf<-glmmTMB(T_loc ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = QTimeW, ziformula = ~., weights = Total, family = binomial)
summary(QLocZeroInf)
check_model(QLocZeroInf, n_bins=10)
emmeans(QLocZeroInf, ~Model, type = "response")
emmeans(QLocZeroInf, ~Model, type = "conditional")
QLocPlot<-plot(emmeans(QLocZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="") + ggtitle("B) The Quoin")
pairs(emmeans(QLocZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(QlocZeroInf, ~Model, type = "response")))
Anova(QLocZeroInf)
QLocPlot

#Cutting data back into photos and applying vigilance into an ordered factor
QTime_by_frame <-
  QTime %>%
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
QTimeW_by_frame<-QTime_by_frame[which(QTime_by_frame$Species == "Wallaby"),]
QVigModelW<-clmm(vigilance_level~Model+n_overlaps+Moon_Phase+Position+(1|Site)+Wind_Speed+Min_Temp+Rainfall, data=QTimeW_by_frame)
QVigModelW<-clmm(vigilance_level~Model+(1|Site)+Moon_Phase+n_overlaps, data=QTimeW_by_frame)
QVigModelW
summary(QVigModelW)
Anova(QVigModelW)

#Outputting modeled effect of statue type on vigilance with P-values
emmeans(QVigModelW,~Model, type='response')
summary(emmeans(QVigModelW, ~Model), type = 'response')
pairs(emmeans(QVigModelW, ~Model))
QVigLevelPlot<-plot(emmeans(QVigModelW, ~Model, type="response", mode="exc.prob"), xlab="", ylab="Treatment Type") + ggtitle("B) The Quoin")
plot(pairs(emmeans(QVigModelW, ~Model)), comparisons = TRUE)  

plot(emmeans(QVigModelW,~Model | vigilance_level, mode="prob", type='response'), horizontal=TRUE, xlab="Likelihood of Vigilance Level", ylab="Treatment Type")
QVigLevelFullPlot<-plot(emmeans(QVigModelW,~Model | vigilance_level, mode="prob", type='response'), horizontal=TRUE, xlab="", ylab="Treatment Type") + ggtitle("B) The Quoin")

#Defining preferred distance as a vector
QTime$Pref_Dist<-factor(QTime$Pref_Dist, levels=c("<5m", "5-10m", "10-15m", ">15m"), ordered=TRUE)
QTime$Pref_Dist
QTimeW<-QTime[which(QTime$Species=="Wallaby"),]


#Creating a model for preferred distance
QDistModelW<-clmm(Pref_Dist~Model+(1|Site)+Moon_Phase+n_overlaps, data=QTimeW)
QDistModelW
summary(QDistModelW)
check_model(QDistModelW, residual_type="normal")
model_performance(QDistModelW)
Anova(QDistModelW)

#Outputting and displaying modeled effect of statue on preferred distance
emmeans(QDistModelW,~Model, mode="exc.prob", type='response')
emmeans(QDistModelW, ~n_overlaps, mode="exc.prob", type="response")
summary(emmeans(QDistModelW, ~Model, mode="exc.prob", type = 'response'))
pairs(emmeans(QDistModelW, ~Model, mode="exc.prob", type = 'response'))
plot(emmeans(QDistModelW, ~Model, mode="exc.prob", type = 'response'))
plot(pairs(emmeans(QDistModelW, ~Model, mode="exc.prob", type = 'response')))
plot(pairs(emmeans(QDistModelW, ~Model, mode="exc.prob", type = 'response')), comparisons = TRUE)

pulkrob.chisq(QDistModelW, c("Model", "Site"))
lipsitz.test(QDistModelW)
autoplot.clm(QDistModelW, what = c("qq", "fitted", "covariate"))
plot(n_overlaps ~ Moon_Phase, data=QTimeW)
with(QTimeW, cor(n_overlaps, Moon_Phase, method = 'kendall'))

plot(emmeans(QDistModelW,~Model | Pref_Dist, mode="prob", type='response'), horizontal=TRUE, xlab="EMM of probability of selection", ylab="Treatment Type", main="B) The Quoin")
QFullDistPlot<-plot(emmeans(QDistModelW,~Model | Pref_Dist, mode="prob", type='response'), horizontal=TRUE, xlab="", ylab="Treatment Type", main="B) The Quoin") +ggtitle("B) The Quoin")
QFullDistPlot
QDistPlot<-plot(emmeans(QDistModelW, ~Model, mode="exc.prob", type = 'response'), ylab="Treatment Type", xlab="") + ggtitle("B) The Quoin")


pulkrob.chisq(QDistModelW, c("Model", "Site"))
lipsitz.test(QDistModelW)
autoplot.clm(QDistModelW, what = c("qq", "fitted", "covariate"))
plot(n_overlaps ~ Moon_Phase, data=QTimeW)
with(QTimeW, cor(n_overlaps, Moon_Phase, method = 'kendall'))

emmeans(QDistModelW,~Position, type='response')
summary(emmeans(QDistModelW, ~Position), type = 'response')
pairs(emmeans(QDistModelW, ~Position))
plot(pairs(emmeans(QDistModelW, ~Position)))
plot(pairs(emmeans(QDistModelW, ~Position)), comparisons = TRUE)

#NEW idea for binary regression of distance
QBinModelW<-glmmTMB(Pref_Dist2 ~Model+(1|Site)+Moon_Phase+Position, family="binomial", data = QTimeW)
QBinModelW
summary(QBinModelW)
emmeans(QBinModelW, ~Model, type="response")
pairs(emmeans(QBinModelW, ~Model, type="response"))
check_model(QBinModelW, n_bins=10)
model_performance(QBinModelW)
Anova(QBinModelW)
plot(emmeans(QBinModelW, ~Model, type="response"))
