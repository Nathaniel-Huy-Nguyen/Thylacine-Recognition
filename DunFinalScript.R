Dun<-read.csv("DunTime3.csv")
Dup<-Dun[which(Dun$Species != "Dog"), ]
DupW<-Dup[which(Dup$Species == "Wallaby"),]
DupW
nrow(DupW)/nrow(Dun)

aggregate(DupW$Model, by=list(DupW$Model), FUN=length)
aggregate(Dun$Species, by=list(Dun$Species), FUN=length)

mean(DupWT$Total)
mean(DupWC$Total)

DunPie<-c(410, 204,  1, 3, 382, 12, 11, 47, 28, 74)
names(DunPie)<-c("Bennett's Wallaby","Swamp Wallaby","Common Ringtail Possum","Short-Beaked Echidna","Feral Rabbit", "Feral Red Deer", "Feral Cat", "Feral Red Fox", "Feral Dog", "Unknown")
DunSpeciesPlot<-barchart(DunPie, col=c("palegreen", "palegreen3", "paleturquoise3", "paleturquoise4", "papayawhip", 'peachpuff', "palevioletred", "palevioletred2", "violetred2", "white"), ylab="Species", xlab="Number of Observations", main="C) Mt Sturgeon")
names(DunPie)<-c("Bennett's Wallaby", "Swamp Wallaby", "", "", "Feral Rabbit")
pie(DunPie, col=c("palegreen", "palegreen3", "paleturquoise3", "paleturquoise4", "papayawhip", 'peachpuff', "palevioletred", "palevioletred2", "violetred2", "white"))
colours()

table(DupW$Model, DupW$Pref_Dist)
chisq.test(table(DupW$Model, DupW$Pref_Dist))

Dupn_by_model<-aggregate(DupW$Model, by=list(Model=DupW$Model, Site=DupW$Site), FUN=length, drop=FALSE)
Dupn_by_model
Dupn_by_model$x <- ifelse(is.na(Dupn_by_model$x),0,Dupn_by_model$x)
DunVisitModel<-glmmTMB(x ~ Model + (1|Site), data = Dupn_by_model, family = "nbinom1")
check_model(DunVisitModel)
summary(DunVisitModel)
emmeans(DunVisitModel, ~ Model, type = "response")
pairs(emmeans(DunVisitModel, ~ Model, type = "response"))
DunTotalPlot<-plot(emmeans(DunVisitModel, ~ Model, type = "response"), ylab="Treatment Type", xlab="Expected Number of Wallabies Per Site") + ggtitle("C) Mt Sturgeon")
DunTotalPlot

#Poisson GLM for total visits
Dupn_by_model<-aggregate(DupW$Model, by=list(Model=DupW$Model, Site=DupW$Site), FUN=length)
DupVisitModel<-glmer(x ~ Model + (1|Site), data = Dupn_by_model, family = "poisson")
check_model(DupVisitModel)
summary(DupVisitModel)
emmeans(DupVisitModel, ~ Model, type = "response")
pairs(emmeans(DupVisitModel, ~ Model, type = "response"))
plot(emmeans(DupVisitModel, ~ Model, type = "response"))
Anova(DupVisitModel)

## use DHARMa functions instead:
Dunsimres <- simulateResiduals(DupVisitModel)

## standard DHARMa diagnostic plots
plot(Dunsimres)

#Converting times to POSIX and creating Duration
DunTime <- mutate(Dup,
                  ## create POSIXct objects for the first & last times
                  First_datetime = strptime(paste(Date,First),format = '%d/%m/%Y %H:%M:%S'),
                  Last_datetime = strptime(paste(Date,Last),format = '%d/%m/%Y %H:%M:%S'),
                  ## adjust if they span the date boundary
                  Last_datetime = Last_datetime + ifelse(Last_datetime >= First_datetime, 0, 60*60*24),
                  Duration = as.numeric(difftime(Last_datetime, First_datetime, units = 'secs')))

#GROUP SIZE STUFF
install.packages("ivs")
library(ivs)


Sites <- unique(DunTime$Site)
site <- Sites[1]
## create a new variable for the number of overlaps
DunTime$n_overlaps <- NA
for(site in Sites){
  irows <- which(DunTime$Site == site)
  ## create an "interval vector" (which contains the start/end times) for this site
  IV <- with(DunTime[irows,],
             ## if the start-time and end-time are identical, add 1 second to the end-time
             iv(First_datetime, Last_datetime + ifelse(First_datetime == Last_datetime, 1, 0)))
  ## Number of overlapping overservations *for this site*
  ##
  ## The '-1' excludes overlap with "self", so it is the count of
  ## the number of other observations that overlap with it
  DunTime$n_overlaps[irows] <- iv_count_overlaps(IV,IV) - 1
  ## to look at what the overlaps are (if you want to check):
  ## subset(iv_locate_overlaps(IV,IV), needles != haystack)
  ## the "subset( ..., needles != haystack)" code filters out self-overlaps
}
DunTimeW<-DunTime[which(DunTime$Species == "Wallaby"),]

mean(DunTimeW$n_overlaps)
pairwise.t.test((DunTimeW$n_overlaps), DunTimeW$Model, na.rm=TRUE)

DunOverlapModel<-glm(n_overlaps ~ Model + Site + Position + Moon_Phase, data=DunTimeW)
summary(DunOverlapModel)
anova(DunOverlapModel)
plot(predictorEffect("Model", DunOverlapModel))
emmeans(DunOverlapModel, ~Model)
pairs(emmeans(DunOverlapModel, ~Model))

#Zero-inflated poisson model for number of overlaps
DunOverZeroInf<-glmmTMB(n_overlaps ~ Model + (1|Site) + Moon_Phase, data = DunTimeW, ziformula = ~ Model + (1|Site) + Position + Moon_Phase, family=poisson)
DunOverZeroInf
summary(DunOverZeroInf)
check_model(DunOverZeroInf, n_bins=10)
emmeans(DunOverZeroInf, ~Model, type = "response")
emmeans(DunOverZeroInf, ~Model, type = "conditional")
DunOverPlot<-plot(emmeans(DunOverZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="EMM of Number of Overlaps")  + ggtitle("C) Mt Sturgeon")
pairs(emmeans(DunOverZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(DunOverZeroInf, ~Model, type = "response")))

DunDurationModel<-lmer(Duration ~ Model + (1|Site) + Position + Moon_Phase, data=DunTimeW)
summary(DunDurationModel)
check_model(DunDurationModel)
anova(DunDurationModel)
DunDurationPlot<-plot(emmeans(DunDurationModel, ~Model), ylab="Treatment Type", xlab="EMM of Known Duration")  + ggtitle("C) Mt Sturgeon")
emmeans(DunDurationModel, ~Model)
pairs(emmeans(DunDurationModel, ~Model))

#Zero inflated weighted binomial for grazing
DunGrazeZeroInf<-glmmTMB(Grazing ~ Model + (1|Site) + Position + Moon_Phase+n_overlaps, data = DunTimeW, ziformula = ~., weights = Total, family = binomial)
summary(DunGrazeZeroInf)
check_model(DunGrazeZeroInf, n_bins=10)
emmeans(DunGrazeZeroInf, ~Model, type = "response")
emmeans(DunGrazeZeroInf, ~Model, type = "response") 
DunGrazePlot<-plot(emmeans(DunGrazeZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="Grazing Probability") + ggtitle("C) Mt Sturgeon")
pairs(emmeans(DunGrazeZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(DunGrazeZeroInf, ~Model, type = "response")))
Anova(DunGrazeZeroInf)

#Zero-Inflated weighted binomial for investigating
DunInvZeroInf<-glmmTMB(Investigating ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = DunTimeW, ziformula = ~., weights = Total, family = binomial)
summary(DunInvZeroInf)
check_model(DunInvZeroInf, n_bins=10)
emmeans(DunInvZeroInf, ~Model, type = "response")
emmeans(DunInvZeroInf, ~Model, type = "conditional")
DunInvPlot<-plot(emmeans(DunInvZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="Investigating Probability") + ggtitle("C) Mt Sturgeon")
pairs(emmeans(DunInvZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(DunInvZeroInf, ~Model, type = "response")))
Anova(DunInvZeroInf)

#Zero-Inflated weighted binomial for total vigilance
DunTimeW$T_vig<-(DunTimeW$V_low+DunTimeW$V_mid+DunTimeW$V_high)
DunVigZeroInf<-glmmTMB(T_vig ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = DunTimeW, ziformula = ~., weights = Total, family = binomial)
summary(DunVigZeroInf)
check_model(DunVigZeroInf, n_bins=10)
emmeans(DunVigZeroInf, ~Model, type = "response")
emmeans(DunVigZeroInf, ~Model, type = "conditional")
DunVigPlot<-plot(emmeans(DunVigZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="Vigilance Probability", xlim=0.55) + ggtitle("C) Mt Sturgeon")
pairs(emmeans(DunVigZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(DunVigZeroInf, ~Model, type = "response")))
Anova(DunVigZeroInf)

#Zero-Inflated weighted binomial for total locomotion
DunTimeW$T_loc<-(DunTimeW$L_towards+DunTimeW$L_away+DunTimeW$L_across)
DunLocZeroInf<-glmmTMB(T_loc ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = DunTimeW, ziformula = ~., weights = Total, family = binomial)
summary(DunLocZeroInf)
check_model(DunLocZeroInf, n_bins=10)
emmeans(DunLocZeroInf, ~Model, type = "response")
emmeans(DunLocZeroInf, ~Model, type = "conditional")
DunLocPlot<-plot(emmeans(DunLocZeroInf, ~Model, type = "response"), ylab="Treatment Type", xlab="Locomotion Probability") + ggtitle("C) Mt Sturgeon")
pairs(emmeans(DunLocZeroInf, ~Model, type = "response"))
plot(pairs(emmeans(DunLocZeroInf, ~Model, type = "response")))
Anova(DunLocZeroInf)
DunLocPlot

DunTimeW$V_direction<-((DunTimeW$V_toward+1)/(DunTimeW$V_away+1))
DunVigDir<-glmmTMB(V_direction ~ Model + (1|Site) + Position + Moon_Phase + n_overlaps, data = DunTimeW)
summary(DunVigDir)
check_model(DunVigDir, n_bins=10)
emmeans(DunVigDir, ~Model, type="response")
pairs(emmeans(DunVigDir, ~Model, type="response"))

#Cutting data back into photos and applying vigilance into an ordered factor
DunTime_by_frame <-
  DunTime %>%
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
DunTimeW_by_frame<-DunTime_by_frame[which(DunTime_by_frame$Species == "Wallaby"),]
DunVigModelW<-clmm(vigilance_level~Model+n_overlaps+Moon_Phase+Position+(1|Site)+Wind_Speed+Min_Temp+Rainfall, data=DunTimeW_by_frame)
DunVigModelW<-clmm(vigilance_level~Model+(1|Site)+Moon_Phase+n_overlaps, data=DunTimeW_by_frame)
DunVigModelW
summary(DunVigModelW)
check_mode(DunVigModelW)
Anova(DunVigModelW)


#Outputting modeled effect of statue type on vigilance with P-values
emmeans(DunVigModelW,~Model, type='response')
summary(emmeans(DunVigModelW, ~Model), mode = "exc.prob", type = 'response')
pairs(emmeans(DunVigModelW, ~Model))
DunVigLevelPlot<-plot(emmeans(DunVigModelW, ~Model, type="response", mode="exc.prob"), ylab="Treatment Type", xlab="EMM of Exceedance Probability") +ggtitle("C) Mt Sturgeon")
plot(pairs(emmeans(DunVigModelW, ~Model)), comparisons = TRUE)  

plot(emmeans(DunVigModelW,~Model | vigilance_level, mode="prob", type='response'), horizontal=TRUE, xlab="Likelihood of Vigilance Level", ylab="Treatment Type")
DunVigLevelFullPlot<-plot(emmeans(DunVigModelW,~Model | vigilance_level, mode="prob", type='response'), horizontal=TRUE, xlab="EMM of Probability of Vigilance Level", ylab="Treatment Type") + ggtitle("A) Mt Sturgeon")

#Defining preferred distance as a vector
DunTime$Pref_Dist<-factor(DunTime$Pref_Dist, levels=c("<5m", "5-10m", "10-15m", ">15m"), ordered=TRUE)
DunTime$Pref_Dist
DunTimeW<-DunTime[which(DunTime$Species=="Wallaby"),]


#Creating a model for preferred distance
DunDistModelW<-clmm(Pref_Dist~Model+(1|Site)+Moon_Phase+n_overlaps, data=DunTimeW)
DunDistModelW
summary(DunDistModelW)
check_model(DunDistModelW)
model_performance(DunDistModelW)
Anova(DunDistModelW)

#Outputting and displaying modeled effect of statue on preferred distance
emmeans(DunDistModelW,~Model, mode="exc.prob", type='response')
plot(emmeans(DunDistModelW,~Model, mode="exc.prob", type='response'))
summary(emmeans(DunDistModelW, ~Model, mode="exc.prob", type = 'response'))
pairs(emmeans(DunDistModelW, ~Model, mode="exc.prob", type = 'response'))
plot(pairs(emmeans(DunDistModelW, ~Model, mode="exc.prob", type = 'response')))
plot(pairs(emmeans(DunDistModelW, ~Model, mode="exc.prob", type = 'response')), comparisons = TRUE)

plot(emmeans(DunDistModelW,~Model | Pref_Dist, mode="prob", type='response'), horizontal=TRUE, xlab="EMM of Probability of Selection", ylab="Treatment Type", main="C) Mount Sturgeon")
DunFullDistPlot<-plot(emmeans(DunDistModelW,~Model | Pref_Dist, mode="prob", type='response'), horizontal=TRUE, xlab="EMM of Probability of Selection", ylab="Treatment Type", main="C) Mount Sturgeon") + ggtitle("C) Mount Sturgeon")
DunDistPlot<-plot(emmeans(DunDistModelW, ~Model, mode="exc.prob", type = 'response'), ylab="Treatment Type", xlab="EMM of Exceedance Probability") + ggtitle("C) Mt Sturgeon")


emmeans(DunDistModelW,~Position, type='response')
summary(emmeans(DunDistModelW, ~Position), type = 'response')
pairs(emmeans(DunDistModelW, ~Position))
plot(pairs(emmeans(DunDistModelW, ~Position)))
plot(pairs(emmeans(DunDistModelW, ~Position)), comparisons = TRUE)

#NEW idea for binary regression of distance
DunBinModelW<-glmmTMB(Pref_Dist2 ~Model+(1|Site)+Moon_Phase+Position, family="binomial", data = DunTimeW)
DunBinModelW
summary(DunBinModelW)
emmeans(DunBinModelW, ~Model, type="response")
pairs(emmeans(DunBinModelW, ~Model, type="response"))
check_model(DunBinModelW, n_bins=10)
model_performance(DunBinModelW)
Anova(DunBinModelW)
plot(emmeans(DunBinModelW, ~Model, type="response"))
