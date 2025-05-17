
#######################
### Amazon bubbles###
#######################

#The code contains a number of analyses, tests and calculations not included in the actual paper. 
#If you need help to find the code for a specific graph or analysis, please contact Stina Edelfeldt, nyaka.epost@gmail.com.



library(mgcv)
library(lme4)
library(car)
#library('lmerTest')
library(multcomp)
library(ggplot2)
library(ggsignif)

library(nlme)
library(mvtnorm)
library(survival)
library(TH.data)
library(MASS)
library(lattice)
library(Formula)
library(gamm4)

library(drc) 
library (Metrics)
library("readxl")
library(moments)
library(pscl)
library(boot)
library(nnet)
library("AER")
library(broom)
library(Matrix)

library(tidyverse)
library(rstatix)
library(ggpubr)
library(lmtest)
library(sandwich)
library(boot)
library(Rmisc)
library(dplyr)
library(ggpmisc) 
library(gamlss)
library(forcats)
library(Hmisc)
library(ggbreak) 
library(DescTools)


library(ltm)

library("lsmeans")
library("BWStest")

library(rcompanion)
library(emmeans)

library(FSA)

library(performance)




#Rafaela libraries
library(rmarkdown)
library(ggthemes)
library(stringi)
library(stringr)
library(knitr)
library(ggsci)
#library(ggsflabel)
library(conflicted)
#library(scico)
library(colorspace)
library(viridis)

#Tonya libraries
library(datasets)
library(multcompView)

#######
#Fit a model
#fit <- lm(Fertility ~ . , data = swiss)

# Function for Root Mean Squared Error
#RMSE <- function(error) { sqrt(mean(error^2)) }
#RMSE(fit$residuals)

# Function for Mean Absolute Error
#mae <- function(error) { mean(abs(error)) }
#mae(fit$residuals)


#Just one care you should take, if there are NAs in the data, use na.rm=T in the functions.


citation("pscl")


######



datBeforeMerge <- read.csv("C:/Users/stied67/OneDrive - Linköpings universitet/Amazonas/Track results/AmazonBubblesCSV.csv")

#Summary TSc data to count file - only one time
datToCountMeanTSc <- datBeforeMerge %>%
  group_by(ORIG_FID_Section) %>% 
  dplyr::summarise(across(where(is.numeric),
                   list(mean = ~ mean(MeanTSc, na.rm = TRUE), count = ~ sum(!is.na(WaterDepth)))))
  write.csv(datToCountMeanTSc, "datToCountMeanTSc.csv")

datCountBeforeMerge <- read.csv("C:/Users/stied67/OneDrive - Linköpings universitet/Amazonas/Track results/AmazonBubblesCount.csv")
datTotPingBeforeMerge <- read.csv("C:/Users/stied67/OneDrive - Linköpings universitet/Amazonas/Track results/AmazonBubblesTotPing.csv")
datBathBeforeMerge <- read.csv("C:/Users/stied67/OneDrive - Linköpings universitet/Amazonas/Track results/Bathymetri_Points_Depth_CSV.csv")

datWaterBody <- read.csv("C:/Users/stied67/OneDrive - Linköpings universitet/Amazonas/Track results/AmazonBubblesWaterBody.csv")

dat <- merge(x = datBeforeMerge, y = datWaterBody, by = "Echogram", all.x = TRUE)
datCount <- merge(x = datCountBeforeMerge, y = datWaterBody, by = "Echogram", all.x = TRUE)
datTotPing <- merge(x = datTotPingBeforeMerge, y = datWaterBody, by = "Echogram", all.x = TRUE)
datBath <- merge(x = datBathBeforeMerge, y = datWaterBody, by = "Echogram", all.x = TRUE)

dat$LocationSplit <-paste(dat$Location, dat$WaterBody, sep = "")
datCount$LocationSplit <-paste(datCount$Location, datCount$WaterBody, sep = "")
datTotPing$LocationSplit <-paste(datTotPing$Location, datTotPing$WaterBody, sep = "")
datBath$LocationSplit <-paste(datBath$Location, datBath$WaterBody, sep = "")

print (dat %>% group_by(LocationSplit) %>% tally(), n=50)
print (datCount %>% group_by(LocationSplit) %>% tally(), n=50)
print (datTotPing %>% group_by(LocationSplit) %>% tally(), n=50)
print (datBath %>% group_by(LocationSplit) %>% tally(), n=50)



dat$DepthPercent <- dat$TargetDepth/dat$WaterDepth
dat$DepthCategory <- as.factor(ifelse(dat$DepthPercent <=0.2, '0-20%', 
                                    ifelse(dat$DepthPercent <=0.4, '21-40%',
                                           ifelse(dat$DepthPercent <=0.6, '41-60%',
                                                  ifelse(dat$DepthPercent <=0.8, '61-80%',
                                                         ifelse(dat$DepthPercent <=1.0, '81-100%',
                                                         ('Missing')))))))


dat$BottomToTarget <- dat$WaterDepth-dat$TargetDepth
dat$TotPing <- dat$ToPing - dat$FromPing - 1

dat$LocationSplit <- as.factor(dat$LocationSplit)
dat$Echogram <- as.factor(dat$Echogram)


dat$WaterTypeAS <- substr(dat$LocationSplit, 1, 1)
datCount$WaterTypeAS <- substr(datCount$LocationSplit, 1, 1)
datTotPing$WaterTypeAS <- substr(datTotPing$LocationSplit, 1, 1)
datBath$WaterTypeAS <- substr(datBath$LocationSplit, 1, 1)


#Make normalized TS
dat$AdjustedTSa <- ((14.9*0.368) -57.7)/((dat$TargetDepth *0.368) -57.7)*dat$MeanTSc

dat$AdjustedTSb <- ((14.9*-0.289) -49.5)/((dat$BottomToTarget *-0.289) -49.5)*dat$MeanTSc

dat$AdjustedTS <- ((14.9*0.368) -57.7)/((dat$TargetDepth *0.368) -57.7)*((14.9*-0.289) -49.5)/((dat$BottomToTarget *-0.289) -49.5)*dat$MeanTSc

#Calculate volume adn Flux from Tonyas Size does matter 2015
#This is mean Flux per bubble spot
#Vz (m/s) x Vb (ml/m3) = Flux (ml/m2/d) (omräkning ingår)
dat$Mean.Vz.day <- dat$Mean.Vz.*86400
dat$σbs <- 10^(dat$MeanTSc/10)
dat$Vb_ml <- 846095*dat$σbs^1.344
dat$diam_mm <- (6*dat$Vb/pi)^(1/3)*10
dat$Flux <- dat$Mean.Vz.day*dat$Vb_ml*-1
dat <- dat[dat$Flux >= 0, ]

#Densitet Ch4 0,657 kg/m³ (dock svårt veta hur mycket Ch4 i bubbla)

aggregate(x = dat$Vb_ml,                # Specify data column
          by = list(dat$LocationSplit),              # Specify group indicator
          FUN = sum) 



#Setorderoffactorsinriverorder-i.e.numberorder
datCount$LocationSplit = factor(datCount$LocationSplit, levels = c("S1L","S1M","S2L","N3M","N3R","N4M","N5R","N6L","N6M","N7M","A7aR","A7bM","A7cM","A8L","A9L","A9M","T10L","T10M","T11F","T11M","T11R"))
dat$LocationSplit = factor(dat$LocationSplit, levels = c("S1L","S1M","S2L","N3M","N3R","N4M","N5R","N6L","N6M","N7M","A7aR","A7bM","A7cM","A8L","A9L","A9M","T10L","T10M","T11F","T11M","T11R"))
datBath$LocationSplit = factor(datBath$LocationSplit, levels = c("S1L","S1M","S2L","N3M","N3R","N4M","N5R","N6L","N6M","N7M","A7aR","A7bM","A7cM","A8L","A9L","A9M","T10L","T10M","T11F","T11M","T11R"))
datTotPing$LocationSplit = factor(datTotPing$LocationSplit, levels = c("S1L","S1M","S2L","N3M","N3R","N4M","N5R","N6L","N6M","N7M","A7aR","A7bM","A7cM","A8L","A9L","A9M","T10L","T10M","T11F","T11M","T11R"))

datCount$WaterTypeAS = factor(datCount$WaterTypeAS, levels=c("S", "N", "A", "T"))
dat$WaterTypeAS = factor(dat$WaterTypeAS, levels=c("S", "N", "A", "T"))
datBath$WaterTypeAS = factor(datBath$WaterTypeAS, levels=c("S", "N", "A", "T"))
datTotPing$WaterTypeAS = factor(datTotPing$WaterTypeAS, levels=c("S", "N", "A", "T"))

datCount$WaterBody = factor(datCount$WaterBody, levels=c("F", "L", "M", "R"))
dat$WaterBody = factor(dat$WaterBody, levels=c("F", "L", "M", "R"))
datBath$WaterBody = factor(datBath$WaterBody, levels=c("F", "L", "M", "R"))
datTotPing$WaterBody = factor(datTotPing$WaterBody, levels=c("F", "L", "M", "R"))


#Calculate volume and area density - rounded to integer
datCount$Section_Volume <- tan(DegToRad(7/2))*100*(datCount$DepthMean^2-2.0^2)
datCount$Volume_Density <- datCount$Join_Count/datCount$Section_Volume
datCount$Volume_DensityInt <- as.integer(datCount$Join_Count/datCount$Section_Volume*10000)
datCount$Area_Density_Decimal <- (datCount$DepthMean-2.0)*datCount$Volume_Density
datCount$Area_Density <- as.integer((datCount$DepthMean-2.0)*datCount$Volume_Density*10000)


#Make new datasets with only N4
datN4 <-subset(dat, LocationSplit %in% c("N4M"))
datCountN4 <-subset(datCount, LocationSplit %in% c("N4M"))
datBathN4 <-subset(datBath, LocationSplit %in% c("N4M"))

#Make new datasets with no N4
datNoN4 <-subset(dat, LocationSplit!="N4M")
datCountNoN4 <-subset(datCount, LocationSplit!="N4M")
datBathNoN4 <-subset(datBath, LocationSplit!="N4M")
datTotPingNoN4 <- subset(datTotPing, LocationSplit!="N4M")
datCountZerosNoN4 <-subset(datCountNoN4, Join_Count==0)
datCountNoZerosNoN4 <-subset(datCountNoN4, Join_Count>0)


levels(datNoN4$WaterTypeAS)[levels(datNoN4$WaterTypeAS)=='N'] <- 'NNoN4'
levels(datCountNoN4$WaterTypeAS)[levels(datCountNoN4$WaterTypeAS)=='N'] <- 'NNoN4'
levels(datBathNoN4$WaterTypeAS)[levels(datBathNoN4$WaterTypeAS)=='N'] <- 'NNoN4'
levels(datTotPingNoN4$WaterTypeAS)[levels(datTotPingNoN4$WaterTypeAS)=='N'] <- 'NNoN4'
levels(datCountZerosNoN4$WaterTypeAS)[levels(datCountZerosNoN4$WaterTypeAS)=='N'] <- 'NNoN4'

levels(datNoN4$WaterBody)[levels(datNoN4$WaterBody)=='M'] <- 'MNoN4'
levels(datCountNoN4$WaterBody)[levels(datCountNoN4$WaterBody)=='M'] <- 'MNoN4'
levels(datBathNoN4$WaterBody)[levels(datBathNoN4$WaterBody)=='M'] <- 'MNoN4'
levels(datTotPingNoN4$WaterBody)[levels(datTotPingNoN4$WaterBody)=='M'] <- 'MNoN4'
levels(datCountZerosNoN4$WaterBody)[levels(datCountZerosNoN4$WaterBody)=='M'] <- 'MNoN4'

#Make new datasets with no T11Triburaty
datNoT11Tributary <-subset(dat,  LocationSplit!="T11R")
datCountNoT11Tributary <-subset(datCount, LocationSplit!="T11R")
datBathNoT11Tributary <-subset(datBath, LocationSplit!="T11R")
datTotPingNoT11Tributary <- subset(datTotPing, LocationSplit!="T11R")
datCountZerosNoT11Tributary <-subset(datCountNoT11Tributary, Join_Count==0)
datCountNoZerosNoT11Tributary <-subset(datCountNoT11Tributary, Join_Count>0)


levels(datNoT11Tributary$WaterTypeAS)[levels(datNoT11Tributary$WaterTypeAS)=='T'] <- 'TNoT11R'
levels(datCountNoT11Tributary$WaterTypeAS)[levels(datCountNoT11Tributary$WaterTypeAS)=='T'] <- 'TNoT11R'
levels(datBathNoT11Tributary$WaterTypeAS)[levels(datBathNoT11Tributary$WaterTypeAS)=='T'] <- 'TNoT11R'
levels(datTotPingNoT11Tributary$WaterTypeAS)[levels(datTotPingNoT11Tributary$WaterTypeAS)=='T'] <- 'TNoT11R'
levels(datCountZerosNoT11Tributary$WaterTypeAS)[levels(datCountZerosNoT11Tributary$WaterTypeAS)=='T'] <- 'TNoT11R'

levels(datNoT11Tributary$WaterBody)[levels(datNoT11Tributary$WaterBody)=='R'] <- 'RNoT11R'
levels(datCountNoT11Tributary$WaterBody)[levels(datCountNoT11Tributary$WaterBody)=='R'] <- 'RNoT11R'
levels(datBathNoT11Tributary$WaterBody)[levels(datBathNoT11Tributary$WaterBody)=='R'] <- 'RNoT11R'
levels(datTotPingNoT11Tributary$WaterBody)[levels(datTotPingNoT11Tributary$WaterBody)=='R'] <- 'RNoT11R'
levels(datCountZerosNoT11Tributary$WaterBody)[levels(datCountZerosNoT11Tributary$WaterBody)=='R'] <- 'RNoT11R'


#Make new datasets with no zeros
datCountNoZero <-subset(datCount, Join_Count>0)
datCountNoN4NoZero <-subset(datCountNoN4, Join_Count>0)
datCountNoT11TributaryNoZero <-subset(datCountNoT11Tributary, Join_Count>0)

#Make new subset with only 0-sections
datCountZeros <-subset(datCount, Join_Count==0)
datCountNoN4Zeros <-subset(datCountNoN4, Join_Count==0)
datCountNoT11TributaryZeros <-subset(datCountNoT11Tributary, Join_Count==0)





#Calculate bubbles per area - rounded to integer _ NOT USE INSTEAD AREA DENSITY
#datCount$Section_Area <- 100 * (datCount$DepthMean-2.0)
#datCount$BubblesPerAreaDecimal <- datCount$Join_Count/datCount$Section_Area
#datCount$BubblesPerArea <- as.integer(datCount$Join_Count/datCount$Section_Area*10000)










#Create equation for conversion from depth till volume scannad - or just create formula for compensating volume bigger further down
#depth <- c(1.6, 2, 2.5, 3, 4, 5, 7, 10, 13, 15, 18, 20, 23, 25, 28, 30, 33)

#volume <- c(0.002, 0.016, 0.041, 0.079, 0.204, 0.409, 1.141, 3.350, 7.373, 11.332, 19.590, 26.877, 40.882, 52.505, 73.770, 90.736, 120.773)

#sonar.df <- data.frame(depth, volume)

#fit1 <- lm(log(volume) ~ log(depth), data = sonar.df)
#
#you were missing a `*` and had bad starting values
#m <- nls(volume ~ a*(depth^b), sonar.df, 
#         start = list(a = exp(coef(fit1)[1]), b = coef(fit1)[2])) # power formula: y = a*x^b
#summary(m)
#plot(volume ~ depth)
#curve(predict(m, newdata = data.frame(depth = x)), add = TRUE)
#Equation volume (m^3 per ping)= 0.003352*depth^3.001 (Vol. SED-beam pr. ping(m^3))










summary(datNoN4)

summary(dat)
summary(dat$LocationSplit)
summary(datCount)
summary(datBath)
summary(datCountSize)
summary(dat$LocationSplit)
summary(datCountNoZero)

lapply(dat, class)
lapply(datCount, class)


#head(datCountN4NoZero)
head(datCount)
head(datCountNoN4)
options(max.print = 1000)



#alt theme in ggplot
#theme_bw() + 




write.csv(dat, "datcsv.csv")
write.csv(datCount, "datCountcsv.csv")
write.csv(datCountNoZero, "datCountNoZerocsv.csv")
write.csv(datCountNoN4, "datCountNoN4csv.csv")
write.csv(datCountNoT11Tributary, "datCountNoT11Tributarycsv.csv")







###SUMMARY DATA###

#Number of bubbles per location
print(dat %>% group_by(LocationSplit) %>% tally(), n=50) 

#Number of pings per location
aggregate(x = datTotPing$Totping,                # Specify data column
          by = list(datTotPing$LocationSplit),              # Specify group indicator
          FUN = sum) 

#Seconds per location (time bars)
aggregate(x = datTotPing$Totsek,                # Specify data column
          by = list(datTotPing$LocationSplit),              # Specify group indicator
          FUN = sum) 

#Number of sections per location
print (datCount %>% group_by(LocationSplit) %>% tally(), n=50)

#Number of 0-sections per location
print (datCountZeros %>% group_by(LocationSplit) %>% tally(), n=50)


#Number of bubbles per WaterType
dat %>% group_by(WaterTypeAS) %>% tally()
datNoN4 %>% group_by(WaterTypeAS) %>% tally() 
datNoT11Tributary %>% group_by(WaterTypeAS) %>% tally()

#Number of pings per Color WaterType
aggregate(x = datTotPing$Totping,                # Specify data column
          by = list(datTotPing$WaterTypeAS),              # Specify group indicator
          FUN = sum) 
aggregate(x = datTotPingNoN4$Totping,                # Specify data column
          by = list(datTotPingNoN4$WaterTypeAS),              # Specify group indicator
          FUN = sum) 
aggregate(x = datTotPingNoT11Tributary$Totping,                # Specify data column
          by = list(datTotPingNoT11Tributary$WaterTypeAS),              # Specify group indicator
          FUN = sum) 

#Seconds per Color WaterType (time bar)
aggregate(x = datTotPing$Totsek,                # Specify data column
          by = list(datTotPing$WaterTypeAS),              # Specify group indicator
          FUN = sum) 
aggregate(x = datTotPingNoN4$Totsek,                # Specify data column
          by = list(datTotPingNoN4$WaterTypeAS),              # Specify group indicator
          FUN = sum) 
aggregate(x = datTotPingNoT11Tributary$Totsek,                # Specify data column
          by = list(datTotPingNoT11Tributary$WaterTypeAS),              # Specify group indicator
          FUN = sum) 

#Number of sections per WaterType
datCount %>% group_by(WaterTypeAS) %>% tally()  
datCountNoN4 %>% group_by(WaterTypeAS) %>% tally() 
datCountNoT11Tributary %>% group_by(WaterTypeAS) %>% tally()  

#Number of 0-sections per WaterType
datCountZeros %>% group_by(WaterTypeAS) %>% tally()
datCountZerosNoN4 %>% group_by(WaterTypeAS) %>% tally()
datCountZerosNoT11Tributary %>% group_by(WaterTypeAS) %>% tally()



#Number of bubbles per waterbody
dat %>% group_by(WaterBody) %>% tally()  
datNoN4 %>% group_by(WaterBody) %>% tally()  
datNoT11Tributary %>% group_by(WaterBody) %>% tally()

#Number of pings per waterbody
aggregate(x = datTotPing$Totping,                # Specify data column
          by = list(datTotPing$WaterBody),              # Specify group indicator
          FUN = sum) 
aggregate(x = datTotPingNoN4$Totping,                # Specify data column
          by = list(datTotPingNoN4$WaterBody),              # Specify group indicator
          FUN = sum) 
aggregate(x = datTotPingNoT11Tributary$Totping,                # Specify data column
          by = list(datTotPingNoT11Tributary$WaterBody),              # Specify group indicator
          FUN = sum) 

#Seconds per waterbody (time bar)
aggregate(x = datTotPing$Totsek,                # Specify data column
          by = list(datTotPing$WaterBody),              # Specify group indicator
          FUN = sum) 
aggregate(x = datTotPingNoN4$Totsek,                # Specify data column
          by = list(datTotPingNoN4$WaterBody),              # Specify group indicator
          FUN = sum) 
aggregate(x = datTotPingNoT11Tributary$Totsek,                # Specify data column
          by = list(datTotPingNoT11Tributary$WaterBody),              # Specify group indicator
          FUN = sum) 

#Number of sections per waterbody
datCount %>% group_by(WaterBody) %>% tally()  
datCountNoN4 %>% group_by(WaterBody) %>% tally() 
datCountNoT11Tributary %>% group_by(WaterBody) %>% tally()  

#Number of 0-sections per waterbody
datCountZeros %>% group_by(WaterBody) %>% tally()
datCountZerosNoN4 %>% group_by(WaterBody) %>% tally()
datCountZerosNoT11Tributary %>% group_by(WaterBody) %>% tally()


#Mean depth per LocationSplit, WaterType and WaterBody
aggregate(x = datBath$Depth,                # Specify data column
          by = list(datBath$LocationSplit),              # Specify group indicator
          FUN = mean) 


aggregate(x = datBath$Depth,                # Specify data column
          by = list(datBath$WaterType),              # Specify group indicator
          FUN = mean) 

aggregate(x = datBathNoN4$Depth,                # Specify data column
          by = list(datBathNoN4$WaterType),              # Specify group indicator
          FUN = mean) 

aggregate(x = datBathNoT11Tributary$Depth,                # Specify data column
          by = list(datBathNoT11Tributary$WaterType),              # Specify group indicator
          FUN = mean) 


aggregate(x = datBath$Depth,                # Specify data column
          by = list(datBath$WaterBody),              # Specify group indicator
          FUN = mean) 

aggregate(x = datBathNoN4$Depth,                # Specify data column
          by = list(datBathNoN4$WaterBody),              # Specify group indicator
          FUN = mean) 

aggregate(x = datBathNoT11Tributary$Depth,                # Specify data column
          by = list(datBathNoT11Tributary$WaterBody),              # Specify group indicator
          FUN = mean) 


#sd per LocationSplit, WaterType and WaterBody
aggregate(x = datBath$Depth,                # Specify data column
          by = list(datBath$LocationSplit),              # Specify group indicator
          FUN = sd) 

aggregate(x = datBath$Depth,                # Specify data column
          by = list(datBath$WaterType),              # Specify group indicator
          FUN = sd) 

aggregate(x = datBathNoN4$Depth,                # Specify data column
          by = list(datBathNoN4$WaterType),              # Specify group indicator
          FUN = sd) 

aggregate(x = datBathNoT11Tributary$Depth,                # Specify data column
          by = list(datBathNoT11Tributary$WaterType),              # Specify group indicator
          FUN = sd) 

aggregate(x = datBath$Depth,                # Specify data column
          by = list(datBath$WaterBody),              # Specify group indicator
          FUN = sd) 

aggregate(x = datBathNoN4$Depth,                # Specify data column
          by = list(datBathNoN4$WaterBody),              # Specify group indicator
          FUN = sd) 





###DEPTH###


summary(dat$WaterDepth)
summary(datBath$Depth)
head(datBath)


#Exclude Forest
datNoF <- subset (dat, !(WaterBody =="F"))
datNoN4NoF <- subset (datNoN4, !(WaterBody =="F"))
datNoT11TributaryNoF <- subset (datNoT11Tributary, !(WaterBody =="F"))

datBathNoN20 <-subset(datBathNoN4, !(Depth < 2.0))
datBathNoMissing <-subset(datBath, !(WaterTypeAS %in% c("Missing")))
datBathNoMissing20 <-subset(datBathNoMissing,  !(Depth < 2.0))
datBathNoMissing20NoF <-subset(datBathNoMissing20, !(WaterBody =="F"))

datBathNoMissingNoN4 <-subset(datBathNoN4, !(WaterTypeAS %in% c("Missing")))
datBathNoMissing20NoN4 <-subset(datBathNoMissingNoN4,  !(Depth < 2.0))
datBathNoMissing20NoN4NoF <-subset(datBathNoMissing20NoN4, !(WaterBody =="F"))

datBathNoMissingNoT11Tributary <-subset(datBathNoT11Tributary, !(WaterTypeAS %in% c("Missing")))
datBathNoMissing20NoT11Tributary <-subset(datBathNoMissingNoT11Tributary,  !(Depth < 2.0))
datBathNoMissing20NoT11TributaryNoF <-subset(datBathNoMissing20NoT11Tributary, !(WaterBody =="F"))



#General test if water depth normally distributed (neccessary?)
hist(dat$WaterDepth, breaks = 100,main="Histogram of water depth", xlab="Water depth (m)")
ggqqplot(dat$WaterDepth)
qqnorm(dat$WaterDepth)
qqline(dat$WaterDepth)
#shapiro.test(dat$AdjustedTS) funkar inte - för många


#KS-test "too sensistive" for comparing distributions?
#ks.test(datBathAmazonMain$Depth, datAmazonMain$WaterDepth)

#Something someone recommended - dont work right now for some reason
#pv.obs = t.test(datBathAmazonMain$Depth,datAmazonMain$WaterDepth)$pval; pv.obs
#pv = replicate(10^5, t.test(datBathAmazonMain$Depth, sample(datAmazonMain$WaterDepth)))
#mean(pv <= pv.obs)


#Wilcoxon rank sum test differences between Bath and bubble water depth one by one
datAmazonMain <- subset(dat,WaterTypeAS %in% c("A"))
datBathAmazonMain <- subset(datBathNoMissing20,WaterTypeAS %in% c("A"))
wilcox.test(datBathAmazonMain$Depth, datAmazonMain$WaterDepth, paired = FALSE)

datBlack <- subset(dat,WaterTypeAS %in% c("N"))
datBathBlack <- subset(datBathNoMissing20,WaterTypeAS %in% c("N"))
wilcox.test(datBathBlack$Depth, datBlack$WaterDepth, paired = FALSE)

datClear <- subset(dat,WaterTypeAS %in% c("T"))
datBathClear <- subset(datBathNoMissing20,WaterTypeAS %in% c("T"))
wilcox.test(datBathClear$Depth, datClear$WaterDepth, paired = FALSE)

datWhite <- subset(dat,WaterTypeAS %in% c("S"))
datBathWhite <- subset(datBathNoMissing20,WaterTypeAS %in% c("S"))
wilcox.test(datBathWhite$Depth, datWhite$WaterDepth, paired = FALSE)

datBlackNoN4 <- subset(datNoN4,WaterTypeAS %in% c("NNoN4"))
datBathBlackNoN4 <- subset(datBathNoMissing20NoN4,WaterTypeAS %in% c("NNoN4"))
wilcox.test(datBathBlackNoN4$Depth, datBlackNoN4$WaterDepth, paired = FALSE)

datClearNoT11Tributary <- subset(datNoT11Tributary,WaterTypeAS %in% c("TNoT11R"))
datBathClearNoT11Tributary <- subset(datBathNoMissing20NoT11Tributary,WaterTypeAS %in% c("TNoT11R"))
wilcox.test(datBathClearNoT11Tributary$Depth, datClearNoT11Tributary$WaterDepth, paired = FALSE)


#Mann-Whitney-Wilcoxon Test  Wilcoxon test differences between Bath and bubble water depth one by one
datForest <- subset(dat,WaterBody %in% c("F"))
datBathForest <- subset(datBathNoMissing20,WaterBody %in% c("F"))
wilcox.test(datBathForest$Depth, datForest$WaterDepth, paired = FALSE)

datLake <- subset(dat,WaterBody %in% c("L"))
datBathLake <- subset(datBathNoMissing20,WaterBody %in% c("L"))
wilcox.test(datBathLake$Depth, datLake$WaterDepth, paired = FALSE)

datMain <- subset(dat,WaterBody %in% c("M"))
datBathMain <- subset(datBathNoMissing20,WaterBody %in% c("M"))
wilcox.test(datBathMain$Depth, datMain$WaterDepth, paired = FALSE)

datTrib <- subset(dat,WaterBody %in% c("R"))
datBathTrib <- subset(datBathNoMissing20,WaterBody %in% c("R"))
wilcox.test(datBathTrib$Depth, datTrib$WaterDepth, paired = FALSE)

datMainNoN4 <- subset(datNoN4,WaterBody %in% c("MNoN4"))
datBathMainNoN4 <- subset(datBathNoMissing20NoN4,WaterBody %in% c("MNoN4"))
wilcox.test(datBathMainNoN4$Depth, datMainNoN4$WaterDepth, paired = FALSE)

datTribNoT11R <- subset(datNoT11Tributary,WaterBody %in% c("RNoT11R"))
datBathTribNoT11R<- subset(datBathNoMissing20NoT11Tributary,WaterBody %in% c("RNoT11R"))
wilcox.test(datBathTribNoT11R$Depth, datTribNoT11R$WaterDepth, paired = FALSE) 






#Mann-Whitney-Wilcoxon Test  Wilcoxon test differences between Bath and bubble water depth one by one
datS1L <- subset(dat,LocationSplit %in% c("S1L"))
datBathS1L <- subset(datBathNoMissing20,LocationSplit %in% c("S1L"))
wilcox.test(datBathS1L$Depth, datS1L$WaterDepth, paired = FALSE)

datS1M <- subset(dat,LocationSplit %in% c("S1M"))
datBathS1M <- subset(datBathNoMissing20,LocationSplit %in% c("S1M"))
wilcox.test(datBathS1M$Depth, datS1M$WaterDepth, paired = FALSE)

datS2L <- subset(dat,LocationSplit %in% c("S2L"))
datBathS2L <- subset(datBathNoMissing20,LocationSplit %in% c("S2L"))
wilcox.test(datBathS2L$Depth, datS2L$WaterDepth, paired = FALSE)

datN3M <- subset(dat,LocationSplit %in% c("N3M"))
datBathN3M <- subset(datBathNoMissing20,LocationSplit %in% c("N3M"))
wilcox.test(datBathN3M$Depth, datN3M$WaterDepth, paired = FALSE)

datN3R <- subset(dat,LocationSplit %in% c("N3R"))
datBathN3R <- subset(datBathNoMissing20,LocationSplit %in% c("N3R"))
wilcox.test(datBathN3R$Depth, datN3R$WaterDepth, paired = FALSE)

datN4M <- subset(dat,LocationSplit %in% c("N4M"))
datBathN4M <- subset(datBathNoMissing20,LocationSplit %in% c("N4M"))
wilcox.test(datBathN4M$Depth, datN4M$WaterDepth, paired = FALSE)

datN5R <- subset(dat,LocationSplit %in% c("N5R"))
datBathN5R <- subset(datBathNoMissing20,LocationSplit %in% c("N5R"))
wilcox.test(datBathN5R$Depth, datN5R$WaterDepth, paired = FALSE)

datN6L <- subset(dat,LocationSplit %in% c("N6L"))
datBathN6L <- subset(datBathNoMissing20,LocationSplit %in% c("N6L"))
wilcox.test(datBathN6L$Depth, datN6L$WaterDepth, paired = FALSE)

datN6M <- subset(dat,LocationSplit %in% c("N6M"))
datBathN6M <- subset(datBathNoMissing20,LocationSplit %in% c("N6M"))
wilcox.test(datBathN6M$Depth, datN6M$WaterDepth, paired = FALSE)

datA7aR <- subset(dat,LocationSplit %in% c("A7aR"))
datBathA7aR <- subset(datBathNoMissing20,LocationSplit %in% c("A7aR"))
wilcox.test(datBathA7aR$Depth, datA7aR$WaterDepth, paired = FALSE)

datA7bM <- subset(dat,LocationSplit %in% c("A7bM"))
datBathA7bM  <- subset(datBathNoMissing20,LocationSplit %in% c("A7bM"))
wilcox.test(datBathA7bM$Depth, datA7bM$WaterDepth, paired = FALSE)

datA7cM <- subset(dat,LocationSplit %in% c("A7cM"))
datBathA7cM <- subset(datBathNoMissing20,LocationSplit %in% c("A7cM"))
wilcox.test(datBathA7cM$Depth, datA7cM$WaterDepth, paired = FALSE)

datA8L <- subset(dat,LocationSplit %in% c("A8L"))
datBathA8L <- subset(datBathNoMissing20,LocationSplit %in% c("A8L"))
wilcox.test(datBathA8L$Depth, datA8L$WaterDepth, paired = FALSE)

datA9L <- subset(dat,LocationSplit %in% c("A9L"))
datBathA9L <- subset(datBathNoMissing20,LocationSplit %in% c("A9L"))
wilcox.test(datBathA9L$Depth, datA9L$WaterDepth, paired = FALSE)

datA9M <- subset(dat,LocationSplit %in% c("A9M"))
datBathA9M <- subset(datBathNoMissing20,LocationSplit %in% c("A9M"))
wilcox.test(datBathA9M$Depth, datA9M$WaterDepth, paired = FALSE)

datT10L <- subset(dat,LocationSplit %in% c("T10L"))
datBathT10L <- subset(datBathNoMissing20,LocationSplit %in% c("T10L"))
wilcox.test(datBathT10L$Depth, datT10L$WaterDepth, paired = FALSE)

datT10M <- subset(dat,LocationSplit %in% c("T10M"))
datBathT10M <- subset(datBathNoMissing20,LocationSplit %in% c("T10M"))
wilcox.test(datBathT10M$Depth, datT10M$WaterDepth, paired = FALSE)

datT11F <- subset(dat,LocationSplit %in% c("T11F"))
datBathT11F <- subset(datBathNoMissing20,LocationSplit %in% c("T11F"))
wilcox.test(datBathT11F$Depth, datT11F$WaterDepth, paired = FALSE)


datT11M <- subset(dat,LocationSplit %in% c("T11M"))
datBathT11M <- subset(datBathNoMissing20,LocationSplit %in% c("T11M"))
wilcox.test(datBathT11M$Depth, datT11M$WaterDepth, paired = FALSE)

datT11R <- subset(dat,LocationSplit %in% c("T11R"))
datBathT11R <- subset(datBathNoMissing20,LocationSplit %in% c("T11R"))
wilcox.test(datBathT11R$Depth, datT11R$WaterDepth, paired = FALSE)



#This is the Water depth Water type region graph we use
ggplot(dat, aes(x=LocationSplit, y=WaterDepth)) + 
  geom_boxplot(width=0.2, position=position_dodge(),fill="grey")+ 
  geom_boxplot(data = datBathNoMissing20, aes(x=LocationSplit, y=Depth),width=0.2,position = position_nudge(0.3) )+
  scale_x_discrete("Locations") +
  xlab("River systems") + ylab("Water depth (m)")+
  theme_hc() +
  theme(axis.text = element_text(size = 24)) +
  theme(axis.title = element_text(size = 28)) +
  theme(plot.title = element_text(size = 28)) +
  scale_y_continuous(breaks = seq(0, 55, by=10), limits=c(0,55)) +
  coord_cartesian(ylim = c(50, 0))+
  stat_summary(size= 6, data=dat, fun.y = min, fun.ymin = length,
               geom = "text", aes(label = ..ymin..), vjust = -0.5) +
  stat_summary(size= 6, data=datBathNoMissing20, fun.y = min, fun.ymin = length,
               geom = "text", aes(x= LocationSplit, y= Depth, label = ..ymin..), vjust = -1.5,position = position_nudge(0.3)) +
  scale_x_discrete(labels = c("     S1L", "     S1M*", "     S2L", "     N3M*", "     N3R*", "     N4M*", "     N5R", "     N6L", "     N6M*", "     N7M", "     A7aR*", "     A7bM", "     A7cM*", "     A8L*", "     A9L*", "     A9M", "     T10L", "     T10M", "     T11F*", "     T11M*", "     T11R*"))
 
  # if I run it once and then and then add an + and mark both this and the one below I get both beside one another...
  


  annotations <- data.frame(
    X = c(Inf),
    Y =  c(-Inf),
    x_adjust = c(1),
    y_adjust = c(2) 
    )


#This is the Water depth Water type region graph we use
ggplot(dat, aes(x=WaterTypeAS, y=WaterDepth)) + 
  geom_boxplot(width=0.3, position=position_dodge(),fill="grey")+ 
  geom_boxplot(data = datBathNoMissing20, aes(x=WaterTypeAS, y=Depth),width=0.3,position = position_nudge(0.4) )+
  geom_boxplot(data= datNoN4,width=0.3, position=position_dodge(),fill="grey")+ 
  geom_boxplot(data = datBathNoMissing20NoN4, aes(x=WaterTypeAS, y=Depth),width=0.3,position = position_nudge(0.4) )+
  geom_boxplot(data= datNoT11Tributary,width=0.3, position=position_dodge(),fill="grey")+ 
  geom_boxplot(data = datBathNoMissing20NoT11Tributary, aes(x=WaterTypeAS, y=Depth),width=0.3,position = position_nudge(0.4) )+
  scale_x_discrete("River systems") +
  xlab("River systems") + ylab("Water depth (m)")+
  theme_hc() +
  theme(axis.text = element_text(size = 24)) +
  theme(axis.title = element_text(size = 28)) +
  theme(plot.title = element_text(size = 28)) +
  scale_y_continuous(breaks = seq(0, 55, by=10), limits=c(0,55))+
  coord_cartesian(ylim = c(50, -5))+
  stat_summary(size= 6, data=dat, fun.y = min, fun.ymin = length,
             geom = "text", aes(label = ..ymin..), vjust = -1) + 
  stat_summary(size= 6, data=datNoN4, fun.y = min, fun.ymin = length,
             geom = "text", aes(label = ..ymin..), vjust = -1) +
  stat_summary(size= 6, data=datNoT11Tributary, fun.y = min, fun.ymin = length,
               geom = "text", aes(label = ..ymin..), vjust = -1) +
  stat_summary(size= 6, data=datBathNoMissing20, fun.y = min, fun.ymin = length,
               geom = "text", aes(x= WaterTypeAS, y= Depth, label = ..ymin..), vjust = -1,position = position_nudge(0.4)) +
  stat_summary(size= 6, data=datBathNoMissing20NoN4, fun.y = min, fun.ymin = length,
               geom = "text", aes(x= WaterTypeAS, y= Depth, label = ..ymin..), vjust = -1,position = position_nudge(0.4) ) +
  stat_summary(size= 6, data=datBathNoMissing20NoT11Tributary, fun.y = min, fun.ymin = length,
               geom = "text", aes(x= WaterTypeAS, y= Depth, label = ..ymin..), vjust = -1,position = position_nudge(0.4) ) +
  geom_text(size = 10, data=annotations, aes(
    x=X,y=Y,hjust=x_adjust,vjust=y_adjust, label=c("(a)")))+
  scale_x_discrete(labels = c("     A*", "     N*", "     NNoN4M*", "     T", "     TNoT11R", "     S")) +

  # if I run it once and then and then add an x and mark both this and the one below I get both beside one another...


#This is the Water depth Waterbody graph we use
ggplot(datNoF, aes(x=WaterBody, y=WaterDepth)) + 
  geom_boxplot(width=0.3, position=position_dodge(),fill="grey")+ 
  geom_boxplot(data = datBathNoMissing20NoF, aes(x=WaterBody, y=Depth),width=0.3,position = position_nudge(0.4) )+
  geom_boxplot(data= datNoN4NoF,width=0.3, position=position_dodge(),fill="grey")+ 
  geom_boxplot(data = datBathNoMissing20NoN4NoF, aes(x=WaterBody, y=Depth),width=0.3,position = position_nudge(0.4) )+
  geom_boxplot(data= datNoT11TributaryNoF,width=0.3, position=position_dodge(),fill="grey")+ 
  geom_boxplot(data = datBathNoMissing20NoT11TributaryNoF, aes(x=WaterBody, y=Depth),width=0.3,position = position_nudge(0.4) )+
  scale_x_discrete("Waterbody type") +
  xlab("Waterbody type") + ylab("Water depth (m)")+
  theme_hc() +
  theme(axis.text = element_text(size = 24)) +
  theme(axis.title = element_text(size = 28)) +
  theme(plot.title = element_text(size = 28)) +
  scale_y_continuous(breaks = seq(0, 55, by=10), limits=c(0,55))+
  coord_cartesian(ylim = c(50, -5))+
  stat_summary(size= 6, data=datNoF, fun.y = min, fun.ymin = length,
               geom = "text", aes(label = ..ymin..), vjust = -1) +   
  stat_summary(size= 6, data=datNoN4NoF, fun.y = min, fun.ymin = length,
               geom = "text", aes(label = ..ymin..), vjust = -1) +
  stat_summary(size= 6, data=datNoT11TributaryNoF, fun.y = min, fun.ymin = length,
               geom = "text", aes(label = ..ymin..), vjust = -1) +
  stat_summary(size= 6, data=datBathNoMissing20NoF, fun.y = min, fun.ymin = length,
               geom = "text", aes(x= WaterBody, y= Depth, label = ..ymin..), vjust = -1,position = position_nudge(0.4)) +
  stat_summary(size= 6, data=datBathNoMissing20NoN4NoF, fun.y = min, fun.ymin = length,
               geom = "text", aes(x= WaterBody, y= Depth, label = ..ymin..), vjust = -1,position = position_nudge(0.4) ) +
  stat_summary(size= 6, data=datBathNoMissing20NoT11TributaryNoF, fun.y = min, fun.ymin = length,
               geom = "text", aes(x= WaterBody, y= Depth, label = ..ymin..), vjust = -1,position = position_nudge(0.4) ) +
  geom_text(size = 10, data=annotations, aes(
    x=X,y=Y,hjust=x_adjust,vjust=y_adjust, label=c("(b)"))) +
  scale_x_discrete(labels = c("      L*", "      M*", "      MNoN4M*", "      R*", "      RNoT11R*")) +
  theme(axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.y=element_blank(),legend.position="none")



  
################################
  
###NUMBER OF BUBBLES###

################################



###Section groups - Rafaela variant for bubbles - simplified by me after comments###

conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::summarize)  

#Rearranging datCountNoZero in sections - two variants (one in #) LocationSplit or Water type region
datCountNoZeroSections <- datCountNoZero %>%
  mutate(count_group = cut(Join_Count,
                           breaks = c(0, 10, 25, 50, 100, 150, 300, 600)
  )) 

datCountNoZeroSum <- datCountNoZeroSections %>%
  group_by(count_group) %>%
  summarize(n = n(), sum = sum(Join_Count)) %>%
  ungroup()

#Have this instead - only in SI
plotSectionBars <- ggplot() +
  geom_bar(
    data = datCountNoZeroSections,
    mapping =
      aes(
        x = count_group, y = Join_Count
      ),
    linewidth = 0.02,
    stat = "identity", ,color="dark gray",fill="dark grey", width = 0.5
  ) +
  geom_text(
    data = datCountNoZeroSum,
    mapping = aes(
      label = n,
      x = count_group, y = sum
    ), vjust = -1.2, size = 6
  ) +
  labs(
    x = "Section group (based on number of bubbles per section)",
    y = "Total number of bubbles in sections",
  ) +
  scale_x_discrete(labels = c("1-10", "11-25", "25-50", "51-100", "101-150", "151-300", "301-600")) +
  scale_y_continuous(breaks=seq(0,1700,100)) +
  theme_few() +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 18)) +
  theme(plot.title = element_text(size = 24)) +
  theme(legend.title=element_text(size=18)) + 
  theme(legend.text=element_text(size=14)) +
  coord_cartesian(ylim = c(0, 1700))
plot(plotSectionBars)



#Rearranging datCountNoZero in sections - two variants (one in #) LocationSplit or Water type region
datCountNoZeroSections <- datCount %>%
  mutate(count_group = cut(Area_Density,
                           breaks = c(-1, 1, 100, 250, 500, 1000, 1500, 3000, 6000)
  )) 

table(datCountNoZeroSections$count_group)

options(max.print=999999)
table(datCount$Join_Count,datCount$LocationSplit)



datCountNoZeroSum <- datCountNoZeroSections %>%
  group_by(count_group) %>%
  summarize(n = n(), sum = sum(Join_Count)) %>%
  ungroup()

#Have this instead - only in SI
plotSectionBars <- ggplot() +
  geom_bar(
    data = datCountNoZeroSections,
    mapping =
      aes(
        x = count_group, y = Area_Density
      ),
    linewidth = 0.02,
    stat = "identity", ,color="dark gray",fill="dark grey", width = 0.5
  ) +
  geom_text(
    data = datCountNoZeroSum,
    mapping = aes(
      label = n,
      x = count_group, y = sum
    ), vjust = -1.2, size = 6
  ) +
  labs(
    x = "Section group (based on number of bubbles per section)",
    y = "Total number of bubbles in sections",
  ) +
  scale_x_discrete(labels = c("1-10", "11-25", "25-50", "51-100", "101-150", "151-300", "301-600")) +
  scale_y_continuous(breaks=seq(0,1700,100)) +
  theme_few() +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 18)) +
  theme(plot.title = element_text(size = 24)) +
  theme(legend.title=element_text(size=18)) + 
  theme(legend.text=element_text(size=14)) +
  coord_cartesian(ylim = c(0, 1700))
plot(plotSectionBars)



#Number of bubbles model with LocationSplit and depth

#This does not converge...
BubblesModelFinalLoc <- zeroinfl(Area_Density ~0 + LocationSplit*DepthMean|0 + LocationSplit*DepthMean,
                              data = datCount, dist = "negbin",  na.action = na.omit)

#Run only location
BubblesModelFinalLoc <- zeroinfl(Area_Density ~LocationSplit|LocationSplit,
                              data = datCount, dist = "negbin",  na.action = na.omit)

#Run only depth
BubblesModelFinalLoc <- zeroinfl(Area_Density ~DepthMean|DepthMean,
                              data = datCount, dist = "negbin",  na.action = na.omit)

#Run depth fo each location
####



summary(BubblesModelFinalLoc)


#Analysis of Deviance Table (Type II tests)
Anova(BubblesModelFinalLoc,type="II",test="Chisq")

#Using emeans to get comparisons
marginal = emmeans(BubblesModelFinalLoc, ~ LocationSplit)

pairs(marginal, adjust="tukey")

cld(marginal,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")   ### Tukey adjustment for multiple comparisons - blir "sidak"

#Converting the results from "marginal" into a data frame (copy-paste via Excel)
datMarginalLoc <- data.frame(LocationSplit=c("S1L(ab)","S1M(ab)","S2L(b)","N3M(abc)","N3R(ab)","N4M(d)","N5R(ab)","N6L(a)","N6M(ab)","N7M(ab)","A7aR(ab)","A7bM(ab)","A7cM(ab)","A8L(abc)","A9L(abc)","A9M(ab)","T10L(abc)","T10M(ab)","T11F(abc)","T11M(ab)","T11R(cd)"), 
                 Area_Density=c(22.2,54.06,40.46,128.65,19.03,881.4,24.42,5.06,22.89,21.39,69.84,28.37,17.67,74.34,223.14,24.9,164.42,14.52,172.35,70.25,610.46), 
                 se=c(15.78,22.47,7.89,58.35,8.22,75.61,7.2,3.28,8.04,17.86,18.02,15.13,18,47.27,99.02,10.73,165.11,6.02,164.1,24.61,145.78),
                 Sections = c(34,46,414,108,129,512,415,163,101,57,202,48,30,40,24,110,15,157,6,88,39))

#Plots of emmeans and SE and significance group
ggplot(datMarginalLoc, aes(x=LocationSplit, y=Area_Density)) + 
  geom_point(size = 5)+ 
   geom_errorbar(aes(ymin=Area_Density-se, ymax=Area_Density+se), width=0,lty=1, size=0.3) +
  scale_x_discrete("Location") +
  xlab("Location") + ylab("Areal Density (bubbles per ha)")+
  theme_hc() +
  theme(axis.text = element_text(size = 24)) +
  theme(axis.title = element_text(size = 28)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  coord_cartesian(ylim = c(0, 100))+
  geom_text(aes(y = Area_Density + se, label = Sections), 
    size = 8, vjust = 0, nudge_y = 20)



#Number of bubbles model with WaterType (River system) and depth

#This does work, but as not anyone of the other do, I dont use it for consistency
BubblesModelFinalRS <- zeroinfl(Area_Density ~0 + WaterTypeAS*DepthMean|0 + WaterTypeAS*DepthMean,
                              data = datCount, dist = "negbin",  na.action = na.omit)

BubblesModelFinalRS <- zeroinfl(Area_Density ~0 + WaterTypeAS|0 + WaterTypeAS,
                              data = datCount, dist = "negbin",  na.action = na.omit)

BubblesModelFinalRSNoN4 <- zeroinfl(Area_Density ~WaterTypeAS|WaterTypeAS,
                                data = datCountNoN4, dist = "negbin",  na.action = na.omit)

BubblesModelFinalRSNoT11R <- zeroinfl(Area_Density ~WaterTypeAS|WaterTypeAS,
                                data = datCountNoT11Tributary, dist = "negbin",  na.action = na.omit)

summary(BubblesModelFinalRS)
summary(BubblesModelFinalRSNoN4)
summary(BubblesModelFinalRSNoT11R)

#Analysis of Deviance Table (Type II tests)
Anova(BubblesModelFinalRS,type="II",test="Chisq")
Anova(BubblesModelFinalRSNoN4,type="II",test="Chisq")
Anova(BubblesModelFinalRSNoT11R,type="II",test="Chisq")

#Using emeans to get comparisons
marginalRS = emmeans(BubblesModelFinalRS, ~ WaterTypeAS)
marginalRSNoN4 = emmeans(BubblesModelFinalRSNoN4, ~ WaterTypeAS)
marginalRSNoT11R = emmeans(BubblesModelFinalRSNoT11R, ~ WaterTypeAS)

pairs(marginalRS, adjust="tukey")
pairs(marginalRSNoN4, adjust="tukey")
pairs(marginalRSNoT11R, adjust="tukey")

cld(marginalRS,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")   ### Tukey adjustment for multiple comparisons - blir "sidak"

cld(marginalRSNoN4,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")   ### Tukey adjustment for multiple comparisons - blir "sidak"

cld(marginalRSNoT11R,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")   ### Tukey adjustment for multiple comparisons - blir "sidak"


#Converting the results from "marginal" into a data frame (copy-paste via Excel)
datMarginalRS <- data.frame(WaterTypeAS=c("A(ab)","N(c)","S(a)","T(b)"), 
                             Area_Density=c(59.6,324.7,40.5,117.3), 
                             se=c(11.21,27.91,7.38,22.94),
                             Sections = c(454, 1485, 494, 305 ))

datMarginalRSNoN4 <- data.frame(WaterTypeAS=c("NNoN4(a)"), 
                            Area_Density=c(31.7), 
                            se=c(4.96),
                            Sections = c(973))

datMarginalRSNoT11R <- data.frame(WaterTypeAS=c("TNoT11R(a)"), 
                            Area_Density=c(45), 
                            se=c(12.07),
                            Sections = c(266))




#Number of bubbles model with WaterBody type and depth

#This does not work
BubblesModelFinalWB <- zeroinfl(Area_Density ~0 + WaterBody*DepthMean|0 + WaterBody*DepthMean,
                                data = datCount, dist = "negbin",  na.action = na.omit)

#Use this
BubblesModelFinalWB <- zeroinfl(Area_Density ~WaterBody|WaterBody,
                                data = datCount, dist = "negbin",  na.action = na.omit)

BubblesModelFinalWBNoN4 <- zeroinfl(Area_Density ~WaterBody|WaterBody,
                                data = datCountNoN4, dist = "negbin",  na.action = na.omit)

BubblesModelFinalWBNoT11R <- zeroinfl(Area_Density ~WaterBody|WaterBody,
                                data = datCountNoT11Tributary, dist = "negbin",  na.action = na.omit)

summary(BubblesModelFinalWB)
summary(BubblesModelFinalWBNoN4)
summary(BubblesModelFinalWBNoT11R)

#Analysis of Deviance Table (Type II tests)
Anova(BubblesModelFinalWB,type="II",test="Chisq")
Anova(BubblesModelFinalWBNoN4,type="II",test="Chisq")
Anova(BubblesModelFinalWBNoT11R,type="II",test="Chisq")

#Using emeans to get comparisons
marginalWB = emmeans(BubblesModelFinalWB, ~ WaterBody)
marginalWBNoN4 = emmeans(BubblesModelFinalWBNoN4, ~ WaterBody)
marginalWBNoT11R = emmeans(BubblesModelFinalWBNoT11R, ~ WaterBody)

pairs(marginalWB, adjust="tukey")
pairs(marginalWBNoN4, adjust="tukey")
pairs(marginalWBNoT11R, adjust="tukey")

cld(marginalWB,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")   ### Tukey adjustment for multiple comparisons - blir "sidak"

cld(marginalWBNoN4,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")   ### Tukey adjustment for multiple comparisons - blir "sidak"

cld(marginalWBNoT11R,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")   ### Tukey adjustment for multiple comparisons - blir "sidak"

#Converting the results from "marginal" into a data frame (copy-paste via Excel)
datMarginalWB <- data.frame(WaterBody=c("F(ab)","L(a)","M(b)","R(a)"), 
                            Area_Density=c(172.3,42.2,385.3,64.3), 
                            se=c(178.95,7.36,32.36,10.39),
                            Sections = c(6, 690, 1257, 785 ))

datMarginalWBNoN4 <- data.frame(WaterBody=c("MNoN4(a)"), 
                                Area_Density=c(44.3), 
                                se=c(6.1),
                                Sections = c(745))

datMarginalWBNoT11R <- data.frame(WaterBody=c("RNoT11R(a)"), 
                                  Area_Density=c(35.8), 
                                  se=c(6.91),
                                  Sections = c(746))

#Plotting River System and Waterbody Type Areal Density - emmeans and SE and significance group
#River system 
ggplot(datMarginalRS, aes(x=WaterTypeAS, y=Area_Density)) + 
  geom_point(size = 5)+ 
  geom_errorbar(aes(ymin=Area_Density-se, ymax=Area_Density+se), width=0,lty=1, size=0.3) +
  geom_point(data= datMarginalRSNoN4, size = 5)+ 
  geom_errorbar(data= datMarginalRSNoN4, aes(ymin=Area_Density-se, ymax=Area_Density+se), width=0,lty=1, size=0.3) +
  geom_point(data= datMarginalRSNoT11R, size = 5)+ 
  geom_errorbar(data= datMarginalRSNoT11R, aes(ymin=Area_Density-se, ymax=Area_Density+se), width=0,lty=1, size=0.3) +
  scale_x_discrete("River System") +
  xlab("River System") + ylab("Areal Density (bubbles per ha)")+
  theme_hc() +
  theme(axis.text = element_text(size = 24)) +
  theme(axis.title = element_text(size = 28)) +
  coord_cartesian(ylim = c(0, 420))+
  geom_text(aes(y = Area_Density + se, label = Sections), 
            size = 8, vjust = 0, nudge_y = 10)+
  geom_text(data= datMarginalRSNoN4, aes(y = Area_Density + se, label = Sections), 
            size = 8, vjust = 0, nudge_y = 10)+
  geom_text(data= datMarginalRSNoT11R, aes(y = Area_Density + se, label = Sections), 
            size = 8, vjust = 0, nudge_y = 10)+
  geom_text(size = 9, aes(x=c(Inf),y=c(Inf),hjust=c(2),vjust=c(1), label=c("(a)")))+

#Waterbody type
ggplot(datMarginalWB, aes(x=WaterBody, y=Area_Density)) + 
  geom_point(size = 5)+ 
  geom_errorbar(aes(ymin=Area_Density-se, ymax=Area_Density+se), width=0,lty=1, size=0.3) +
  geom_point(data= datMarginalWBNoN4, size = 5)+ 
  geom_errorbar(data= datMarginalWBNoN4, aes(ymin=Area_Density-se, ymax=Area_Density+se), width=0,lty=1, size=0.3) +
  geom_point(data= datMarginalWBNoT11R, size = 5)+ 
  geom_errorbar(data= datMarginalWBNoT11R, aes(ymin=Area_Density-se, ymax=Area_Density+se), width=0,lty=1, size=0.3) +
  scale_x_discrete("Waterbody type") +
  xlab("Waterbody type") + ylab("Areal Density (bubbles per ha)")+
  theme_hc() +
  theme(axis.text = element_text(size = 24)) +
  theme(axis.title = element_text(size = 28)) +
  coord_cartesian(ylim = c(0, 420))+
  geom_text(aes(y = Area_Density + se, label = Sections), 
            size = 8, vjust = 0, nudge_y = 10)+
  geom_text(data= datMarginalWBNoN4, aes(y = Area_Density + se, label = Sections), 
            size = 8, vjust = 0, nudge_y = 10)+
  geom_text(data= datMarginalWBNoT11R, aes(y = Area_Density + se, label = Sections), 
            size = 8, vjust = 0, nudge_y = 10)+
  geom_text(size = 9, aes(x=c(Inf),y=c(Inf),hjust=c(2),vjust=c(1), label=c("(b)")))+
  theme(axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.y=element_blank(),legend.position="none")


#Influence of depth 



#Influence of depth for each instance
datCountS <- subset(datCount, WaterTypeAS=="S")
datCountN <- subset(datCount, WaterTypeAS=="N")
datCountA <- subset(datCount, WaterTypeAS=="A")
datCountT <- subset(datCount, WaterTypeAS=="T")
datCountNNoN4 <- subset(datCountNoN4, WaterTypeAS=="NNoN4")
datCountTNoT11R <- subset(datCountNoT11Tributary, WaterTypeAS=="TNoT11R")

datCountF <- subset(datCount, WaterBody=="F")
datCountL <- subset(datCount, WaterBody=="L")
datCountM <- subset(datCount, WaterBody=="M")
datCountR <- subset(datCount, WaterBody=="R")
datCountMNoN4 <- subset(datCountNoN4, WaterBody=="MNoN4")
datCountRNoT11R <- subset(datCountNoT11Tributary, WaterBody=="RNoT11R")

BubblesModelFinalAll <- zeroinfl(Area_Density ~WaterTypeAS * DepthMean,
                               data = datCount, dist = "negbin",  na.action = na.omit)

summary(BubblesModelFinalAll)
datCount$PredictAreaDensity <- predict(BubblesModelFinalAll, datCount)


BubblesModelFinalAllNoN4 <- zeroinfl(Area_Density ~WaterTypeAS * DepthMean,
                                 data = datCountNoN4, dist = "negbin",  na.action = na.omit)

summary(BubblesModelFinalAllNoN4)
datCountNoN4$PredictAreaDensity <- predict(BubblesModelFinalAllNoN4, datCountNoN4)


BubblesModelFinalAllNoT11 <- zeroinfl(Area_Density ~WaterTypeAS * DepthMean,
                                 data = datCountNoT11Tributary, dist = "negbin",  na.action = na.omit)

summary(BubblesModelFinalAllNoT11)
datCountNoT11Tributary$PredictAreaDensity <- predict(BubblesModelFinalAllNoT11, datCountNoT11Tributary)





#Area density - River system Depth models
BubblesModelFinalS <- zeroinfl(Area_Density ~DepthMean,
                              data = datCountS, dist = "negbin",  na.action = na.omit)

BubblesModelFinalN <- zeroinfl(Area_Density ~DepthMean,
                               data = datCountN, dist = "negbin",  na.action = na.omit)

BubblesModelFinalA <- zeroinfl(Area_Density ~DepthMean,
                               data = datCountA, dist = "negbin",  na.action = na.omit)

BubblesModelFinalT <- zeroinfl(Area_Density ~DepthMean,
                               data = datCountT, dist = "negbin",  na.action = na.omit)

BubblesModelFinalNNoN4 <- zeroinfl(Area_Density ~DepthMean,
                               data = datCountNNoN4, dist = "negbin",  na.action = na.omit)

BubblesModelFinalTNoT11R <- zeroinfl(Area_Density ~DepthMean,
                                   data = datCountTNoT11R, dist = "negbin",  na.action = na.omit)



summary(BubblesModelFinalS)
summary(BubblesModelFinalN)
summary(BubblesModelFinalA)
summary(BubblesModelFinalT)
summary(BubblesModelFinalNNoN4)
summary(BubblesModelFinalTNoT11R)

r2_zeroinflated(BubblesModelFinalS, method = c("correlation"))
r2_zeroinflated(BubblesModelFinalN, method = c("correlation"))
r2_zeroinflated(BubblesModelFinalA, method = c("correlation"))
r2_zeroinflated(BubblesModelFinalT, method = c("correlation"))
r2_zeroinflated(BubblesModelFinalNNoN4, method = c("correlation"))
r2_zeroinflated(BubblesModelFinalTNoT11R, method = c("correlation"))

datCountNoF <- subset(datCount, datCount$WaterBody!= "F")
datCountNoT11TributaryNoF <- subset(datCountNoT11Tributary, datCountNoT11Tributary$WaterBody!= "F")


#Area density - Waterbody Depth models ALL
BubblesModelFinalAllWB <- zeroinfl(Area_Density ~WaterBody * DepthMean,
                                 data = datCountNoF, dist = "negbin",  na.action = na.omit)

summary(BubblesModelFinalAllWB)
datCountNoF$PredictAreaDensityWB <- predict(BubblesModelFinalAllWB, datCountNoF)

BubblesModelFinalAllNoT11WB <- zeroinfl(Area_Density ~WaterBody * DepthMean,
                                      data = datCountNoT11TributaryNoF, dist = "negbin",  na.action = na.omit)

summary(BubblesModelFinalAllNoT11WB)
datCountNoT11TributaryNoF$PredictAreaDensityWB <- predict(BubblesModelFinalAllNoT11WB, datCountNoT11TributaryNoF)


#Area density - Waterbody Depth models 
BubblesModelFinalF <- zeroinfl(Area_Density ~DepthMean,
                               data = datCountF, dist = "negbin",  na.action = na.omit)

BubblesModelFinalL <- zeroinfl(Area_Density ~DepthMean,
                               data = datCountL, dist = "negbin",  na.action = na.omit)

BubblesModelFinalM <- zeroinfl(Area_Density ~DepthMean,
                               data = datCountM, dist = "negbin",  na.action = na.omit)

BubblesModelFinalR <- zeroinfl(Area_Density ~DepthMean,
                               data = datCountR, dist = "negbin",  na.action = na.omit)

BubblesModelFinalMNoN4 <- zeroinfl(Area_Density ~DepthMean,
                                   data = datCountMNoN4, dist = "negbin",  na.action = na.omit)

BubblesModelFinalRNoT11R <- zeroinfl(Area_Density ~DepthMean,
                                     data = datCountRNoT11R, dist = "negbin",  na.action = na.omit)


summary(BubblesModelFinalF)
summary(BubblesModelFinalL)
summary(BubblesModelFinalM)
summary(BubblesModelFinalR)
summary(BubblesModelFinalMNoN4)
summary(BubblesModelFinalRNoT11R)

r2_zeroinflated(BubblesModelFinalF, method = c("correlation"))
r2_zeroinflated(BubblesModelFinalL, method = c("correlation"))
r2_zeroinflated(BubblesModelFinalM, method = c("correlation"))
r2_zeroinflated(BubblesModelFinalR, method = c("correlation"))
r2_zeroinflated(BubblesModelFinalMNoN4, method = c("correlation"))
r2_zeroinflated(BubblesModelFinalRNoT11R, method = c("correlation"))


#Test about significances etc when removia a few values
datTest <- subset(datCountNNoN4,datCountNNoN4$Area_Density<2000 )
BubblesModelTest<- zeroinfl(Area_Density ~DepthMean,
                            data = datTest, dist = "negbin",  na.action = na.omit)
summary(BubblesModelTest)

datTest2 <- subset(datCountMNoN4,datCountMNoN4$Area_Density<2000 )
BubblesModelTest2<- zeroinfl(Area_Density ~DepthMean,
                            data = datTest2, dist = "negbin",  na.action = na.omit)
summary(BubblesModelTest2)





#Plotting figures
ndatCount <- rbind(datCount, datCountNoN4, datCountNoT11Tributary)
ndatCount$WaterTypeAS = factor(ndatCount$WaterTypeAS, levels=c("S", "N", "NNoN4","A", "T", "TNoT11R"))
ndatCount$WaterBody = factor(ndatCount$WaterBody, levels=c("F", "L", "M", "MNoN4","R", "RNoT11R"))

ndatCountZeros <- rbind(datCountZeros, datCountNoN4Zeros, datCountNoT11TributaryZeros)
ndatCountZeros$WaterTypeAS = factor(ndatCountZeros$WaterTypeAS, levels=c("S", "N", "NNoN4","A", "T", "TNoT11R"))
ndatCountZeros$WaterBody = factor(ndatCountZeros$WaterBody, levels=c("F", "L", "M", "MNoN4","R", "RNoT11R"))

ndatCountNoZero <- rbind(datCountNoZero, datCountNoN4NoZero, datCountNoT11TributaryNoZero)
ndatCountNoZero$WaterTypeAS = factor(ndatCountNoZero$WaterTypeAS, levels=c("S", "N", "NNoN4","A", "T", "TNoT11R"))
ndatCountNoZero$WaterBody = factor(ndatCountNoZero$WaterBody, levels=c("F", "L", "M", "MNoN4","R", "RNoT11R"))




R2valuesWT <- data.frame(WaterTypeAS = c("N","S","A","T","NNoN4", "TNoT11R"), Lab = c("R^2=0.027","R^2=0.016","R^2=0.060","R^2=0.062","", "R^2=0.047"))
R2valuesWS <- data.frame(WaterBody = c("F","L","M","R","MNoN4", "RNoT11R"), Lab = c("","R^2=0.070","R^2=0.022","R^2=0.068","", "R^2=0.021"))

ggplot() +
  theme_bw() +
  geom_point(
    data = ndatCount,
    aes(x = Area_Density, y = DepthMean), color = "red", size = 1
  ) +
  geom_smooth(
    data = datCount,
    aes(x = PredictAreaDensity, y = DepthMean), color = "black", linewidth = 1,
  ) +
  geom_smooth(
    data = datCountNoT11Tributary,
    aes(x = PredictAreaDensity, y = DepthMean), color = "black", linewidth = 1,
  ) +

    facet_wrap(vars(WaterTypeAS), 
             labeller = labeller(WaterTypeAS = 
                                   c("N" ="                     N*                (a)",
                                     "S" = "S*",
                                     "A" = "A*", 
                                     "T" = "T*", 
                                     "NNoN4" = "NNoN4",
                                     "TNoT11R" = "TNoT11R*")), ncol = 2) +
  labs(x = "Area density (bubbles per ha)", y = "Mean section depth (m)") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 24),
    strip.text = element_text(size = 22)
  ) +
  geom_text(aes(x = 8000, y = 35, label = Lab), size = 6, data = R2valuesWT) +
  coord_cartesian(xlim = c(0, 10000),ylim = c(40, 0) )+

  
ggplot(ndatCount, aes(x = Area_Density, y = DepthMean)) +
  geom_point(
   color = "red", size = 1
  ) +
  labs(x = "Area density (bubbles per ha)", y = "") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 24),
    strip.text = element_text(size = 22)
  ) +
  geom_smooth(
    data = datCountNoF,
    aes(x = PredictAreaDensityWB, y = DepthMean), color = "black", linewidth = 1,
  ) +
  geom_smooth(
    data = datCountNoT11TributaryNoF,
    aes(x = PredictAreaDensityWB, y = DepthMean), color = "black", linewidth = 1,
  ) +
  facet_wrap(vars(WaterBody), 
             labeller = labeller(WaterBody = 
                      c("F" = "F",
                        "L" ="                      L*                (b)",
                        "M" = "M*", 
                        "R" = "R*", 
                        "MNoN4" = "MNoN4",
                        "RNoT11R" = "RNoT11R*")), ncol = 2) +
  geom_text(aes(x = 8000, y = 35, label = Lab), size = 6, data = R2valuesWS) +
  coord_cartesian(xlim = c(0, 10000),ylim = c(40, 0) )
# + xlim(c(0, 2500)) +
#ylim(c(40, 0))


#Checking locations - just to see
ggplot() +
  theme_bw() +
  geom_point(
    data = ndatCountNoZero,
    aes(x = Area_Density, y = DepthMean), color = "red", size = 1
  ) +
  facet_wrap(vars(LocationSplit), ncol = 3) +
  labs(x = "Area density (bubbles per ha)", y = "Mean section depth (m)") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  coord_cartesian(xlim = c(0, 10000),ylim = c(20, 0) )

#Checking locations - just to see
ggplot() +
  theme_bw() +
  geom_point(
    data = ndatCount,
    aes(x = Area_Density, y = DepthMean), color = "red", size = 1
  ) +
  facet_wrap(vars(LocationSplit), ncol = 3) +
  labs(x = "Area density (bubbles per ha)", y = "Mean section depth (m)") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18)
  ) +
  coord_cartesian(xlim = c(0, 10000),ylim = c(20, 0) )


#############################

###SIZE###

#############################


summary(dat)


#MeanTSc vs TargetDepth. Used for the normalizing/adjusting equation (for depth)
ggplot(dat,aes(TargetDepth,MeanTSc))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("Bubble depth (m)") + ylab("TS (dB)") +
  stat_poly_eq(use_label(c("eq", "R2")), size = 6) +
  theme_hc()+
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) +
  theme(plot.title = element_text(size = 26)) +
  geom_text(size = 8, aes(x=c(Inf),y=c(Inf),hjust=c(2),vjust=c(1), label=c("(a)")))+
  coord_cartesian(ylim = c(-70, -20)) +
#If run with + after at the same time as graph below print beside one another
 
#MeanTSc vs BottomToDepth. Used for the normalizing/adjusting equation (for depth)
ggplot(dat,aes(BottomToTarget,MeanTSc))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab("Bubble distance from bottom (m)") + 
  stat_poly_eq(use_label(c("eq", "R2")), size = 6) +
  theme_hc()+
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) +
  theme(plot.title = element_text(size = 26))+
  geom_text(size = 8, aes(x=c(Inf),y=c(Inf),hjust=c(2),vjust=c(1), label=c("(b)")))+
  coord_cartesian(ylim = c(-70, -20)) +
  theme(axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.y=element_blank(),legend.position="none")

# stat_fit_glance(method = 'lm',
#                method.args = list(formula = my.formula),
#               geom = 'text',
# aes(label = paste("P-value = ", signif(after_stat(p.value), digits = 3), sep = "")),
# label.x = 'center', label.y = 0.35, size = 5)



#Histogram over Mean TSc
hist(dat$MeanTSc, xlim=c(-80, -20), ylim=c(0, 350),breaks = 100,
     xlab="Mean TS (dB)", main=NULL,
     cex.lab=1.5, cex.axis=1.5)

#Histogram over Adjusted TSc  
hist(dat$AdjustedTS, xlim=c(-80, -20), ylim=c(0, 350),breaks = 100,
     xlab="Adjusted Mean TS (dB)", main=NULL, 
     cex.lab=1.5, cex.axis=1.5)

#Histogram over water depth in general
hist(dat$WaterDepth, breaks = 100,xlab="Mean TSc (dB)")



#Correlation water depth adjustedTS
cor(dat$WaterDepth, dat$AdjustedTS,  method = "pearson", use = "complete.obs")

#NOrmalization check
qqnorm(dat$AdjustedTS)
qqline(dat$AdjustedTS)
ggqqplot(dat$AdjustedTS)
qqnorm(dat$MeanTSc)
qqline(dat$MeanTSc)
skewness(dat$AdjustedTS)
head(dat)

#Wilcoxon test differences between TS between River systems(if I write pairwise before - same result- should work even on multiple comparisons)
#All methods with 
#Also tested Dunn.test - gives very similar results
#People recommened Holm (and actually "none", but that was discussed...)
dat %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "bonferroni")
dat %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "holm")
dat %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "hochberg")
dat %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "hommel")
dat %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BH")
dat %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BY")
dat %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "none")

#ex dunntest
dunnTest(AdjustedTS ~ WaterTypeAS,data=dat,method="none")
dunnTest(AdjustedTS ~ WaterTypeAS,data=dat,method="holm")

dat  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "bonferroni")
dat  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "holm")
dat  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BY")
dat  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "none")

datNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "holm")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "hochberg")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "hommel")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BH")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BY")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "none")

datNoN4  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoN4  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "holm")
datNoN4  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BY")
datNoN4  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "none")

datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "holm")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "hochberg")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "hommel")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BH")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BY")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "none")

datNoT11Tributary  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoT11Tributary  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "holm")
datNoT11Tributary  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BY")
datNoT11Tributary  %>% wilcox_effsize(AdjustedTS ~ WaterTypeAS, p.adjust.method = "none")


datNoT11TributaryNoN4 <-subset(datNoT11Tributary, LocationSplit!="N4")
levels(datNoT11TributaryNoN4$WaterTypeAS)[levels(datNoT11TributaryNoN4$WaterTypeAS)=='N'] <- 'NNoN4'

datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "holm")
datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "hochberg")
datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "hommel")
datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BH")
datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterTypeAS, p.adjust.method = "BY")


#Wilcoxon test differences between TS WaterBody (if I write pairwise before - same result- should work even on multiple comparisons)
#All methods with 
#Also tested Dunn.test - gives very similar results
#People recommened Holm (and actually "none", but that was discussed...)
dat %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "bonferroni")
dat %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "holm")
dat %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "hochberg")
dat %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "hommel")
dat %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "BH")
dat %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "BY")
dat %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "none")

#ex dunntest
dunnTest(AdjustedTS ~ WaterBody,data=dat,method="none")
dunnTest(AdjustedTS ~ WaterBody,data=dat,method="holm")


datNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "bonferroni")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "holm")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "hochberg")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "hommel")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "BH")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "BY")
datNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "none")

datNoN4  %>% wilcox_effsize(AdjustedTS ~ WaterBody, p.adjust.method = "bonferroni")
datNoN4  %>% wilcox_effsize(AdjustedTS ~ WaterBody, p.adjust.method = "holm")
datNoN4  %>% wilcox_effsize(AdjustedTS ~ WaterBody, p.adjust.method = "BY")
datNoN4  %>% wilcox_effsize(AdjustedTS ~ WaterBody, p.adjust.method = "none")

datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "bonferroni")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "holm")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "hochberg")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "hommel")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "BH")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "BY")
datNoT11Tributary %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "none")

datNoT11Tributary  %>% wilcox_effsize(AdjustedTS ~ WaterBody, p.adjust.method = "bonferroni")
datNoT11Tributary  %>% wilcox_effsize(AdjustedTS ~ WaterBody, p.adjust.method = "holm")
datNoT11Tributary  %>% wilcox_effsize(AdjustedTS ~ WaterBody, p.adjust.method = "BY")
datNoT11Tributary  %>% wilcox_effsize(AdjustedTS ~ WaterBody, p.adjust.method = "none")


datNoT11TributaryNoN4 <-subset(datNoT11Tributary, LocationSplit!="N4")
levels(datNoT11TributaryNoN4$WaterBody)[levels(datNoT11TributaryNoN4$WaterBody)=='M'] <- 'MNoN4'

datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "bonferroni")
datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "holm")
datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "hochberg")
datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "hommel")
datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "BH")
datNoT11TributaryNoN4 %>% wilcox_test(AdjustedTS ~ WaterBody, p.adjust.method = "BY")






#This is the Water depth River System graph we use
ggplot(dat, aes(x=WaterTypeAS, y=AdjustedTS)) + 
  geom_boxplot(width=0.3, position=position_dodge(),fill="lightgrey")+ 
  geom_boxplot(data= datNoN4,width=0.3, position=position_dodge(),fill="lightgrey")+ 
  geom_boxplot(data= datNoT11Tributary,width=0.3, position=position_dodge(),fill="lightgrey")+ 
  scale_x_discrete("River system") +
  xlab("River system") + ylab("Adjusted TS (dB)")+
  stat_summary(fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoN4,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoT11Tributary,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  theme_hc() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) +
  theme(plot.title = element_text(size = 28)) +
  coord_cartesian(ylim = c(-70, -20))+
  stat_summary(size= 6, data=dat, fun = max, fun.max = length,
               geom = "text", aes(label = after_stat(ymax)), vjust = -1) +   
  stat_summary(size= 6, data=datNoN4, fun = max, fun.max = length,
               geom = "text", aes(label = after_stat(ymax)), vjust = -1) +
  stat_summary(size= 6, data=datNoT11Tributary, fun = max, fun.max = length,
               geom = "text", aes(label = after_stat(ymax)), vjust = -1) +
  geom_text(size = 8, aes(x=c(Inf),y=c(Inf),hjust=c(2),vjust=c(1), label=c("(a)")))+
  scale_x_discrete(labels = c("A(bc)", "N(a)", "NNoN4(c)", "S(ab)", "T(ab)", "TNoT11R(d)")) +
  # if I run it once and then and then mark both this and the one below I get both beside one another...
  
  
  #This is the Water depth Water body graph we use
  ggplot(dat, aes(x=WaterBody, y=AdjustedTS)) + 
  geom_boxplot(width=0.3, position=position_dodge(),fill="lightgrey")+ 
  geom_boxplot(data= datNoN4,width=0.3, position=position_dodge(),fill="lightgrey")+ 
  geom_boxplot(data= datNoT11Tributary,width=0.3, position=position_dodge(),fill="lightgrey")+ 
  scale_x_discrete("Waterbody type") +
  xlab("Waterbody type") + ylab("Adjusted TS (dB)")+
  stat_summary(fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoN4,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoT11Tributary,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  theme_hc() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) +
  theme(plot.title = element_text(size = 28)) +
  coord_cartesian(ylim = c(-70, -20))+
  stat_summary(size= 6, data=dat, fun = max, fun.max = length,
               geom = "text", aes(label = after_stat(ymax)), vjust = -1) +   
  stat_summary(size= 6, data=datNoN4, fun = max, fun.max = length,
               geom = "text", aes(label = after_stat(ymax)), vjust = -1) +
  stat_summary(size= 6, data=datNoT11Tributary, fun = max, fun.max = length,
               geom = "text", aes(label = after_stat(ymax)), vjust = -1) +
  geom_text(size = 8, aes(x=c(Inf),y=c(Inf),hjust=c(2),vjust=c(1), label=c("(b)")))+
  scale_x_discrete(labels = c("F(abc)", "L(a)", "M(a)", "MNoN4(c)", "R(a)", "RNoT11R(b)")) + 
  theme(axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.y=element_blank(),legend.position="none")
  





#This is the Water depth River System graph we use
ggplot(dat, aes(x=LocationSplit, y=AdjustedTS)) + 
  geom_boxplot(width=0.3, position=position_dodge(),fill="lightgrey")+ 
  geom_boxplot(data= datNoN4,width=0.3, position=position_dodge(),fill="lightgrey")+ 
  geom_boxplot(data= datNoT11Tributary,width=0.3, position=position_dodge(),fill="lightgrey")+ 
  scale_x_discrete("Location") +
  xlab("Location") + ylab("Adjusted TS (dB)")+
  stat_summary(fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoN4,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoT11Tributary,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  theme_hc() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) +
  theme(plot.title = element_text(size = 28)) +
  coord_cartesian(ylim = c(-70, -20))+
  stat_summary(size= 6, data=dat, fun = max, fun.max = length,
               geom = "text", aes(label = after_stat(ymax)), vjust = -1) +   
  stat_summary(size= 6, data=datNoN4, fun = max, fun.max = length,
               geom = "text", aes(label = after_stat(ymax)), vjust = -1) +
  stat_summary(size= 6, data=datNoT11Tributary, fun = max, fun.max = length,
               geom = "text", aes(label = after_stat(ymax)), vjust = -1) +
  #scale_x_discrete(labels = c("A(bc)", "N(a)", "NNoN4(c)", "S(ab)", "T(ab)", "TNoT11R(d)")) 
  # if I run it once and then and then mark both this and the one below I get both beside one another...
  
options(max.print=999999)

gla <- dat %>% wilcox_test(AdjustedTS ~ LocationSplit, p.adjust.method = "bonferroni")
gla <- dat %>% wilcox_test(AdjustedTS ~ LocationSplit, p.adjust.method = "holm")
dat %>% wilcox_test(AdjustedTS ~ LocationSplit, p.adjust.method = "hochberg")
dat %>% wilcox_test(AdjustedTS ~ LocationSplit, p.adjust.method = "hommel")
dat %>% wilcox_test(AdjustedTS ~ LocationSplit, p.adjust.method = "BH")
dat %>% wilcox_test(AdjustedTS ~ LocationSplit, p.adjust.method = "BY")
dat %>% wilcox_test(AdjustedTS ~ LocationSplit, p.adjust.method = "none")

print(gla, n=300)



#Means...
dat %>% 
  group_by(WaterTypeAS) %>% 
  dplyr::summarise(AdjustedTS_means = mean(AdjustedTS, na.rm = TRUE)) 
dat %>% 
  group_by(WaterTypeAS) %>% 
  dplyr::summarise(AdjustedTS_means = median(AdjustedTS, na.rm = TRUE)) 

dat %>% 
  group_by(WaterBody) %>% 
  dplyr::summarise(AdjustedTS_means = mean(AdjustedTS, na.rm = TRUE)) 
dat %>% 
  group_by(WaterBody) %>% 
  dplyr::summarise(AdjustedTS_means = median(AdjustedTS, na.rm = TRUE))  


#If compare means also... WatertypeAS
anova <- aov(AdjustedTS ~ WaterTypeAS, data = dat)
anovaNoNA <- aov(AdjustedTS ~ WaterTypeAS, data = datNoN4)
anovaNoT11Tributary<- aov(AdjustedTS ~ WaterTypeAS, data = datNoT11Tributary)
anovaNoT11TributaryNoN4<- aov(AdjustedTS ~ WaterTypeAS, data = datNoT11TributaryNoN4)

tukey <- TukeyHSD(anova)
tukeyNoNA <- TukeyHSD(anovaNoNA)
tukeyNoT11Tributary <- TukeyHSD(anovaNoT11Tributary)
tukeyNoT11TributaryNoN4 <- TukeyHSD(anovaNoT11TributaryNoN4)

cld <- multcompLetters4(anova, tukey)
cldNoNA <- multcompLetters4(anovaNoNA,tukeyNoNA)
cldNoT11Tributary <- multcompLetters4(anovaNoT11Tributary, tukeyNoT11Tributary)
cldNoT11TributaryNoN4 <- multcompLetters4(anovaNoT11TributaryNoN4, tukeyNoT11TributaryNoN4)

print(cld)
print(cldNoNA)
print(cldNoT11Tributary)
print(cldNoT11TributaryNoN4)



datNoT11TributaryNoN4 <-subset(datNoT11Tributary, LocationSplit!="N4")
levels(datNoT11TributaryNoN4$WaterTypeAS)[levels(datNoT11TributaryNoN4$WaterTypeAS)=='N'] <- 'NNoN4'

datNoT11TributaryNoN4 <-subset(datNoT11Tributary, LocationSplit!="N4")
levels(datNoT11TributaryNoN4$WaterBody)[levels(datNoT11TributaryNoN4$WaterBody)=='M'] <- 'MNoN4'

#If compare means also... WaterBody
anova <- aov(AdjustedTS ~ WaterBody, data = dat)
anovaNoNA <- aov(AdjustedTS ~ WaterBody, data = datNoN4)
anovaNoT11Tributary<- aov(AdjustedTS ~ WaterBody, data = datNoT11Tributary)
anovaNoT11TributaryNoN4<- aov(AdjustedTS ~ WaterBody, data = datNoT11TributaryNoN4)

tukey <- TukeyHSD(anova)
tukeyNoNA <- TukeyHSD(anovaNoNA)
tukeyNoT11Tributary <- TukeyHSD(anovaNoT11Tributary)
tukeyNoT11TributaryNoN4 <- TukeyHSD(anovaNoT11TributaryNoN4)

cld <- multcompLetters4(anova, tukey)
cldNoNA <- multcompLetters4(anovaNoNA,tukeyNoNA)
cldNoT11Tributary <- multcompLetters4(anovaNoT11Tributary, tukeyNoT11Tributary)
cldNoT11TributaryNoN4 <- multcompLetters4(anovaNoT11TributaryNoN4, tukeyNoT11TributaryNoN4)

print(cld)
print(cldNoNA)
print(cldNoT11Tributary)
print(cldNoT11TributaryNoN4)





#TS vs Depth including BlackNoN4 and ClearNoT11Tributary


anova(lm(AdjustedTS ~ WaterTypeAS*WaterDepth,
         data = dat))

anova(lm(AdjustedTS ~ WaterBody*WaterDepth,
         data = dat))



#Subsets with only one group
datAmazonMain <- subset(dat,WaterTypeAS %in% c("A"))
datBlack <- subset(dat,WaterTypeAS %in% c("N"))
datClear <- subset(dat,WaterTypeAS %in% c("T"))
datWhite <- subset(dat,WaterTypeAS %in% c("S"))
datBlackNoN4 <- subset(datNoN4,WaterTypeAS %in% c("NNoN4"))
datClearNoT11Tributary <- subset(datNoT11Tributary,WaterTypeAS %in% c("TNoT11R"))
datLake <- subset(dat,WaterBody %in% c("L"))
datMain <- subset(dat,WaterBody %in% c("M"))
datTrib <- subset(dat,WaterBody %in% c("R"))
datMainNoN4 <- subset(datNoN4,WaterBody %in% c("MNoN4"))
datTribNoT11R <- subset(datNoT11Tributary,WaterBody %in% c("RNoT11R"))


#For river system
SizeModelA <- lm(AdjustedTS ~ WaterDepth, # regression formula
                    data=datAmazonMain, na.action = na.omit) # data set

summary(SizeModelA)
sd(datAmazonMain$AdjustedTS)

SizeModelN <- lm(AdjustedTS ~ WaterDepth, # regression formula
                 data=datBlack, na.action = na.omit) # data set

summary(SizeModelN)
sd(datBlack$AdjustedTS)

SizeModelT <- lm(AdjustedTS ~ WaterDepth, # regression formula
                 data=datClear, na.action = na.omit) # data set

summary(SizeModelT)
sd(datClear$AdjustedTS)

SizeModelS <- lm(AdjustedTS ~ WaterDepth, # regression formula
                 data=datWhite, na.action = na.omit) # data set

summary(SizeModelS)
sd(datWhite$AdjustedTS)

SizeModelNNoN4 <- lm(AdjustedTS ~ WaterDepth, # regression formula
                 data=datBlackNoN4, na.action = na.omit) # data set

summary(SizeModelNNoN4)
sd(datBlackNoN4$AdjustedTS)

SizeModelTNoT11R <- lm(AdjustedTS ~ WaterDepth, # regression formula
                     data=datClearNoT11Tributary, na.action = na.omit) # data set

summary(SizeModelTNoT11R)
sd(datClearNoT11Tributary$AdjustedTS)


#For water body
SizeModelLake <- lm(AdjustedTS ~ WaterDepth, # regression formula
                      data=datLake, na.action = na.omit) # data set

summary(SizeModelLake)
sd(datLake$AdjustedTS)

SizeModelMain <- lm(AdjustedTS ~ WaterDepth, # regression formula
                    data=datMain, na.action = na.omit) # data set

summary(SizeModelMain)
sd(datMain$AdjustedTS)

SizeModelTrib <- lm(AdjustedTS ~ WaterDepth, # regression formula
                    data=datTrib, na.action = na.omit) # data set

summary(SizeModelTrib)
sd(datTrib$AdjustedTS)

SizeModelMainNoN4 <- lm(AdjustedTS ~ WaterDepth, # regression formula
                    data=datMainNoN4, na.action = na.omit) # data set

summary(SizeModelMainNoN4)
sd(datMainNoN4$AdjustedTS)

SizeModelTribNoT11R <- lm(AdjustedTS ~ WaterDepth, # regression formula
                        data=datTribNoT11R, na.action = na.omit) # data set

summary(SizeModelTribNoT11R)
sd(datTribNoT11R$AdjustedTS)








ndat <- rbind(dat, datNoN4, datNoT11Tributary)


#Relationsship between depth and nTS - by WaterTypeAS
ggplot(ndat, aes(x = AdjustedTS, y = WaterDepth)) +
  geom_point(color = "red", size = 1)+
  geom_smooth(method="lm", se=TRUE, data = subset(ndat, WaterTypeAS %in% c('S', 'N', 'A', 'T')), color = "black", linewidth = 1) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label), "*\", \"*", 
                                  after_stat(p.value.label), "*\".\"",
                                  sep = "")),
               label.x= "right", label.y= "bottom",
               formula = y~x, parse = TRUE, size = 6)+  
  facet_wrap(vars(WaterTypeAS), ncol = 2,
             labeller = labeller(WaterTypeAS = 
                                   c("N" ="                N               (a)",
                                     "S" = "S",
                                     "A" = "A", 
                                     "T" = "T", 
                                     "NNoN4" = "NNoN4",
                                     "TNoT11R" = "TNoT11R"))) +  
  coord_cartesian(ylim = c(40, 0)) +
  labs(x = "AdjustedTS (dB)", y = "Water depth (m)") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 22),
    strip.text = element_text(size = 24),
  )+

#Relationsship between depth and nTS - by WaterBody
ggplot(ndat, aes(x = AdjustedTS, y = WaterDepth)) +
  geom_point(color = "red", size = 1)+
  geom_smooth(method="lm", data = subset(ndat, WaterBody %in% c('L', 'M', 'R', 'MNoN4', 'RNoT11R')),
              color = "black", linewidth = 1) + 
  stat_poly_eq(data = subset(ndat, WaterBody %in% c('L', 'M', 'R', 'MNoN4', 'RNoT11R')), aes(label =  paste(after_stat(rr.label), "*\", \"*", 
                                  after_stat(p.value.label), "*\".\"",
                                  sep = "")),
               label.x= "right", label.y= "bottom",
               formula = y~x, parse = TRUE, size = 6)+  
  facet_wrap(vars(WaterBody), ncol = 2,
             labeller = labeller(WaterBody = 
                                   c("L" ="                L               (b)",
                                     "F" = "F",
                                     "M" = "M", 
                                     "R" = "R", 
                                     "MNoN4" = "MNoN4",
                                     "RNoT11R" = "RNoT11R"))) +  
  coord_cartesian(ylim = c(40, 0)) +
  labs(x = "AdjustedTS (dB)", y = "") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 22),
    strip.text = element_text(size = 24),
  )












################################
###   FLUX   ###################
###############################

#Histogram mm over volume and Flux
hist(dat$Vb_ml, breaks = 100,main="Histogram of bubble volume", xlab="volume(ml)")
hist(dat$Flux, breaks = 100,main="Histogram of bubble Flux", xlab="Flux")
datTSc30 <- subset(dat, MeanTSc<=-30)
hist(datTSc30$Vb_ml, breaks = 100,main="Histogram of bubble volume", xlab="volume(ml)")
hist(datTSc30$Flux, breaks = 100,main="Histogram of bubble Flux", xlab="Flux")

aggregate(x = dat$Vb_ml,                # Specify data column
          by = list(dat$LocationSplit),              # Specify group indicator
          FUN = mean) 

aggregate(x = datTSc30$Vb_ml,                # Specify data column
          by = list(datTSc30$LocationSplit),              # Specify group indicator
          FUN = mean) 


#NOrmalization check
qqnorm(dat$Flux)
qqline(dat$Flux)
ggqqplot(dat$Flux)
skewness(dat$Flux)
head(dat)


datNoT11TributaryNoN4 <-subset(datNoT11Tributary, LocationSplit!="N4")
levels(datNoT11TributaryNoN4$WaterTypeAS)[levels(datNoT11TributaryNoN4$WaterTypeAS)=='Black'] <- 'BlackNoN4'


#Log10 Fluxdata
dat$Flux10 <- log(dat$Flux)
datNoN4$Flux10 <- log(datNoN4$Flux)
datNoT11Tributary$Flux10 <- log(datNoT11Tributary$Flux)
datNoT11TributaryNoN4$Flux10 <- log(datNoT11TributaryNoN4$Flux)

#NOrmalization check logFlux
qqnorm(dat$Flux10)
qqline(dat$Flux10)
ggqqplot(dat$Flux10)
skewness(dat$Flux10)
hist(dat$Flux10, breaks = 100,main="Histogram of bubble Flux", xlab="Flux")


#Wilcoxon test differences between TS (pairwize - should work even on multiple comparisons)
dat %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoN4%>% wilcox_test(Flux ~ WaterTypeAS)
datNoT11Tributary%>% wilcox_test(Flux ~ WaterTypeAS)
datNoT11TributaryNoN4%>% wilcox_test(Flux ~ WaterTypeAS)







#Wilcoxon test differences between Flux between River systems(if I write pairwise before - same result- should work even on multiple comparisons)
#All methods with 
#Also tested Dunn.test - gives very similar results
#People recommened Holm (and actually "none", but that was discussed...)
dat %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "bonferroni")
dat %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "holm")
dat %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "hochberg")
dat %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "hommel")
dat %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "BH")
dat %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "BY")
dat %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "none")

dat %>% wilcox_test(Flux10 ~ WaterTypeAS, p.adjust.method = "bonferroni")
dat %>% wilcox_test(Flux10 ~ WaterTypeAS, p.adjust.method = "holm")
dat %>% wilcox_test(Flux10 ~ WaterTypeAS, p.adjust.method = "hochberg")
dat %>% wilcox_test(Flux10 ~ WaterTypeAS, p.adjust.method = "hommel")
dat %>% wilcox_test(Flux10 ~ WaterTypeAS, p.adjust.method = "BH")
dat %>% wilcox_test(Flux10 ~ WaterTypeAS, p.adjust.method = "BY")
dat %>% wilcox_test(Flux10 ~ WaterTypeAS, p.adjust.method = "none")

#ex dunntest
dunnTest(Flux ~ WaterTypeAS,data=dat,method="none")
dunnTest(Flux ~ WaterTypeAS,data=dat,method="holm")

dat  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "bonferroni")
dat  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "holm")
dat  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "BY")
dat  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "none")

datNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "holm")
datNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "hochberg")
datNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "hommel")
datNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "BH")
datNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "BY")
datNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "none")

datNoN4  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoN4  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "holm")
datNoN4  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "BY")
datNoN4  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "none")

datNoT11Tributary %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "holm")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "hochberg")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "hommel")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "BH")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "BY")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "none")

datNoT11Tributary  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoT11Tributary  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "holm")
datNoT11Tributary  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "BY")
datNoT11Tributary  %>% wilcox_effsize(Flux ~ WaterTypeAS, p.adjust.method = "none")


datNoT11TributaryNoN4 <-subset(datNoT11Tributary, LocationSplit!="N4")
levels(datNoT11TributaryNoN4$WaterTypeAS)[levels(datNoT11TributaryNoN4$WaterTypeAS)=='N'] <- 'NNoN4'

datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "bonferroni")
datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "holm")
datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "hochberg")
datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "hommel")
datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "BH")
datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterTypeAS, p.adjust.method = "BY")


#Wilcoxon test differences between Flux WaterBody (if I write pairwise before - same result- should work even on multiple comparisons)
#All methods with 
#Also tested Dunn.test - gives very similar results
#People recommened Holm (and actually "none", but that was discussed...)
dat %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "bonferroni")
dat %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "holm")
dat %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "hochberg")
dat %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "hommel")
dat %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "BH")
dat %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "BY")
dat %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "none")

#ex dunntest
dunnTest(Flux ~ WaterBody,data=dat,method="none")
dunnTest(Flux ~ WaterBody,data=dat,method="holm")


datNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "bonferroni")
datNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "holm")
datNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "hochberg")
datNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "hommel")
datNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "BH")
datNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "BY")
datNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "none")

datNoN4  %>% wilcox_effsize(Flux ~ WaterBody, p.adjust.method = "bonferroni")
datNoN4  %>% wilcox_effsize(Flux ~ WaterBody, p.adjust.method = "holm")
datNoN4  %>% wilcox_effsize(Flux ~ WaterBody, p.adjust.method = "BY")
datNoN4  %>% wilcox_effsize(Flux ~ WaterBody, p.adjust.method = "none")

datNoT11Tributary %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "bonferroni")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "holm")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "hochberg")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "hommel")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "BH")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "BY")
datNoT11Tributary %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "none")

datNoT11Tributary  %>% wilcox_effsize(Flux ~ WaterBody, p.adjust.method = "bonferroni")
datNoT11Tributary  %>% wilcox_effsize(Flux ~ WaterBody, p.adjust.method = "holm")
datNoT11Tributary  %>% wilcox_effsize(Flux ~ WaterBody, p.adjust.method = "BY")
datNoT11Tributary  %>% wilcox_effsize(Flux ~ WaterBody, p.adjust.method = "none")


datNoT11TributaryNoN4 <-subset(datNoT11Tributary, LocationSplit!="N4")
levels(datNoT11TributaryNoN4$WaterBody)[levels(datNoT11TributaryNoN4$WaterBody)=='M'] <- 'MNoN4'

datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "bonferroni")
datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "holm")
datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "hochberg")
datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "hommel")
datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "BH")
datNoT11TributaryNoN4 %>% wilcox_test(Flux ~ WaterBody, p.adjust.method = "BY")



give.n <- function(x){
  return(c(y = max(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

#Flux (per bubble spot) per Water Type including BlackNoN4 and ClearNoT11Tributary
ggplot(dat, aes(x=WaterTypeAS, y=Flux)) +   
  geom_boxplot(width=0.4, position=position_dodge(),fill="lightgrey")+
  geom_boxplot(data=datNoN4, aes(x=WaterTypeAS, y=Flux), width=0.4, position=position_dodge(),fill="lightgrey")+
  geom_boxplot(data=datNoT11Tributary, aes(x=WaterTypeAS, y=Flux), width=0.4, position=position_dodge(),fill="lightgrey")+
  xlab("River system") + ylab("Flux (ml m^-2 d^-1)")+
  scale_y_log10() +
  stat_summary(fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoN4,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoT11Tributary,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(size= 6, fun.data = give.n, geom = "text") +
  stat_summary(data= datNoN4,size= 6, fun.data = give.n, geom = "text") +
  stat_summary(data= datNoT11Tributary,size= 6, fun.data = give.n, geom = "text") +
  theme_hc() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) +
  theme(plot.title = element_text(size = 28)) +
  geom_text(size = 8, aes(x=c(Inf),y=c(Inf),hjust=c(2),vjust=c(1), label=c("(c)")))+
  #coord_cartesian(ylim = c(0, 50000))+
  scale_x_discrete(labels = c("A(bc)", "N(a)", "NNoN4(c)", "S(ab)", "T(abc)", "TNoT11T(d)")) +

#Flux (per bubble spot) per Water Type including BlackNoN4 and ClearNoT11Tributary
ggplot(dat, aes(x=WaterBody, y=Flux)) +   
  geom_boxplot(width=0.4, position=position_dodge(),fill="lightgrey")+
  geom_boxplot(data=datNoN4, aes(x=WaterBody, y=Flux), width=0.4, position=position_dodge(),fill="lightgrey")+
  geom_boxplot(data=datNoT11Tributary, aes(x=WaterBody, y=Flux), width=0.4, position=position_dodge(),fill="lightgrey")+
  xlab("Waterbody type") + ylab("")+
  scale_y_log10() +
  stat_summary(fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoN4,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoT11Tributary,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(size= 6, fun.data = give.n, geom = "text") +
  stat_summary(data= datNoN4,size= 6, fun.data = give.n, geom = "text") +
  stat_summary(data= datNoT11Tributary,size= 6, fun.data = give.n, geom = "text") +
  theme_hc() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) +
  theme(plot.title = element_text(size = 28)) +
  geom_text(size = 8, aes(x=c(Inf),y=c(Inf),hjust=c(2),vjust=c(1), label=c("(d)")))+
  #coord_cartesian(ylim = c(0, 50000))+
  scale_x_discrete(labels = c("F(abc)", "L(a)", "M(a)", "MNoN4(b)", "R(a)", "RNoT11R(c)"))




#Flux (per bubble spot) per Water Type including BlackNoN4 and ClearNoT11Tributary
ggplot(dat, aes(x=LocationSplit, y=Flux)) +   
  geom_boxplot(width=0.4, position=position_dodge(),fill="lightgrey")+
  geom_boxplot(data=datNoN4, aes(x=LocationSplit, y=Flux), width=0.4, position=position_dodge(),fill="lightgrey")+
  geom_boxplot(data=datNoT11Tributary, aes(x=LocationSplit, y=Flux), width=0.4, position=position_dodge(),fill="lightgrey")+
  xlab("Location") + ylab("Flux (ml m^-2 d^-1)")+
  scale_y_log10() +
  stat_summary(fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoN4,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(data= datNoT11Tributary,fun = mean, color = "steelblue", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 5,
               show.legend = FALSE) +
  stat_summary(size= 6, fun.data = give.n, geom = "text") +
  stat_summary(data= datNoN4,size= 6, fun.data = give.n, geom = "text") +
  stat_summary(data= datNoT11Tributary,size= 6, fun.data = give.n, geom = "text") +
  theme_hc() +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) +
  theme(plot.title = element_text(size = 28))




#Means...
dat %>% 
  group_by(WaterTypeAS) %>% 
  dplyr::summarise(Flux10_means = mean(Flux10, na.rm = TRUE)) 
dat %>% 
  group_by(WaterTypeAS) %>% 
  dplyr::summarise(Flux10_means = median(Flux10, na.rm = TRUE)) 

dat %>% 
  group_by(WaterBody) %>% 
  dplyr::summarise(Flux10_means = mean(Flux10, na.rm = TRUE)) 
dat %>% 
  group_by(WaterBody) %>% 
  dplyr::summarise(Flux10_means = median(Flux10, na.rm = TRUE))  


#If compare means also... WatertypeAS
anova <- aov(Flux10 ~ WaterTypeAS, data = dat)
anovaNoNA <- aov(Flux10 ~ WaterTypeAS, data = datNoN4)
anovaNoT11Tributary<- aov(Flux10 ~ WaterTypeAS, data = datNoT11Tributary)
anovaNoT11TributaryNoN4<- aov(Flux10 ~ WaterTypeAS, data = datNoT11TributaryNoN4)

tukey <- TukeyHSD(anova)
tukeyNoNA <- TukeyHSD(anovaNoNA)
tukeyNoT11Tributary <- TukeyHSD(anovaNoT11Tributary)
tukeyNoT11TributaryNoN4 <- TukeyHSD(anovaNoT11TributaryNoN4)

cld <- multcompLetters4(anova, tukey)
cldNoNA <- multcompLetters4(anovaNoNA,tukeyNoNA)
cldNoT11Tributary <- multcompLetters4(anovaNoT11Tributary, tukeyNoT11Tributary)
cldNoT11TributaryNoN4 <- multcompLetters4(anovaNoT11TributaryNoN4, tukeyNoT11TributaryNoN4)

print(cld)
print(cldNoNA)
print(cldNoT11Tributary)
print(cldNoT11TributaryNoN4)

datNoT11TributaryNoN4 <-subset(datNoT11Tributary, LocationSplit!="N4")
levels(datNoT11TributaryNoN4$WaterTypeAS)[levels(datNoT11TributaryNoN4$WaterTypeAS)=='N'] <- 'NNoN4'

datNoT11TributaryNoN4 <-subset(datNoT11Tributary, LocationSplit!="N4")
levels(datNoT11TributaryNoN4$WaterBody)[levels(datNoT11TributaryNoN4$WaterBody)=='M'] <- 'MNoN4'

#If compare means also... WaterBody
anova <- aov(Flux10 ~ WaterBody, data = dat)
anovaNoNA <- aov(Flux10 ~ WaterBody, data = datNoN4)
anovaNoT11Tributary<- aov(Flux10 ~ WaterBody, data = datNoT11Tributary)
anovaNoT11TributaryNoN4<- aov(Flux10 ~ WaterBody, data = datNoT11TributaryNoN4)

tukey <- TukeyHSD(anova)
tukeyNoNA <- TukeyHSD(anovaNoNA)
tukeyNoT11Tributary <- TukeyHSD(anovaNoT11Tributary)
tukeyNoT11TributaryNoN4 <- TukeyHSD(anovaNoT11TributaryNoN4)

cld <- multcompLetters4(anova, tukey)
cldNoNA <- multcompLetters4(anovaNoNA,tukeyNoNA)
cldNoT11Tributary <- multcompLetters4(anovaNoT11Tributary, tukeyNoT11Tributary)
cldNoT11TributaryNoN4 <- multcompLetters4(anovaNoT11TributaryNoN4, tukeyNoT11TributaryNoN4)

print(cld)
print(cldNoNA)
print(cldNoT11Tributary)
print(cldNoT11TributaryNoN4)




#Flux vs Depth including BlackNoN4 and ClearNoT11Tributary

anova(lm(Flux10 ~ WaterTypeAS*WaterDepth,
         data = dat))

anova(lm(Flux10 ~ WaterBody*WaterDepth,
         data = dat))



#Subsets with only one group
datAmazonMain <- subset(dat,WaterTypeAS %in% c("A"))
datBlack <- subset(dat,WaterTypeAS %in% c("N"))
datClear <- subset(dat,WaterTypeAS %in% c("T"))
datWhite <- subset(dat,WaterTypeAS %in% c("S"))
datBlackNoN4 <- subset(datNoN4,WaterTypeAS %in% c("NNoN4"))
datClearNoT11Tributary <- subset(datNoT11Tributary,WaterTypeAS %in% c("TNoT11R"))
datLake <- subset(dat,WaterBody %in% c("L"))
datMain <- subset(dat,WaterBody %in% c("M"))
datTrib <- subset(dat,WaterBody %in% c("R"))
datMainNoN4 <- subset(datNoN4,WaterBody %in% c("MNoN4"))
datTribNoT11R <- subset(datNoT11Tributary,WaterBody %in% c("RNoT11R"))


#For river system
SizeModelA <- lm(Flux10 ~ WaterDepth, # regression formula
                 data=datAmazonMain, na.action = na.omit) # data set

summary(SizeModelA)

SizeModelN <- lm(Flux10 ~ WaterDepth, # regression formula
                 data=datBlack, na.action = na.omit) # data set

summary(SizeModelN)

SizeModelT <- lm(Flux10 ~ WaterDepth, # regression formula
                 data=datClear, na.action = na.omit) # data set

summary(SizeModelT)

SizeModelS <- lm(Flux10 ~ WaterDepth, # regression formula
                 data=datWhite, na.action = na.omit) # data set

summary(SizeModelS)

SizeModelNNoN4 <- lm(Flux10 ~ WaterDepth, # regression formula
                     data=datBlackNoN4, na.action = na.omit) # data set

summary(SizeModelNNoN4)

SizeModelTNoT11R <- lm(Flux10 ~ WaterDepth, # regression formula
                       data=datClearNoT11Tributary, na.action = na.omit) # data set

summary(SizeModelTNoT11R)


#For water body
SizeModelLake <- lm(Flux10 ~ WaterDepth, # regression formula
                    data=datLake, na.action = na.omit) # data set

summary(SizeModelLake)

SizeModelMain <- lm(Flux10 ~ WaterDepth, # regression formula
                    data=datMain, na.action = na.omit) # data set

summary(SizeModelMain)

SizeModelTrib <- lm(Flux10 ~ WaterDepth, # regression formula
                    data=datTrib, na.action = na.omit) # data set

summary(SizeModelTrib)

SizeModelMainNoN4 <- lm(Flux10 ~ WaterDepth, # regression formula
                        data=datMainNoN4, na.action = na.omit) # data set

summary(SizeModelMainNoN4)

SizeModelTribNoT11R <- lm(Flux10 ~ WaterDepth, # regression formula
                          data=datTribNoT11R, na.action = na.omit) # data set

summary(SizeModelTribNoT11R)



ndat <- rbind(dat, datNoN4, datNoT11Tributary)


#Relationsship between depth and nTS - by WaterTypeAS
ggplot(ndat, aes(x = Flux, y = WaterDepth)) +
  geom_point(color = "red", size = 1)+
  geom_smooth(method="lm", , data = subset(ndat, WaterTypeAS %in% c('S', 'N', 'NNoN4', 'A', 'T')), color = "black", linewidth = 1) + 
  stat_poly_eq(aes(label =  paste(after_stat(rr.label), "*\", \"*", 
                                  after_stat(p.value.label), "*\".\"",
                                  sep = "")),
               label.x= "right", label.y= "bottom",
               formula = y~x, parse = TRUE, size = 5)+  
  facet_wrap(vars(WaterTypeAS), ncol = 2,
             labeller = labeller(WaterTypeAS = 
                                   c("N" ="                N               (c)",
                                     "S" = "S",
                                     "A" = "A", 
                                     "T" = "T", 
                                     "NNoN4" = "NNoN4",
                                     "TNoT11R" = "TNoT11R"))) +  
  scale_x_log10() +
  coord_cartesian(ylim = c(40, 0)) +
  labs(x = "Flux (dB)", y = "Water depth (m)") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 22),
    strip.text = element_text(size = 24),
  )+
  
  
  #Relationsship between depth and nTS - by WaterBody
  ggplot(ndat, aes(x = Flux, y = WaterDepth)) +
  geom_point(color = "red", size = 1)+
  geom_smooth(method="lm", data = subset(ndat, WaterBody %in% c('L', 'M', 'R', 'MNoN4', 'RNoT11R')),
              color = "black", linewidth = 1) + 
  stat_poly_eq(data = subset(ndat, WaterBody %in% c('L', 'M', 'R', 'MNoN4', 'RNoT11R')), aes(label =  paste(after_stat(rr.label), "*\", \"*", 
                                                                                                            after_stat(p.value.label), "*\".\"",
                                                                                                            sep = "")),
               label.x= "right", label.y= "bottom",
               formula = y~x, parse = TRUE, size = 5)+  
  facet_wrap(vars(WaterBody), ncol = 2,
             labeller = labeller(WaterBody = 
                                   c("L" ="                L               (d)",
                                     "F" = "F",
                                     "M" = "M", 
                                     "R" = "R", 
                                     "MNoN4" = "MNoN4",
                                     "RNoT11R" = "RNoT11R"))) +  
  scale_x_log10() +
  coord_cartesian(ylim = c(40, 0)) +
  labs(x = "Flux (dB)", y = "") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 22),
    strip.text = element_text(size = 24),
  )












