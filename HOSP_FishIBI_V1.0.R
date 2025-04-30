install.packages("installr")
library(installr)
updateR()

                                                        ##HOSP Fish IBI Analysis##

#After installing packages on your PC once, you do not need to run the packages again, but you do need to run
#the libraries every time you reopen RStudio.
install.packages('vegan')
library(vegan)
install.packages("tidyverse")
library(tidyverse)                   
                                                  #####################################
                                                  ##Start here after loading packages##
                                                  #####################################

                                                          #######STEP 1#######

#Copy the file address for the exported Database spreadsheets here. Make sure they are saved as .csv files.
Fish<-read.csv("HOSP_qsel_Export_FishCommSppIndiv_EF.csv")  #This is the raw fish species data.
Effort<-read.csv("HOSP_qsel_Export_FishCommSampling_EF.csv") #This contains sampling effort.
FGR<-read.csv("tlu_TaxaSpecies.csv") #This contains FG and reproductive classes.

#Before we get started, let's change any -999 values or blanks to NA's in our dataframe.
Fish[Fish == "-999"] <- as.numeric("NA")

Fish[Fish==""]<-NA

#The next several lines add a Tag column to both dataframes which will allow stats to be organized and averaged 
#according to specific reaches across years.
new<-data.frame(cbind(Fish$Season, Fish$LocationID))
new$Tag<-apply(new,1,paste,collapse="")
Fish<-mutate(Fish, Tag = new$Tag)
new2<-data.frame(cbind(Effort$Season, Effort$LocationID))
new2$Tag<-apply(new2,1,paste,collapse="")
Effort<-mutate(Effort, Tag = new2$Tag)

#Now we will remove unnecessary columns from the FGR and Effort dataframes so we can bind just the relevant columns to the main Fish dataframe.
FGR <- FGR[,c(1,3,6,7,8)]
Effort <- Effort[,c(6,12)]

#Now we'll combine the total effort in seconds for all sampling bouts throughout each reach and convert the total seconds to minutes.
Effort <- Effort %>% 
  group_by(Tag) %>% 
  summarise(minutes = sum(SamplingEffort_sec)/60)

#Now we'll merge the FGR dataframe.
Fish = merge(Fish, FGR, all=T)

#Create a TotalFish column and merge it as well.
TotFish <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Tag) %>% 
  summarise(TotalFish = sum(NumObs))

#Let's condense the Fish dataframe into just the columns we need to produce IBI scores now.
Fish <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Season, LocationID, Tag, TaxonCode, FamilyName, TaxaFishNumber, AnomalyID, NumObs, ToleranceCode, ReproductiveClassification, 
           TrophicClassification) %>%
  summarise()

#Time for the initial summaries of the 8 metrics in this IBI.
#Metric 1 Species Richness
Summary1.1 <- Fish %>% 
  group_by(Tag) %>% 
  summarise(TotalSpp = n_distinct(TaxonCode))

#Metric 2 Number of Cyprinid (minnow) Species 
Fish2 <- Fish %>%
  group_by(Tag, TaxonCode, FamilyName, ToleranceCode) %>%
  summarise()
Summary1.2 <- Fish2 %>% 
  group_by(Tag) %>% 
  summarise(CyprinidSpp = length(which(FamilyName=='Cyprinidae')))


#Metric 3 Number of Sensitive Species
Summary1.3 <- Fish2 %>% 
  group_by(Tag) %>% 
  summarise(TCI = length(which(ToleranceCode=='I')))

#Metric 4 Percent Individuals as Green Sunfish
CompGS = aggregate(NumObs~TaxonCode + Tag, subset(Fish, TaxonCode=='lepcya'), sum)
CompGS = merge(TotFish, CompGS, all=T)
CompGS$CompGS = ((CompGS$NumObs/CompGS$TotalFish)*100)
Summary1.4 <- CompGS %>% 
  filter(!is.na(Tag)) %>%
  group_by(Tag, CompGS) %>% 
  summarise()


#Metric 5 Ratio of Generalist to Specialists  

#Numerator calculates number of Generalists, ie TrophicClassification of DAI, DAIP, and AIP
Numerator = aggregate(NumObs~TrophicClassification + Tag, subset(Fish, TrophicClassification=='DAI' | TrophicClassification=='DAIP' |
                                                                 TrophicClassification=='AIP'), sum)
Numerator <- Numerator %>% 
  group_by(Tag) %>% 
  summarise(NumObsN = sum(NumObs))

#Denominator calculates number of specialists, ie TrophicClassification NOT DAI, DAIP, AIP, U(unknown)
Denominator = aggregate(NumObs~TrophicClassification + Tag, subset(Fish, !(TrophicClassification %in% c('DAI', 'DAIP', 'AIP','U'))), sum)

Denominator <- Denominator %>% 
  group_by(Tag) %>% 
  summarise(NumObsD = sum(NumObs))

Ratio = merge(Numerator, Denominator, all=T)
Ratio[is.na(Ratio)] <- 0

Summary1.5 <- Ratio %>% 
  group_by(Tag) %>% 
  summarise(FGratio = NumObsN/NumObsD)

#Metric 6  Percent Individuals as Top Carnivores
CompTP = aggregate(NumObs~TrophicClassification + Tag, subset(Fish, TrophicClassification=='TP'), sum)
CompTP = merge(TotFish, CompTP, all=T)
CompTP[is.na(CompTP)] <- 0
CompTP$CompTP = ((CompTP$NumObs/CompTP$TotalFish)*100)

Summary1.6 <- CompTP %>% 
  filter(!is.na(Tag)) %>%
  group_by(Tag, CompTP) %>% 
  summarise()

#Metric 7 Catch per Unit Effort (Minute) of sampling time
Fish3 = merge(TotFish, Effort, all=T)
Summary1.7 <- Fish3 %>% 
  group_by(Tag) %>% 
  summarise(CPUE = TotalFish/minutes)

#Metric 8 Percent Individuals with Disease or Anomaly
CompAnom = aggregate(NumObs~AnomalyID + Tag, subset(Fish, AnomalyID!='N'), sum)
CompAnom = merge(TotFish, CompAnom, all=T)
Summary1.8 <- CompAnom %>% 
  group_by(Tag, TotalFish) %>% 
  summarise(NumObs = sum(NumObs))
Summary1.8[is.na(Summary1.8)] <- 0
Summary1.8$CompAnom = ((Summary1.8$NumObs/Summary1.8$TotalFish)*100)
Summary1.8[,2:3] <- NULL

#The following code creates a Raw Metric summary for export.
MRSummary <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Season, LocationID, Tag) %>%
  summarise()

#MRSumary has the raw metric values for the 8 metrics for each LocationID and Season (year). THESE RAW METRIC VALUES SHOULD BE INCLUDED IN SUMMARY REPORTS ALONG WITH IBI SCORE AND IBI RATING
MRSummary <- cbind(MRSummary, NumSpp = Summary1.1$TotalSpp, NumCyprinSp = Summary1.2$CyprinidSpp, NumSensSp = Summary1.3$TCI, PerGreenSun = Summary1.4$CompGS, 
                   GeneralToSpecial = Summary1.5$FGratio, PerTopCarn = Summary1.6$CompTP, CPUE = Summary1.7$CPUE, PerDisease = Summary1.8$CompAnom)
MRSummary[,3] <- NULL

                                                          #######STEP 2#######

#Now we'll calculate the final Metric Scores to produce an IBI score.

#Metric 1 Score for Total number of species
Summary1.1 <- Summary1.1 %>% 
  mutate(UT = ifelse(TotalSpp>=13, 5, 0),
         MT = ifelse(TotalSpp<13 & TotalSpp>6, 3, 0),
         LT = ifelse(TotalSpp<=6, 1, 0))

Summary1.1 <- Summary1.1 %>%
  group_by(Tag) %>%
  summarise(MS1 = sum(UT,MT,LT))

#Metric 2 Score for Number of Cyprinidae species
Summary1.2 <- Summary1.2 %>% 
  mutate(UT = ifelse(CyprinidSpp>=4, 5, 0),
         MT = ifelse(CyprinidSpp<4 & CyprinidSpp>1, 3, 0),
         LT = ifelse(CyprinidSpp<=1, 1, 0))

Summary1.2 <- Summary1.2 %>%
  group_by(Tag) %>%
  summarise(MS2 = sum(UT,MT,LT))

#Metric 3 Score for Number of Sensitive Species
Summary1.3 <- Summary1.3 %>% 
  mutate(UT = ifelse(TCI>=5, 5, 0),
         MT = ifelse(TCI<5 & TCI>2, 3, 0),
         LT = ifelse(TCI<=2, 1, 0))
Summary1.3 <- Summary1.3 %>%
  group_by(Tag) %>%
  summarise(MS3 = sum(UT,MT,LT))

#Metric 4 Score for Percent Individuals as Green Sunfish
Summary1.4 <- Summary1.4 %>% 
  mutate(UT = ifelse(CompGS<=2.1, 5, 0),
         MT = ifelse(CompGS>2.1 & CompGS<4.4, 3, 0),
         LT = ifelse(CompGS>=4.4, 1, 0))
Summary1.4 <- Summary1.4 %>%
  group_by(Tag) %>%
  summarise(MS4 = sum(UT,MT,LT))

#Metric 5 Score for Ratio of Generalist to Specialist feeders
Summary1.5 <- Summary1.5 %>% 
  mutate(UT = ifelse(FGratio<3, 5, 0),
         MT = ifelse(FGratio>=3 & FGratio<=6, 3, 0),
         LT = ifelse(FGratio>6, 1, 0))
Summary1.5 <- Summary1.5 %>%
  group_by(Tag) %>%
  summarise(MS5 = sum(UT,MT,LT))

#Metric 6 Score for Percent Individuals as Top Carnivores (piscivores)
Summary1.6 <- Summary1.6 %>% 
  mutate(UT = ifelse(CompTP>=2.3, 5, 0),
         MT = ifelse(CompTP<2.3 & CompTP>1.1, 3, 0),
         LT = ifelse(CompTP<=1.1, 1, 0))
Summary1.6 <- Summary1.6 %>%
  group_by(Tag) %>%
  summarise(MS6 = sum(UT,MT,LT))

#Metric 7 Score for Catch per Unit Effort(ie Minute)
Summary1.7 <- Summary1.7 %>% 
  mutate(UT = ifelse(CPUE>20, 3, 0),
         MT = ifelse(CPUE<=20 & CPUE>=6, 5, 0),
         LT = ifelse(CPUE<6, 1, 0))
Summary1.7 <- Summary1.7 %>%
  group_by(Tag) %>%
  summarise(MS7 = sum(UT,MT,LT))

#Metric 8 Score for Percent of Individuals with Disease or Anomaly
Summary1.8 <- Summary1.8 %>% 
  mutate(UT = ifelse(CompAnom<=1, 5, 0),
         MT = ifelse(CompAnom>1 & CompAnom<=3.8, 3, 0),
         LT = ifelse(CompAnom>3.8, 1, 0))
Summary1.8 <- Summary1.8 %>%
  group_by(Tag) %>%
  summarise(MS8 = sum(UT,MT,LT))

                                                          #######STEP 3#######

#Finally, we'll bind all of the metric scores to the same dataframe and sum them to come up with a final IBI score.
Fish <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Season, LocationID, Tag) %>%
  summarise()

Fish <- cbind(Fish, ScoreNumSpp = Summary1.1$MS1, ScoreNumCypSp = Summary1.2$MS2, ScoreNumSensSp = Summary1.3$MS3, ScorePerGreenSun = Summary1.4$MS4, ScoreGeneralToSpecial = Summary1.5$MS5, 
             ScorePerTopCarn = Summary1.6$MS6, ScoreCPUE = Summary1.7$MS7, ScorePerDisease = Summary1.8$MS8)

MSSummary <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Season, LocationID, ScoreNumSpp, ScoreNumCypSp, ScoreNumSensSp, ScorePerGreenSun, ScoreGeneralToSpecial, ScorePerTopCarn, ScoreCPUE, ScorePerDisease) %>%
  summarise()



FinalSummary <- Fish %>%
  group_by(Season, LocationID) %>%
  summarise(IBI = (ScoreNumSpp+ScoreNumCypSp+ScoreNumSensSp+ScorePerGreenSun+ScoreGeneralToSpecial+ScorePerTopCarn+ScoreCPUE+ScorePerDisease))


#Gives the IBI score a rating from Excellent to Poor

FinalSummary <- FinalSummary %>% 
    mutate(Excellent = ifelse(IBI>=33,'Excellent', ''),
         Good = ifelse(IBI>=25 & IBI<33,'Good',''),
         Fair = ifelse(IBI>=17 & IBI<25,'Fair', ''),
         Poor = ifelse(IBI<17,'Poor', ''))
FinalSummary$IBIrating <-paste(FinalSummary$Excellent, FinalSummary$Good, FinalSummary$Fair, FinalSummary$Poor)

IBIScoreAndRating <-FinalSummary[,c(1,2,3,8)]
  

