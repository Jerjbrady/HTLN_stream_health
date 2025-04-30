install.packages("installr")
library(installr)
updateR()

                                                          ##EFMO Fish IBI Analysis##

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
Fish<-read.csv("C:/HopesFiles/HTLN/AquaticShares/FishRscripts/EFMOoutput/EFMO_FishCommSppIndiv_EF.csv")  #This is the raw fish species data.
FGR<-read.csv("C:/HopesFiles/HTLN/AquaticShares/FishRscripts/tlu_TaxaSpecies.csv") #This contains FG and reproductive classes.

#Before we get started, let's change any -999 values or blanks to NA's in our dataframe.
Fish[Fish == "-999"] <- as.numeric("NA")

Fish[Fish==""]<-NA

#The next several lines add a Tag column to both dataframes which will allow stats to be organized and averaged 
#according to specific reaches across years.
new<-data.frame(cbind(Fish$Season, Fish$LocationID))
new$Tag<-apply(new,1,paste,collapse="")
Fish<-mutate(Fish, Tag = new$Tag)

#Now we will remove unnecessary columns from the FGR and Effort dataframes so we can bind just the relevant columns to the main Fish dataframe.
FGR <- FGR[,c(1,3,6,8)]

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
  group_by(Season, LocationID, Tag, TaxonCode, FamilyName, TaxaFishNumber, NumObs, ToleranceCode, TrophicClassification) %>%
  summarise()


#Calculating the raw values of the 5 metrics in this IBI.

#Metric 1 Percent Intolerant Species
Fish2 <- Fish %>%
  group_by(Tag, TaxonCode, FamilyName, ToleranceCode) %>%
  summarise()
Summary1.1 <- Fish2 %>% 
  group_by(Tag) %>% 
  summarise(TCI = length(which(ToleranceCode=='I')))

#Metric 2 Percent Tolerant Individuals
CompT = aggregate(NumObs~ToleranceCode + Tag, subset(Fish, ToleranceCode=='T'), sum)
CompT = merge(TotFish, CompT, all=T)
Summary1.2 <- CompT %>% 
  group_by(Tag) %>% 
  summarise(CompT = (NumObs/TotalFish)*100)

#Metric 3 Percent Top Carnivore/Predator Individuals
CompTP = aggregate(NumObs~TrophicClassification + Tag, subset(Fish, TrophicClassification=='TP'), sum)
CompTP = merge(TotFish, CompTP, all=T)
Summary1.3 <- CompTP %>% 
  group_by(Tag) %>% 
  summarise(CompTP = (NumObs/TotalFish)*100)
Summary1.3[is.na(Summary1.3)] <- 0

#Metric 4 Percent Cool or Coldwater Individuals
CompSC = aggregate(NumObs~TaxonCode + Tag, subset(Fish, TaxonCode=='ichfos' | TaxonCode=='ichgag' | TaxonCode=='letapp' | TaxonCode=='clielo' |
                                                  TaxonCode=='hybhan' | TaxonCode=='chreos' | TaxonCode=='chrneo' | TaxonCode=='esomas' | 
                                                  TaxonCode=='lotlot' | TaxonCode=='culinc' | TaxonCode=='cotbai' | TaxonCode=='saltru' | 
                                                  TaxonCode=='oncmyk'), sum)
CompSC <- CompSC %>% 
group_by(Tag) %>% 
summarise(NumObs = sum(NumObs))
CompSC = merge(TotFish, CompSC, all=T)
CompSC$CompSC = ((CompSC$NumObs/CompSC$TotalFish)*100)
Summary1.4 <- CompSC %>% 
  group_by(Tag, CompSC) %>% 
  summarise()

#Metric 5 Percent Salmonids that are Brook Trout
CompSd = aggregate(NumObs~TaxonCode + Tag, subset(Fish, TaxonCode=='salfon' | TaxonCode=='saltru' | TaxonCode=='oncmyk'), sum)
CompSd <- CompSd %>% 
  group_by(Tag) %>% 
  summarise(CompSd = sum(NumObs))
CompSn <- Fish %>%
  mutate(CompSn = ifelse(TaxonCode=='salfon', sum(NumObs), 0))
CompSn <- CompSn %>% 
  group_by(Tag) %>% 
  summarise(CompSn = sum(CompSn))
CompS = merge(CompSn, CompSd, all=T)
CompS[is.na(CompS)] <- 0
CompS$CompS = ((CompS$CompSn/CompS$CompSd)*100)
Summary1.5 <- CompS %>% 
  group_by(Tag, CompS) %>% 
  summarise()
Summary1.5[is.na(Summary1.5)] <- 0

#The following code creates a Raw Metric summary for export.
MRSummary <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Season, LocationID, Tag) %>%
  summarise()

MRSummary <- cbind(MRSummary, NumIntolSp = Summary1.1$TCI, PerTol = Summary1.2$CompT, PerTopCarn = Summary1.3$CompTP, PerColdwater = Summary1.4$CompSC, 
                   PerBrookTrout = Summary1.5$CompS)
MRSummary[,3] <- NULL

write.csv(MRSummary, "C:/HopesFiles/HTLN/AquaticShares/FishRscripts/EFMOoutput/EFMO_fish_RawMetricSummary.csv") #You will need to change the 
                                                                                                               #name and address accordingly.

                                                            #######STEP 2#######

#Calculation of metric scores from the raw metric values to produce an IBI score.


#Metric Score 1 Number of Intolerant Species
Summary1.1 <- Summary1.1 %>% 
  mutate(UT = ifelse(TCI>=2, 20, 0),
         MT = if_else(TCI==1, 10, 0),
         LT = if_else(TCI==0, 0, 0))
Summary1.1 <- Summary1.1 %>%
  group_by(Tag) %>%
  summarise(ScoreNumIntolSp = sum(UT,MT,LT))

#Metric 2 Percent Tolerant Individuals
Summary1.2 <- Summary1.2 %>% 
  mutate(UT = ifelse(CompT<=5, 20, 0),
         MT = ifelse(CompT>5 & CompT<23, 10, 0),
         LT = ifelse(CompT>=23, 0, 0))
Summary1.2 <- Summary1.2 %>%
  group_by(Tag) %>%
  summarise(ScorePerTol = sum(UT,MT,LT))

#Metric 3 Percent Top Carnivores
Summary1.3 <- Summary1.3 %>% 
  mutate(UT = ifelse(CompTP>=46, 20, 0),
         MT = ifelse(CompTP<46 & CompTP>14, 10, 0),
         LT = ifelse(CompTP<=14, 0, 0))
Summary1.3 <- Summary1.3 %>%
  group_by(Tag) %>%
  summarise(ScorePerTopCarn = sum(UT,MT,LT))

#Metric 4 Percent Cool and Coldwater Individuals
Summary1.4 <- Summary1.4 %>% 
  mutate(UT = ifelse(CompSC>=86, 20, 0),
         MT = ifelse(CompSC<86 & CompSC>42, 10, 0),
         LT = ifelse(CompSC<=42, 0, 0))
Summary1.4 <- Summary1.4 %>%
  group_by(Tag) %>%
  summarise(ScorePerColdwater = sum(UT,MT,LT))

#Metric 5 Percent of Salmonids that are Brook Trout
Summary1.5 <- Summary1.5 %>% 
  mutate(UT = ifelse(CompS>=96, 20, 0),
         MT = ifelse(CompS<96 & CompS>4, 10, 0),
         LT = ifelse(CompS<=4, 0, 0))
Summary1.5 <- Summary1.5 %>%
  group_by(Tag) %>%
  summarise(ScorePerBrownTrout = sum(UT,MT,LT))

                                                            #######STEP 3#######

#Bind all of the metric scores to the same dataframe and sum them to come up with a final IBI score.
Fish <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Season, LocationID, Tag) %>%
  summarise()


Fish <- cbind(Fish, ScoreNumIntolSp = Summary1.1$ScoreNumIntolSp, ScorePerTol = Summary1.2$ScorePerTol, ScorePerTopCarn = Summary1.3$ScorePerTopCarn, ScorePerColdwater = Summary1.4$ScorePerColdwater, ScorePerBrownTrout = Summary1.5$ScorePerBrownTrout)

MSSummary <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Season, LocationID, ScoreNumIntolSp, ScorePerTol, ScorePerTopCarn, ScorePerColdwater, ScorePerBrownTrout) %>%
  summarise()

write.csv(MSSummary, "C:/HopesFiles/HTLN/AquaticShares/FishRscripts/EFMOoutput/EFMO_fish_MetricScoreSummary.csv") #You will need to change the 
                                                                                                               #name and address accordingly.

FinalSummary <- Fish %>%
  group_by(Season, LocationID) %>%
  summarise(IBI = (ScoreNumIntolSp+ScorePerTol+ScorePerTopCarn+ScorePerColdwater+ScorePerBrownTrout))

#Lastly we'll save this summary as .csv file so that we can load them into a new R script that we'll use to make various plots.
write.csv(FinalSummary, "C:/HopesFiles/HTLN/AquaticShares/FishRscripts/EFMOoutput/EFMO_fish_IBI.csv") #You will need to change the 
                                                                                                            #name and address accordingly.
