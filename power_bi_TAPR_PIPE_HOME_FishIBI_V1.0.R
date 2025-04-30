#install.packages("installr")
#library(installr)
#updateR()

                                                  ##TAPR, PIPE and HOME Fish IBI Analysis##

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
fish_csvs <- list(
HOME = "HOME_FishCommSppIndiv_SE.csv",   #This is the raw fish species data.
#TAPR = "TAPR_FishCommSppIndiv_SE.csv",
PIPE = "PIPE_FishCommSppIndiv_SE.csv")


Reach_csvs <- list(
  HOME = "HOME_FishCommSampling_SE.csv",
  #TAPR = "TAPR_FishCommSampling_SE.csv", 
  PIPE = "PIPE_FishCommSampling_SE.csv") #This contains the lengths of reaches sampled.

FGR<-read.csv("tlu_TaxaSpecies.csv") #This contains FG, TC and reproductive classes.
#Now we will remove unnecessary columns from the FGR dataframe so we can bind just the relevant columns to the main Fish dataframe.

#Keeps TaxonCode, FamilyName, ScientificName, ToleranceCode, TrophicClassification
FGR <- FGR[,c(1,3,4,6,8)]

all_MRSummary <- list()
all_IBI <- list()

for (site in names(fish_csvs)){
tryCatch({
  message("Processing site: ", site)

Fish <- read.csv(fish_csvs[[site]])
ReachLength <- read.csv(Reach_csvs[[site]])
#Before we get started, let's change any -999 values or blanks to NA's in our dataframe.
Fish[Fish == "-999"] <- as.numeric("NA")

Fish[Fish==""]<-NA

#The next several lines add a Tag column to both dataframes which will allow stats to be organized and averaged 
#according to specific reaches across years.
new<-data.frame(cbind(Fish$Season, Fish$LocationID))
new$Tag<-apply(new,1,paste,collapse="")
Fish<-mutate(Fish, Tag = new$Tag)
new2<-data.frame(cbind(ReachLength$Season, ReachLength$LocationID))
new2$Tag<-apply(new2,1,paste,collapse="")
ReachLength<-mutate(ReachLength, Tag = new2$Tag)


#Let's get the total length sampled for all reaches.
ReachLength <- ReachLength %>%
  filter(!is.na(Tag)) %>%
  group_by(Tag) %>%
  summarise(Totlength = sum(SeineReachLength_m))


#Create a TotalFish dataframe to merge later.
TotFish <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Tag) %>% 
  summarise(TotalFish = sum(NumObs))

#Now we'll merge the FGR dataframe.
Fish = merge(Fish, FGR, all=T)

#Let's condense the Fish dataframe back to the orignial number of fish records
Fish <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Season, LocationID, Tag, TaxonCode, SiteNumber, AnomalyID, FamilyName, ScientificName, TaxaFishNumber, NumObs, ToleranceCode, TrophicClassification) %>%
  summarise()

#Deletes all fish not identified to species (i.e. unknowns) and deletes hybrids to get accurate species richness
Fish_MinusUnknowns = aggregate(NumObs~TaxonCode + Tag, subset(Fish, !(TaxonCode %in% c('cotspp', 'cypspp', 'cypunk', 'ethspp', 'ichamm',
                                                                                        'lapamm', 'lepcxh','lepcxm', 'lephxm', 'lepixm',
                                                                                       'lepmec','lepmxc','lepmxi','lepmxm', 'lepspp',
                                                                                       'luxspp','moxspp','nofish','notspp','ntuspp','pimpspp'))), sum)


#Time for the initial summaries of the 12 metrics in this IBI.
#Metric 1 Species Richness
Summary1.1 <- Fish_MinusUnknowns %>% 
  group_by(Tag) %>% 
  summarise(TotSp = n_distinct(TaxonCode))

#Metric 2 Number of Darter species
Fish2 <- Fish %>%
  group_by(Tag, TaxonCode, FamilyName, ToleranceCode) %>%
  summarise()
Summary1.2 <- Fish2 %>% 
  group_by(Tag) %>% 
  summarise(NumEP = length(which(TaxonCode=='ethcra' | TaxonCode=='etheux' | TaxonCode=='ethaut' | TaxonCode=='ethzon' | TaxonCode=='permac' |
                                 TaxonCode=='ethuni' | TaxonCode=='etheri' | TaxonCode=='ethfla' | TaxonCode=='perevi' | TaxonCode=='ethble' | 
                                 TaxonCode=='ethexi' | TaxonCode=='ethnig' | TaxonCode=='percap' | TaxonCode=='ethrad' | TaxonCode=='ethspe' |
                                 TaxonCode=='ethcae' | TaxonCode=='ethmih' | TaxonCode=='perfla' | TaxonCode=='ethjul')))

#Metric 3 Number of Sunfish Species
Summary1.3 <- Fish2 %>% 
  group_by(Tag) %>% 
  summarise(NumL = length(which(TaxonCode=='lepcya' | TaxonCode=='lepgul' | TaxonCode=='lephum' | TaxonCode=='lepmac' | TaxonCode=='lepmeg' |
                                   TaxonCode=='lepmic' | TaxonCode=='lepmin')))

#Metric 4 Number of Sucker species
Summary1.4 <- Fish2 %>% 
  group_by(Tag) %>% 
  summarise(NumCat = length(which(FamilyName=='Catostomidae')))

#Metric 5 Number of Intolerant species
Summary1.5 <- Fish2 %>% 
  group_by(Tag) %>% 
  summarise(TCI = length(which(ToleranceCode=='I')))

#Metric 6 Proportion of individuals as Green Sunfish
CompGS <- Fish %>%
  filter(TaxonCode=='lepcya') %>%
  group_by(Tag) %>% 
  summarise(NumGS =  sum(NumObs, na.rm = TRUE))

CompGS = merge(TotFish, CompGS, all=T)

Summary1.6 <- CompGS %>% 
  group_by(Tag) %>% 
  summarise(CompGS = (NumGS/TotalFish)*100)
Summary1.6[is.na(Summary1.6)] <- 0

#Metric 7 Proportion of individuals as omnivores
CompO <- Fish %>%
  filter(TrophicClassification=='DAI' | TrophicClassification=='DAIP' |
           TrophicClassification=='AIP' ) %>%
  group_by(Tag) %>% 
  summarise(NumO =  sum(NumObs, na.rm = TRUE))

CompO = merge(TotFish, CompO, all=T)

Summary1.7 <- CompO %>% 
  group_by(Tag) %>% 
  summarise(CompO = (NumO/TotalFish)*100)
Summary1.7[is.na(Summary1.7)] <- 0

#Metric 8 Proportion of individuals as insectivorous minnows
CompCI <- Fish %>%
  filter(FamilyName=='Cyprinidae' & TrophicClassification=='I') %>%
  group_by(Tag) %>% 
  summarise(NumCI =  sum(NumObs, na.rm = TRUE))

CompCI = merge(TotFish, CompCI, all=T)

Summary1.8 <- CompCI %>% 
  group_by(Tag) %>% 
  summarise(CompCI = (NumCI/TotalFish)*100)
Summary1.8[is.na(Summary1.8)] <- 0

#Metric 9 Proportion of individuals as top carnivores/piscivores
CompTP <- Fish %>%
  filter(TrophicClassification=='TP') %>%
  group_by(Tag) %>% 
  summarise(NumTP =  sum(NumObs, na.rm = TRUE))

CompTP = merge(TotFish, CompTP, all=T)
Summary1.9 <- CompTP %>% 
  group_by(Tag) %>% 
  summarise(CompTP = (NumTP/TotalFish)*100)
Summary1.9[is.na(Summary1.9)] <- 0

#Metric 10   Number of individuals in the sample - Catch per 50ft seine haul
#converts reach length from meters to feet and divides by 50 feet
ReachLengthFt <- ReachLength %>%  
  group_by(Tag) %>%
  summarise(TotLengthPer50ft = (Totlength*3.2808399)/50) 

CatchL = merge(TotFish, ReachLengthFt, all=T)

Summary1.10 <- CatchL %>% 
  group_by(Tag) %>% 
  summarise(CatchPer50Ft = (TotalFish/TotLengthPer50ft))
Summary1.10[is.na(Summary1.10)] <- 0

#Metric 11 proportion of individuals as hybrids
CompH <- Fish %>%
  filter(TaxonCode=='lepmxc' | TaxonCode=='lepmxi' | TaxonCode=='lepcxh' | TaxonCode=='lepmxm' |
           TaxonCode=='lepmec' | TaxonCode=='lephxm' | TaxonCode=='lepcxm' | TaxonCode=='lepixm') %>%
  group_by(Tag) %>% 
  summarise(NumH =  sum(NumObs, na.rm = TRUE))


CompH = merge(TotFish, CompH, all=T)

Summary1.11 <- CompH %>% 
  group_by(Tag) %>% 
  summarise(CompH = (NumH/TotalFish)*100)
Summary1.11[is.na(Summary1.11)] <- 0

#Metric 12 Proportion of individuals with anomalies
CompAnom <- Fish %>%
  filter(AnomalyID != 'N') %>%
  group_by(Tag) %>% 
  summarise(NumAnom =  sum(NumObs, na.rm = TRUE))


CompAnom = merge(TotFish, CompAnom, all=T)

Summary1.12 <- CompAnom %>% 
  group_by(Tag) %>% 
  summarise(CompAnom = (NumAnom/TotalFish)*100)
Summary1.12[is.na(Summary1.12)] <- 0



#The following code creates a Raw Metric summary for export.
MRSummary <- Fish %>% 
  filter(!is.na(Tag)) %>%
  group_by(Season, LocationID, Tag) %>%
  summarise()

#MR1= Species Richness excluding unknowns/hybrids, MR2 = Number Darter Spp, MR3= Number of Sunfish Spp, MR4=Number of Sucker spp, MR5= Number Intolerant Spp
#MR6= Proportion Green sunfish, MR7=Proportion of Omnivores, MR8=Proportion Insectivorous minnows, MR9= Proportion Top Piscivores
#MR10=Catch per 50ft seine haul, MR11= Proportion of Hybrids, MR12 Proportion with Anomalies
MRSummary <- cbind(MRSummary, MR1 = Summary1.1$TotSp, MR2 = Summary1.2$NumEP, MR3 = Summary1.3$NumL, MR4 = Summary1.4$NumCat, MR5 = Summary1.5$TCI,
                   MR6 = Summary1.6$CompGS, MR7 = Summary1.7$CompO, MR8 = Summary1.8$CompCI , MR9 = Summary1.9$CompTP , MR10 = Summary1.10$CatchPer50Ft,
                   MR11 = Summary1.11$CompH, MR12 = Summary1.12$CompAnom)

MSSummary <- MRSummary #Metric Raw Values (MRSummary) are placed in dataframe for Metric Score Summary (MSSummary)
MRSummary[,3] <- NULL

MRSummary$ParkCode <- site

all_MRSummary[[site]] <- MRSummary
                                                            #######STEP 2#######

#Now we'll calculate the Metric Scores to come up with a final IBI.
#Metrics 1-12

#Assigns stream order to each sample reach at HOME, PIPE, and TAPR
MSSummaryWithSO <- MSSummary %>% 
mutate(SO4 = ifelse(LocationID=='HOMEShiner01lower' | LocationID=='TAPRShiner35lower', 4, 0),
       SO3 = ifelse(LocationID=='PIPEShiner01above' | LocationID=='PIPEShiner01lower' | LocationID=='PIPEShiner01middl' |
                   LocationID=='PIPEShiner01upper' | LocationID=='TAPRShiner22lower' | LocationID=='TAPRShiner24lower' |
                   LocationID=='TAPRShiner36middl' | LocationID=='TAPRShiner36upper', 3, 0),
       SO2 = ifelse(LocationID=='TAPRShiner01lower' | LocationID=='TAPRShiner01middl' | LocationID=='TAPRShiner10middl' |
                   LocationID=='TAPRShiner12middl' | LocationID=='TAPRShiner34lower' | LocationID=='TAPRShiner22ltupp' |
                   LocationID=='TAPRShiner22rtupp' | LocationID=='TAPRShiner24middl', 2, 0),
       SO1 = ifelse(LocationID=='TAPRShiner02lower' | LocationID=='TAPRShiner04middl' | LocationID=='TAPRShiner17upper' |
                   LocationID=='TAPRShiner23middl' | LocationID=='TAPRShiner01upper' | LocationID=='TAPRShiner02middl' |
                   LocationID=='TAPRShiner05middl' | LocationID=='TAPRShiner10upper' | LocationID=='TAPRShiner17lower' |
                                        LocationID=='TAPRShiner18middl', 1, 0))

#Puts the stream order for each reach into variable called SO
MSSummaryWithSO <- MSSummaryWithSO %>% 
  mutate(SO = SO4 + SO3 + SO2 + SO1)



#Assigns a score of 5, 3, or 1 based on upper and lower thresholds for each of the 12 metrics
MSSummaryWithSO <- MSSummaryWithSO %>% 
  group_by(Tag, Season, LocationID) %>%
  summarise(MS1a = ifelse(SO==4 & MR1>12, 5, 0),
            MS1b = ifelse(SO==4 & MR1>=7 & MR1<=12, 3, 0),
            MS1c = ifelse(SO==4 & MR1<7, 1, 0),
            MS1d = ifelse(SO==3 & MR1>9, 5, 0),
            MS1e = ifelse(SO==3 & MR1>=5 & MR1<=9, 3, 0),
            MS1f = ifelse(SO==3 & MR1<5, 1, 0),
            MS1g = ifelse(SO<=2 & MR1>7, 5, 0),
            MS1h = ifelse(SO<=2 & MR1>=4 & MR1<=7, 3, 0),
            MS1i = ifelse(SO<=2 & MR1<4, 1, 0),
            MS2a = ifelse(SO==4 & MR2>2, 5, 0),
            MS2b = ifelse(SO==4 & MR2==2, 3, 0),
            MS2c = ifelse(SO==4 & MR2<2, 1, 0),
            MS2d = ifelse(SO<=3 & MR2>1, 5, 0),
            MS2e = ifelse(SO<=3 & MR2==1, 3, 0),
            MS2f = ifelse(SO<=3 & MR2<1, 1, 0),
            MS3a = ifelse(SO==4 & MR3>2, 5, 0),
            MS3b = ifelse(SO==4 & MR3==2, 3, 0),
            MS3c = ifelse(SO==4 & MR3<2, 1, 0),
            MS3d = ifelse(SO<=3 & MR3>1, 5, 0),
            MS3e = ifelse(SO<=3 & MR3==1, 3, 0),
            MS3f = ifelse(SO<=3 & MR3<1, 1, 0),
            MS4a = ifelse(SO==4 & MR4>2, 5, 0),
            MS4b = ifelse(SO==4 & MR4==2, 3, 0),
            MS4c = ifelse(SO==4 & MR4<2, 1, 0),
            MS4d = ifelse(SO<=3 & MR4>1, 5, 0),
            MS4e = ifelse(SO<=3 & MR4==1, 3, 0),
            MS4f = ifelse(SO<=3 & MR4<1, 1, 0),
            MS5a = ifelse(SO==4 & MR5>3, 5, 0),
            MS5b = ifelse(SO==4 & MR5>=2 & MR5<=3, 3, 0),
            MS5c = ifelse(SO==4 & MR5<2, 1, 0),
            MS5d = ifelse(SO==3 & MR5>2, 5, 0),
            MS5e = ifelse(SO==3 & MR5==2, 3, 0),
            MS5f = ifelse(SO==3 & MR5<2, 1, 0),
            MS5g = ifelse(SO<=2 & MR5>1, 5, 0),
            MS5h = ifelse(SO<=2 & MR5==1, 3, 0),
            MS5i = ifelse(SO<=2 & MR5<1, 1, 0),            
            MS6a = ifelse(MR6<5, 5, 0),
            MS6b = ifelse(MR6>=5 & MR6<=20, 3, 0),
            MS6c = ifelse(MR6>20, 1, 0),             
            MS7a = ifelse(MR7<20, 5, 0),
            MS7b = ifelse(MR7>=20 & MR7<=45, 3, 0),
            MS7c = ifelse(MR7>45, 1, 0),
            MS8a = ifelse(MR8>45, 5, 0),
            MS8b = ifelse(MR8<=45 & MR8>=20, 3, 0),
            MS8c = ifelse(MR8<20, 1, 0), 
            MS9a = ifelse(MR9>5, 5, 0),
            MS9b = ifelse(MR9<=5 & MR9>=1, 3, 0),
            MS9c = ifelse(MR9<1, 1, 0),
            MS10a = ifelse(MR10>80, 5, 0),
            MS10b = ifelse(MR10<=80 & MR10>=40, 3, 0),
            MS10c = ifelse(MR10<40, 1, 0),
            MS11a = ifelse(MR11==0, 5, 0),
            MS11b = ifelse(MR11<=1 & MR11>0, 3, 0),
            MS11c = ifelse(MR11>1, 1, 0),
            MS12a = ifelse(MR12==0, 5, 0),
            MS12b = ifelse(MR12<=1 & MR12>0, 3, 0),
            MS12c = ifelse(MR12>1, 1, 0))


MSSummaryWithSO <- MSSummaryWithSO %>% 
  group_by(Season, LocationID, Tag) %>%
  summarise(MS1 = (MS1a + MS1b + MS1c + MS1d + MS1e + MS1f + MS1g + MS1h + MS1i),
            MS2 = (MS2a + MS2b + MS2c + MS2d + MS2e + MS2f),
            MS3 = (MS3a + MS3b + MS3c + MS3d + MS3e + MS3f),
            MS4 = (MS4a + MS4b + MS4c + MS4d + MS4e + MS4f),
            MS5 = (MS5a + MS5b + MS5c + MS5d + MS5e + MS5f + MS5g + MS5h + MS5i),
            MS6 = (MS6a + MS6b + MS6c),
            MS7 = (MS7a + MS7b + MS7c),
            MS8 = (MS8a + MS8b + MS8c),
            MS9 = (MS9a + MS9b + MS9c),
            MS10 = (MS10a + MS10b + MS10c),
            MS11 = (MS11a + MS11b + MS11c),
            MS12 = (MS12a + MS12b + MS12c))


MSSummaryWithSO[,3] <- NULL


                                                            #######STEP 3#######

#Now let's calculate the IBI.
#Metrics 1-11

IBISummary <- MSSummaryWithSO %>%
  group_by(Season, LocationID) %>%
  summarise(IBI = (MS1 + MS2 + MS3 + MS4 + MS5 + MS6 + MS7 + MS8 + MS9 + MS10 + MS11 + MS12))

FinalSummary <- IBISummary %>% 
  mutate (Excellent = ifelse(IBI>=51,'Excellent', ''),
          Good = ifelse(IBI>=41 & IBI<=50,'Good',''),
          Fair = ifelse(IBI>=31 & IBI<=40,'Fair', ''),
          Poor = ifelse(IBI>=21 & IBI<=30,'Poor', ''),
          Very_Poor = ifelse(IBI<=20, 'Very Poor',''))
FinalSummary$IBIrating <-paste(FinalSummary$Excellent, FinalSummary$Good, FinalSummary$Fair, FinalSummary$Poor, FinalSummary$Very_Poor)

IBIScoreAndRating <-FinalSummary[,c(1,2,3,9)]

IBIScoreAndRating$ParkCode <- site

all_IBI[[site]] <- IBIScoreAndRating
}, error = function(e){
  message("Error in site", site,":", conditionMessage(e))
})
}

powerBI_data_MR <- do.call(rbind, all_MRSummary)
powerBI_data_IBI <- do.call(rbind, all_IBI)
                                                                                                          #name and address accordingly.
