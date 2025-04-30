#install.packages("installr")
#library(installr)
#updateR()

                                                    ##EF Fish Community Analysis##

#After installing packages on your PC once, you do not need to run the packages again, but you do need to run
#the libraries every time you reopen RStudio.
#install.packages('vegan')
library(vegan)
#install.packages("tidyverse")
library(tidyverse)                   
                                                #####################################
                                                ##Start here after loading packages##
                                                #####################################

#Copy the file address for the exported Database spreadsheets here. Make sure they are saved as .csv files.
Fish<-read.csv("C:/Users/ccheri/Documents/R/DB Exports/Fish DB Exports/OZAR Fish/FishCommSppIndiv_EF.csv")  #This is the raw fish species data.
Effort<-read.csv("C:/Users/ccheri/Documents/R/DB Exports/Fish DB Exports/OZAR Fish/FishCommSampling_EF.csv") #This contains sampling effort.
SiteWidth<-read.csv("C:/Users/ccheri/Documents/R/DB Exports/Fish DB Exports/OZAR Fish/CrossSectionalInfo_EF.csv") #This contains reach width.

#--------------------------------
#First section of code preps the data and creates functions for metric calculations

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
new3<-data.frame(cbind(Fish$Season, Fish$LocationID, Fish$SamplingGear, Fish$ChannelType, Fish$TaxonCode, Fish$TaxaFishNumber, Fish$BatchID))
new3$Tag<-apply(new3,1,paste,collapse="")
Fish<-mutate(Fish, Tag2 = new3$Tag)
new4<-data.frame(cbind(SiteWidth$Season, SiteWidth$LocationID))
new4$Tag<-apply(new4,1,paste,collapse="")
SiteWidth<-mutate(SiteWidth, Tag = new4$Tag)

#We need to find the average width of every reach. First we'll remove rows without width measurements and then average by reach.
SiteWidth <- SiteWidth %>% filter(Width_m != 'NA')
SiteWidth <- SiteWidth %>% 
  group_by(Tag) %>% 
  summarise(Reachwidth = mean((Width_m)))

#Simplifying the Effort dataframe by reducing its columns to just the Tag, Reach length and effort in seconds.
Effort <- Effort[,c(6,9,12)]

#Now we'll combine the total effort in seconds for all sampling bouts throughout each reach and convert the total seconds to minutes.
Effort <- Effort %>% 
  group_by(Tag) %>% 
  summarise(minutes = sum(SamplingEffort_sec)/60,
            ReachLength = ReachLength_m)

#Condense further to one total sampling effort for the entire LocationID and Year (i.e. Tag)
Effort <- Effort %>% 
  group_by(Tag, minutes, ReachLength) %>% 
  summarise()

#Simplifying the Fish dataframe by reducing its columns to just those variables needed; takes out FederalProtection; IMDProtection; 
#FlaggedData, OpsReviewFlaggedData_Comment and Timestamp.
FishSimplified <- Fish[,c(1:7,10:19,23:24)]

#Now we'll add the effort in minutes to our Simplified Fish dataframe.
FishWithEffort = merge(FishSimplified, Effort, all=T) 

#Next we'll add the average width of each reach so we can calculate area. 
#We'll call this dataframe FishComplete since it has all the raw data needed to calculate metrics.
FishComplete = merge(FishWithEffort, SiteWidth, all=T)
FishComplete$Area <- c(FishComplete$ReachLength*FishComplete$Reachwidth)

#Before we begin calculating the first metrics, let's create a dataframe that will keep track of the scientific and common names 
#that are associated with the taxon codes which we will be using to keep track of species-specific metrics. We'll bind this 
#dataframe to our final output dataframe with species-specific metrics at the end.
fishnames <- Fish[,c(6,7,10)]

fishnames <- fishnames %>%
  group_by(TaxonCode, CommonName, ScientificName) %>%      
  summarize()

#We'll need to create a simpson diversity formula that differs from the default Vegan package formula.
mydiversity <- function(x){
  N <- sum(x)
  1 - sum(((x^2-x)/(N^2-N)))
}

#We need to create formula for standard error
se <- function(x) sqrt(var(x) / length(x)) # or se <- function(x) sd(x)/sqrt(length(x))

#We also need to create a function to calculate 95% confidence intervals.
z = qnorm(0.025)
conf <- function(x) abs(z * se(x))


#-----------------------------------------------
#This section calculates Site Metrics, i.e. the metrics for each Location and Year: 
#Species Richness, Total Catch, Total CPUE, and Diversity

#The following several lines are connected and run together to produce a summary dataframe with Species Richness, 
#Total number of fish caught, and total CPUE at a location and year.
Site_summary <- FishComplete %>%
  filter(!is.na(Tag)) %>%
  group_by(LocationID, Season, Tag) %>%      
  summarize(Taxa_Richness = n_distinct(TaxonCode),
            TotalFish = sum(NumObs), 
            TotalCPUE = TotalFish/minutes)

Site_summary <- Site_summary %>% 
  group_by(LocationID, Season, Tag, Taxa_Richness, TotalFish, TotalCPUE) %>% 
  summarise()

#We need the total number of individuals per species. (i.e. species counts)
SpTot <- FishComplete %>%
  filter(!is.na(Tag)) %>%
  group_by(LocationID, Season, Tag, TaxonCode) %>%      
  summarize(SpTot = sum(NumObs))

#Use number of individuals per species to calculate Simpsons Diversity
SimpDiv <- SpTot %>%
  filter(!is.na(Tag)) %>%
  group_by(LocationID, Season, Tag) %>%      
  summarize(Diversity = mydiversity(SpTot))

#Add Simpsons Diversity to the other site metrics in Site_Summary and call it FishSiteMetrics
#this dataframe is printed to a csv file at the end of this R script
SimpDiv [,1:2] <- NULL
FishSiteMetrics = merge(Site_summary, SimpDiv, all=T)


#-----------------------------------------------------
#Calculate Species Metrics for each Location and Year: Species CPUE, Ave Length, Ave Weight, Composition by Biomass, Composition By Abundance

#Bind Species counts (i.e number individuals per species) to the Fish dataframe that has the Effort and Area included in dataframe.
SpTot[,1:2] <- NULL
FishComplete = merge(FishComplete, SpTot, all=T)

#Now we need to replace 'NA' Weights with 0's and combine them into one column.
FishWeightsNA <- FishComplete %>% replace_na(list(Weight_g=0))
FishWeightsNA <- FishWeightsNA %>% replace_na(list(BatchWT_g=0))

#We'll create a new column that combines all weights.
FishWeightsNA$WeightNew = (FishWeightsNA$Weight_g + FishWeightsNA$BatchWT_g)

#Now we'll create a new column in a new dataframe that we'll merge to our original which will give every Batch a "Batchsum" which
#will act as a flagging column for individual fishes in a batch that do or don't have batch weights associated with them.
FishBatchSum <- FishWeightsNA %>%
  filter(!is.na(Tag)) %>%
  group_by(Season, LocationID, Tag, TaxonCode, BatchID) %>%      
  summarize(Batchsum = sum(BatchWT_g),
            Tag2 = Tag2)

FishBatchSum <- FishBatchSum[,c(6,7)]
FishWT = merge(FishWeightsNA, FishBatchSum, all=T)

#Add an index so that each row can be accessed later.
FishWT$index <- 1:nrow(FishWT)

#Use group_by() to order some important columns. 
Keep1 <- FishWT %>% 
  group_by(LocationID, Season, Tag, TaxonCode, TaxaFishNumber, NumObs, BatchID, WeightNew, SpTot) %>% 
  summarise(Keep1 = n())

#Merge each grouped df with a new df using all = TRUE to have a columns with a new value for each row.
#These values are titled as 'flags' and have a 'yes' or 'no' value depending on whether the fish in each row should be included 
#in weight estimates as long as it doesn't match the preceeding 'ifelse' statement.

#gives a no (flags record) if number observed is not 1, meaning it's a count or no fish were collected
df1 <- merge(FishWT, Keep1, all = TRUE) %>% 
  mutate(flag1 = ifelse(NumObs!='1', "no", "yes"))

#gives a no if individual weights are missing
df2 <- merge(df1, Keep1, all = TRUE) %>% 
  mutate(flag2 = ifelse(NumObs=='1' & BatchID=='0' & Weight_g=='0', "no", "yes"))

#gives a no if batch weight is missing
df3 <- merge(df2, Keep1, all = TRUE) %>% 
  mutate(flag3 = ifelse(BatchID!='0' & Weight_g=='0' & Batchsum=='0', "no", "yes"))

#Remove flagged rows 
FishMinusFlagged <- df3 %>% filter(flag1 == 'yes') %>% 
  filter(flag2 == "yes") %>% 
  filter(flag3 == "yes")

#we can average the weight values now that those flagged records above are removed from the dataframe.
Weight_summary <- FishMinusFlagged %>%
  filter(!is.na(Tag)) %>%
  group_by(LocationID, Season, Tag, TaxonCode) %>%      
  summarize(Weights = mean(WeightNew))

#Bind dataframe with average weights to the species count dataframe to calculate average species biomass below
#This also allows you to keep species that don't have any weights associated with them; i.e. no individual weights or batch wts
SpBiomass = merge(SpTot, Weight_summary, all=T)

#Calculates average biomass for each species at a site and year.
AveBiomass_summary <- SpBiomass %>%
  filter(!is.na(Tag)) %>%
  group_by(LocationID, Season, Tag, TaxonCode) %>%      
  summarize(Biomass = Weights*SpTot)

#Gets the site area for each site and year to be merged with the AveBiomass data below
SiteArea <- FishComplete %>% 
  group_by(LocationID, Season, Tag, Area) %>% 
  summarize()

#Creates dataframe that has the average species biomass with the site area
AveBiomassWithArea = merge(AveBiomass_summary, SiteArea, all=T)

#Calculates the species Biomass by area.
BiomassPerArea <- AveBiomassWithArea %>%
  filter(!is.na(Tag)) %>%
  group_by(LocationID, Season, Tag, TaxonCode) %>%   
  summarize(BiomassA = Biomass/Area)

#Calculate Total Fish Biomass for a site and year
TotalBiomassPerArea <- BiomassPerArea %>%
  filter(!is.na(Tag)) %>%
  group_by(Tag) %>%      
  summarize(BiomassAsum = sum(BiomassA, na.rm = TRUE))

#Adds Total fish biomass for a site and year to the Species Biomass per Area dataframe
BiomassPerArea = merge(BiomassPerArea, TotalBiomassPerArea, all=T)

#Finally we'll calculate Composition Biomass by area for each species.
PerComp_Biomass <- BiomassPerArea %>%
  filter(!is.na(Tag)) %>%
  group_by(LocationID, Season, Tag, TaxonCode) %>%      
  summarize(CompBio = (BiomassA/BiomassAsum)*100)

#Now we'll produce another summary dataframe to calculate mean species lengths.
SpAveLenth <- FishComplete %>%
  filter(!is.na(Tag)) %>%
  group_by(LocationID, Season, Tag, TaxonCode) %>%      
  summarize(Length = mean(na.omit(TotalLength_mm)),
            Length_SD = sd(na.omit(TotalLength_mm)),
            Length_SE = se(na.omit(TotalLength_mm)),
            Length_CI = conf(na.omit(TotalLength_mm)))

#Need to merge FishSiteMetrics which has total fish count with FishComplete to calculate Composition by Abundance below
FishComplete = merge(FishComplete, FishSiteMetrics)

#This summary will give us percent composition by abundance and species specific CPUE.
CompAbund_SpCPUE <- FishComplete %>%
  filter(!is.na(Tag)) %>%
  group_by(LocationID, Season, Tag, TaxonCode) %>%
  summarize(CompAbu = (SpTot/TotalFish)*100,
            SpCPUE = (SpTot/minutes))

#condense the dataframe to one record per species at a site and year
CompAbund_SpCPUE <- CompAbund_SpCPUE %>% 
  group_by(LocationID, Season, Tag, TaxonCode, CompAbu, SpCPUE) %>% 
  summarise()

#Now we combine all dataframes that have species metrics
#Using the all.x=T for the second merge below allows you to do a left join, keeping LocationID and Season
SpMetrics = merge(SpBiomass,PerComp_Biomass, all=T)
SpMetrics [,3:4] <- NULL
SpMetrics2 = merge(SpAveLenth,SpMetrics, all=T)
SpMetrics3 = merge(SpMetrics2,CompAbund_SpCPUE, all=T)

#Now we'll merge that fishnames dataframe to the species-specific metrics dataframe.
SpMetricsWithNames = merge(SpMetrics3, fishnames,  all=T)

#Here we'll reorder columns and use just the columns we need.
SpMetricsWithNames <- SpMetricsWithNames %>%
  group_by(Season, LocationID, TaxonCode, ScientificName, CommonName, SpTot, SpCPUE, Length, Length_SD, Length_SE, Length_CI, Weights, 
           CompAbu, CompBio) %>%      
  summarize()

FishSiteMetrics <- FishSiteMetrics %>%
  group_by(Season, LocationID, Taxa_Richness, Diversity, TotalFish, TotalCPUE) %>%      
  summarize()

#Lastly we'll save these two summaries as .csv files so that we can load them into a new R script that we'll use to make various plots.
write.csv(FishSiteMetrics, "C:/Users/ccheri/Documents/R/Summaries/Park Fish Summaries/OZAR/OZAR_fish_sitemetrics.csv") #You will need to
write.csv(SpMetricsWithNames, "C:/Users/ccheri/Documents/R/Summaries/Park Fish Summaries/OZAR/OZAR_fish_speciesmetrics.csv") #change the name and 
                                                                                                                             #address accordingly.
