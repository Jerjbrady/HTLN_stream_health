#install.packages("installr")
#library(installr)
#updateR()

                                                           ##VERSION 1.0##


                                                        ##Invert Analysis##

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
Inverts<-read.csv("C:/Users/ccheri/Documents/R/DB Exports/Invert DB Exports/Inverts/BUFF/AllCountData.csv")  #This is the raw invert replicates data.

#Before we get started, let's change any -999 values to NA's in our dataframe.
Inverts[Inverts == "-999"] <- as.numeric("NA")

#We'll also remove any cases of 'unknowns' in the data.
Inverts<-subset(Inverts, TaxonCode!="UNKN")

#The next six lines add three columns to the dataframe; one is a TV*rep count column, and the other two are Tag 
#columns so that stats can be organized and calculated according to specific replicates and riffles across years.
Inverts<-mutate(Inverts, TVxcount = Inverts$RepCount_Calc * Inverts$ToleranceValue)
rep<-data.frame(cbind(Inverts$LocationID, Inverts$Season, Inverts$RiffleNo, Inverts$Replicate))
rep$Tag<-apply(rep,1,paste,collapse="")
rif<-data.frame(cbind(Inverts$LocationID, Inverts$Season, Inverts$RiffleNo))
rif$Tag2<-apply(rif,1,paste,collapse="")
Inverts<-mutate(Inverts, Tag = rep$Tag, Tag2 = rif$Tag2)

Inverts[Inverts==""]<-NA

#The code below is creating new copies of the dataframe above, but without rows that violate the lowest taxon rule 
#for TR, and Diversity estimates.

#The following lines will keep higher taxa like Oligocheates by adding a genus and/or order value for them in the
#dataframe.
Inverts$Genus <- ifelse(Inverts$Order != "", Inverts$Genus, "blank")
Inverts$Order <- ifelse(Inverts$TaxonCode != "OLIG", Inverts$Order, "OLIG")     #will be moved below once assigned order.
Inverts$Order <- ifelse(Inverts$TaxonCode != "NEMA", Inverts$Order, "NEMA")
Inverts$Order <- ifelse(Inverts$TaxonCode != "NEMP", Inverts$Order, "NEMP")     #will be moved below once assigned order.
Inverts$Order <- ifelse(Inverts$TaxonCode != "HIR", Inverts$Order, "HIR")       #will be moved below once assigned order.
Inverts$Order <- ifelse(Inverts$TaxonCode != "BIVAL", Inverts$Order, "BIVAL")
Inverts$Order <- ifelse(Inverts$Class != "Gastropoda", Inverts$Order, "GAST")   #most Gastropods lack assigned orders.

#The following several lines of code will filter out duplicates where there is a lower classification

#Add an index so that each row can be accessed later
Inverts$index <- 1:nrow(Inverts)

#Use group_by() to find all subfamily duplicates within a date, site, riffle, replicate 
subfamilyduplicates <- Inverts %>% 
  group_by(Season, LocationID, RiffleNo, Replicate, EventID, Order, Family, Subfamily) %>% 
  summarise(subfamilyduplicates = n())

#Use group_by() to find all family duplicates within a date, site, riffle, replicate 
familyduplicates <- Inverts %>% 
  group_by(Season, LocationID, RiffleNo, Replicate, EventID, Order, Family) %>% 
  summarise(familyduplicates = n())

#Use group_by() to find all order duplicates within a date, site, riffle, replicate 
orderduplicates <- Inverts %>% 
  group_by(Season, LocationID, RiffleNo, Replicate, EventID, Order) %>% 
  summarise(orderduplicates = n())

#Merge each grouped df with original df using all = TRUE to have a column with number of duplicates,
#Then flag if there are no lower classifications
df1 <- merge(Inverts, subfamilyduplicates, all = TRUE) %>% 
  mutate(flag1 = ifelse(subfamilyduplicates > 1 & is.na(Genus), "yes", "no"))

df2 <- merge(df1, familyduplicates, all = TRUE) %>% 
  mutate(flag2 = ifelse(familyduplicates > 1 & is.na(Genus) & is.na(Subfamily), "yes", "no"))

df3 <- merge(df2, orderduplicates, all = TRUE) %>% 
  mutate(flag3 = ifelse(orderduplicates > 1 & is.na(Genus) & is.na(Family) & is.na(Subfamily), "yes", "no"))

#Remove flagged rows 
Inverts2 <- df3 %>% filter(flag1 == "no") %>% 
  filter(flag2 == "no") %>% 
  filter(flag3 == "no") 

Inverts2

#The following six lines are connected and run together to produce a summary dataframe of TR and Diversity 
#stats for every replicate.
Replicate_summary <- Inverts2 %>%
  filter(!is.na(Inverts2$Tag)) %>%
  group_by(Inverts2$LocationID, Inverts2$Season, Inverts2$Tag, Inverts2$Tag2) %>%      
  summarize(Taxa_Richness = n(),
            Shannon = diversity(RepCount_Calc, index = "shannon"))

#The following lines of code are messy but the only way I've yet figured out how to accurately calculate EPT ratio
#and bind it to the dataframe we produced above to create our complete summary of replicate stats.
Denominator_Taxa <- Inverts %>%
  filter(!is.na(Inverts$Tag)) %>%
  summarize(EPT2 = aggregate(RepCount_Calc~Order + Tag, subset(Inverts, Order=='Ephemeroptera' | Order=='Plecoptera' | Order=='Trichoptera' | Family=='Chironomidae'), sum))

EPTd <- aggregate(Denominator_Taxa$EPT2$RepCount_Calc~Denominator_Taxa$EPT2$Tag, Denominator_Taxa$EPT2, sum)

Numerator_Taxa <- subset(Denominator_Taxa, EPT2$Order=='Ephemeroptera' | EPT2$Order=='Plecoptera' | EPT2$Order=='Trichoptera')
EPTn <- aggregate(Numerator_Taxa$EPT2$RepCount_Calc~Numerator_Taxa$EPT2$Tag, Numerator_Taxa, sum)

names(EPTd) <- c('Tag', 'Repcount_Calcd')
names(EPTn) <- c('Tag', 'Repcount_Calcn')

EPTA = merge(EPTd,EPTn,all=T) #This step is absolutely necessary in cases where there are more observations in the 
                              #denominator than the numerator.  
EPTr <- as.data.frame(EPTA[,3]/EPTA[,2])

#Sometimes there are no EPT taxa in the sample, so we produce NA's. We'll convert those to 0 since there's no ratio.
EPTr[is.na(EPTr)] <- 0

#Now let's find EPT Richness for each replicate.
EPTTR <- Inverts2 %>% group_by(Inverts2$Tag) %>% tally(Order=='Ephemeroptera' | Order=='Plecoptera' | Order=='Trichoptera')

#Then we'll find Shannon Evenness for each replicate.
HE <- Replicate_summary %>% 
  filter(!is.na(`Inverts2$Tag`)) %>% 
  group_by(`Inverts2$Tag`) %>% summarise((Shannon)/(log(Taxa_Richness)))

#And now we'll calculate HBI which first requires us to create a small dataframe from our original data that has rows
#without Tolerance Values removed.
A <- as.data.frame(cbind(Inverts$LocationID, Inverts$Season, Inverts$Tag, Inverts$Tag2)) #R fussy, won't bind all 
B <- cbind(Inverts$RepCount_Calc, Inverts$TVxcount)                                      #col w/o making them chars
InvertsHBI <- cbind(A, B)
names(InvertsHBI) <- c('Location', 'Season','Tag', 'Tag2', 'RepCount_Calc', 'TVxcount')
InvertsHBI <- na.omit(InvertsHBI)

HBI_summary <- InvertsHBI %>%
  filter(!is.na(InvertsHBI$Tag)) %>%
  group_by(InvertsHBI$Location, InvertsHBI$Season, InvertsHBI$Tag, InvertsHBI$Tag2) %>%      
  summarize(HBI = sum(TVxcount)/sum(RepCount_Calc))

#Finally, we can bind the last the replicate stats that we calculated to our Replicate summary dataframe.
Final_Rep_summary1 <- cbind(Replicate_summary, HBI_summary$HBI, EPTTR$n, HE[,2])

#EPTr might have a different number of rows if a sample doesn't have any EPT or chironomids so we'll bind it separately.
newdf = cbind(EPTA,EPTr)
names(newdf) <- c('Inverts2$Tag','d','n','EPTr')
newdf$d <- NULL
newdf$n <- NULL
Final_Rep_summary2 <- merge(Final_Rep_summary1,newdf,all=T)

#These two lines of code change the column names and spit out the stats in the console.
names(Final_Rep_summary2) <- c('Tag', 'Location', 'Season', 'Tag2', 'TR', 'Shannon', 'HBI', 'EPTTR', 'HE', 'EPTr')  #note for Cameron >> !!!!!Check Col names when using with other data (i.e. habitat) if changing Tag column name!!!!!!
#Replace any EPTr NA's with 0.
Final_Rep_summary2 <- replace(Final_Rep_summary2, is.na(Final_Rep_summary2), 0)
Final_Rep_summary2

                                                    #########################
                                                    ##RIFFLE STATS ANALYSIS##
                                                    #########################

#Now we'll remove the Tag column because we are done with replicate stats and moving on to riffle stats.
Final_Rep_summary2$Tag <- NULL

#Similar to before, the following code will filter the replicate stats and produce a riffle stats summary.
Riffle_summary <- Final_Rep_summary2 %>%
  filter(!is.na(Tag2)) %>%
  group_by(Location, Season, Tag2) %>%      
  summarise(TR = mean(TR),
            Shannon = mean(Shannon),
            HE = mean(HE),
            HBI = mean(HBI),
            EPTr = mean(EPTr),
            EPTTR = mean(EPTTR))

#Now all that's left to do is find the final stats (averages, standard error, variance, etc.) for the each season.
#As before, we'll start by removing the Tag2 column, because we're summarizing all riffle values within a season.
Riffle_summary$Tag2 <- NULL

#We need to create formula for standard error
se <- function(x) sqrt(var(x) / length(x)) # or se <- function(x) sd(x)/sqrt(length(x))

#We also need to create a function to calculate 95% confidence intervals.
z = qnorm(0.025)
conf <- function(x) abs(z * se(x))

#The last several lines of script below give us our final summary statistics for all Invert count data from our site
#across all years.
Invert_summary <- Riffle_summary %>%
  group_by(Location, Season) %>%
  summarize(TR_Avg = mean(TR), TR_Error = se(TR), TR_Stdev = sd(TR), TR_Min = min(TR), TR_Max = max(TR), TR_conf = conf(TR),
            Shannon_Avg = mean(Shannon), Shannon_Error = se(Shannon), Shannon_Stdev = sd(Shannon), Shannon_Min = min(Shannon), Shannon_Max = max(Shannon), Shannon_conf = conf(Shannon),
            HE_Avg = mean(HE), HE_Error = se(HE), HE_Stdev = sd(HE), HE_Min = min(HE), HE_Max = max(HE), HE_conf = conf(HE),
            HBI_Avg = mean(HBI), HBI_Error = se(HBI), HBI_Stdev = sd(HBI), HBI_Min = min(HBI), HBI_Max = max(HBI), HBI_conf = conf(HBI),
            EPTr_Avg = mean(EPTr), EPTr_Error = se(EPTr), EPTr_Stdev = sd(EPTr), EPTr_Min = min(EPTr), EPTr_Max = max(EPTr), EPTr_conf = conf(EPTr),
            EPTTR_Avg = mean(EPTTR), EPTTR_Error = se(EPTTR), EPTTR_Stdev = sd(EPTTR), EPTTR_Min = min(EPTTR), EPTTR_Max = max(EPTTR), EPTTR_conf = conf(EPTTR))

#Let's see those final stats!
Invert_summary

#Lastly we'll save the summary as a .csv file so that we can load it into a new R script that we'll use to make
#various plots.

write.csv(Invert_summary, "C:/Users/ccheri/Documents/R/Summaries/Park Invert Summaries/BUFF_InvertsummaryALL.csv") #You will need to change the 
#name and address accordingly.

#############################################################################################################################

#The following code creates a dataframe with the total counts for all taxa at each site per year and the percent composition of 
#each taxon at each site per year.You need the "Inverts2" dataframe from the stats summary code above to continue.

#Here we create numerator values to find the ratio of each taxon at a site.
Inverts2numerators <- Inverts2 %>%
  filter(!is.na(Tag)) %>%
  reframe(count = aggregate(RepCount_Calc~LocationID + Season + TaxonCode, subset(Inverts2, TaxonCode!=''), sum))
Inverts2numerators <- Inverts2numerators %>%
  group_by(count$LocationID, count$Season) %>%
  reframe(count$Season, count$TaxonCode, count$RepCount_Calc)
names(Inverts2numerators) <- c('Location','Season','TaxonCode','numcount')

#Next we create a denominator value of the total site count.
Inverts2denominators <- Inverts2 %>%
  reframe(count = aggregate(RepCount_Calc~LocationID + Season, subset(Inverts2, TaxonCode!=''), sum))
Inverts2denominators <- Inverts2denominators %>%
  group_by(count$LocationID, count$Season) %>%
  reframe(count$Season, count$RepCount_Calc)
names(Inverts2denominators) <- c('Location','Season','denomcount')

#Merge the numerator and denominator counts into one dataframe.
Inverts2ratio = merge(Inverts2denominators,Inverts2numerators,all=T)

#Divide the numerator value by the denominator value for each taxon per site per year.
Inverts2ratio2 <- as.data.frame(Inverts2ratio[,5]/Inverts2ratio[,3])

#Bind this ratio to the dataframe with numerator values which are also the total individual taxon counts that we also want.
Inverts2ratiofinal = cbind(Inverts2ratio,Inverts2ratio2)

#Remove denominator counts
Inverts2ratiofinal = Inverts2ratiofinal[,c(1:2,4:6)]
names(Inverts2ratiofinal) <- c('Location','Season','TaxonCode','Taxoncount','Ratio')

#Express the ratio as a percentage.
UniqueTaxaStats <- Inverts2ratiofinal %>%
  group_by(Location, Season, TaxonCode, Taxoncount) %>%
  reframe(Percent = (Ratio*100))

#The following code can be used to quickly look at percent composition for any particular Taxon in the UniqueTaxaStats dataframe.
#####################################
#    Taxon <- UniqueTaxaStats %>%
#     reframe(count = aggregate(Percent~Location + Season + TaxonCode, subset(UniqueTaxaStats, TaxonCode=='CHIR'), sum))

#    Taxon <- Taxon %>%
#     group_by(count$Location, count$Season, count$TaxonCode) %>%
#     reframe(count$Location, count$Season, count$TaxonCode, count$Percent)
#    names(Taxon) <- c('Location','Season','TaxonCode','Percent')
#####################################

#Let's see those final stats!
UniqueTaxaStats

#Lastly we'll save the Taxon counts and ratios as a .csv file so that we can load it into a new R script that we'll use to 
#make various plots.

write.csv(UniqueTaxaStats, "C:/Users/ccheri/Documents/R/Park Invert Summaries/BUFF_UniqueTaxaStats.csv") #You will need to change the 
#name and address accordingly.
