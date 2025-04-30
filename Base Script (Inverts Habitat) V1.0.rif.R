#install.packages("installr")
#library(installr)
#updateR()

                                                 ##VERSION 1.0##


                                              ##Habitat Analysis##

#After installing packages on your PC once, you do not need to run the packages again, but you do need to run
#the libraries every time you reopen RStudio.
#install.packages('vegan')
library(vegan)
#install.packages("dplyr")
library(dplyr)                       #####################################
                                     ##Start here after loading packages##
                                     #####################################

                                       #INDIVIDUAL SAMPLE HABITAT DATA#

#Copy the file address for the exported Database spreadsheets here. Make sure they are saved as .csv files.
Hab<-read.csv("C:/Users/ccheri/Documents/R/DB Exports/Invert DB Exports/Habitat/BUFF/AllExportSamplesData.csv")  #This is the raw riffle habitat data.

#Before we get started, let's change any -999 values to NA's in our dataframe.
Hab[Hab == "-999"] <- as.numeric("NA")

#The next three lines add a Tag column to the dataframe which will allow stats to be organized and averaged 
#according to specific riffles across years.
new<-data.frame(cbind(Hab$Season, Hab$RiffleNo))
new$Tag<-apply(new,1,paste,collapse="")
Hab<-mutate(Hab, Tag = new$Tag)

#The following code will filter the replicate stats and produce a riffle mean stats summary.
Rif_summary <- Hab %>%
  filter(!is.na(Tag)) %>%
  group_by(Tag, Season, LocationID) %>%      
  summarize(Depth = mean(Depth_cm),
            Velocity = mean(Velocity_ms),
            Substrate = mean(Substrate),
            Sub_size = mean(Substrate_mm),
            Embeddedness = mean(Embeddedness),
            Embed_cover = mean(Embeddedness.),
            Vegetation = mean(Vegetation),
            Veg_cover= mean(Vegetation.),
            Periphyton= mean(Periphyton),
            Peri_cover= mean(Periphyton.))

#Now all that's left to do is find the final stats (averages, standard error, variance, etc.) for the each season.
#We'll start by removing the Tag column, because we're summarizing all riffle values within a season.
Rif_summary$Tag <- NULL

#We need to create formula for standard error
se <- function(x) sqrt(var(x) / length(x)) # or se <- function(x) sd(x)/sqrt(length(x))

#We also need to create a function to calculate 95% confidence intervals.
z = qnorm(0.025)
conf <- function(x) abs(z * se(x))

#The following many lines of script below give us our final summary statistics for all Riffle Habitat data from our 
#sites across all years.
Hab_summary <- Rif_summary %>%
  filter(!is.na(Season)) %>%
  group_by(LocationID, Season) %>%
  summarize(Depth_Avg = mean(Depth), Depth_Error = se(Depth), Depth_Stdev = sd(Depth), Depth_Min = min(Depth), 
            Depth_Max = max(Depth), Depth_conf = conf(Depth), #Depth_count = n(Depth),
            Velocity_Avg = mean(Velocity), Velocity_Error = se(Velocity), Velocity_Stdev = sd(Velocity), 
            Velocity_Min = min(Velocity), Velocity_Max = max(Velocity), Velocity_conf = conf(Velocity), #Velocity_count = n(Velocity),
            Substrate_Avg = mean(Substrate), Substrate_Error = se(Substrate), Substrate_Stdev = sd(Substrate), 
            Substrate_Min = min(Substrate), Substrate_Max = max(Substrate), Substrate_conf = conf(Substrate), #Substrate_count = n(Substrate),
            Sub_size_Avg = mean(Sub_size), Sub_size_Error = se(Sub_size), Sub_size_Stdev = sd(Sub_size), 
            Sub_size_Min = min(Sub_size), Sub_size_Max = max(Sub_size), Sub_size_conf = conf(Sub_size), #Sub_size_count = n(Sub_size),
            Embeddedness_Avg = mean( Embeddedness),  Embeddedness_Error = se( Embeddedness),  
            Embeddedness_Stdev = sd( Embeddedness),  Embeddedness_Min = min( Embeddedness),  
            Embeddedness_Max = max( Embeddedness), Embeddedness_conf = conf(Embeddedness), #Embeddedness_count = n(Embeddedness),
            Embed_cover_Avg = mean(Embed_cover), Embed_cover_Error = se(Embed_cover), 
            Embed_cover_Stdev = sd(Embed_cover), Embed_cover_Min = min(Embed_cover), 
            Embed_cover_Max = max(Embed_cover), Embed_cover_conf = conf(Embed_cover), #EPTTR_count = n(EPTTR),
            Vegetation_Avg = mean(Vegetation), Vegetation_Error = se(Vegetation), Vegetation_Stdev = sd(Vegetation), 
            Vegetation_Min = min(Vegetation), Vegetation_Max = max(Vegetation), Vegetation_conf = conf(Vegetation), #Vegetation_count = n(Vegetation),
            Veg_cover_Avg = mean(Veg_cover), Veg_cover_Error = se(Veg_cover), Veg_cover_Stdev = sd(Veg_cover), 
            Veg_cover_Min = min(Veg_cover), Veg_cover_Max = max(Veg_cover), Veg_cover_conf = conf(Veg_cover), #Veg_cover_count = n(Veg_cover),
            Periphyton_Avg = mean(Periphyton), Periphyton_Error = se(Periphyton), Periphyton_Stdev = sd(Periphyton), 
            Periphyton_Min = min(Periphyton), Periphyton_Max = max(Periphyton), Periphyton_conf = conf(Periphyton), #Periphyton_count = n(Periphyton),
            Peri_cover_Avg = mean(Peri_cover), Peri_cover_Error = se(Peri_cover), Peri_cover_Stdev = sd(Peri_cover), 
            Peri_cover_Min = min(Peri_cover), Peri_cover_Max = max(Peri_cover),Peri_cover_conf = conf(Peri_cover)) #Peri_cover_count = n(Peri_cover))

#Let's see those final stats!
Hab_summary

#Lastly we'll save the summary as a .csv file so that we can load it into a new R script that we'll use to make
#various plots.

write.csv(Hab_summary, "C:/Users/ccheri/Documents/R/Summaries/Park Habitat Summaries/BUFF_HabsummaryALL.csv") #You will need to change the name and
#address accordingly.


                                                  #RIFFLE HABITAT DATA#

#Copy the file address for the exported Database spreadsheets here. Make sure they are saved as .csv files.
Rif<-read.csv("C:/Users/ccheri/Documents/R/DB Exports/Invert DB Exports/Habitat/BUFF/AllExportRiffleData.csv")  #This is the raw riffle habitat data.

#Before we get started, let's change any -999 values to NA's in our dataframe.
Rif[Rif == "-999"] <- as.numeric("NA")

#We need to create formula for standard error
se <- function(x) sqrt(var(x) / length(x)) # or se <- function(x) sd(x)/sqrt(length(x))

#We also need to create a function to calculate 95% confidence intervals.
z = qnorm(0.025)
conf <- function(x) abs(z * se(x))

#The following many lines of script below give us our final summary statistics for all Riffle Habitat data from our 
#sites across all years.
Riff_summary <- Rif %>%
  group_by(LocationID, Season) %>%      
  summarize(Temp_Avg = mean(Temperature_C), Temp_Error = se(Temperature_C), Temp_Stdev = sd(Temperature_C), 
            Temp_min = min(Temperature_C), Temp_max = max(Temperature_C), Temp_conf = conf(Temperature_C),
            Cond_Avg = mean(SpCond_microS), Cond_Error = se(SpCond_microS), Cond_Stdev = sd(SpCond_microS), 
            Cond_min = min(SpCond_microS), Cond_max = max(SpCond_microS), Cond_conf = conf(SpCond_microS),
            pH_Avg = mean(pH), pH_Error = se(pH), pH_Stdev = sd(pH), pH_min = min(pH), pH_max = max(pH), pH_conf = conf(pH),
            DO_Avg = mean(DO_mgperL), DO_Error = se(DO_mgperL), DO_Stdev = sd(DO_mgperL),
            DO_min = min(DO_mgperL), DO_max = max(DO_mgperL), DO_conf = conf(DO_mgperL))

#Let's see those final stats!
Riff_summary

#Lastly we'll save the summary as a .csv file so that we can load it into a new R script that we'll use to make
#various plots.

write.csv(Riff_summary, "C:/Users/ccheri/Documents/R/Summaries/Park Habitat Summaries/BUFF_RifsummaryALL.csv") #You will need to change the name and
#address accordingly.






#Additionally we can create smaller, specific tables for easier summary stat investigation.

#Smaller summary for depth and velocity.
Flow_summary <- Rif_summary %>%
  filter(!is.na(Season)) %>%
  group_by(Season) %>%
  summarize(Depth_Avg = mean(Depth), Depth_Error = se(Depth), Depth_Stdev = sd(Depth), Depth_Min = min(Depth), 
            Depth_Max = max(Depth), #Depth_count = n(Depth),
            Velocity_Avg = mean(Velocity), Velocity_Error = se(Velocity), Velocity_Stdev = sd(Velocity), 
            Velocity_Min = min(Velocity), Velocity_Max = max(Velocity)) #Velocity_count = n(Velocity),

#Smaller summary for substrate.
Substrate_summary <- Rif_summary %>%
  filter(!is.na(Season)) %>%
  group_by(Season) %>%
  summarize(Substrate_Avg = mean(Substrate), Substrate_Error = se(Substrate), Substrate_Stdev = sd(Substrate), 
            Substrate_Min = min(Substrate), Substrate_Max = max(Substrate), #Substrate_count = n(Substrate),
            Sub_size_Avg = mean(Sub_size), Sub_size_Error = se(Sub_size), Sub_size_Stdev = sd(Sub_size), 
            Sub_size_Min = min(Sub_size), Sub_size_Max = max(Sub_size)) #Sub_size_count = n(Sub_size),
            
#Smaller summary for embeddedness.
Embeddedness_summary <- Rif_summary %>%
  filter(!is.na(Season)) %>%
  group_by(Season) %>%
  summarize(Embeddedness_Avg = mean( Embeddedness),  Embeddedness_Error = se( Embeddedness),  
            Embeddedness_Stdev = sd( Embeddedness),  Embeddedness_Min = min( Embeddedness),  
            Embeddedness_Max = max( Embeddedness), # Embeddedness_count = n( Embeddedness)
            Embed_cover_Avg = mean(Embed_cover), Embed_cover_Error = se(Embed_cover), 
            Embed_cover_Stdev = sd(Embed_cover), Embed_cover_Min = min(Embed_cover), 
            Embed_cover_Max = max(Embed_cover)) #EPTTR_count = n(EPTTR)
            
#Smaller summary for vegetation.
Veg_summary <- Rif_summary %>%
  filter(!is.na(Season)) %>%
  group_by(Season) %>%
  summarize(Vegetation_Avg = mean(Vegetation), Vegetation_Error = se(Vegetation), Vegetation_Stdev = sd(Vegetation), 
            Vegetation_Min = min(Vegetation), Vegetation_Max = max(Vegetation), #Vegetation_count = n(Vegetation),
            Veg_cover_Avg = mean(Veg_cover), Veg_cover_Error = se(Veg_cover), Veg_cover_Stdev = sd(Veg_cover), 
            Veg_cover_Min = min(Veg_cover), Veg_cover_Max = max(Veg_cover)) #Veg_cover_count = n(Veg_cover),
            
#Smaller summary for periphyton.
Periphyton_summary <- Rif_summary %>%
  filter(!is.na(Season)) %>%
  group_by(Season) %>%
  summarize(Periphyton_Avg = mean(Periphyton), Periphyton_Error = se(Periphyton), Periphyton_Stdev = sd(Periphyton), 
            Periphyton_Min = min(Periphyton), Periphyton_Max = max(Periphyton), #Periphyton_count = n(Periphyton),
            Peri_cover_Avg = mean(Peri_cover), Peri_cover_Error = se(Peri_cover), Peri_cover_Stdev = sd(Peri_cover), 
            Peri_cover_Min = min(Peri_cover), Peri_cover_Max = max(Peri_cover)) #Peri_cover_count = n(Peri_cover))
            














